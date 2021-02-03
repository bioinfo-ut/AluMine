
#!/bin/bash

# Description: This is a pipeline for discovery of REF-plus polymorphic ALU
# element insertions from human reference genome.
# Authors: Tarmo Puurand and Maido Remm, 2018-2019

# Dependencies
# blastall from NCBI legacy BLAST package ftp://ftp.ncbi.nlm.nih.gov/blast/executables/legacy.NOTSUPPORTED/
# glistmaker and glistquery from GenomeTester4 package https://github.com/bioinfo-ut/GenomeTester4/

# Define paths
human_chr_path="/storage9/db/human_37/data/chr"
chimp_chr_path="/storage9/db/chimp_214/data/chr"
chimp_list_path="/storage9/db/kmer_lists/" # Must be writable
chimp_list_filename="pan_troglodytes_32.list"

# Create 32-mer list from the chimpanzee genome, if it does not yet exist.
# NB! The list file size is ca 30 GB
chimp_list="${chimp_list_path}/${chimp_list_filename}"
if [ ! -f $chimp_list ];then
	if [ ! -x glistmaker ];then
		cp ../bin/glistmaker .
		chmod 755 glistmaker
	fi
	echo "Chimp list does not exist. Generating..."
	./glistmaker "${chimp_chr_path}/*fa" -w 32 -o pan_troglodytes
fi

if [ ! -x glistquery ];then
	cp ../bin/glistquery .
	chmod 755 glistquery
fi

# Create initial kmer database and FASTA file with candidate Alu sequences.
# This script generates two files: Alu.pre-blast.fas and Alu.pre-blast.kmer.db
echo ""
echo "Extracting candidates from the reference genome"
./find_ref_plus_candidates.pl $human_chr_path
echo -n "Number of REF-plus candidates with TSD sequence: "
cat Alu.pre-blast.kmer.db | wc -l

# Run BLAST detected candidates and keep only the candidates that have >100 bits of homology with known Alu elements
blastall -i Alu.pre-blast.fas -d Alu.fas -p blastn -G 1 -E 1 -F F -b 1 -v 1 -m 8 | awk '{if($12>100)print}' | cut -f 1 | sort -n | uniq > Alu.REF-plus.blast
./extract_kmers_by_name.pl Alu.REF-plus.blast Alu.pre-blast.kmer.db > Alu.after-blast.kmer.db
./extract_fasta_by_name.pl Alu.REF-plus.blast Alu.pre-blast.fas > Alu.REF-plus.fas
echo -n "Number of candidates after BLAST homology search: "
cat Alu.after-blast.kmer.db | wc -l

# Check whether the 25+7 bp sequence seen in human genome is also present in chimp genome, allowing 2 mismatches. Remove those that are present in both species.
cut -f 3 Alu.after-blast.kmer.db > Alu.32-mer.txt
./glistquery $chimp_list -mm 2 -f Alu.32-mer.txt | awk '{if ($2<1)print}' > Alu.32-mer.not_in_chimp_genome.txt
./extract_kmers_by_32mer.pl Alu.32-mer.not_in_chimp_genome.txt Alu.after-blast.kmer.db > Alu.not_in_chimp_genome.kmer.db
cut -f 1 Alu.not_in_chimp_genome.kmer.db > Alu.not_in_chimp_genome.names.txt
./extract_fasta_by_name.pl Alu.not_in_chimp_genome.names.txt Alu.pre-blast.fas > Alu.REF-plus_not_in_chimp_genome.fas
echo -n "Number of candidates after checking the chimp genome: "
cat Alu.not_in_chimp_genome.kmer.db | wc -l

# Remove candidates that have identical k-mers. gmer_counter can't handle duplicate k-mers
cat Alu.not_in_chimp_genome.kmer.db | ./sort_kmer_db.sh | ./remove_all_duplicate_kmers.pl > REF-plus.kmer.db
echo -n "Number of candidates after removing duplicate k-mers: "
cat REF-plus.kmer.db | wc -l

# Clean up
#rm -f Alu.REF-plus.blast
rm -f Alu.32-mer.txt
rm -f Alu.not_in_chimp_genome.txt
#rm -f Alu.not_in_chimp_genome.kmer.db
#rm -f Alu.pre-blast.kmer.db
#rm -f Alu.pre-blast.fas
echo "REF-plus database created"
echo ""
