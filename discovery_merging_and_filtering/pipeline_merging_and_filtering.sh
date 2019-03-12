#!/bin/bash

# Description: This is a collection of scripts for merging and filtering
# kmer database after the REF-plus and REF-minus discovery steps
# Authors: Tarmo Puurand and Maido Remm, 2018-2019

# Download some examples for testing:
if [ ! -f SNV32.kmer.db ];then
  echo "Downloading SNV database..."
  wget http://bioinfo.ut.ee/FastGT/downloads/SNV32.kmer.db
fi

# Merge REF-minus, REF-plus and 32-mer SNV database. SNV database helps to estimate correct model for calling genotypes
echo "Merging REF-minus, REF-plus and SNV databases..."
cat ../discovery_REF-minus/REF-minus.kmer.db ../discovery_REF-plus/REF-plus.kmer.db SNV32.kmer.db > tmp.kmer.db

# Remove all k-mers that are located within 25 bp from each other
# Remove gmer_counter needs unique kmers only
echo "Filtering databases..."
cat tmp.kmer.db | remove_duplicate_kmers.pl > ALL.kmer.db
rm -f tmp.kmer.db
echo -n "Finished filtering databases. "
echo -n "$(cat ALL.kmer.db | wc -l)"
echo " kmer pairs are in the file ALL.kmer.db"

########## Use ALL.kmer.db to genotype some real samples #######################
# Convert ALL.kmer.db to binary database:
### gmer_counter -db ALL.kmer.db -w ALL.kmer.dbb (ca 5 minutes)

# Run gmer_counter on real individuals (ca 1 hour per individual)
### for ind in all_individuals; do
### gmer_counter -dbb ALL.kmer.dbb ind.fastq > ind.counts
### done

# Run gmer_caller on real individuals (ca 5 minutes per individual)
### for ind in all_individuals; do
### gmer_caller ind.counts > ind.calls
### done

# Script filter_gt.pl goes through all .calls files and
### tags markers that have potential problems (ca 15 seconds per individual)
### perl filter_gt.pl

# Prints final ALU k-mer database with only markers that behave properly. Output goes to ALU_v1.kmer.db
### perl extract_kmers_from_db.pl

# Prints final stats: number of detected AB or BB genotypes for REF-minus and REF-plus elements separately
### perl filter_gt_final_stats.pl

# ==============================================================================
#Total number of files 2241 filter_gt_files.txt
#Total number of candidates: 28962 filter_gt_markers.txt
#Remaining markers: 10436
#Remaining REF-plus: 904 + 10845 = 11749
#Remaining REF-minus: 7439
#Remaining total: 19188
#
#Removed for EGT: 7409
#Removed for HWE: 53
#Removed for MFB: 48
#Removed for ZRM: 2264
#Not Removed, ZRP: 10845
#Not Removed OK: 8343
