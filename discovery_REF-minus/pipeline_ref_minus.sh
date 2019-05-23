#! /bin/bash

# Description: This is a pipeline for discovery of REF-minus polymorphic ALU
# element insertions from BAM or FASTQ files.
# Authors: Tarmo Puurand and Maido Remm, 2018-2019

# Dependencies:
# samtools, if you want to take sequences from BAM files, not from FASTQ
# List/Util.pm (https://perldoc.perl.org/List/Util.html) used in ref_minus_post_gtester.pl

############ Download some examples for testing: ###############################
# if [ ! -s NA12891_S1.bam ];then
#   wget http://bioinfo.ut.ee/FastGT/downloads/NA12891_S1.bam
# fi
# if [ ! -s NA12892_S1.bam ];then
#   wget http://bioinfo.ut.ee/FastGT/downloads/NA12892_S1.bam
# fi
# sample_path="."
# samples=(NA12891_S1 NA12892_S1)

############# Alternatively define path and ID-s of your own samples: ##########
# These must be matched with BAM or FASTQ file names
sample_path="/storage9/db/Illumina/Platinum/ERP001960/bam"
samples=(NA12877_S1)
# sample_path="/ctg"
# samples=(V19411-clone V26498-clone)
########## Define human reference genome path for creating gtester indices #####
human_chr_path="/storage9/db/human_37/data/chr" # Each chr in separate file, cannot use multifasta files
gtester_path="/storage9/db/kmer_indexes"
gtester_prefix="${gtester_path}/human_37"
gtester_index="${gtester_prefix}_25.index"
gtester_names="${gtester_prefix}.names"
############# Done with path settings #############################################


#################### Generating gtester index file ###############################
if [ ! -s $gtester_index ] || [ ! -s $gtester_names ];then
  if [ ! -x gindexer ];then
    cp ../bin/gindexer .
    chmod 755 gindexer
  fi
  echo "gtester 25-mer index does not exist. Generating... This will take an hour or more"
  ./gindexer -n 25 -o ${gtester_prefix} -i \
"${human_chr_path}/1.fa" \
"${human_chr_path}/2.fa" \
"${human_chr_path}/3.fa" \
"${human_chr_path}/4.fa" \
"${human_chr_path}/5.fa" \
"${human_chr_path}/6.fa" \
"${human_chr_path}/7.fa" \
"${human_chr_path}/8.fa" \
"${human_chr_path}/9.fa" \
"${human_chr_path}/10.fa" \
"${human_chr_path}/11.fa" \
"${human_chr_path}/12.fa" \
"${human_chr_path}/13.fa" \
"${human_chr_path}/14.fa" \
"${human_chr_path}/15.fa" \
"${human_chr_path}/16.fa" \
"${human_chr_path}/17.fa" \
"${human_chr_path}/18.fa" \
"${human_chr_path}/19.fa" \
"${human_chr_path}/20.fa" \
"${human_chr_path}/21.fa" \
"${human_chr_path}/22.fa" \
"${human_chr_path}/X.fa" \
"${human_chr_path}/Y.fa" \
"${human_chr_path}/MT.fa"
fi

################## Discovery of potential REF-minus Alus ######################
echo "Starting discovery of REF-minus elements. This takes ca 2 hours per sample."
for id in ${samples[@]}
do
   echo "Starting REF-minus discovery for $id"
   if [ ! -s $id.gtester.input.txt ];then
     # Define path to BAM files
     infile="${sample_path}/${id}.bam"
     if [ -s $infile ];then
       samtools view $infile | perl find_ref_minus_candidates_bam.pl | sort  | uniq  >  $id.gtester.input.txt
       # If the infile(s) are in FASTQ format use another script:
       # cat $infile.fastq | perl find_ref_minus_candidates_fastq.pl | sort  | uniq  >  $id.gtester.input.txt
     else
       echo "$infile not found."
     fi
   fi
done
echo "Finished discovery of REF-minus elements."

################### Searching for locations with gtester ######################
# It is significantly faster to run all gtester jobs together
echo "Starting gtester runs. This takes 2-6 hours for the first sample and few minutes for each subsequent sample."
for id in ${samples[@]}
do
  # Copy gtester binary here if necessary
  if [ ! -x gtester ];then
    cp ../bin/gtester .
    chmod 775 gtester
  fi

  # Find locations and generate 32-mer database
  if [ -s $id.gtester.input.txt ] && [ ! -s $id.gtester.output.txt ];then
    # Run gtester for each sample
    ./gtester -i $gtester_index -g $gtester_names -3p 10 -f $id.gtester.input.txt > $id.gtester.output.txt
    # post-process the output, generating kmer database
    cat $id.gtester.output.txt | perl ref_minus_post_gtester.pl > $id.REF-minus.kmer.db
    # Clean up
    rm -f $id.gtester.input.txt
    rm -f $id.gtester.output.txt
    echo "gtester done for $id"
  fi

done
echo "Finished gtester runs."

############ Merging REF-minuses from all samples into one database ############
echo "Merging REF-minus elements from all samples..."
# cd /storage7/analyysid/alu_insetrion_minus_181008/tester/
# for i in V*db
# do
#   echo "$i "
#   cat $i | ~/sort_and_filter_kmer_db.sh > ~/filtered_REF_minus/${i}
#   cat $i | ~/sort_and_filter_kmer_db.sh | ~/count_GC_kmer_db.sh > ~/filtered_REF_minus2/${i}
# done

rm -f tmp.kmer.db
touch tmp.kmer.db
for id in ${samples[@]}
do
  cat $id.REF-minus.kmer.db >> tmp.kmer.db
done

################ Filtering the final database #######################
# Remove closely located candidates (within 25 bp from each other)
# and those with GC% >= 30/32 or GC% <= 2/32
# and those that have identical k-mer in the database
# (gmer_counter is confused by k-mers with identical sequences)
cat tmp.kmer.db | ./sort_kmer_db.sh | uniq | ./remove_closely_located_and_GC_rich_kmers.pl | ./remove_all_duplicate_kmers.pl > REF-minus.kmer.db
rm -f tmp.kmer.db
echo "Finished merging REF-minus elements. The results are in the file REF-minus.kmer.db"
