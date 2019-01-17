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
cat tmp.kmer.db | remove_duplicate_kmers.pl > ALL.db
rm -f tmp.kmer.db
echo -n "Finished filtering databases. "
echo -n "$(cat ALL.db | wc -l)"
echo " kmer pairs are in the file ALL.db"

########## Use ALL.db to genotype some real samples ############################
# Convert ALL.db to binary database:
# gmer_counter -db ALL.db -w ALL.dbb (ca 5 minutes)

# Run gmer_counter on real individuals (ca 1 hour per individual)
# gmer_counter -dbb ALL.dbb sample.fastq > sample.counts

# Run gmer_caller on real individuals (ca 5 minutes per individual)
# gmer_caller sample.counts > sample.calls

# Script filter_gt.pl goes through all .calls files and
# tags markers that have potential problems (ca 15 seconds per individual)
# perl filter_gt.pl

# Prints final k-mer database with only markers that behave properly
# perl extract_kmers_from_db.pl > ALU_v1.kmer.db

#Total number of candidates: 32786 filter_gt_markers.txt
#Remaining markers: 10436
#Remaining REF-plus: 1825
#Remaining REF-minus: 8611
#
#Removed for EGT: 8985
#Removed for HWE: 140
#Removed for MFB: 46
#Removed for ZRM: 1808
#Removed for ZRP: 11371
