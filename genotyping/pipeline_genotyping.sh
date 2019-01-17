#!/bin/bash

# Download pre-made database of k-mers
if [ ! -s ALU_v1.kmer.db ];then
  wget http://bioinfo.ut.ee/FastGT/downloads/ALU_v1.kmer.db
fi

# Alternatively, use the ALU.db created during the discovery phase
# cp ../discovery_merging_and_filtering/ALU.db .

# Download some example FASTQ files with sequencing reads of the individual NA12877 (coded as ERR194146) and unpack.
# You can use other FASTQ files if you already have them on your server.
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR194/ERR194146/ERR194146.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR194/ERR194146/ERR194146_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR194/ERR194146/ERR194146_2.fastq.gz
gunzip ERR194146*.fastq

# Alternatively create FASTQ files from BAM
# samtools fastq ERR194146.bam > ERR194146.fastq

if [ ! -x gmer_counter ];then
  cp ../bin/gmer_counter .
  chmod 755 gmer_counter
fi

if [ ! -x gmer_caller ];then
  cp ../bin/gmer_caller .
  chmod 755 gmer_caller
fi

# Genotyping
gmer_counter -db ALU_v1.kmer.db ERR194146*.fastq > ERR194146.counts
gmer_caller --info ERR194146.counts > ERR194146.calls
