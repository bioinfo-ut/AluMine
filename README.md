# AluMine

This GitHub repository stores various scripts required to discover and genotype polymorphic Alu element insertions. There are four different workflows: REF-plus discovery, REF-minus discovery, merging and filtering workflow and genotyping workflow. REF-plus and REF-minus workflows generate text file containing 32-mer pairs that can be subsequently used for genotyping using FastGT package. The scripts are written in PERL and bash.

Downloading to local server:  
```
git clone https://github.com/bioinfo-ut/AluMine  
```
### REF-plus discovery scripts  
The key steps in REF-plus discovery pipeline are:
* Search for all potential full-length Alu elements with the script find_ref_plus_candidates.pl. This script searches the reference genome for 10bp Alu element signatures (with 1 mismatch) and for Target Site Duplication sequences within 270-350 bp.
* BLAST search that checks whether detected candidate elements are homologous to known Alu elements.
* Search against chimpanzee genome using 25-mer lists. This step removes older elements that are likely to be fixed in both species.

To run the scripts yourself open pipeline_ref_plus.sh in text editor and define paths to FASTA files of the reference genome and chimpanzee genome.
```
cd ~/AluMine/discovery_REF-plus
bash pipeline_ref_plus.sh  
```  

### REF-minus discovery scripts
The key elements in REF-minus discovery pipeline are: find_ref_minus_candidates_bam.pl, gtester, and ref_minus_post_gtester.pl.  

* Scripts find_ref_minus_candidates_bam.pl or find_ref_minus_candidates_fastq.pl search for 10bp Alu signature sequences from BAM or FASTQ files and
write out all potential signatures together with 25bp flanking sequence.
* Then 'gtester' is used to localize the 25bp flanking sequences in the
reference genome.
* ref_minus_post_gtester.pl removes fixed Alu elements (present in the human reference genome) and prints out pair of 32-mers
for those that are not present in reference genome.

Before running the REF-minus discovery pipeline, open file pipeline_ref_minus.sh in text editor and define paths to sample files and path to human chromosome files.  
```
cd ~/AluMine/discovery_REF-minus
bash pipeline_ref_minus.sh
```
pipeline_ref_minus.sh supports input in BAM, single FASTQ, multiple FASTQ and gnuzipped FASTQ formats. 

### Merging and filtering the k-mer databases
The following steps are required for additional filtering using genotype data from real individuals
* Alu-element candidates with identical k-mers were removed.
* Alu-element candidates that are located within 25bp of each other were removed.
* REF-minus and REF-plus k-mers need to be merged with 30M SNV k-mers. SNVs help to build more accurate model for genotype calling.

Remaining Alu-elements were genotyped on 2200 individuals and additional filtering steps were used:
* Alu elements with >10% calls with unexpected ploidy were removed.
* Alu elements with Hardy-Weinberg Equilibrium P-value<1.5E-6 were removed.

These steps are performed separately, their description can be seen in file filter_gt.pl.
```
cd ~/AluMine/discovery_merging_and_filtering
more pipeline_merging_and_filtering.sh
```

### Genotyping
It is possible to skip the discovery phase and use our database of known Alu insertion polymorphisms (32,786 candidate polymorphisms).  
The k-mer database for genotyping (ALU_v1.kmer.db) is available at [FastGT webpage](http://bioinfo.ut.ee/FastGT/index.php?r=site/page&view=kmers).

```
cd ~/AluMine/genotyping
bash pipeline_genotyping.sh
```
### Citing
Please cite: Puurand T, Kukuškina V, Pajuste F-D, Remm M. (2019). AluMine: alignment-free method for the discovery of polymorphic Alu element insertions. Mobile DNA 10:31.
[https://doi.org/10.1186/s13100-019-0174-3](doi: https://doi.org/10.1186/s13100-019-0174-3).

### Additional data
Additional material can be downloaded from our webpage at [http://bioinfo.ut.ee/?page_id=167&lang=en](http://bioinfo.ut.ee/?page_id=167&lang=en). Pre-compiled human 25-mer index for REF-minus discovery (57GB) and pre-compiled chimp 32-mer list for REF-plus discovery (27GB) can be downloaded from [http://bioinfo.ut.ee/AluMine/](http://bioinfo.ut.ee/AluMine/).
