# Number of nucleotides in the reference genome was determined as following: 
grep -v \> human_37_all.fa | tr -d '\n' | wc 

# Number of unique k-mers was determined as following: 
# Make list of all k-mers
glistmaker human_37_all.fa -w 25 -o human37

# Overall number of 25-mers without N: 
glistquery /storage9/db/kmer_lists/human37_25.list -stat

# Number of 25-mers with unique location: 
glistquery /storage9/db/kmer_lists/human37_25.list -distribution 1

