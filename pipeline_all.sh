#!/bin/bash

# Warning!
# Before you run this script:
# - edit pipeline_ref_plus.sh and set proper paths (human and chimp genomes)
# - edit pipeline_ref_minus.sh and set proper paths (human reference) and sample names (FASTQ or BAM)


# REF-plus discovery
cd discovery_REF-plus
pipeline_ref_plus.sh
cd ..

# REF-minus discovery
cd discovery_REF-minus
pipeline_ref_minus.sh
cd ..

# Merging databases
cd discovery_merging_and_filtering
pipeline_merging_and_filtering.sh
cd ..

if [ -s discovery_merging_and_filtering/ALL.db ];then
  echo "Database created successfully"
else
  echo "Database was not created. Check for errors"
fi
