#!/bin/bash

# Add output to the file sims.txt and then run figure_sims.pl to generate the table
for k in 1 2 3 4 5 6 7 8 9 10
do
	a="/smallea/tmp/maido/sim"
	dir=$a$k
	for j in M F
	do
		for i in 0 1 2 3 4 5 6 7 8 9 10
		do
			cd $dir
			echo "Mutation rate $i in $j" 
			echo -n "N-rich regions: "
			grep -c NNN inserted_alu_elements_${i}_${j}.txt
			grep -v NNN inserted_alu_elements_${i}_${j}.txt | cut -f 1 > in_${i}_${j}
			echo -n "True positives: "
			grep -c -f in_${i}_${j} gtester.${i}.${j}.REF-minus.kmer.db
			rm -f in_${i}_${j}
			echo -n "Total number detected: "
			wc -l gtester.${i}.${j}.REF-minus.kmer.db
		done
        done
done
