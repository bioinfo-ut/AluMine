#!/bin/bash

for i in 0 1 2 3 4 5 6 7 8 9 10
do
   for sex in M F
   do
      /storage5/lauris/mapper/gtester -i /bigea/test_directory/maido/human_37_25.index -g /storage7/ctg/uued/tester/human_37.names -3p 10 -f detected_${i}_${sex}.txt | perl /storage10/tarmo/Alu_artikli_materjalid/tarmo_ref_minus_post_gtester.pl | sed -e 's/^X:/23:/' -e 's/^Y:/24:/' | sort -n | sed -e 's/^23:/X:/' -e 's/^24:/Y:/' > gtester.${i}.${sex}.REF-minus.kmer.db
      echo "Done with $i"
   done
done
