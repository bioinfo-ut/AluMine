@cov = (5, 10, 15, 20, 25, 30, 35, 40, 45, 50);

foreach $i (@cov){
   $outfile = "NA12878_" . $i . "x.pre-gtester.txt";
   print STDERR "PRE-TESTER COVERAGE=$i\n";
   system ("cat /storage9/db/Illumina/Platinum/ERP001960/fastq/ERR194146_*.fastq | /mambakodu/mremm/depth_of_coverage_dilutions.pl $i | perl tarmo_ref_minus_fastq.pl | sort  | uniq > $outfile");
}

foreach $i (@cov){
   $infile = "NA12878_" . $i . "x.pre-gtester.txt";
   $outfile = "NA12878_" . $i . "x.post-gtester.txt";
   print STDERR "TESTER COVERAGE=$i\n";
   system ("/storage5/lauris/mapper/gtester -i /bigea/test_directory/maido/human_37_25.index  -g /storage7/ctg/uued/tester/human_37.names -3p 10 -f $infile > $outfile");
}

foreach $i (@cov){
   $infile = "NA12878_" . $i . "x.post-gtester.txt";
   $outfile = "NA12878_" . $i . "x.REF-minus.kmer.db";
   print STDERR "POST-TESTER COVERAGE=$i\n";
   system ("cat $infile | perl /storage10/tarmo/Alu_artikli_materjalid/tarmo_ref_minus_post_gtester.pl > $outfile");
}
