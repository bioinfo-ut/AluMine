foreach $SEX ("men", "women"){
$ref_genome = "/storage9/db/human_37/data/chr/human_37_" . $SEX . ".fa";
$alt_genome = "human_37_alt_" . $SEX . ".fa";
$S = "M" if ($SEX eq "men");
$S = "F" if ($SEX eq "women");

$wgsim = "/usr/local/bin/wgsim";
$ref_minus_script = "/storage10/tarmo/Alu_artikli_materjalid/tarmo_ref_minus_fastq.pl";
$nr_lines = `wc -l $ref_genome`;
$nr_insertions = 1000;

$alu_seq = "GGTGGCTCACGCCTGTAATCCCAGCACTTTGGGAGGCCGAGGCGGGTGGATCATGAGGTCAGGAGATCGAGACCATCCTGGCTAACAAGGTGAAACCCCGTCTCTACTAAAAATACAAAAAATTAGCCGGGCGCGGTGGCGGGCGCCTGTAGTCCCAGCTACTCGGGAGGCTGAGGCAGGAGAATGGCGTGAACCCGGGAAGCGGAGCTTGCAGTGAGCCGAGATTGCGCCACTGCAGTCCGCAGTCCGACCTGGGCGACAGAGCGAGACTCCGTCTC";

%sig_frequencies = (
'GGCCGGGCGC',	'457',
'GGCCGGGCGA',	'4',
'GGCCGGGCGG',	'2',
'GGCCGGGCGT',	'49',
'GGCCGGGCAC',	'96',
'GGCCGGGCCC',	'2',
'GGCCGGGCTC',	'5',
'GGCCGGGAGC',	'8',
'GGCCGGGGGC',	'3',
'GGCCGGGTGC',	'78',
'GGCCGGACGC',	'5',
'GGCCGGCCGC',	'1',
'GGCCGGTCGC',	'1',
'GGCCGAGCGC',	'4',
'GGCCGCGCGC',	'1',
'GGCCGTGCGC',	'1',
'GGCCAGGCGC',	'107',
'GGCCCGGCGC',	'4',
'GGCCTGGCGC',	'5',
'GGCAGGGCGC',	'6',
'GGCGGGGCGC',	'4',
'GGCTGGGCGC',	'100',
'GGACGGGCGC',	'2',
'GGGCGGGCGC',	'1',
'GGTCGGGCGC',	'5',
'GACCGGGCGC',	'5',
'GCCCGGGCGC',	'2',
'GTCCGGGCGC',	'2',
'AGCCGGGCGC',	'23',
'CGCCGGGCGC',	'6',
'TGCCGGGCGC',	'11');

# Generate a collection of Alu signatures according to frequencies given in sig_frequencies
@sig = ();
foreach $s (keys %sig_frequencies){
   $n = $sig_frequencies{$s};
   for (1..$n){
      push (@sig,$s);
   }
}

$i = 0;
foreach $i (0..10){ # $different error rates
   $output_file = "detected_" . $i . "_" . $S . ".txt";
   $mut_file = "mutations_" . $i . "_" . $S . ".txt";
   $stderr_file = "wgsim_seed_" . $i . "_" . $S . ".txt";
   $log_file = "inserted_alu_elements_" . $i . "_" . $S . ".txt";

   #### Generate random line numbers where Alu elements should be inserted
   @alu_locations = @sorted_alu_locations = ();
   for ($j=0;$j<1000;$j++){
       push @alu_locations, int(rand()*$nr_lines);
   }
   @sorted_alu_locations = sort { $a <=> $b } @alu_locations;

   #### Generate random signatures for inserted Alu elements
   @sorted_alu_signatures = ();
   for ($j=0;$j<1000;$j++){
      $random_int = int(rand()*1000);
      push @sorted_alu_signatures, $sig[$random_int];
   }

   #### Write a genome sequence with Alu insertions #####
   open REF, "$ref_genome" or die;
   open ALU,">$alt_genome" or die;
   open LOG, ">$log_file" or die;
   $k = 0;$chr="";$pos=0;
   while(<REF>){
      print ALU "$_";
      chomp;
      if (/\>/){
       	 @tmp = split(/\>|\s+/);
       	 $chr = $tmp[1];
       	 $pos = 0;
      }
      else{
         $pos += length($_);
      }

      if($. == $sorted_alu_locations[$k]){
         $signature = $sorted_alu_signatures[$k];
         $tsd_seq = substr($_,-15);
         $polyA = "A" x (int(rand()*75) + 5);
         if(rand() < 0.5){
           print ALU $signature . $alu_seq . $polyA . $tsd_seq . "\n";
           printf(LOG "%s:%09d:fw:REF-minus:%s\t%s\n", $chr, $pos, $signature, $tsd_seq);
         }
         else{
           print ALU revcom($signature . $alu_seq . $polyA) . $tsd_seq . "\n";
           $posrv = $pos - length($tsd_seq);
           printf(LOG "%s:%09d:rv:REF-minus:%s\t%s\n", $chr, $posrv, $signature, $tsd_seq);
         }
         $k++;
      }
   }
   close REF;
   close ALU;
   close LOG;
   print "Genome generation finished: ";
   print `date`;

   #### SIMULATIONS #####
   $fq1 = "reads_1" . $S . ".fq";
   $fq2 = "reads_2" . $S . ".fq";
   # 30x nucleotide coverage
   # default mutation_rate r=0 (default 0.0015)
   # fraction of indels R=0 (default 0.15)
   # fraction of ambiguous nucleotide cutoff A=1.00 (default 0.05) all-N sequences allowed
   $r = $i * 0.0001; # up to 0.010 or 1% per HAPLOID genome. This is 2% per diploid genome.
   system("$wgsim -h -e 0.005 -1 151 -2 151 -r $r -A 1.0 -d 500 -N 306000000 $alt_genome $fq1 $fq2 1>> $mut_file 2>> $stderr_file");
   print "Simulations finished: ";
   print `date`;

   #### Detect REF-minus elements #####
   system("cat $fq1 $fq2 | $ref_minus_script > $output_file");
   system("rm -f $fq1 $fq2 $alt_genome");
   print "Finished processing $S with error rate $r:";
   print `date`;
}
}
sub revcom{
   my $seqref = shift;
   my $rc = reverse $seqref;
   $rc =~ tr [AaCcGgTt] [TtGgCcAa];
   return $rc;
}
