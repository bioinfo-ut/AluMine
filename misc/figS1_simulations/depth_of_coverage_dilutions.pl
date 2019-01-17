#!/bin/perl

$requested_depth_of_coverage = $ARGV[0];
$reads = (3252722312 + 3252722312 + 34630608) / 4; # NA12878 DOC = 55.0456332026667
$read_length = 101;
$total_nucleotides = $reads * $read_length;
$depth_of_coverage = $total_nucleotides / 3000000000;
$resolution = 1000;
$dilution = int(0.5 + $resolution * $requested_depth_of_coverage / $depth_of_coverage);
$dilution = $resolution if ($dilution > $resolution);

#print STDERR "DOC = $depth_of_coverage\n";
#print STDERR "RDOC= $requested_depth_of_coverage\n";
#print STDERR "DIL = $dilution / $resolution\n";

$i = 0;
while(<STDIN>){
   if ($i % 4*$resolution == 0){
       for($j=0;$j<4*$resolution;$j++){
          print if ($j <= 4*$dilution);
          $_ = <STDIN>; $i++;
       }
   }
}

