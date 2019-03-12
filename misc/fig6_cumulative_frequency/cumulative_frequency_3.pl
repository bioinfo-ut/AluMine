#!/bin/perl

$workingdir="/mambakodu/mremm/AluMine/misc/fig6_cumulative_frequency/";
$datadir="/storage7/analyysid/alu_insetrion_minus_181008/tester/";
#$datadir="/mambakodu/mremm/puhastatud_REF-minus/";
#$namefile="/mambakodu/mremm/nimed_2221.txt";
$namefile="/storage7/analyysid/alu_insetrion_minus_181008/tester/nimed_kmer_db.txt";
$resultfile="fig6.incremental.txt";
$final_db="/storage10/tarmo/Alu_artikli_materjalid/ALU2.kmer.db";
%seen = ();

open OUT, ">$resultfile" or die;
open KUMU, $namefile or die;
open KMER, $final_db or die;
while (<KMER>){
  @tmp=split(/\t/);
  $exists{$tmp[0]} = 1;
}
close KMER;

$kum_seen = 0;
while(<KUMU>){
   chomp;
   $id = $_;
   #$file = $datadir . $id . ".kmer.db";
   $file = $datadir . $id;
   $i = 0;
   open IND, $file or die;
   while(<IND>){
      chomp;
      $i++;
      @tmp = split(/\t/);
      @tmp_name = split(/:/,$tmp[0]);
      if ($tmp_name[2] eq "fw" and $tmp_name[3] eq "REF-minus"){
         $tmp_name[1] =~ s/^0+//; # convert to number
         $pos = $tmp_name[1] + 25;
         $tmp_name[1] = sprintf("%09d", $pos);
      }
      $tmp[0] = join (":", @tmp_name);
      if ($exists{$tmp[0]}){
        $seen{$tmp[0]} = 1;
      }
   }
   close IND;
   $cumulatively_seen = 0;
   for $k (keys %seen){
     $cumulatively_seen++ if $seen{$k};
   }
   print OUT "$id\t$i\t$cumulatively_seen\n";
}
close KUMU;

for $k (keys %exists){
   print STDERR "$k\n" if($k =~ /minus/ and not $seen{$k});
}
