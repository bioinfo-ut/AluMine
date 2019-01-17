#!/bin/perl -w

$diff = 1000000;
$prev_diff = 1000000;
$pos = 0;
$prev_pos = 0;
$prev_chr = "";
$prev_line = "";
$prev_gc1 = 16; 
$prev_gc2 = 16;


while (<>){
   $line = $_;
   chomp;
   @tmp = split(/:/);
   @tab = split(/\t/);
   $chr = $tmp[0];
   $pos = $tmp[1];
   $pos =~ s/^0+//; # convert to number

   if ($chr eq $prev_chr){
      $diff = $pos - $prev_pos;
   }
   else{
      $diff = 1000000;
   }
   $gc1 = ($tab[2] =~ tr/[C|G]/N/); # GC count in kmer_A
   $gc2 = ($tab[3] =~ tr/[C|G]/N/); # kmer_B

   if ($diff > 25 and $prev_diff > 25 and $prev_gc1 < 30 and $prev_gc2 < 30 and $prev_gc1 > 2 and $prev_gc2 > 2){
      print $prev_line;
   }
   $prev_gc1 = $gc1;
   $prev_gc2 = $gc2;
   $prev_line = $line;
   $prev_chr = $chr;
   $prev_pos = $pos;
   $prev_diff = $diff;
}
 
