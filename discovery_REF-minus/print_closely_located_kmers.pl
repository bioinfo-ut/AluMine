#!/bin/perl -w

$diff = 1000000;
$prev_diff = 1000000;
$pos = 0;
$prev_pos = 0;
$prev_chr = "";
$prev_line = "";

while (<>){
   $line = $_;
   @tmp = split(/:/);
   $chr = $tmp[0];
   $pos = $tmp[1];
   $pos =~ s/^0+//; # convert to number
   if ($chr eq $prev_chr){
      $diff = $pos - $prev_pos;
   }
   else {
      $diff = 1000000;
   }

   if ($diff < 25 or $prev_diff < 25){
      print $prev_line;
   }
   $prev_line = $line;
   $prev_chr = $chr;
   $prev_pos = $pos;
   $prev_diff = $diff;
}
 
