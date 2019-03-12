#!/bin/perl

%tag = ();
open OUT, ">/mambakodu/mremm/ALU_v1.kmer.db" or die;
open F, "/mambakodu/mremm/filter_gt_markers.txt" or die;
while (<F>){
   @tmp = split(/\t/);
   $tag{$tmp[1]} = "OK" if ($tmp[0] eq "OK " or $tmp[0] eq "ZRP");
}

open DB, "/mambakodu/mremm/ALL.kmer.db" or die;
while (<DB>){
   @tmp = split(/\t/);
   print OUT if ($tag{$tmp[0]} eq "OK" or $tmp[0] !~ /REF/);
}
close DB;
