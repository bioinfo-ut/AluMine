#!/bin/perl

%tag = ();
open F, "filter_gt_markers.txt" or die;
while (<F>){
   @tmp = split(/\t/);
   $tag{$tmp[1]} = "OK" if ($tmp[0] eq "OK " or $tmp[0] eq "ZRP");
}

open DB, "ALL.db" or die;
while (<DB>){
   @tmp = split(/\t/);
   print if ($tag{$tmp[0]} eq "OK" or $tmp[0] !~ /REF/);
}
close DB;
