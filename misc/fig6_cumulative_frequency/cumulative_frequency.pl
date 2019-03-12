#!/bin/perl

open KUMU,"/storage7/analyysid/alu_insetrion_minus_181008/tester/nimed_kmer_db.txt" or die;
$kum= 0;
while(<KUMU>){
   chomp;
   open IND,"$_" or die;
   $tk = 0;
   $nimi = substr($_,0,6);
   while(<IND>){
      chomp;
      @tmp = split(/\t/);
      $kum++ unless exists $koos{$tmp[0]};
      $koos{$tmp[0]} = 1;
      $tk++;
   }
   close IND;
   print"$nimi\t$tk\t$kum\n";
}
close KUMU;
