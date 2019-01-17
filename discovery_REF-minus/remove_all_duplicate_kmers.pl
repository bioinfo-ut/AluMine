#!/bin/perl

# Removes all copies of kmers that have duplicate in k-mer database.
while (<>){
   $line[$.] = $_;
   next if /^#/;
   chomp;
   @tmp = split(/\t/);
   $id = $tmp[0];
   $n = $tmp[1];
   for $i (1..$n){
      $fw_mer = $tmp[1+$i];
      $rv_mer = &rc(\$fw_mer);
      $kmer = (($fw_mer lt $rv_mer)?$fw_mer:$rv_mer);
      $count{$kmer} += 1;
   }
}

for $i (1..$.){
   $duplicate = 0;
   $_ = $line[$i];
   chomp;
   @tmp = split(/\t/);
   $id = $tmp[0];
   $n = $tmp[1];
   for $i (1..$n){
      $fw_mer = $tmp[1+$i];
      $rv_mer = &rc(\$fw_mer);
      $kmer = (($fw_mer lt $rv_mer)?$fw_mer:$rv_mer);
      if ($count{$kmer} > 1){
         $duplicate = 1;
         $duplicate_id{$kmer} = $id;
      }
   }
   print $_ . "\n" unless $duplicate;
}


sub rc {
  my $seqref = shift;
  (my $revcomp = reverse $$seqref) =~ tr [AaCcGgTt] [TtGgCcAa];
  return $revcomp;
}
