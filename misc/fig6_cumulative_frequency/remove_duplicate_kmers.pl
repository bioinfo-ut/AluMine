#!/bin/perl

# Removes duplicates from the gmer_counter database
# One copy remains in the database
$number_of_duplicates = 0;

while (<>){
   $duplicate = 0;
   chomp;
   @tmp = split(/\t/);
   $id = $tmp[0];
   $n = $tmp[1];
   for $i (1..$n){
      $fw_mer = $tmp[1+$i];
      $rv_mer = &rc(\$fw_mer);
      $kmer = (($fw_mer lt $rv_mer)?$fw_mer:$rv_mer);
      # print STDERR "$fw_mer\t$rv_mer\t$kmer\n";
      $count{$kmer} += 1;
      if ($count{$kmer} > 1){
         $duplicate_id{$kmer} = $id;
         $duplicate = 1;
         $number_of_duplicates += 1;
      }
   }
   print $_ . "\n" unless $duplicate;
}

#for $k (keys %duplicate_id){
#      print STDERR "$count{$k}\t$id{$k}\t$k\n" if ($duplicate_id{$k});
#}

#print STDERR "$number_of_duplicates duplicates\n";



sub rc {
  my $seqref = shift;
  (my $revcomp = reverse $$seqref) =~ tr [AaCcGgTt] [TtGgCcAa];
  return $revcomp;
}
