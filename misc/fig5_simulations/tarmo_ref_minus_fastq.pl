#!/bin/perl

# cat FASTQ_FILE.fastq | perl CURRENT_SCRIPT.pl > POTENTIAL_REF-minus_Alus.txt
@alus = (
'GGCCGGGCGC',
'GGCCGGGCGA','GGCCGGGCGG','GGCCGGGCGT',
'GGCCGGGCAC','GGCCGGGCCC','GGCCGGGCTC',
'GGCCGGGAGC','GGCCGGGGGC','GGCCGGGTGC',
'GGCCGGACGC','GGCCGGCCGC','GGCCGGTCGC',
'GGCCGAGCGC','GGCCGCGCGC','GGCCGTGCGC',
'GGCCAGGCGC','GGCCCGGCGC','GGCCTGGCGC',
'GGCAGGGCGC','GGCGGGGCGC','GGCTGGGCGC',
'GGACGGGCGC','GGGCGGGCGC','GGTCGGGCGC',
'GACCGGGCGC','GCCCGGGCGC','GTCCGGGCGC',
'AGCCGGGCGC','CGCCGGGCGC','TGCCGGGCGC'
);

while(<>){
   next if (($. % 4) != 2);
   next unless (/GGCCG/ or /GGCGC/ or /CGGCC/ or /GCGCC/); # Two common 5-mers and their reverse complements
   chomp;
   $fw_read = $_;
   $rv_read = &rev_compl($fw_read); 
   foreach $alu (@alus){
      $pos = index($fw_read,$alu,25);
      if ($pos >= 25){
         $seq = substr($fw_read,$pos-25,35);
         $sum{$seq}++;
      }
      $pos = index($rv_read,$alu,25);
      if ($pos >= 25){
       	 $seq = substr($rv_read,$pos-25,35);
         $sum{$seq}++;
      }
   }
}

foreach $kmer (keys %sum){
   print "$kmer\n" if ($kmer !~ /N/ and $sum{$kmer} > 4 and $sum{$kmer} < 101);
}

sub rev_compl {
   my $seqref = shift;
   my $revcomp = reverse $seqref;
   $revcomp =~ tr [AaCcGgTt] [TtGgCcAa];
   return $revcomp;
}
