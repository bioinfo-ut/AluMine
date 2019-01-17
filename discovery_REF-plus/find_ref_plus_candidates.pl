#!/bin/perl

$chr_path = $ARGV[0]; #"/storage9/db/human_37/data/chr";
@chrs = ('1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22','X','Y');
#@chrs = ('Y');
@alus = (
'GGCCGGGCGC',
'AGCCGGGCGC','CGCCGGGCGC','TGCCGGGCGC',
'GACCGGGCGC','GCCCGGGCGC','GTCCGGGCGC',
'GGACGGGCGC','GGGCGGGCGC','GGTCGGGCGC',
'GGCAGGGCGC','GGCGGGGCGC','GGCTGGGCGC',
'GGCCAGGCGC','GGCCCGGCGC','GGCCTGGCGC',
'GGCCGAGCGC','GGCCGCGCGC','GGCCGTGCGC',
'GGCCGGACGC','GGCCGGCCGC','GGCCGGTCGC',
'GGCCGGGAGC','GGCCGGGGGC','GGCCGGGTGC',
'GGCCGGGCAC','GGCCGGGCCC','GGCCGGGCTC',
'GGCCGGGCGA','GGCCGGGCGG','GGCCGGGCGT'
);

$TSD_LENGTH = 5;
$MIN_LENGTH = 270;
$MAX_LENGTH = 350;

foreach $alu (@alus){
  $with_tsd{$alu}=0;
  $without_tsd{$alu}=0;
}

#$fastafilename = "Alu.TSD" . $TSD_LENGTH . "." . $MIN_LENGTH . "-" . $MAX_LENGTH . ".fas";
#$kmerfilename = "Alu.TSD" . $TSD_LENGTH . "." . $MIN_LENGTH . "-" . $MAX_LENGTH . ".pre-blast.kmer.db";
$fastafilename = "Alu.pre-blast.fas";
$kmerfilename = "Alu.pre-blast.kmer.db";

open (FAS, ">", $fastafilename) or die;
open (KMER, ">", $kmerfilename) or die;
# load chromosome sequence in
foreach $chr (@chrs) {
  $chr_seq = "";
  $chr_file = $chr_path . "/" . $chr . ".fa";
  open (CHR,"$chr_file") or die "Could not open the chr file $chr_file\n";
  while(<CHR>){
    chomp;
    next if /^\s*$/;
    next if /^\s*>/;
    $chr_seq .= $_;
  }
  close(CHR);
  foreach $alu (@alus){
    $pos = 0;
    while(($pos = index($chr_seq, $alu, $pos+1)) != -1) {
      $kmer_A = $kmer_B = $tsd = "";

      $kmer_A = substr($chr_seq, $pos-25,32);
      $kmer_B = substr($chr_seq, $pos-25,25);
      next if ($kmer_A =~ /N/ or $kmer_B =~ /N/); # N-s won't work in gmer_counter
      $tsd = substr($chr_seq, $pos-$TSD_LENGTH,$TSD_LENGTH);
      for($end_pos = $pos+$MIN_LENGTH; $end_pos < $pos+$MAX_LENGTH; $end_pos++){
        if($tsd eq substr($chr_seq, $end_pos,$TSD_LENGTH)){
          $alu_seq = substr($chr_seq,$pos,$end_pos-$pos+$TSD_LENGTH+7);
          $kmer_B .= substr($chr_seq,$end_pos+$TSD_LENGTH,7); # this situation might happen several times
          last if (length($kmer_B) == 32);
        }
      }
      $alu_len = $end_pos - $pos;
      if (length($kmer_B) == 32){
        printf (KMER "%s:%09d:fw:REF-plus:%s:%d\t2\t%s\t%s\n", $chr, $pos, $alu, $alu_len, $kmer_A, $kmer_B);
        printf (FAS  ">%s:%09d:fw:REF-plus:%s:%d\n%s\n", $chr, $pos, $alu, $alu_len, $alu_seq);
        $with_tsd{$alu} += 1;
      }
      else{
        $without_tsd{$alu} += 1;
      }
    }
  }
  ######### Reverse strand ###########
  $chr_rev = reverse_complement(\$chr_seq);
  foreach $alu (@alus){
    $pos = 0;
    while(($pos = index($chr_rev, $alu, $pos+1)) != -1) {
      $kmer_A = $kmer_B = $tsd = "";

      $kmer_A = substr($chr_rev, $pos-25,32);
      $kmer_B = substr($chr_rev, $pos-25,25);
      next if ($kmer_A =~ /N/ or $kmer_B =~ /N/);
      $tsd = substr($chr_rev, $pos-$TSD_LENGTH, $TSD_LENGTH);
      for($end_pos = $pos+$MIN_LENGTH; $end_pos < $pos+$MAX_LENGTH; $end_pos++){
        if($tsd eq substr($chr_rev, $end_pos, $TSD_LENGTH)){
          $alu_seq = substr($chr_rev,$pos,$end_pos-$pos+$TSD_LENGTH+7);
          $kmer_B .= substr($chr_rev,$end_pos+$TSD_LENGTH,7); # this situation might happen several times
          last if (length($kmer_B) == 32);
        }
      }
      $alu_len = $end_pos - $pos;
      if (length($kmer_B) == 32){
        $rev_pos = length($chr_rev) - $pos;
        printf (KMER "%s:%09d:rv:REF-plus:%s:%d\t2\t%s\t%s\n", $chr, $rev_pos, $alu, $alu_len, $kmer_A, $kmer_B);
        printf (FAS  ">%s:%09d:rv:REF-plus:%s:%d\n%s\n", $chr, $rev_pos, $alu, $alu_len, $alu_seq);
        $with_tsd{$alu} += 1;
      }
      else{
        $without_tsd{$alu} += 1;
      }
    }
  }
}
close $fastafile;
print STDERR "SIGNATURE\tALL\tWITH_TSD\n";
$total = 0;
foreach $alu (@alus){
  $all_alu = $with_tsd{$alu} + $without_tsd{$alu};
  $total += $all_alu;
  print STDERR "$alu\t$all_alu\t$with_tsd{$alu}\n";
}
print STDERR "Number of REF-plus candidates: $total\n";

sub reverse_complement {
  my $seqref = shift;
  (my $revcomp = reverse $$seqref) =~ tr [AaCcGgTt] [TtGgCcAa];
  return $revcomp;
}
