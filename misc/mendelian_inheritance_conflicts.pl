#!/bin/perl

@cephs = ('NA12877_S1.REF-minus.kmer.db','NA12878_S1.REF-minus.kmer.db','NA12879_S1.REF-minus.kmer.db','NA12880_S1.REF-minus.kmer.db','NA12881_S1.REF-minus.kmer.db','NA12882_S1.REF-minus.kmer.db','NA12883_S1.REF-minus.kmer.db','NA12884_S1.REF-minus.kmer.db','NA12885_S1.REF-minus.kmer.db','NA12886_S1.REF-minus.kmer.db','NA12887_S1.REF-minus.kmer.db','NA12888_S1.REF-minus.kmer.db','NA12889_S1.REF-minus.kmer.db','NA12890_S1.REF-minus.kmer.db','NA12891_S1.REF-minus.kmer.db','NA12892_S1.REF-minus.kmer.db','NA12893_S1.REF-minus.kmer.db');
@kombin = ([13,14,1],[15,16,2]);
#@kombin = ([1,2,3],[1,2,4],[1,2,5],[1,2,6],[1,2,7],[1,2,8],[1,2,9],[1,2,10],[1,2,11],[1,2,12],[1,2,17],[13,14,1],[15,16,2]);
chdir ("/storage7/analyysid/ceph_alu/minus2");
$nr = 0;

open F, "REF-minus.kmer.db" or die;
%passed_marker = ();
while(<F>){
  chomp;
  @tmp = split(/:/);
  $id = $tmp[0] . ":" . $tmp[1]; 
  $passed_marker{$id} = 1;
}

# Initialize table
%tabel = ();
for $i (keys(%passed_marker)){
   for $j (0..scalar(@cephs)){
      $tabel{$i}[$j] = 0;
   }
}

# Fill table with discoveries
foreach $ceph (@cephs){
   open SISSE, "$ceph" or die;
   while(<SISSE>){
      chomp;
      @tmp = split(/\t/);
      @tmp2 = split(/:/,$tmp[0]);
      $id = $tmp2[0] . ":" . $tmp2[1];
      next unless $passed_marker{$id};
      $tabel{$id}[0] = $tmp[0];
      $tabel{$id}[$nr+1] = 1;
   }
   close SISSE;
   $nr++;
}

$trios = 0;
$trios_with_refminus_in_child = 0;
$conflicts = 0;
$alus = 0;
$conflicting_alus = 0;
$alus_with_ref_minus = 0;

foreach $alu (keys %tabel){
   ++$alus;
   $conflicting_alu = 0;
   $alu_in_child = 0;
   print "$tabel{$alu}[0]";
   for($r = 0; $r < scalar(@kombin); $r++){ # kombin
       $pattern = $tabel{$alu}[$kombin[$r][0]] . $tabel{$alu}[$kombin[$r][1]] . $tabel{$alu}[$kombin[$r][2]];
       $alu_in_child = 1 if (substr($pattern,2,1) eq "1");
       $conflicting_alu = 1 if ($pattern eq "001");
       ++$trios;
       ++$trios_with_refminus_in_child if $alu_in_child;
       ++$conflicts if $conflicting_alu;
       print "\t$pattern";
   }
   ++$alus_with_ref_minus if $alu_in_child;
   ++$conflicting_alus if $conflicting_alu;    
   print "\n";
}
print "Analyzed $alus Alu elements and $trios trios.\n";
printf ("REF-minus element was discovered in child in $trios_with_refminus_in_child trios, $conflicts with conflicts (%.2f%).\n", 100* $conflicts / $trios_with_refminus_in_child);
printf ("REF-minus element was discovered in child in $alus_with_ref_minus elements, $conflicting_alus with conflicts (%.2f%).\n", 100*$conflicting_alus/$alus_with_ref_minus); 
