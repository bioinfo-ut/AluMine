#! /bin/perl
# Finds lines from file, given names in reference file

$namefile = "$ARGV[0]";
$infile= "$ARGV[1]";

%names = ();
$i = $j = 0;

open (REF, "$namefile") or die "$namefile not found.
Usage:
$0 <File_with_names> <File to be filtered>";
while (<REF>){
    next if (/^\s+$/);
    chomp;
    @tmp=split(/\t/);
    $names{$tmp[0]} = 1;
    $i++;
}
close REF;
# print STDERR "$i names \n";

# Now go through the kmer file
open IN, $infile or die "$infile not found.
Usage:
$0 <File_with_names> <File to be filtered>";
while(<IN>){
  @tmp = split(/\t/);
  $kmer = $tmp[2];
  if($names{$kmer}){
    print;
    $j++;
  }
}
close IN;
# print STDERR "Done, $j lines detected\n";
