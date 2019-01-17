#! /usr/bin/perl
# Finds sequences from fasta file all.dat if given names in reference file
# file name from command-line!

$reffile= "$ARGV[0]";
$infile= "$ARGV[1]"; # "all.dat";
$i = 0;

open (REF, "$reffile") or die "$reffile not found
Usage: $0 <Reference file with file names> <multifasta with sequences>\n";
while (<REF>){
    next if (/^\s+$/);
    s/\>//;
    s/\*//g;
    chomp;
    $seen{$_} = 1;
    ++$i;
}
close REF;
#print STDERR "$i names \n";

# Now go through the fasta file
$name = "unknown";
open IN, $infile or die "$infile not found";
while(<IN>){
   chomp;
   if(/\>/){
      chomp;
      s/\>//;
      $name = $_;
#      print STDERR $name . "\n";
      print ">$name\n" if ($seen{$name});
   }
   else {
      print $_ . "\n" if ($seen{$name});
   }
}
close IN;
