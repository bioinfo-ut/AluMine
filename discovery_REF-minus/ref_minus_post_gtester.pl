#!/bin/perl
# Reads gtester output and discards non-locatable, non-unique and fixed Alu candidates

@tmp = ();
$i1 = $i2 = $i3 = $i4 = 0;
while(<>){
   chomp;
   if (/^Query/) { # New query
      if ($. > 1){ # process the previous query
         if ($locations_in_genome == 0){
            $i1++; # 25-mer not detected in the genome
         }
         elsif ($locations_in_genome > 1){
            $i2++; # 25-mer not unique in the genome
         }
         elsif ($locations_in_genome == 1){ # 25-mer unique in the genome
            if($alu_signature eq $ref_sequence){ # 100% identical strings
               $i3++; # 25-mer unique in the genome, but Alu is fixed
            }
            else{ # Calculate Levenshtein distance
               $dist = &levenshtein($alu_signature,$ref_sequence);
               if (($dist <= 2) or ($kmer_A eq $kmer_B)){
                  $i3++; # 25-mer unique in the genome, but Alu is fixed
               }
               else { # Edit distance between two 10 bp strings > 2
                  $i4++;
                  if ($count{$kmer_A} != 1 and $count{$kmer_B} != 1 and $count{&rc(\$kmer_A)} != 1 and $count{&rc(\$kmer_B)} != 1){
                     $i5++; # 25-mer unique in the genome and Alu is NOT fixed
                     printf("%s:%09d:%s:REF-minus:%s\t2\t%s\t%s\n", $chr, $pos, $strand, $alu_signature, $kmer_A, $kmer_B);
                     $count{$kmer_A} = 1;
                     $count{$kmer_B} = 1;
                  }
               }
            }
         }
         else {
            print STDERR "Unexpected number of locations: $locations_in_genome\n";
         }
      }
      # reset variables
      $locations_in_genome = 0;
      $chr = "";
      $pos = 0;
      $strand = "";
      $alu_signature = $ref_sequence = "";
      $kmer_A = $kmer_B = "";

      @tmp = split(/\s+/);
      $alu_signature = substr($tmp[2],25,10);
      $kmer_B = substr($tmp[2],0,32);
   }
   elsif($locations_in_genome == 0){
      @tmp = split(/\s+/);
      $chr = $tmp[1];
      if ($tmp[2] eq '+'){
          $strand = "fw";
          $pos = $tmp[3] + 25; # location of Alu insertion is at the END of 25-mer
      }
      if ($tmp[2] eq '-'){
        $strand = "rv";
        $pos = $tmp[3];
      }

      $locations_in_genome += 1;
      $_ = <>; # next line contains ref sequence 25+10
      $ref_sequence = substr($_,25,10);
      $kmer_A = substr($_,0,32);
   }
   else{
      $locations_in_genome += 1;
      $_ = <>;
   }
}

# Process the last query:
if ($locations_in_genome == 0){
   $i1++; # 25-mer not detected in the genome
}
elsif ($locations_in_genome > 1){
   $i2++; # 25-mer not unique in the genome
}
elsif ($locations_in_genome == 1){ # 25-mer unique in the genome
   if($alu_signature eq $ref_sequence){ # 100% identical strings
      $i3++; # 25-mer unique in the genome, but Alu is fixed
   }
   else{ # Calculate Levenshtein distance
      $dist = &levenshtein($alu_signature,$ref_sequence);
      if ($dist <= 2  or ($kmer_A eq $kmer_B)){
         $i3++; # 25-mer unique in the genome, but Alu is fixed
      }
      else { # Edit distance between two 10 bp strings > 2
         $i4++;
         if ($count{$kmer_A} != 1 and $count{$kmer_B} != 1 and $count{&rc(\$kmer_A)} != 1 and $count{&rc(\$kmer_B)}){
            $i5++; # 25-mer unique in the genome and Alu is NOT fixed
            printf("%s:%09d:%s:REF-minus:%s\t2\t%s\t%s\n", $chr, $pos, $strand, $alu_signature, $kmer_A, $kmer_B);
         }
      }
   }
}
else {
   print STDERR "Unexpected number of locations: $locations_in_genome\n";
}

printf (STDERR "%d 25-mers tested\n", $i1 + $i2 + $i3 + $i4);
printf (STDERR "%d 25-mers with location\n", $i2 + $i3 + $i4);
printf (STDERR "%d 25-mers with unique location\n", $i3 + $i4);
printf (STDERR "%d 25-mers with unique location and Alu is not fixed\n", $i4);
printf (STDERR "%d 25-mers with unique location and Alu is not fixed and non-redundant k-mer\n", $i5);

################## FUNCTIONS ##################################################################################################
# https://en.wikibooks.org/wiki/Algorithm_Implementation/Strings/Levenshtein_distance
use List::Util qw(min);
sub levenshtein {
    my ($str1, $str2) = @_;
    my @ar1 = split //, $str1;
    my @ar2 = split //, $str2;

    my @dist;
    $dist[$_][0] = $_ foreach (0 .. @ar1);
    $dist[0][$_] = $_ foreach (0 .. @ar2);

    foreach my $i (1 .. @ar1){
        foreach my $j (1 .. @ar2){
            my $cost = $ar1[$i - 1] eq $ar2[$j - 1] ? 0 : 1;
            $dist[$i][$j] = min(
                        $dist[$i - 1][$j] + 1,
                        $dist[$i][$j - 1] + 1,
                        $dist[$i - 1][$j - 1] + $cost );
        }
    }
    return $dist[@ar1][@ar2];
}

sub rc {
  my $seqref = shift;
  (my $revcomp = reverse $$seqref) =~ tr [AaCcGgTt] [TtGgCcAa];
  return $revcomp;
}
