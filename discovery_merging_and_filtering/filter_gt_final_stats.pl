#!/bin/perl

# Filter markers based on their predicted genotypes in 2280 individuals
# remove individuals with low proportion of expected genotypes (expected are diploid genotypes for autosomes, haploid for X and Y in men )
# mark EGT low proportion of expected genotypes (expected are diploid genotypes for autosomes, haploid for X and Y in men )
# mark HWE low probability for Hardy-Weinberg Equilibrium
# mark MFB B-allele significantly different in men and women
# mark ZRM zero freq REF-minus
# mark ZRP zero freq REF-plus (OK to keep)
# mark the rest as OK

$aludbfile = "/mambakodu/mremm/ALU_v1.kmer.db";
$outfile2 = "/mambakodu/mremm/filter_gt_stats.txt";
chdir "/storage7/ctg/uued/calling_gmer/calls/ALU_SNP/";

$critical_proportion_of_expected_gt_for_individual = 0.50;
$critical_likelihood_of_the_genotype = 0.9; # don't use any genotypes that have lower likelihood than this
$critical_total_kmer_count = 5;             # don't use any genotypes that are based on k-mer frequency lower than 5 (median diploid depth of coverage in this set is 22.8)

opendir(my $dh, ".") || die "Can't opendir: $!";
my @files = readdir($dh);
closedir $dh;

%tag = ();
open(DB, '<', $aludbfile) or die;
while (<DB>){
    @tmp = split(/\t/);
    $tag{$tmp[0]} = "OK" if ($tmp[0] =~ /REF/);
}
close DB;

open(F2, '>', $outfile2) or die;
foreach $f (@files){
    next if (-s $f < 1000000000 or $f !~ /V538/); # V538* for testing
    $sex = ""; $gt = "";
    $expected_genotype_proportion = 0;
    $expected_genotype = $unexpected_genotype = 0;
    $alt_minus_genotypes = $alt_plus_genotypes = 0;
    $alt_minus_het_genotypes = $alt_minus_hom_genotypes = 0;
    $alt_plus_het_genotypes = $alt_plus_hom_genotypes = 0;
    $i = 0;
    open F1, $f or die;
    while (<F1>){
        if (/^#/){
            chomp;
            ($dummy,$sex) = split (/\s+/) if (/Sex/);
            ($dummy,$cov) = split (/\s+/) if (/EstimatedCoverage/);
            next;
        }
        next if ($_ !~ /REF/);
        @tmp = split(/\t/);
        next if ($tag{$tmp[0]} ne "OK"); # Count genotypes for final set of markers only
        ++$i;
        $gt = $tmp[1];

        # Additional quality control of genotypes
        if ($tmp[1] ne "NC"){
            $gt = "NN" if (($tmp[2] < $critical_likelihood_of_the_genotype) or (($tmp[3] + $tmp[4]) < $critical_total_kmer_count));
            $gt = "NN" if (($tmp[1] eq "AB") and (($tmp[3] < $critical_total_kmer_count) or ($tmp[4] < $critical_total_kmer_count)));
        }

        ($chr,$pos,$strand,$type, $dummy) = split(/:/,$tmp[0]);
        if (($sex eq "M") and ($chr eq "X" or $chr eq "Y") and ($gt eq "A" or $gt eq "B")) {
            $expected_genotype += 1;
        }
        elsif ($gt eq "AA" or $gt eq "AB" or $gt eq "BB"){
            $expected_genotype += 1;
        }
        else{
            $unexpected_genotype += 1;
        }
        # Count the number of genotypes that are different from reference
        if (($sex eq "M") and ($chr eq "X" or $chr eq "Y") and ($gt eq "B")) {
            $alt_minus_genotypes += 1 if ($type eq "REF-minus");
            $alt_plus_genotypes += 1 if ($type eq "REF-plus");
        }
        elsif ($gt eq "AB" or $gt eq "BB"){
            if ($type eq "REF-minus"){
                $alt_minus_genotypes += 1;
                $alt_minus_het_genotypes += 1 if (($chr ne "X" and $chr ne "Y") and ($gt eq "AB"));
                $alt_minus_hom_genotypes += 1 if (($chr ne "X" and $chr ne "Y") and ($gt eq "BB"));
            }
            if ($type eq "REF-plus"){
                $alt_plus_genotypes += 1;
                $alt_plus_het_genotypes += 1 if (($chr ne "X" and $chr ne "Y") and ($gt eq "AB"));
                $alt_plus_hom_genotypes += 1 if (($chr ne "X" and $chr ne "Y") and ($gt eq "BB"));
            }
        }
    }
    close F1;
    $expected_genotype_proportion = $expected_genotype / ($expected_genotype + $unexpected_genotype);
    # Don't use individuals with poor overall expected genotype rate
    next if ($expected_genotype_proportion < $critical_proportion_of_expected_gt_for_individual);

    printf (F2 "%s Processed %d lines: SEX: %s COV: %.3f EGP: %.3f EXP: %d ALT_PLUS: %d ALT_MINUS: %d AB_PLUS: %d BB_PLUS: %d AB_MINUS: %d BB_MINUS: %d\n",
      $f, $i, $sex, $cov, $expected_genotype_proportion, $expected_genotype, $alt_plus_genotypes, $alt_minus_genotypes,
      $alt_plus_het_genotypes, $alt_plus_hom_genotypes, $alt_minus_het_genotypes, $alt_minus_hom_genotypes);

}
close F2;
