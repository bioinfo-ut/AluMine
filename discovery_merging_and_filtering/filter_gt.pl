#!/bin/perl -w

# Filter markers based on their predicted genotypes in 2280 individuals
# remove individuals with low proportion of expected genotypes (expected are diploid genotypes for autosomes, haploid for X and Y in men )
# mark EGT low proportion of expected genotypes (expected are diploid genotypes for autosomes, haploid for X and Y in men )
# mark HWE low probability for Hardy-Weinberg Equilibrium
# mark MFB B-allele significantly different in men and women
# mark ZRM zero freq REF-minus
# mark ZRP zero freq REF-plus (OK to keep)
# mark the rest as OK


chdir "/storage7/ctg/uued/calling_gmer/calls/ALU_SNP/";
$critical_proportion_of_expected_gt_for_individual = 0.50;
$critical_proportion_of_expected_gt_for_marker = 0.90;
$critical_hwe_chisquare = 25.97828096; # P<0.01, with multiple testing: P=0.01/28,962 => 3.45E-7 => chidist(25.978,1)
$critical_mfb_chisquare = 6.634896601; # P<0.01, without multiple testing: difference between men BAF and women BAF is suspicious
$critical_likelihood_of_the_genotype = 0.9; # don't use any genotypes that have lower likelihood than this
$critical_total_kmer_count = 5;             # don't use any genotypes that are based on k-mer frequency lower than 5 (median diploid depth of coverage in this set is 22.8)


$outfile2 = "/mambakodu/mremm/filter_gt_files.txt";
$outfile3 = "/mambakodu/mremm/filter_gt_markers.txt";
# Frequency of autosomal genotypes
%A_marker_count = ();
%FA_AA = (); %MA_AA = (); # female count, male count
%FA_AB = (); %MA_AB = ();
%FA_BB = (); %MA_BB = ();
%FA_NN = (); %MA_NN = ();

# Frequency of chrX genotypes
%X_marker_count = ();
%FX_AA = ();
%FX_AB = ();
%FX_BB = ();
%FX_NN = ();
%MX_A = ();
%MX_B = ();
%MX_N = ();

# Frequency of chrY genotypes
%Y_marker_count = ();
%MY_A = ();
%MY_B = ();
%MY_N = ();

opendir(my $dh, ".") || die "Can't opendir: $!";
my @files = readdir($dh);
closedir $dh;

open(F2, '>', $outfile2) or die;
foreach $f (@files){
    next if (-s $f < 1000000000 or $f !~ /calls/); # V538* for testing
    $sex = ""; $gt = ""; $id = ""; %gt = ();
    $expected_genotype_proportion = 0;
    $expected_genotype = $unexpected_genotype = 0;
    $alt_minus_genotypes = $alt_plus_genotypes = 0;
    $alt_minus_het_genotypes = $alt_minus_hom_genotypes = 0;
    $alt_plus_het_genotypes = $alt_plus_hom_genotypes = 0;
    $i = -1;
    # Calculate the fraction of expected genotypes for given individual
    open F1, $f or die;
    while (<F1>){
        if (/^#/){
            chomp;
            ($dummy,$sex) = split (/\s+/) if (/Sex/);
            ($dummy,$cov) = split (/\s+/) if (/EstimatedCoverage/);
            next;
        }
        next if ($_ !~ /REF/);
        ++$i;
        @tmp = split(/\t/);
        $gt = $tmp[1];

        # Additional quality control of genotypes
        if ($tmp[1] ne "NC"){
            $gt = "NN" if (($tmp[2] < $critical_likelihood_of_the_genotype) or (($tmp[3] + $tmp[4]) < $critical_total_kmer_count));
            $gt = "NN" if (($tmp[1] eq "AB") and (($tmp[3] < $critical_total_kmer_count) or ($tmp[4] < $critical_total_kmer_count)));
        }

        $gt{$tmp[0]} = $gt;
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
    printf (F2 "%s Processed %d lines: SEX: %s COV: %.3f EGP: %.3f EXP: %d ALT_PLUS: %d ALT_MINUS: %d AB_PLUS: %d BB_PLUS: %d AB_MINUS: %d BB_MINUS: %d\n",
      $f, $i, $sex, $cov, $expected_genotype_proportion, $expected_genotype, $alt_plus_genotypes, $alt_minus_genotypes,
      $alt_plus_het_genotypes, $alt_plus_hom_genotypes, $alt_minus_het_genotypes, $alt_minus_hom_genotypes);

    # Don't use individuals with poor overall expected genotype rate
    next if ($expected_genotype_proportion < $critical_proportion_of_expected_gt_for_individual);

    # Process all genotypes of given individual
    for $id (keys %gt){
        $gt = $gt{$id};
        ($chr,$dummy) = split(/:/,$id);
        if ($chr eq "X" and $sex eq "F"){
            $X_marker_count{$id} += 1;
            if    ($gt eq "AA"){$FX_AA{$id} += 1;$FX_AB{$id} += 0;$FX_BB{$id} += 0;$FX_NN{$id} += 0;$ec+=1;}
            elsif ($gt eq "AB"){$FX_AA{$id} += 0;$FX_AB{$id} += 1;$FX_BB{$id} += 0;$FX_NN{$id} += 0;$ec+=1;}
            elsif ($gt eq "BB"){$FX_AA{$id} += 0;$FX_AB{$id} += 0;$FX_BB{$id} += 1;$FX_NN{$id} += 0;$ec+=1;}
            else               {$FX_AA{$id} += 0;$FX_AB{$id} += 0;$FX_BB{$id} += 0;$FX_NN{$id} += 1;$uc+=1;}
        }
        elsif ($chr eq "X" and $sex eq "M"){
            $X_marker_count{$id} += 1;
            if    ($gt eq "A"){$MX_A{$id} += 1;$MX_B{$id} += 0;$MX_N{$id} += 0;$ec+=1;}
            elsif ($gt eq "B"){$MX_A{$id} += 0;$MX_B{$id} += 1;$MX_N{$id} += 0;$ec+=1;}
            else              {$MX_A{$id} += 0;$MX_B{$id} += 0;$MX_N{$id} += 1;$uc+=1;}
        }
        elsif ($chr eq "Y" and $sex eq "M"){
            $Y_marker_count{$id} += 1;
            if    ($gt eq "A"){$MY_A{$id} += 1;$MY_B{$id} += 0;$MY_N{$id} += 0;$ec+=1;}
            elsif ($gt eq "B"){$MY_A{$id} += 0;$MY_B{$id} += 1;$MY_N{$id} += 0;$ec+=1;}
            else              {$MY_A{$id} += 0;$MY_B{$id} += 0;$MY_N{$id} += 1;$uc+=1;}
        }
        elsif ($sex eq "M"){ # male autosomes
            $A_marker_count{$id} += 1;
            if    ($gt eq "AA"){$MA_AA{$id} += 1;$MA_AB{$id} += 0;$MA_BB{$id} += 0;$MA_NN{$id} += 0;$ec+=1;}
            elsif ($gt eq "AB"){$MA_AA{$id} += 0;$MA_AB{$id} += 1;$MA_BB{$id} += 0;$MA_NN{$id} += 0;$ec+=1;}
            elsif ($gt eq "BB"){$MA_AA{$id} += 0;$MA_AB{$id} += 0;$MA_BB{$id} += 1;$MA_NN{$id} += 0;$ec+=1;}
            else               {$MA_AA{$id} += 0;$MA_AB{$id} += 0;$MA_BB{$id} += 0;$MA_NN{$id} += 1;$uc+=1;}
        }
        elsif ($sex eq "F") { # female autosomes
            $A_marker_count{$id} += 1;
            if    ($gt eq "AA"){$FA_AA{$id} += 1;$FA_AB{$id} += 0;$FA_BB{$id} += 0;$FA_NN{$id} += 0;$ec+=1;}
            elsif ($gt eq "AB"){$FA_AA{$id} += 0;$FA_AB{$id} += 1;$FA_BB{$id} += 0;$FA_NN{$id} += 0;$ec+=1;}
            elsif ($gt eq "BB"){$FA_AA{$id} += 0;$FA_AB{$id} += 0;$FA_BB{$id} += 1;$FA_NN{$id} += 0;$ec+=1;}
            else               {$FA_AA{$id} += 0;$FA_AB{$id} += 0;$FA_BB{$id} += 0;$FA_NN{$id} += 1;$uc+=1;}
        }
        else{
            print STDERR "Unclassified marker $id\n";
        }
    }
}
close F2;

# Process all Alu markers and divide them into categories EGR, HWE, BAF, ZRM, ZRP, RRP, OK
open(F3, '>', $outfile3) or die;
print F3 "STATUS\tNAME\tBAF\tFBAF\tMBAF\tAA\tAB\tBB\tNN\tA\tB\tN\n";
foreach $k (sort keys %A_marker_count){
    ($chr,$loc,$strand,$type,$dummy) = split(/:/,$k);
    $AA = $AB = $BB = $NN = 0;
    $AA = $FA_AA{$k} + $MA_AA{$k};
    $AB = $FA_AB{$k} + $MA_AB{$k};
    $BB = $FA_BB{$k} + $MA_BB{$k};
    $NN = $FA_NN{$k} + $MA_NN{$k};
    $expected_calls = $AA + $AB + $BB;
    $all_calls = $expected_calls + $NN;
    $proportion_of_expected_gt = $expected_calls / $all_calls;
    $HWE = &calc_HWE($AA,$AB,$BB); # chisquare statistic
    $BAF = &calc_pBA($AA,$AB,$BB); # b-allele frequency
    $MBAF = &calc_pBA($MA_AA{$k},$MA_AB{$k},$MA_BB{$k});
    $FBAF = &calc_pBA($FA_AA{$k},$FA_AB{$k},$FA_BB{$k});
    $MFB = &calc_chisquare(2*$MA_AA{$k}+$MA_AB{$k},2*$MA_BB{$k}+$MA_AB{$k},2*$FA_AA{$k}+$FA_AB{$k}, 2*$FA_BB{$k}+$FA_AB{$k}); # chisquare statistic

    if ($proportion_of_expected_gt < $critical_proportion_of_expected_gt_for_marker){
        print F3 "EGT\t";
    }
    elsif ($HWE ne "NA" and $HWE > $critical_hwe_chisquare){
        print F3 "HWE\t";
    }
    elsif ($MFB ne "NA" and $MFB > $critical_mfb_chisquare){
        print F3 "MFB\t";
    }
    elsif ($type eq "REF-minus" and ($AB + $BB) == 0){
        print F3 "ZRM\t";
    }
    elsif ($type eq "REF-plus" and ($AB + $BB) == 0){
        print F3 "ZRP\t";
    }
    else{
        print F3 "OK \t";
    }
    printf (F3 "%s\t%s\t%s\t%s\t%d\t%d\t%d\t%d\t-\t-\t-\n", $k,$BAF,$FBAF,$MBAF,$AA,$AB,$BB,$NN);
}

foreach $k (sort keys %X_marker_count){
    ($chr,$loc,$strand,$type,$dummy) = split(/:/,$k);
    $expected_calls = $FX_AA{$k} + $FX_AB{$k} + $FX_BB{$k} + $MX_A{$k} + $MX_B{$k};
    $all_calls = $expected_calls + $FX_NN{$k} + $MX_N{$k};
    $proportion_of_expected_gt = $expected_calls / $all_calls;
    $HWE = &calc_HWE($FX_AA{$k},$FX_AB{$k},$FX_BB{$k}); # chisquare statistic
    $BAF = &calc_pBX($FX_AA{$k},$FX_AB{$k},$FX_BB{$k},$MX_A{$k},$MX_B{$k}); # b-allele frequency
    $MBAF = &calc_pBY($MX_A{$k},$MX_B{$k});
    $FBAF = &calc_pBA($FX_AA{$k},$FX_AB{$k},$FX_BB{$k});
    $MFB = &calc_chisquare($MX_A{$k},$MX_B{$k},2*$FX_AA{$k}+$FX_AB{$k}, 2*$FX_BB{$k}+$FX_AB{$k}); # chisquare statistic

    if ($proportion_of_expected_gt < $critical_proportion_of_expected_gt_for_marker){
        print F3 "EGT\t";
    }
    elsif ($HWE ne "NA" and $HWE > $critical_hwe_chisquare){
        print F3 "HWE\t";
    }
    elsif ($MFB ne "NA" and $MFB > $critical_mfb_chisquare){
        print F3 "MFB\t";
    }
    elsif ($type eq "REF-minus" and ($FX_AB{$k} + $FX_BB{$k} + $MX_B{$k}) == 0){
        print F3 "ZRM\t";
    }
    elsif($type eq "REF-plus" and ($FX_AB{$k} + $FX_BB{$k} + $MX_B{$k}) == 0){
        print F3 "ZRP\t";
    }
    else{
        print F3 "OK \t";
    }
    printf (F3 "%s\t%s\t%s\t%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n", $k,$BAF,$FBAF,$MBAF,$FX_AA{$k},$FX_AB{$k},$FX_BB{$k},$FX_NN{$k},$MX_A{$k},$MX_B{$k},$MX_N{$k});
}

foreach $k (sort keys %Y_marker_count){
    ($chr,$loc,$strand,$type,$dummy) = split(/:/,$k);
    $expected_calls = $MY_A{$k} + $MY_B{$k};
    $all_calls = $expected_calls + $MY_N{$k};
    $proportion_of_expected_gt = $expected_calls / $all_calls;
    $BAF = &calc_pBY($MY_A{$k},$MY_B{$k});
    $MBAF = &calc_pBY($MY_A{$k},$MY_B{$k});
    if ($proportion_of_expected_gt < $critical_proportion_of_expected_gt_for_marker){
        print F3 "EGT\t";
    }
    elsif($type eq "REF-minus" and $MY_B{$k} == 0){
        print F3 "ZRM\t";
    }
    elsif($type eq "REF-plus" and $MY_B{$k} == 0){
        print F3 "ZRP\t";
    }
    else{
        print F3 "OK \t";
    }
    printf (F3 "%s\t%s\t-\t%s\t-\t-\t-\t-\t%d\t%d\t%d\n", $k,$BAF,$MBAF,$MY_A{$k},$MY_B{$k},$MY_N{$k});
}
close F3;

sub calc_HWE(){ # returns string or chisquare
    my ($oaa, $oab, $obb) = @_; # observed counts
    my $total = $oaa + $oab + $obb;
    return "NA" unless ($total);
    my $pA = (2.0 * $oaa + $oab) / (2 * $total);
    my $pB = 1.0 - $pA;
    return "NA" unless ($pA * $pB);
    my $eaa = $pA * $pA * $total; # expected counts
    my $eab = 2.0 * $pA * $pB * $total;
    my $ebb = $pB * $pB * $total;
    my $chisquare = ($oaa-$eaa)*($oaa-$eaa)/$eaa + ($oab-$eab)*($oab-$eab)/$eab + ($obb-$ebb)*($obb-$ebb)/$ebb;
    return $chisquare;
}

sub calc_pBA(){ # returns string
    my ($oaa, $oab, $obb) = @_;
    my $total = 2 * ($oaa + $oab + $obb);
    return "NA" unless ($total);
    return sprintf("%.6f", ((2 * $obb + $oab) / $total));
}
sub calc_pBY(){
    my ($oa, $ob) = @_;
    my $total = $oa + $ob;
    return "NA" unless ($total);
    return sprintf("%.6f", ($ob / $total));
}
sub calc_pBX(){
    my ($oaa, $oab, $obb, $oa, $ob) = @_;
    my $total = 2 * ($oaa + $oab + $obb) + $oa + $ob;
    return "NA" unless ($total);
    return sprintf("%.6f", ((2 * $obb + $oab + $ob) / $total));
}
sub calc_chisquare(){ # input four integers
    my ($oa, $ob, $oc, $od) = @_;
    my $sum_ab=$oa+$ob;
    my $sum_cd=$oc+$od;
    my $sum_ac=$oa+$oc;
    my $sum_bd=$ob+$od;
    my $total = $oa + $ob + $oc + $od;
    return "NA" unless ($total);
    my $ea = $sum_ab * $sum_ac / $total;
    my $eb = $sum_ab * $sum_bd / $total;
    my $ec = $sum_cd * $sum_ac / $total;
    my $ed = $sum_cd * $sum_bd / $total;
    return "NA" unless ($ea*$eb*$ec*$ed);
    my $chisquare = ($oa-$ea)*($oa-$ea)/$ea + ($ob-$eb)*($ob-$eb)/$eb + ($oc-$ec)*($oc-$ec)/$ec + ($od-$ed)*($od-$ed)/$ed;
    return sprintf("%.3f",$chisquare);
}
