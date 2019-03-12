$ind = 0;
while (<>){
	chomp;
	@tmp = split(/\s+/);
	if ($tmp[0] eq "Mutation"){
		$mut = $tmp[2];
		$ind += 1 if $mut == 0;
		$sex[$ind] = $tmp[4];
	}
	$N[$ind][$mut] = $tmp[2] if ($tmp[0] eq "N-rich");
	$TP[$ind][$mut] = $tmp[2] if ($tmp[0] eq "True");
	$T[$ind][$mut] = $tmp[3] if ($tmp[0] eq "Total"); 
}

for $j (-1..10){
	for $i (1..$ind){
		print "\t" unless $i == 1;
		if ($j < 0){ print $i . $sex[$i];}
		else{print $N[$i][$j];}
	}
	print "\n";
}
print "\n";

for $j (-1..10){
        for $i (1..$ind){
                print "\t" unless $i ==	1;
                if ($j < 0){ print $i . $sex[$i];}
                else{print $TP[$i][$j];}
        }
        print "\n";
}
print "\n";

for $j (-1..10){
        for $i (1..$ind){
                print "\t" unless $i ==	1;
                if ($j < 0){ print $i . $sex[$i];}
                else{print $T[$i][$j];}
        }
        print "\n";
}
print "\n";
