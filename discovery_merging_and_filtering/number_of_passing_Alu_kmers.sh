for i in "filter_gt_markers.txt" "filter_gt_markers_OK.txt"
do
	echo -n "Total number of candidates: "
	wc -l $i
	echo ""
	echo -n "Remaining markers: "
        grep -c -E 'OK|RRP' $i
	echo -n "Remaining REF-plus: "
	grep -E 'OK|RRP' $i | grep -c plus 
	echo -n "Remaining REF-minus: "
	grep -E 'OK|RRP' $i | grep -c minus
	echo ""
	echo -n "Removed for EGT: "
	grep -c 'EGT' $i
	echo -n "Removed for HWE: "
	grep -c 'HWE' $i 
	echo -n "Removed for MFB: "
	grep -c 'MFB' $i 
	echo -n "Removed for ZRM: "
	grep -c 'ZRM' $i
	echo -n "Removed for ZRP: "
 	grep -c 'ZRP' $i
	echo ""
done
