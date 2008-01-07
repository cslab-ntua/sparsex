#!/usr/bin/zsh

limit=1.1
for file in $*
do

	# CRSVI
	uvals=$(sed '1d; s/[[:digit:]]\+ [[:digit:]]\+[[:space:]]\+//' <$file | sort -u | wc -l) 
	vals=$(head -1 $file | awk '{ print $3 }')
	if [ $(( ($vals+0.0)/($uvals+0.0) > 2 )) -eq 1 ] ; then
		continue
	fi
	if [ $(( ($vals+0.0)/($uvals+0.0) > $limit )) -eq 1 ] ; then
		echo ./spmv_crsvi $file $uvals
	fi

done | tee ~/spm-results.crsvi.$(date +%Y%m%d.%k:%M:%S).$(hostname)
