#!/usr/bin/zsh

for file in $*
do
	# CRS
	./spmv_crs $file

	# Delta
	SPM_DELTA_SEQ_LIMIT=10000 python2.5 python/delta.py $file
	python2.5 python/delta.py $file

	# CRSVI
	uvals=$(sed '1d; s/[[:digit:]]\+ [[:digit:]]\+[[:space:]]\+//' <$file | sort -u | wc -l) 
	vals=$(head -1 $file | awk '{ print $3 }')
	if [ $(( ($vals+0.0)/($uvals+0.0) > 50 )) -eq 1 ] ; then
		./spmv_crsvi $file $uvals
	fi

done | tee ~/spm-results.$(date +%Y%m%d.%H:%M:%S).$(hostname)
