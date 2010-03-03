#!/bin/bash

#set -x
set -e

if [ -z "$1" ]; then
	#mtxfiles=$(echo tests.mtx/* /s/matrices/sorted/(???~(024|041|042)).*.mtx.sorted)
	mtxfiles=$(echo /s/matrices/sorted/???.*.mtx.sorted)
else
	mtxfiles="$*"
fi

echo "Log started $(iddate) @$(hostname)"
for mtx in $(echo $mtxfiles)
do
	echo "*********************** $mtx **********************"

	#scripts/csrvi_run.sh -c $mtx
	uvals=$(sed '1d; s/[[:digit:]]\+ [[:digit:]]\+[[:space:]]\+//' <$mtx | sort -u | wc -l)
	opts="-c -t -u $uvals"
	for mt_conf in "0" "0,2" "0,4" "0,2,4,6" "0,1,2,3,4,5,6,7"; do
		echo MT_CONF=$mt_conf scripts/csrvi_run.sh $opts $mtx
		MT_CONF=$mt_conf scripts/csrvi_run.sh $opts $mtx | scripts/spmv_avg.py
	done

	echo '----------------------------------------'
	echo

done
echo "Log ended $(iddate) @$(hostname)"
