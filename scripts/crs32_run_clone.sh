#!/bin/zsh

if [ -z "$1" ]; then
	#mtxfiles=$(echo tests.mtx/* /s/matrices/sorted/(???~(024|041|042)).*.mtx.sorted)
	mtxfiles=$(echo /s/matrices/sorted/???.*.mtx.sorted)
else
	mtxfiles="$*"
fi

spmv=./spmv

for mtx in $(echo $mtxfiles)
do
	echo "*********************** $mtx **********************"
	for mt_conf in "0" "0,2"  "0,4"  "0,2,4,6"  "0,1,2,3,4,5,6,7"; do
		echo MT_CONF=$mt_conf $spmv -c -b $mtx spm_crs32_double_mt_multiply
		     MT_CONF=$mt_conf $spmv -c -b $mtx spm_crs32_double_mt_multiply | scripts/spmv_avg.py
	done
	echo '----------------------------------------'
	echo

done
