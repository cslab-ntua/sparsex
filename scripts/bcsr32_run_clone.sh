#!/bin/zsh

set -e

if [ -z "$1" ]; then
	#mtxfiles=$(echo tests.mtx/* /s/matrices/sorted/(???~(024|041|042)).*.mtx.sorted)
	mtxfiles=$(echo /s/matrices/sorted/???.*.mtx.sorted)
else
	mtxfiles="$*"
fi

echo "Log started $(iddate) @$(hostname)"
spmv="./spmv"
for mtx in $(echo $mtxfiles)
do
	echo "*********************** $mtx **********************"

	for mt_conf in "0" "0,2" "0,4" "0,2,4,6" "0,1,2,3,4,5,6,7"; do
		for meth in bcsr_mulv_{1x2,1x3,1x4,2x1,2x2,2x3,2x4,3x1,3x2,4x1,4x2}; do
			echo MT_CONF=$mt_conf $spmv -c -b $mtx $meth
			     MT_CONF=$mt_conf $spmv -c -b $mtx $meth | scripts/spmv_avg.py
		done
	done

	echo '----------------------------------------'
	echo

done
echo "Log ended $(iddate) @$(hostname)"
