#!/bin/zsh

if [ -z "$1" ]; then
	#mtxfiles=$(echo tests.mtx/* /s/matrices/sorted/(???~(024|041|042)).*.mtx.sorted)
	mtxfiles=$(echo /s/matrices/sorted/???.*.mtx.sorted)
else
	mtxfiles="$*"
fi

for mtx in $(echo $mtxfiles)
do
	echo "*********************** $mtx **********************"

	# serial
	echo ./spmv -b $mtx spm_crs32_double_multiply
	./spmv -b $mtx spm_crs32_double_multiply

	echo MT_CONF=0,2 ./spmv -b $mtx spm_crs32_double_mt_multiply
	     MT_CONF=0,2 ./spmv -b $mtx spm_crs32_double_mt_multiply

	echo MT_CONF=0,4 ./spmv -b $mtx spm_crs32_double_mt_multiply
	     MT_CONF=0,4 ./spmv -b $mtx spm_crs32_double_mt_multiply

	echo MT_CONF=0,2,4,6 ./spmv -b $mtx spm_crs32_double_mt_multiply
	     MT_CONF=0,2,4,6 ./spmv -b $mtx spm_crs32_double_mt_multiply

	echo MT_CONF=0,1,2,3,4,5,6,7 ./spmv -b $mtx spm_crs32_double_mt_multiply
	     MT_CONF=0,1,2,3,4,5,6,7 ./spmv -b $mtx spm_crs32_double_mt_multiply

	echo '----------------------------------------'
	echo

done
