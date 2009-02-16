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

	# different cores
	echo MT_CONF=$(seq -s ',' 0 8 8) ./spmv -b $mtx spm_crs32_double_mt_multiply
	     MT_CONF=$(seq -s ',' 0 8 8) ./spmv -b $mtx spm_crs32_double_mt_multiply

	echo MT_CONF=$(seq -s ',' 0 8 24) ./spmv -b $mtx spm_crs32_double_mt_multiply
	     MT_CONF=$(seq -s ',' 0 8 24) ./spmv -b $mtx spm_crs32_double_mt_multiply

	echo MT_CONF=$(seq -s ',' 0 8 56) ./spmv -b $mtx spm_crs32_double_mt_multiply
	     MT_CONF=$(seq -s ',' 0 8 56) ./spmv -b $mtx spm_crs32_double_mt_multiply

	#same core
	echo MT_CONF=$(seq -s ',' 0 1) ./spmv -b $mtx spm_crs32_double_mt_multiply
	     MT_CONF=$(seq -s ',' 0 1) ./spmv -b $mtx spm_crs32_double_mt_multiply

	echo MT_CONF=$(seq -s ',' 0 3) ./spmv -b $mtx spm_crs32_double_mt_multiply
	     MT_CONF=$(seq -s ',' 0 3) ./spmv -b $mtx spm_crs32_double_mt_multiply

	echo MT_CONF=$(seq -s ',' 0 5) ./spmv -b $mtx spm_crs32_double_mt_multiply
	     MT_CONF=$(seq -s ',' 0 5) ./spmv -b $mtx spm_crs32_double_mt_multiply

	echo MT_CONF=$(seq -s ',' 0 7) ./spmv -b $mtx spm_crs32_double_mt_multiply
	     MT_CONF=$(seq -s ',' 0 7) ./spmv -b $mtx spm_crs32_double_mt_multiply

	# all
	echo MT_CONF=$(seq -s ',' 0 15) ./spmv -b $mtx spm_crs32_double_mt_multiply
	     MT_CONF=$(seq -s ',' 0 15) ./spmv -b $mtx spm_crs32_double_mt_multiply

	echo MT_CONF=$(seq -s ',' 0 31) ./spmv -b $mtx spm_crs32_double_mt_multiply
	     MT_CONF=$(seq -s ',' 0 31) ./spmv -b $mtx spm_crs32_double_mt_multiply

	echo MT_CONF=$(seq -s ',' 0 63) ./spmv -b $mtx spm_crs32_double_mt_multiply
	     MT_CONF=$(seq -s ',' 0 63) ./spmv -b $mtx spm_crs32_double_mt_multiply

	echo
	echo
done
