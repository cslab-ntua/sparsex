#!/bin/zsh

if [ -z "$1" ]; then
	#mtxfiles=$(echo tests.mtx/* /s/matrices/sorted/(???~(024|041|042)).*.mtx.sorted)
	mtxfiles=$(echo /s/matrices/sorted/???.*.mtx.sorted)
else
	mtxfiles="$*"
fi

function run_crs32_mt {
	mtx="$1"
	mtcnf="$2"

	if [ -z "$mtcnf" ]; then
		echo "MT_CONF NOT SPECIFIED"
		exit 1
	fi

	echo MT_CONF=$mtcnf ./spmv -b $mtx spm_crs32_double_mt_multiply
	     MT_CONF=$mtcnf ./spmv -b $mtx spm_crs32_double_mt_multiply
}

echo "Log started $(iddate) @$(hostname)"
for mtx in $(echo $mtxfiles)
do
	echo "*********************** $mtx **********************"

	# serial
	echo ./spmv -b $mtx spm_crs32_double_multiply
	./spmv -b $mtx spm_crs32_double_multiply


	#### 2 processors
	# same core, two threads
	run_crs32_mt $mtx "0,1"
	# same die, different cores
	run_crs32_mt $mtx "0,2"
	# different die
	run_crs32_mt $mtx "0,8"

	#### 4 processors
	# same die, different cores
	run_crs32_mt $mtx "0,2,4,6"
	# two dies, two cores each
	run_crs32_mt $mtx "0,4,8,12"

	#### 8 processors
	# same die
	run_crs32_mt $mtx "0,1,2,3,4,5,6,7"
	# different die, all different cores
	run_crs32_mt $mtx "0,2,4,6,8,10,12,14"

	##### 16 processors
	run_crs32_mt $mtx "0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15"
done
echo "Log ended $(iddate) @$(hostname)"
