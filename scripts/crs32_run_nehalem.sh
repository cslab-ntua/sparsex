#!/bin/zsh

set -e

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

	## add more iterations, since nehalem performance has larger variation
	echo MT_CONF=$mtcnf ./spmv -c -b $mtx -L 7 spm_crs32_double_mt_multiply
	     MT_CONF=$mtcnf ./spmv -c -b $mtx -L 7 spm_crs32_double_mt_multiply | scripts/spmv_avg.py
}

echo "Log started $(iddate) @$(hostname)"
for mtx in $(echo $mtxfiles)
do
	echo "*********************** $mtx **********************"

	#### 2 processors
	# 0,1 same core, two threads
	# 0,2 same die, different cores
	# 0,8 different die
	#### 4 processors
	# 0,2,4,6 same die, different cores
	# 0,4,8,12two dies, two cores each
	#### 8 processors
	# 0,1,2,3,4,5,6,7 same die
	# "0,2,4,6,8,10,12,14" different die, all different cores
	##### 16 processors
	#"0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15"
	for mt in \
	"0" "0,1" "0,2" "0,8" \
	"0,2,4,6" "0,4,8,12" "0,1,2,3,4,5,6,7"  \
	"0,2,4,6,8,10,12,14" \
	"0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15"; do
		run_crs32_mt $mtx $mt
	done
done
echo "Log ended $(iddate) @$(hostname)"
