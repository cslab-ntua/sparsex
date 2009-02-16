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
	echo scripts/csrdu_run.sh -j -a $mtx
	scripts/csrdu_run.sh -j -a $mtx

	# different cores
	echo MT_CONF=$(seq -s ',' 0 8 8) scripts/csrdu_run.sh -j -a -t $mtx
	     MT_CONF=$(seq -s ',' 0 8 8) scripts/csrdu_run.sh -j -a -t $mtx

	echo MT_CONF=$(seq -s ',' 0 8 24) scripts/csrdu_run.sh -j -a -t $mtx
	     MT_CONF=$(seq -s ',' 0 8 24) scripts/csrdu_run.sh -j -a -t $mtx

	echo MT_CONF=$(seq -s ',' 0 8 56) scripts/csrdu_run.sh -j -a -t $mtx
	     MT_CONF=$(seq -s ',' 0 8 56) scripts/csrdu_run.sh -j -a -t $mtx

	#same core
	echo MT_CONF=$(seq -s ',' 0 1) scripts/csrdu_run.sh -j -a -t $mtx
	     MT_CONF=$(seq -s ',' 0 1) scripts/csrdu_run.sh -j -a -t $mtx

	echo MT_CONF=$(seq -s ',' 0 3) scripts/csrdu_run.sh -j -a -t $mtx
	     MT_CONF=$(seq -s ',' 0 3) scripts/csrdu_run.sh -j -a -t $mtx

	echo MT_CONF=$(seq -s ',' 0 5) scripts/csrdu_run.sh -j -a -t $mtx
	     MT_CONF=$(seq -s ',' 0 5) scripts/csrdu_run.sh -j -a -t $mtx

	echo MT_CONF=$(seq -s ',' 0 7) scripts/csrdu_run.sh -j -a -t $mtx
	     MT_CONF=$(seq -s ',' 0 7) scripts/csrdu_run.sh -j -a -t $mtx

	# all
	echo MT_CONF=$(seq -s ',' 0 15) scripts/csrdu_run.sh -j -a -t $mtx
	     MT_CONF=$(seq -s ',' 0 15) scripts/csrdu_run.sh -j -a -t $mtx

	echo MT_CONF=$(seq -s ',' 0 31) scripts/csrdu_run.sh -j -a -t $mtx
	     MT_CONF=$(seq -s ',' 0 31) scripts/csrdu_run.sh -j -a -t $mtx

	echo MT_CONF=$(seq -s ',' 0 63) scripts/csrdu_run.sh -j -a -t $mtx
	     MT_CONF=$(seq -s ',' 0 63) scripts/csrdu_run.sh -j -a -t $mtx

	echo
	echo
done
