
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

	for mt_conf in "0" "0,1" "0,2" "0,8" "0,2,4,6" "0,4,8,12" "0,1,2,3,4,5,6,7" "0,2,4,6,8,10,12,14" "0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15"; do
		for seq_opts in "-d 0" "-d 4" "-d 8"; do
			for opts in "-j -a -t" "-j -t" "-a -t" "-t"; do
				     echo MT_CONF=$mt_conf scripts/csrdu_run.sh -N $opts $seq_opts $mtx
				          MT_CONF=$mt_conf scripts/csrdu_run.sh -N $opts $seq_opts $mtx
			done
		done
	done

	echo '----------------------------------------'
	echo

done
echo "Log ended $(iddate) @$(hostname)"
