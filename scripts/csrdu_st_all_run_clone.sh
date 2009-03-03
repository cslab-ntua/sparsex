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

	for seq_opts in "-d 0" "-d 4" "-d 6" "-d 8" "-d 10" "-d 12"; do
		for opts in "-j -a " "-j " "-a " ""; do
			     echo scripts/csrdu_run.sh $opts $seq_opts $mtx
				  scripts/csrdu_run.sh $opts $seq_opts $mtx
		done
	done

	echo '----------------------------------------'
	echo

done
echo "Log ended $(iddate) @$(hostname)"
