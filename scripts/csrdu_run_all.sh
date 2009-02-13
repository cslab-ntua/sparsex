#!/bin/zsh

check_opt=""
while getopts "c" option
do
	case $option in
		c ) check_opt=" -c " ;;
		* ) echo "Unknown option"; exit 1;;
	esac
done

shift $(($OPTIND-1))

if [ -z "$1" ]; then
	#mtxfiles=$(echo tests.mtx/* /s/matrices/sorted/(???~(024|041|042)).*.mtx.sorted)
	mtxfiles=$(echo /s/matrices/sorted/???.*.mtx.sorted)
else
	mtxfiles="$*"
fi

for mtx in $(echo $mtxfiles)
do
	scripts/csrdu_run.sh $check_opt $mtx
	scripts/csrdu_run.sh $check_opt -j $mtx
	scripts/csrdu_run.sh $check_opt -a $mtx
	scripts/csrdu_run.sh $check_opt -j -a $mtx
	echo
done
