#!/bin/zsh

check_opt=""
while getopts "c" option
do
	case $option in
		c ) check_opt="-c" ;;
		* ) echo "Unknown option"; exit 1;;
	esac
done
shift $(($OPTIND-1))

echo $1
if [ -z "$1" ]; then
	mtxfiles=$(echo /s/matrices/sorted/???.*.mtx.sorted)
else
	mtxfiles="$*"
fi

for mtx in $(echo $mtxfiles)
do
	scripts/csrvi_run.sh $check_opt $mtx
	echo
done
