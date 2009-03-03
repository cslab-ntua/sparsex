#!/bin/zsh

usage() {
	echo "Usage: $(basename $0) <options> mmf_file"
	echo "Options:"
	echo " -j ........... use jump"
	echo " -a ........... aligned deltas"
	echo " -t ........... use multithreaded version"
	echo " -d ........... dense min len (0 to disable)"
	echo " -q ........... quiet"
	echo " -c ........... check matrix"
}

jmp=0;
aligned=0
verbose=1
de_minlen=0
mt=0
check_opt=""

while getopts "hjad:tqc" option
do
	case $option in
		h ) usage; exit 0 ;;
		j ) jmp=1 ;;
		a ) aligned=1 ;;
		d ) de_minlen="$OPTARG" ;;
		t ) mt=1 ;;
		q ) verbose=0 ;;
		c ) check_opt="-c" ;;
		* ) echo "Unknown option"; exit 1;;
	esac
done

method="spm_csrdu_double"
[ "$mt" -eq "1" ] && method="${method}_mt"
[ "$aligned" -eq "1" ] && method="${method}_aligned"
[ "$jmp" -eq "1" ] && method="${method}_jmp"
method="${method}_multiply"

shift $(($OPTIND-1))

if [ -z "$1" ]; then
	usage;
	exit 1
fi

set -x
CSRDU_VERBOSE=$verbose CSRDU_DE_MINLEN=$de_minlen CSRDU_JMP=$jmp CSRDU_ALIGNED=$aligned ./spmv -b $check_opt $1 $method
