#!/bin/zsh

ttu_lim=1.5
usage() {
	echo "Usage: $(basename $0) <options> mmf_file"
	echo "Options:"
	echo " -j ........... use jump"
	echo " -a ........... aligned deltas"
	echo " -d ........... dense min len (0 to disable)"
	echo " -q ........... quiet"
	echo " -c ........... check matrix"
	echo " -N ........... use NUMA local allocation"
	echo " -l ........... ttu limit (default:$ttu_lim)"
	echo " -u ........... unique value (default:calculate it)"
	echo " (set MT_CONF for multithreaded version) "
}

jmp=0;
aligned=0
verbose=1
de_minlen=0
mt=0
check_opt=""

while getopts "hjad:qcNl:u:" option
do
	case $option in
		h ) usage; exit 0 ;;
		j ) jmp=1 ;;
		a ) aligned=1 ;;
		d ) de_minlen="$OPTARG" ;;
		q ) verbose=0 ;;
		c ) check_opt="-c" ;;
		N ) numa=1 ;;
		l ) ttu_lim="$OPTARG" ;;
		u ) uvals="$OPTARG" ;;
		* ) echo "Unknown option"; exit 1;;
	esac
done

shift $(($OPTIND-1))

if [ -z "$1" ]; then
	usage;
	exit 1
fi

file="$1"

if [ -z "$uvals" ]; then
	uvals=$(sed '1d; s/[[:digit:]]\+ [[:digit:]]\+[[:space:]]\+//' <$file | sort -u | wc -l)
fi
vals=$(head -1 $file | awk '{ print $3 }')
ttu=$(( ($vals+0.0)/($uvals+0.0) ))

#set -x
set -e
if [ $(( $ttu > $ttu_lim )) -eq 1 ] ; then
	vi_bits=$(scripts/calc_bits $uvals)
	method="spm_csrdu_vi${vi_bits}_double"
	if [ -n "$MT_CONF" ]; then
		method="${method}_mt"
	fi
	[ "$aligned" -eq "1" ] && method="${method}_aligned"
	[ "$jmp" -eq "1" ] && method="${method}_jmp"
	[ "$numa" -eq "1" ] && method="${method}_numa"
	method="${method}_multiply"
	CSRDU_VERBOSE=$verbose CSRDU_DE_MINLEN=$de_minlen CSRDU_JMP=$jmp CSRDU_ALIGNED=$aligned ./spmv -b $check_opt $1 $method
fi
