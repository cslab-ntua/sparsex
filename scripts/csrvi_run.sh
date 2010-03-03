#!/usr/bin/zsh

set -e

ttu_lim=1.5
usage() {
	echo "Usage: $(basename $0) <options> mmf_file"
	echo "Options:"
	echo " -t ........... use multithreaded version"
	echo " -q ........... quiet"
	echo " -c ........... check matrix"
	echo " -N ........... use NUMA local allocation"
	echo " -l ........... ttu limit (default:$ttu_lim)"
	echo " -u ........... unique value (default:calculate it)"
}

while getopts "htqcNl:u:" option
do
	case $option in
		h ) usage; exit 0 ;;
		t ) mt=1 ;;
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
if [ $(( $ttu > $ttu_lim )) -eq 1 ] ; then
	vi_bits=$(scripts/calc_bits $uvals)
	method="spm_crs32_vi${vi_bits}_double"
	[ "$mt" -eq "1" ] && method="${method}_mt"
	[ "$numa" -eq "1" ] && method="${method}_numa"
	method="${method}_multiply"
	./spmv -b $check_opt $file $method
fi
