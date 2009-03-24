#!/usr/bin/zsh

ttu=1.5
usage() {
	echo "Usage: $(basename $0) <options> mmf_file"
	echo "Options:"
	echo " -t ........... use multithreaded version"
	echo " -q ........... quiet"
	echo " -c ........... check matrix"
	echo " -N ........... use NUMA local allocation"
	echo " -l ........... ttu limit (default:$ttu)"
}

while getopts "htqcNl:" option
do
	case $option in
		h ) usage; exit 0 ;;
		t ) mt=1 ;;
		q ) verbose=0 ;;
		c ) check_opt="-c" ;;
		N ) numa=1 ;;
		l ) ttu="$OPTARG" ;;
		* ) echo "Unknown option"; exit 1;;
	esac
done


shift $(($OPTIND-1))
if [ -z "$1" ]; then
	usage;
	exit 1
fi

file="$1"
uvals=$(sed '1d; s/[[:digit:]]\+ [[:digit:]]\+[[:space:]]\+//' <$file | sort -u | wc -l)
vals=$(head -1 $file | awk '{ print $3 }')

set -x
if [ $(( ($vals+0.0)/($uvals+0.0) > $ttu )) -eq 1 ] ; then
	vi_bits=$(scripts/calc_bits $uvals)
	method="spm_crs32_vi${vi_bits}_double"
	[ "$mt" -eq "1" ] && method="${method}_mt"
	method="${method}_multiply"
	./spmv -b $check_opt $file $method
fi
