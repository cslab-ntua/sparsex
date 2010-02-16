#!/bin/bash

program_name=$(basename $0)

function print_usage
{
    cat << _USAGE_EOF
Usage: $program_name MMF_FILE SPM_OUT
_USAGE_EOF
}

if [ $# -lt 2 ]; then
    print_usage
    exit 1
fi

mm_file=$1
spm_out_file=$2

tmp_mm_file=$(mktemp)
tmp_spm_out_file=$(mktemp)

sed -e '1d' $mm_file | \
    awk '{ printf("%ld %ld %.4e\n", $1, $2, $3) }' > $tmp_mm_file
sed -e '1,3 d; s/[[:space:]]*cnt.*//' $spm_out_file | \
    awk '{ printf("%08ld %08ld %.4e\n", $1, $2, $3) }' | \
    sort -n | awk '{ printf("%ld %ld %.4e\n", $1, $2, $3) }' > $tmp_spm_out_file

diff $tmp_mm_file $tmp_spm_out_file

# echo "MMF -> $tmp_mm_file" 1>&2
# echo "OUT -> $tmp_spm_out_file" 1>&2
/bin/rm -f $tmp_mm_file
/bin/rm -f $tmp_spm_out_file
