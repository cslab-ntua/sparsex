#/bin/bash

program_name=$(basename $0)

function print_usage
{
    cat << _USAGE_EOF
Usage: $program_name [MMF_FILE]...
_USAGE_EOF
}

if [ $# -lt 1 ]; then
    print_usage
    exit
fi

for f in $*; do
    diff_out=$(basename $f).diff
    echo -n "Running $f..."
    ./spmv $f | ../scripts/mmf_eq2.sh $f - > $diff_out
    echo " done"
done
