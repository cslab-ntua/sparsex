#!/bin/bash

print_usage()
{
    echo "Usage: $0 <mmf file> [...]"
}

if [ "x$1" == "x" ]; then
	echo "$0: No input" >&2
    print_usage
	exit 1;
fi

tmp_stripped_input=$(mktemp)
tmp_sorted=$(mktemp)

# Strip comments
for f in $*; do
	if [ ! -e $f ]; then
		echo "$0: No such file: $f" >&2
		continue	
	fi

	if [ -d $f ]; then
		echo "$0: $f is a directory" >&2
		continue
	fi

    echo -n "Converting $f ... "
	grep -v -e '^[[:space:]]*%' $f > $tmp_stripped_input
	file_header=$(head -n 1 $tmp_stripped_input)

	# Sort the file and get rid of the last line; it contains the header.
	# Before sorting, however, you need to pad with zeros the two first
	# columns, in order to properly sort them (`sort's -n option refers to the
	# sorting field only, we want numeric sort to second column as well).
	# Finally, discard padding before producing the final output.
	awk '{ printf "%08s %08s %s\n", $1, $2, $3}' < $tmp_stripped_input \
		| sort -nk 1 | awk '{ printf "%ld %ld %s\n", $1, $2, $3}' \
		| head -n -1 > $tmp_sorted
	echo $file_header | cat - $tmp_sorted > "$f".sorted
    echo "done [output in $f.sorted]"
done

# Cleanup
/bin/rm -f $tmp_stripped_input $tmp_sorted
exit 0
