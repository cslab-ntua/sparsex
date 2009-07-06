#!/bin/zsh
if [ -z "$1" ]; then
	echo "Usage: $0 <file>"
	exit 1;
fi

while [ -n "$1" ]; do
	echo "************* Matrix: $1"
	f="/tmp/jit_test-$(basename $1)"
	./jit $1 2>/dev/null 1> $f

	diff -u =(sed -e "1d; s/\([[:digit:]]\+\)[[:space:]]\+\([[:digit:]]\+\).*$/\1 \2/" < $1) =(sort -V $f)

	rm $f
	shift
done
