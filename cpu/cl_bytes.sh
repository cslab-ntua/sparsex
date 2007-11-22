#!/bin/sh

cl=$(sed -n -e 's/cache_alignment\t: \([[:digit:]]\+\)/\1/p' < /proc/cpuinfo | head -1)
if [ -n "$cl" ]; then
	echo $cl
	exit 0
fi

cl=$(sed -n -e 's/clflush size\t: \([[:digit:]]\+\)/\1/p' < /proc/cpuinfo | head -1)
if [ -n "$cl" ]; then
	echo $cl
	exit 0
fi


# give up!
echo "===>ERROR: The $0 script can't find the cacheline size info"
exit 1

