#!/bin/sh

# find cache size info
# Kornilios Kourtis <kkourt@cslab.ece.ntua.gr>

cs=$(sed -r -n '
	s/^cache size[[:space:]]+: /(/
	s/ KB/*1024)/p
	s/ MB/*1024*1024)/p' < /proc/cpuinfo | head -1)

if [ -n "$cs" ]; then
	echo $cs
	exit 0
fi

f="/sys/devices/system/cpu/cpu0/l2_cache_size"
if [ -f $f ]; then
	echo $(cat $f)
	exit 0
fi

echo "===>ERROR: The $0 script can't find the cache size info"
exit 1
