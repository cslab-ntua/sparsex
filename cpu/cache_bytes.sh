#!/bin/sh
cs=$(sed -r -n '
	s/^cache size[[:space:]]+: /(/
	s/ KB/*1024)/p
	s/ MB/*1024*1024)/p' < /proc/cpuinfo | head -1)

if [ -z "$cs" ]; then
	echo "===>ERROR: The $0 script can't find the cache size info"
	exit 1
fi

echo $cs
