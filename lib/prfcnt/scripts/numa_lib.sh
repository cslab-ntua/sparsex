#!/bin/sh
echo "#include <numa.h>" | gcc -E - 2>/dev/null > /dev/null
if [ "$?" -eq "0" ]; then
	if [ -n "$1" ]; then
		echo "$1"
	fi
	exit 0
else
	exit 1
fi
