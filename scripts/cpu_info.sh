#!/bin/bash

# CPUID basic information highest value
# This seems to be safe for now ...
if [ -x "$(dirname $0)/cpuid" ]; then
	hv=$( echo 0 | $(dirname $0)/cpuid | sed -r -e 's/.*eax:([[:alnum:]]+) .*/\1/' )
	if [ "$hv" == "0xa" ]; then
		echo "CORE"
		exit 0
	fi
fi

cb=$(egrep -o '(Xeon|AMD|IA-64|Niagara)' </proc/cpuinfo \
     | head -1 \
     | sed -e '
           y/abcdefghijklmnopqrstuvwxyz/ABCDEFGHIJKLMNOPQRSTUVWXYZ/;
	   s/-/_/')

if [ ! -z "$cb" ]; then
	echo $cb
	exit 0
fi

if [ ! -f $(dirname $0)/cpuid ]; then
	echo "===>WARN: Some bozo (probably you) forgot to compile cpuid" >&2
	echo "===>WARN: Compiling ..." >&2
	make -C $(dirname $0) > /dev/null
fi

if [ -z "$cb" ]; then
	echo "===>ERROR: The $0 script can't find the cpu type"
	exit 1
fi

echo $cb
