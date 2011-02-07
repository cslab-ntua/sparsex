#!/bin/bash

# find CPU type for prfcnt
# --kkourt

cpuid="$(dirname $0)/../cpuid"

# compile cpuid
if [ ! -x "$cpuid"  ]; then
	echo "===>WARN: cpuid is not compiled" >&2
	echo "===>WARN: Compiling ..." >&2
        # set cpu, so that Makefile doesn't call this script again
	CPU="DUMMY" make -C $(dirname $cpuid) cpuid > /dev/null
fi


# CPUID basic information highest value
# This seems to be safe for now ...
if [ -x "$cpuid" ]; then
	hv=$( echo 0 | $cpuid | sed -r -e 's/.*eax:([[:alnum:]]+) .*/\1/' )
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

if [ -z "$cb" ]; then
	echo "===>ERROR: The $0 script can't find the cpu type" >&2
	exit 1
fi

echo $cb
