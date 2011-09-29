#!/bin/bash

#
# Query the cpu family /proc/cpuinfo
#       6: Intel Core
#      15: Intel Netburst
#
cb=$(cat /proc/cpuinfo | grep 'cpu family' | head -n1 | awk -F': ' '{print $2}')
if [ "$cb" = "6" ]; then
	echo "CORE"
	exit 0
elif [ "$cb" = "15" ]; then
	echo "XEON"
	exit 0
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
