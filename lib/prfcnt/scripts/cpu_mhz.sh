#!/bin/bash

mhz=$(sed -r -n 's/^cpu MHz[[:space:]]+: //p' < /proc/cpuinfo | head -1)
if [ -n "$mhz" ]; then
	echo $mhz
	exit 0
fi

f="/sys/devices/system/cpu/cpu0/clock_tick"
if [ -f $f ]; then
	echo "($(cat $f)/(1000*1000))" | bc -l
	exit 0
fi

# failed ...
exit 1
