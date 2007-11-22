#!/bin/sh
sed -r -n 's/^cpu MHz[[:space:]]+: //p' < /proc/cpuinfo | head -1
