#!/bin/bash

spmv_prog="./spmv"
prog_out="dunnington.out"
mt_conf="0,1,2 0,1,2,12,13,14 0,1,2,12,13,14,3,4,5,15,16,17 0,1,2,12,13,14,3,4,5,15,16,17,6,7,8,18,19,20"

matrices=$*

if [ -z "$matrices" ]; then
    echo "$0: too few arguments" 1>&2
    echo "Usage: $0 [MATRIX]..."
fi

for m in $matrices; do
    for c in $mt_conf; do
        echo XFORM=$(seq -s, 1 22) MT_CONF=$c $spmv_prog $m &>> $prog_out
    done
done
