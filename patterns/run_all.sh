#!/bin/bash

spmv_prog="./spmv_clone"
prog_out="clone_wsort.out"
mt_conf="0"

matrices=$*

if [ -z "$matrices" ]; then
    echo "$0: too few arguments" 1>&2
    echo "Usage: $0 [MATRIX]..."
fi

tmp_out=$(mktemp)
for m in $matrices; do
    for c in $mt_conf; do
		echo "MT_CONF=$c" >> $prog_out
        WINDOW_SIZE=255 XFORM_CONF=$(seq -s, 1 22) MT_CONF=$c $spmv_prog $m &> $tmp_out
		cat $tmp_out >> $prog_out
    done 
done
/bin/rm $tmp_out
