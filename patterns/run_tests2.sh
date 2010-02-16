#!/bin/bash

NR_ROW_TYPES=4
NR_COL_TYPES=7
ROW_TYPE_BASE=10
COL_TYPE_BASE=16

SPMV_EXECUTABLE="./spmv"

for ((i = 0; i < NR_ROW_TYPES; i++)); do
    ((xform = ROW_TYPE_BASE + i))
    env MT_CONF=0 XFORM_CONF=$xform $SPMV_EXECUTABLE $* 2>> err.txt
done

for ((i = 0; i < NR_COL_TYPES; i++)); do
    ((xform = COL_TYPE_BASE + i))
    env MT_CONF=0 XFORM_CONF=$xform $SPMV_EXECUTABLE $* 2>> err.txt
done
