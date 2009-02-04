#!/bin/zsh

if [ -z "$1" ]; then
	echo "Usage $0 <mmf_file>"
	exit 1
fi

logf="/tmp/csrdu_test.log"
scripts/mtx_cmp_idx.sh =(./spm_csrdu_test "$1") "$1"  > $logf

if [ ! $? -eq 0 ]; then
	echo "Failed (see diff in $logf)"
	exit 1
fi
