#!/bin/zsh

setopt extended_glob

if [ -z "$1" ]; then
	mtxfiles=$(echo tests.mtx/* /s/matrices/sorted/(???~(024|041|042)).*.mtx.sorted)
else
	mtxfiles="$*"
fi
for mtx in $(echo $mtxfiles)
do
	echo CSRDU_VERBOSE=0 CSRDU_DE_MINLEN=4 CSRDU_JMP=1 scripts/csrdu_test.sh $mtx
	     CSRDU_VERBOSE=0 CSRDU_DE_MINLEN=4 CSRDU_JMP=1 scripts/csrdu_test.sh $mtx
	if [ ! $? -eq 0 ]; then
		exit
	fi

	echo CSRDU_VERBOSE=0 CSRDU_DE_MINLEN=4 CSRDU_JMP=0 scripts/csrdu_test.sh $mtx
	     CSRDU_VERBOSE=0 CSRDU_DE_MINLEN=4 CSRDU_JMP=0 scripts/csrdu_test.sh $mtx
	if [ ! $? -eq 0 ]; then
		exit
	fi

	echo CSRDU_VERBOSE=0 CSRDU_JMP=1 scripts/csrdu_test.sh $mtx
	     CSRDU_VERBOSE=0 CSRDU_JMP=1 scripts/csrdu_test.sh $mtx
	if [ ! $? -eq 0 ]; then
		exit
	fi

	echo CSRDU_VERBOSE=0 CSRDU_JMP=0 scripts/csrdu_test.sh $mtx
	     CSRDU_VERBOSE=0 CSRDU_JMP=1 scripts/csrdu_test.sh $mtx
	if [ ! $? -eq 0 ]; then
		exit
	fi

	echo
done
