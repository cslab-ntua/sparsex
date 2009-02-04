#!/bin/bash
for mtx in tests.mtx/* /s/matrices/sorted/???.*.mtx.sorted
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

done
