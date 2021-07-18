#!/bin/bash
# [UTILITY SCRIPT]
#
FRAGMENTS_LIST_FILE=$1
ALIGNMENTS_DIR=$2


N_POSITIVE=0
N_TOTAL=0
while IFS= read -r INPUT_FILE; do
	BASE_NAME=$(basename ${INPUT_FILE} .txt)
	N_ALIGNMENTS=$(wc -l < ${ALIGNMENTS_DIR}/test-basin-${BASE_NAME}/LAshow.txt)
	if [ $(( ${N_ALIGNMENTS} - 2 )) -lt 4  ]; then
		echo "${BASE_NAME} has just $(( ${N_ALIGNMENTS} - 2 )) alignments"
	else
		N_POSITIVE=$(( ${N_POSITIVE} + 1 ))
	fi
	N_TOTAL=$(( ${N_TOTAL} + 1 ))
done < ${FRAGMENTS_LIST_FILE}
echo "${N_POSITIVE} modules out of ${N_TOTAL} have some fragment-fragment alignment"