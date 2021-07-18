#!/bin/bash
# [UTILITY SCRIPT]
#
DESCRIPTORS_DIR=$1
TMP_FILE="./tmp.txt"

rm -f ${TMP_FILE}
for FILE in $(ls ${DESCRIPTORS_DIR}/basin-*.txt); do 
	HEADER=$(head -n 1 ${FILE})
	BASENAME=$(basename ${FILE})
	echo "${BASENAME},${HEADER}" >> ${TMP_FILE}
done
sort -t , -k 3 -n ${TMP_FILE}
