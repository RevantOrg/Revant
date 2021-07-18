#!/bin/bash
# [UTILITY SCRIPT]
#
DESCRIPTORS_LIST=$1
DESCRIPTORS_DIR=$2
LONG_PERIOD="1000"

for DESCRIPTOR in $(cat ${DESCRIPTORS_LIST}); do
	HEADER=$(head -n 1 ${DESCRIPTORS_DIR}/${DESCRIPTOR})
	TYPE=$(echo "${HEADER}" | cut -d , -f 4)
	if [ ${TYPE} -eq 6 ]; then
		PERIOD=$(echo "${HEADER}" | cut -d , -f 5)
		if [ ${PERIOD} -gt ${LONG_PERIOD} ]; then
			echo "${DESCRIPTOR},LONG"
		else 
			echo "${DESCRIPTOR},short"
		fi
	fi
done