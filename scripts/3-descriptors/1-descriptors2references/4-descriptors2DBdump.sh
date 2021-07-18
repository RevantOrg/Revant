#!/bin/bash
#
# The script builds a list of all the distinct read IDs that appear in a descriptor file.
#
# --------------------------------- CONFIGURATION ----------------------------------------
# This is the only section of the script that needs to be customized.
#
PROJECT_DIR=$1
# ----------------------------------------------------------------------------------------




STEP1_DIR="${PROJECT_DIR}/step1"
STEP4_DIR="${STEP1_DIR}/finalOutput/step4"
STEP5_DIR="${STEP4_DIR}/step5"
DESCRIPTORS_DIR=${STEP4_DIR}
DESCRIPTORS_LIST="${STEP5_DIR}/descriptors.txt"
READS_OUTPUT="${STEP5_DIR}/reads.txt"
TMP_FILE="${STEP5_DIR}/tmp.txt"

rm -rf ${READS_OUTPUT}
rm -f ${TMP_FILE}
while IFS= read -r DESCRIPTOR; do
	READS=$(cut -d ',' -f 1 ${DESCRIPTORS_DIR}/${DESCRIPTOR})
	for R in $(echo ${READS}); do
		echo $(( ${R} + 1 )) >> ${TMP_FILE}  # +1 since our readIDs start from zero
	done
done < ${DESCRIPTORS_LIST}
sort ${TMP_FILE} | uniq > ${READS_OUTPUT}
rm -f ${TMP_FILE}
