#!/bin/bash
#
# Assume that we have a list of read substrings files produced by the previous script,
# and that we want to convert them into strings. The program builds the list of all the
# distinct read IDs that are needed for the conversion. Read IDs in the output file start
# from one.
#
# --------------------------------- CONFIGURATION ----------------------------------------
# This is the only part of the script that needs to be customized.
#
PROJECT_DIR=$1
# ----------------------------------------------------------------------------------------




STEP1_DIR="${PROJECT_DIR}/step1"
STEP5_DIR="${STEP1_DIR}/finalOutput/step4/step5"
FRAGMENTS_DIR="${STEP5_DIR}/fragments"
READS_OUTPUT="${STEP5_DIR}/reads.txt"
TMP_FILE="${STEP5_DIR}/tmp.txt"

rm -f ${READS_OUTPUT}
rm -f ${TMP_FILE}
for FILE in $(find ${FRAGMENTS_DIR} -maxdepth 1 -type f -name "fragments-*"); do
	awk 'BEGIN { FS = "," } { print $1+1 }' ${FILE} >> ${TMP_FILE}
done
sort ${TMP_FILE} | uniq > ${READS_OUTPUT}
rm -f ${TMP_FILE}
