#!/bin/bash
#
# Assume that $DBdump -s$ has been executed on the $reads.txt$ file produced by the
# previous script. This script creates a single file per read sequence from the output of
# $DBdump -s$, naming the files with the read ID provided in an input list.
#
# --------------------------------- CONFIGURATION ----------------------------------------
# This is the only part of the script that needs to be customized.
#
PROJECT_DIR=$1
# ----------------------------------------------------------------------------------------




STEP1_DIR="${PROJECT_DIR}/step1"
STEP5_DIR="${STEP1_DIR}/finalOutput/step4/step5"
DBDUMP_FILE="${STEP5_DIR}/DBdump-output.txt"
READ_IDS_FILE="${STEP5_DIR}/reads.txt"
OUTPUT_DIR="${STEP5_DIR}/reads"
TMP_FILE="${STEP5_DIR}/tmp.txt"

rm -rf ${OUTPUT_DIR}
mkdir ${OUTPUT_DIR}
rm -f ${TMP_FILE}
tail -n +5 ${DBDUMP_FILE} | cut -d ' ' -f 3 | paste -d ',' ${READ_IDS_FILE} - > ${TMP_FILE}
awk -v outputDir="${OUTPUT_DIR}" 'BEGIN { FS="," } { outputFile=sprintf("%s/read%s.txt",outputDir,$1); print $2 > outputFile; close(outputFile); }' ${TMP_FILE}
rm -f ${TMP_FILE}
