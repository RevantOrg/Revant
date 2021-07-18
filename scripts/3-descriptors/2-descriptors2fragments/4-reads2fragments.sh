#!/bin/bash
#
# Given a set of read sequences, and a set of read substrings files, the program extracts
# the corresponding substrings from the reads.
#
# --------------------------------- CONFIGURATION ----------------------------------------
# This is the only part of the script that needs to be customized.
#
PROJECT_DIR=$1
# ----------------------------------------------------------------------------------------




STEP1_DIR="${PROJECT_DIR}/step1"
STEP5_DIR="${STEP1_DIR}/finalOutput/step4/step5"
FRAGMENTS_DIR="${STEP5_DIR}/fragments"
READS_DIR="${STEP5_DIR}/reads"
OUTPUT_DIR="${STEP5_DIR}/fragments-strings"

if [ ! -d ${OUTPUT_DIR} ]; then
	mkdir ${OUTPUT_DIR}
fi
for FILE in $(find ${FRAGMENTS_DIR} -maxdepth 1 -type f -name "*"); do
	ID=$(( 2**14 ))  # Arbitrary, should be bigger than the number of references.
	awk -v id="${ID}" -v outputFile="${OUTPUT_DIR}/$(basename ${FILE})" -v readsDir="${READS_DIR}" 'BEGIN { FS="," } { \
		id=id+1; read=$1+1; start=$2+1; end=$3+1; \
		header=sprintf(">U0/%s/0_%d",id,end-start+1); \
		print header >> outputFile; \
		inputFile=sprintf("%s/read%d.txt",readsDir,read); \
		getline readSequence < inputFile; \
		close(inputFile); \
		print substr(readSequence,start,end-start+1) >> outputFile; \
		close(outputFile); \
	}' ${FILE}
done