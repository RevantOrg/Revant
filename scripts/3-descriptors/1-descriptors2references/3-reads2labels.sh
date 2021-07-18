#!/bin/bash
#
# Given a set of read sequences, and a list of label coordinates files that should be
# converted into strings, the script extracts the corresponding substrings from reads.
#
# Remark: after the program ends, $READS_DIR$ can be deleted.
#
# --------------------------------- CONFIGURATION ----------------------------------------
# This is the only section of the script that needs to be customized.
#
PROJECT_DIR=$1
# ----------------------------------------------------------------------------------------




STEP1_DIR="${PROJECT_DIR}/step1"
STEP5_DIR="${STEP1_DIR}/finalOutput/step4/step5"
JOBS_FILE="${STEP5_DIR}/jobs.txt"
READS_DIR="${STEP5_DIR}/reads"

while IFS= read -r FILE; do
	FILE_LENGTH=${#FILE}
	PREFIX=${FILE%-labelCoordinates-cyclic.txt}
	PREFIX_LENGTH=${#PREFIX}
	if [ ${PREFIX_LENGTH} -lt ${FILE_LENGTH} ]; then
		OUT_FILE="${PREFIX}-labels-cyclic.txt"
	else
		PREFIX=${FILE%-labelCoordinates.txt}
		OUT_FILE="${PREFIX}-labels.txt"
	fi
	rm -f ${OUT_FILE}
	awk -v outputFile="${OUT_FILE}" -v readsDir="${READS_DIR}" 'BEGIN { FS="," } { \
		read=$1+1; start=$2+1; end=$3+1;  \
		inputFile=sprintf("%s/read%d.txt",readsDir,read); \
		getline readSequence < inputFile; \
		close(inputFile); \
		print substr(readSequence,start,end-start+1) >> outputFile;
 	}' ${FILE}
done < ${JOBS_FILE}
