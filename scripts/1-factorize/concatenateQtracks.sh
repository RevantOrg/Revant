#!/bin/bash
#
# [UTILITY SCRIPT] Given a directory that contains the quality tracks (text files) of
# several blocks of a database, the script concatenates all the tracks into a single text
# file.
#
# --------------------------------- CONFIGURATION ----------------------------------------
# This is the only section of the script that needs to be customized.
#
INPUT_DIR=$1
FIRST_FILE=$2
LAST_FILE=$3
INPUT_PREFIX="qv."
INPUT_SUFFIX=".txt"
OUTPUT_FILE="${INPUT_DIR}/DAMAR_Qtrack.txt"
# ----------------------------------------------------------------------------------------




/bin/rm -f ${OUTPUT_FILE}
for FILE in $(seq ${FIRST_FILE} ${LAST_FILE})
do 
   /bin/cat ${INPUT_DIR}/${INPUT_PREFIX}${FILE}${INPUT_SUFFIX} >> ${OUTPUT_FILE}
done
