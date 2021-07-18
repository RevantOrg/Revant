#!/bin/bash
#
# Assume that we have a list of kernel descriptor files, and that we want to convert the
# reference sequence of every descriptor into strings. The script builds a list of all
# the label coordinates files that should be converted into strings, and a list of all
# distinct read IDs that are needed for the conversion (see 
# $IntervalGraphStep3.printBidirectedGraphLabels()$).
#
# --------------------------------- CONFIGURATION ----------------------------------------
# This is the only section of the script that needs to be customized.
#
PROJECT_DIR=$1
# ----------------------------------------------------------------------------------------



STEP1_DIR="${PROJECT_DIR}/step1"
STEP5_DIR="${STEP1_DIR}/finalOutput/step4/step5"
DESCRIPTORS_LIST="${STEP5_DIR}/descriptors.txt"
JOBS_OUTPUT="${STEP5_DIR}/jobs.txt"
READS_OUTPUT="${STEP5_DIR}/reads.txt"
TMP_FILE="${STEP5_DIR}/tmp.txt"

# Computing the necessary $(component,cluster)$ pairs.
rm -f ${TMP_FILE}
cut -d '-' -f 2,3 ${DESCRIPTORS_LIST} | sort | uniq > ${TMP_FILE}

# Listing all label coordinates files of the selected pairs
rm -f ${JOBS_OUTPUT}
while IFS= read -r DESCRIPTOR; do
	COMPONENT_ID=$(echo ${DESCRIPTOR} | cut -d '-' -f 1 -)
	CLUSTER_ID=$(echo ${DESCRIPTOR} | cut -d '-' -f 2 -)
	DIRECTORY=${STEP1_DIR}/${COMPONENT_ID}-clusters
	for FILE in $(find ${DIRECTORY} -maxdepth 1 -type f -name "${CLUSTER_ID}-kernel*-labelCoordinates.txt" 2> /dev/null); do
		echo ${FILE} >> ${JOBS_OUTPUT}
	done
	for FILE in $(find ${DIRECTORY} -maxdepth 1 -type f -name "${CLUSTER_ID}-kernel*-labelCoordinates-cyclic.txt" 2> /dev/null); do
		echo ${FILE} >> ${JOBS_OUTPUT}
	done
done < ${TMP_FILE}
rm -f ${TMP_FILE}

# Listing all the distinct read IDs needed by some label coordinate file
rm -f ${TMP_FILE}
while IFS= read -r FILE; do
	READS=$(cut -d ',' -f 1 ${FILE})
	for R in $(echo ${READS}); do
		echo $(( ${R} + 1 )) >> ${TMP_FILE}  # +1 since our readIDs start from zero
	done
done < ${JOBS_OUTPUT}
sort ${TMP_FILE} | uniq > ${READS_OUTPUT}
rm -f ${TMP_FILE}
