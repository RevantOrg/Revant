#!/bin/bash
#
# Checks fragment-reference alignments for consistency, and builds data structures needed
# by the following script.
#
# Remark: the script assumes that environment variable $REVANT_BINARIES$ is already set
# to the directory that contains REVANT's binaries, i.e. to the directory that contains
# the $de/mpi_cbg/revant/...$ subpath.
#
# --------------------------------- CONFIGURATION ----------------------------------------
# This is the only part of the script that needs to be customized.
#
PROJECT_DIR=$1
MIN_ALIGNMENT_LENGTH="500"  # Value used in fragment-{fragment,reference} alignments.
SMALL_DIFFS_ARE_ERRORS="0"
# REVANT
JAVA_RUNTIME_FLAGS="-Xms2G -Xmx10G"
# ----------------------------------------------------------------------------------------




STEP1_DIR="${PROJECT_DIR}/step1"
STEP5_DIR="${STEP1_DIR}/finalOutput/step4/step5"
ALIGNMENTS_DIR="${STEP5_DIR}/fragments-strings-alignments"
FRAGMENTS_LIST="${ALIGNMENTS_DIR}/list-fragments.txt"
CODE_DIR="/Users/ramseysnow/Dropbox/dropbox private/collaborators/myers/NEW/sw"

while IFS= read -r INPUT_FILE; do
	BASE_NAME=$(basename ${INPUT_FILE} .txt)
	CURRENT_DIR="${ALIGNMENTS_DIR}/test-basin-${BASE_NAME}"
	java ${JAVA_RUNTIME_FLAGS} -classpath "${REVANT_BINARIES}" de.mpi_cbg.revant.factorize.DBdump2ReadLengths ${CURRENT_DIR}/output-DBdump.txt > ${CURRENT_DIR}/reads-lengths.txt
done < ${FRAGMENTS_LIST}
java ${JAVA_RUNTIME_FLAGS} -classpath "${REVANT_BINARIES}" de.mpi_cbg.revant.fragments.FragmentsStep1 ${STEP1_DIR} ${MIN_ALIGNMENT_LENGTH} ${SMALL_DIFFS_ARE_ERRORS}
