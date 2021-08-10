#!/bin/bash
#
# Given a fragment-fragment alignments graph, the script resets the reference strings
# to be the "centers" of the clusters in the graph (see $FragmentsStep2.java$ for
# details). Repeats with no fragments file are not processed.
#
# Remark: at the end of this step, some fragment files might be empty, since it might
# happen that no fragment satisfies some consistency criteria. This is not necessarily an
# error.
#
# --------------------------------- CONFIGURATION ----------------------------------------
# This is the only part of the script that needs to be customized.
#
PROJECT_DIR=$1
MIN_ALIGNMENT_LENGTH="1500"  # Value used in repeat inference
MIN_N_FRAGMENTS="10"
REVANT_BINARIES="/Users/ramseysnow/Dropbox/dropbox private/collaborators/myers/NEW/sw/bin"
JAVA_RUNTIME_FLAGS="-Xms2G -Xmx10G"
# ----------------------------------------------------------------------------------------




STEP1_DIR="${PROJECT_DIR}/step1"
STEP4_DIR="${STEP1_DIR}/finalOutput/step4"
STEP5_DIR="${STEP4_DIR}/step5"
ALIGNMENTS_DIR="${STEP5_DIR}/fragments-strings-alignments"
FRAGMENTS_LIST="${ALIGNMENTS_DIR}/list-fragments.txt"
FRAGMENT_STRINGS_DIR="${STEP5_DIR}/fragments-strings"
REFERENCES_STRINGS_DIR="${STEP5_DIR}/references-strings"
NEW_REFERENCES_DIR="${ALIGNMENTS_DIR}/references-strings-new"
NEW_FRAGMENTS_DIR="${ALIGNMENTS_DIR}/fragments-strings-new"
STATS_FILE="${ALIGNMENTS_DIR}/step1-stats.txt"

rm -rf ${NEW_REFERENCES_DIR} ${NEW_FRAGMENTS_DIR}
mkdir ${NEW_REFERENCES_DIR} ${NEW_FRAGMENTS_DIR}
while IFS= read -r INPUT_FILE; do
	BASE_NAME=$(basename ${INPUT_FILE} .txt)
	ID=${BASE_NAME#fragments-}
	CURRENT_DIR="${ALIGNMENTS_DIR}/test-basin-${BASE_NAME}"
	NREADS=$(wc -l < ${CURRENT_DIR}/reads-lengths.txt)
	NREADS=${NREADS##*( )}
	FRAGMENT_STRINGS="${FRAGMENT_STRINGS_DIR}/${BASE_NAME}.txt"
	REFERENCE_IS_GOOD=$(grep ^${ID}, ${STATS_FILE} | cut -d "," -f 2)
	if [ ${REFERENCE_IS_GOOD} -eq 0 ]; then
		echo "Running FragmentsStep2 on ${BASE_NAME}..."
		java ${JAVA_RUNTIME_FLAGS} -classpath "${REVANT_BINARIES}" de.mpi_cbg.revant.fragments.FragmentsStep2 ${NREADS} ${CURRENT_DIR} ${STEP4_DIR} ${ID} ${MIN_N_FRAGMENTS} ${FRAGMENT_STRINGS} ${MIN_ALIGNMENT_LENGTH} ${NEW_REFERENCES_DIR} ${NEW_FRAGMENTS_DIR}
		EXIT_STATUS=$?
		if [ ${EXIT_STATUS} -eq 255 ]; then
			echo "Keeping the existing reference of ${BASE_NAME}"
			java ${JAVA_RUNTIME_FLAGS} -classpath "${REVANT_BINARIES}" de.mpi_cbg.revant.fragments.FragmentsStep2Trivial ${NREADS} ${CURRENT_DIR} ${ID} ${FRAGMENT_STRINGS} ${NEW_FRAGMENTS_DIR}
			EXIT_STATUS=$?
			if [ ${EXIT_STATUS} -ne 0 ]; then
			    echo "The following command returned error ${EXIT_STATUS}:"
				echo "java FragmentsStep2Trivial ${NREADS} ${CURRENT_DIR} ${ID} ${FRAGMENT_STRINGS} ${NEW_FRAGMENTS_DIR}"
				exit ${EXIT_STATUS}
			fi
			cp ${REFERENCES_STRINGS_DIR}/reference-${ID}.txt ${NEW_REFERENCES_DIR}
			REFERENCE_LENGTH=$(head -n 1 ${REFERENCES_STRINGS_DIR}/reference-${ID}.txt | awk -F "_" '{ print $2 }')
			echo ${REFERENCE_LENGTH} > ${NEW_REFERENCES_DIR}/reference-${ID}-lengths.txt
		elif [ ${EXIT_STATUS} -ne 0 ]; then
		    echo "The following command returned error ${EXIT_STATUS}:"
			echo "java FragmentsStep2 ${NREADS} ${CURRENT_DIR} ${STEP4_DIR} ${ID} ${MIN_N_FRAGMENTS} ${FRAGMENT_STRINGS} ${MIN_ALIGNMENT_LENGTH} ${NEW_REFERENCES_DIR} ${NEW_FRAGMENTS_DIR}"
			exit ${EXIT_STATUS}
		fi
	else 
		echo "Keeping the existing reference of ${BASE_NAME}"
		java ${JAVA_RUNTIME_FLAGS} -classpath "${REVANT_BINARIES}" de.mpi_cbg.revant.fragments.FragmentsStep2Trivial ${NREADS} ${CURRENT_DIR} ${ID} ${FRAGMENT_STRINGS} ${NEW_FRAGMENTS_DIR}
		EXIT_STATUS=$?
		if [ ${EXIT_STATUS} -ne 0 ]; then
		    echo "The following command returned error ${EXIT_STATUS}:"
			echo "java FragmentsStep2Trivial ${NREADS} ${CURRENT_DIR} ${ID} ${FRAGMENT_STRINGS} ${NEW_FRAGMENTS_DIR}"
			exit ${EXIT_STATUS}
		fi
		cp ${REFERENCES_STRINGS_DIR}/reference-${ID}.txt ${NEW_REFERENCES_DIR}
		REFERENCE_LENGTH=$(head -n 1 ${REFERENCES_STRINGS_DIR}/reference-${ID}.txt | awk -F "_" '{ print $2 }')
		echo ${REFERENCE_LENGTH} > ${NEW_REFERENCES_DIR}/reference-${ID}-lengths.txt
	fi
done < ${FRAGMENTS_LIST}