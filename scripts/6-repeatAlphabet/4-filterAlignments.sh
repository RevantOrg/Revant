#!/bin/bash
# 
# Prints a bitvector that tells which read-read alignments to remove because likely
# repeat-induced.
#
# Remark: the script assumes that environment variable $REVANT_BINARIES$ is already set
# to the directory that contains REVANT's binaries, i.e. to the directory that contains
# the $de/mpi_cbg/revant/...$ subpath.
#
# --------------------------------- CONFIGURATION ----------------------------------------
# This is the only section of the script that needs to be customized.
#
INPUT_DIR=$1
BROKEN_READS=$2  # 1=TRUE
ALIGNMENTS_FILE="${INPUT_DIR}/LAshow-reads-reads.txt"
MIN_ALIGNMENT_LENGTH="500"  # In read-read alignments
MAX_K_UNIQUE_INTERVALS="8"  # Same as in $3-getUniqueSubstrings.sh$
FILTERING_MODE="0"  # 0-loose, 1=tight, 2=tight with matching characters.
N_THREADS="4"
# REVANT
JAVA_RUNTIME_FLAGS="-Xms2G -Xmx10G"
# ----------------------------------------------------------------------------------------

READ_LENGTHS_FILE="${INPUT_DIR}/reads-lengths.txt"
READ_IDS_FILE="${INPUT_DIR}/reads-ids.txt"
N_READS=$(wc -l < ${READ_IDS_FILE})
TMPFILE_NAME="filterAlignments-tmp"
TMPFILE_PATH="${INPUT_DIR}/${TMPFILE_NAME}"

echo "Splitting the alignments file..."
rm -f ${TMPFILE_PATH}-1-*
if [ ${BROKEN_READS} -eq 1 ]; then
	# Reusing the chunks of the read-read alignments file that are already there (we
	# assume that they all have the header).
	for FILE in $(ls ${INPUT_DIR}/breakReads-tmp-2-*.txt ); do
		ID=$(basename ${FILE} .txt)
		ID=${ID#breakReads-tmp-2-}
		mv ${INPUT_DIR}/breakReads-tmp-2-${ID}.txt ${TMPFILE_PATH}-1-${ID}.txt
	done
else
	N_ALIGNMENTS=$(( $(wc -l < ${ALIGNMENTS_FILE}) - 2 ))
	LAST_READA_FILE="${INPUT_DIR}/LAshow-reads-reads-lastReadA.txt"
	java ${JAVA_RUNTIME_FLAGS} -classpath "${REVANT_BINARIES}" de.mpi_cbg.revant.factorize.SplitAlignments ${N_ALIGNMENTS} ${N_THREADS} ${ALIGNMENTS_FILE} ${TMPFILE_PATH}-1- ${LAST_READA_FILE}
fi
echo "Alignments filtered and split in ${N_THREADS} parts"

echo "Filtering alignments..."
READS_TRANSLATED_FILE="${INPUT_DIR}/reads-translated-disambiguated.txt"
READS_TRANSLATED_BOUNDARIES="${INPUT_DIR}/reads-translated-boundaries-new.txt"
FULLY_UNIQUE_FILE="${INPUT_DIR}/reads-fullyUnique-new.txt"
N_FULLY_UNIQUE=$(wc -l < ${FULLY_UNIQUE_FILE})
FULLY_CONTAINED_FILE="${INPUT_DIR}/reads-fullyContained-new.txt"
N_FULLY_CONTAINED=$(wc -l < ${FULLY_CONTAINED_FILE})
UNIQUE_INTERVALS_FILE="${INPUT_DIR}/unique-intervals-k1-${MAX_K}.txt"
ALPHABET_FILE="${INPUT_DIR}/alphabet-cleaned.txt"
rm -f ${TMPFILE_PATH}-2-*
rm -f ${TMPFILE_PATH}-3-*

function filterThread() {
	local ALIGNMENTS_FILE_ID=$1
	java ${JAVA_RUNTIME_FLAGS} -classpath "${REVANT_BINARIES}" de.mpi_cbg.revant.apps.FilterAlignments ${TMPFILE_PATH}-1-${ALIGNMENTS_FILE_ID}.txt ${N_READS} ${READ_LENGTHS_FILE} ${READ_IDS_FILE} ${READS_TRANSLATED_FILE} ${READS_TRANSLATED_BOUNDARIES} ${FULLY_UNIQUE_FILE} ${N_FULLY_UNIQUE} ${FULLY_CONTAINED_FILE} ${N_FULLY_CONTAINED} ${INPUT_DIR}/unique-intervals-k1-${MAX_K_UNIQUE_INTERVALS}.txt ${FILTERING_MODE} ${ALPHABET_FILE} ${TMPFILE_PATH}-2-${ALIGNMENTS_FILE_ID} ${MIN_ALIGNMENT_LENGTH}
	if [ ${BROKEN_READS} -eq 1 ]; then
		NEW2OLD_FILE="${INPUT_DIR}/broken2unbroken.txt"
		OLD2NEW_FILE="${INPUT_DIR}/unbroken2broken.txt"
		NREADS_OLD=$(wc -l < "${INPUT_DIR}/reads-ids-unbroken.txt")
		READ_LENGTHS_FILE_OLD="${INPUT_DIR}/reads-lengths-unbroken.txt"
		ALIGNMENTS_FILE_OLD="${INPUT_DIR}/breakReads-tmp-1-${ALIGNMENTS_FILE_ID}.txt"
		java ${JAVA_RUNTIME_FLAGS} -classpath "${REVANT_BINARIES}" de.mpi_cbg.revant.apps.BreakReads5 ${TMPFILE_PATH}-2-${ALIGNMENTS_FILE_ID} ${TMPFILE_PATH}-1-${ALIGNMENTS_FILE_ID}.txt ${NEW2OLD_FILE} ${N_READS} ${OLD2NEW_FILE} ${NREADS_OLD} ${ALIGNMENTS_FILE_OLD} ${READ_LENGTHS_FILE_OLD} ${TMPFILE_PATH}-3-${ALIGNMENTS_FILE_ID}
	fi
}

if [ -e ${TMPFILE_PATH}-1-${N_THREADS}.txt ]; then
	TO=${N_THREADS}
else
	TO=$(( ${N_THREADS} - 1 ))
fi
for THREAD in $(seq 0 ${TO}); do
	filterThread ${THREAD} &
done
wait
echo "Alignments filtered successfully"
OUTPUT_BITVECTOR="${ALIGNMENTS_FILE}.mode${FILTERING_MODE}.bitvector"
rm -f ${OUTPUT_BITVECTOR}
for THREAD in $(seq 0 ${TO}); do
	if [ ${BROKEN_READS} -eq 1 ]; then
		cat ${TMPFILE_PATH}-3-${THREAD} >> ${OUTPUT_BITVECTOR}
	else
		cat ${TMPFILE_PATH}-2-${THREAD} >> ${OUTPUT_BITVECTOR}
	fi
done
