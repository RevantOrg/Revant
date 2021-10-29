#!/bin/bash
# 
# ------->
#
# Remark: the script assumes that environment variable $REVANT_BINARIES$ is already set
# to the directory that contains REVANT's binaries, i.e. to the directory that contains
# the $de/mpi_cbg/revant/...$ subpath.
#
# --------------------------------- CONFIGURATION ----------------------------------------
# This is the only section of the script that needs to be customized.
#
INPUT_DIR=$1
N_HAPLOTYPES="2"
COVERAGE="5"  # Of one haplotype
MAX_K="10"
N_THREADS="4"
# REVANT
JAVA_RUNTIME_FLAGS="-Xms2G -Xmx10G"
# ----------------------------------------------------------------------------------------


READ_LENGTHS_FILE="${INPUT_DIR}/reads-lengths.txt"
READ_IDS_FILE="${INPUT_DIR}/reads-ids.txt"
N_READS=$(wc -l < ${READ_IDS_FILE})
TMPFILE_NAME="getUniqueSubstrings-tmp"
TMPFILE_PATH="${INPUT_DIR}/${TMPFILE_NAME}"
READS_TRANSLATED_FILE="${INPUT_DIR}/reads-translated-new.txt"
READS_TRANSLATED_BOUNDARIES="${INPUT_DIR}/reads-translated-boundaries-new.txt"
ALPHABET_FILE="${INPUT_DIR}/alphabet.txt"
MIN_FREQUENCY_UNIQUE=${COVERAGE}
MAX_FREQUENCY_UNIQUE=$(( ${COVERAGE}*${N_HAPLOTYPES} + (${COVERAGE}-1) ))
# We use open blocks if they match just one character, since they are the only way to
# detect e.g. a transposon that is longer than every read and that occurs just once in
# the genome, or an extremely long satellite that occurs just once in the genome.
UNIQUE_MODE="1"; OPEN_MODE="0"; MULTI_MODE="1"
MAX_HISTOGRAM_COUNT="10000"  # Arbitrary
rm -rf ${TMPFILE_PATH}*


echo "Collecting ${K}-mers in parallel..."
function collectionThread() {
	LOCAL_K=$1
	TRANSLATED_READS_FILE=$2
	OUTPUT_FILE=$3
	java ${JAVA_RUNTIME_FLAGS} -classpath "${REVANT_BINARIES}" de.mpi_cbg.revant.apps.CollectKmers ${LOCAL_K} ${TRANSLATED_READS_FILE} ${ALPHABET_FILE} ${UNIQUE_MODE} ${OPEN_MODE} ${MULTI_MODE} ${OUTPUT_FILE}
}
split -l $(( ${N_READS} / ${N_THREADS} )) ${READS_TRANSLATED_FILE} "${TMPFILE_PATH}-0-"
for K in $(seq 1 ${MAX_K}); do
	SORT_OPTIONS_KMERS=""
	for i in $(seq 1 ${K}); do
		SORT_OPTIONS_KMERS="${SORT_OPTIONS_KMERS} -k ${i},${i}n"
	done
	for FILE in $(find -s ${INPUT_DIR} -name "${TMPFILE_NAME}-0-*" ); do
		THREAD_ID=${FILE#${INPUT_DIR}/${TMPFILE_NAME}-0-}
		collectionThread ${K} ${FILE} ${TMPFILE_PATH}-${K}-${THREAD_ID} &
	done
	wait
	sort --parallel ${N_THREADS} -m -t , ${SORT_OPTIONS_KMERS} ${TMPFILE_PATH}-${K}-* > ${TMPFILE_PATH}-${K}.txt
	OUTPUT_FILE_HISTOGRAM="${INPUT_DIR}/histogram-k${K}.txt"
	UNIQUE_KMERS_FILE="${INPUT_DIR}/unique-k${K}.txt"
	java ${JAVA_RUNTIME_FLAGS} -classpath "${REVANT_BINARIES}" de.mpi_cbg.revant.apps.CompactKmers ${TMPFILE_PATH}-${K}.txt ${K} ${MIN_FREQUENCY_UNIQUE} ${MAX_FREQUENCY_UNIQUE} ${UNIQUE_KMERS_FILE} ${MAX_HISTOGRAM_COUNT} ${OUTPUT_FILE_HISTOGRAM}
done