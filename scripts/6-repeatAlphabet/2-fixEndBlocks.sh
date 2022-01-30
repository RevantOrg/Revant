#!/bin/bash
# 
# For increasing values of k, the script tries to disambiguate the endpoints of reads
# using a context of length k.
#
# Remark: since collecting k-mers is expensive, one might want to try just few values of
# k.
#
# Remark: the script assumes that environment variable $REVANT_BINARIES$ is already set
# to the directory that contains REVANT's binaries, i.e. to the directory that contains
# the $de/mpi_cbg/revant/...$ subpath.
#
# --------------------------------- CONFIGURATION ----------------------------------------
# This is the only section of the script that needs to be customized.
#
INPUT_DIR=$1
HAPLOTYPE_COVERAGE="5"  # Of one haplotype
TIGHT_MODE="0"
MIN_K="2"  # One plus the min length of a context used for disambiguation
MAX_K="5"  # One plus the max length of a context used for disambiguation
N_THREADS="4"
# REVANT
JAVA_RUNTIME_FLAGS="-Xms2G -Xmx10G"
# ----------------------------------------------------------------------------------------


TMPFILE_NAME="fixEndBlocks-tmp"
TMPFILE_PATH="${INPUT_DIR}/${TMPFILE_NAME}"
READ_IDS_FILE="${INPUT_DIR}/reads-ids.txt"
N_READS=$(wc -l < ${READ_IDS_FILE})
READS_TRANSLATED_FILE="${INPUT_DIR}/reads-translated-new.txt"
READS_DISAMBIGUATED_FILE="${INPUT_DIR}/reads-translated-disambiguated.txt"
ALPHABET_FILE="${INPUT_DIR}/alphabet-cleaned.txt"
MIN_FREQUENCY_UNIQUE=${HAPLOTYPE_COVERAGE}
UNIQUE_MODE="1"; OPEN_MODE="1"; MULTI_MODE="0"  # Open blocks are forbidden
rm -f ${TMPFILE_PATH}*

function kmersThread() {
	local LOCAL_K=$1
	local LOCAL_TRANSLATED_READS_FILE=$2
	local LOCAL_KMERS_FILE=$3
	java ${JAVA_RUNTIME_FLAGS} -classpath "${REVANT_BINARIES}" de.mpi_cbg.revant.apps.CollectKmers ${LOCAL_K} ${LOCAL_TRANSLATED_READS_FILE} ${ALPHABET_FILE} ${UNIQUE_MODE} ${OPEN_MODE} ${MULTI_MODE} null ${LOCAL_KMERS_FILE}
}

function fixThread() {
	local LOCAL_OLD_TRANSLATED_FILE=$1
	local LOCAL_KMERS_FILE=$2
	local LOCAL_K=$3
	local LOCAL_NEW_TRANSLATED_FILE=$4
	java ${JAVA_RUNTIME_FLAGS} -classpath "${REVANT_BINARIES}" de.mpi_cbg.revant.apps.FixEndBlocks ${ALPHABET_FILE} ${LOCAL_OLD_TRANSLATED_FILE} ${LOCAL_KMERS_FILE} ${LOCAL_K} ${TIGHT_MODE} ${LOCAL_NEW_TRANSLATED_FILE}
}

rm -f "${TMPFILE_PATH}-1-*"
split -l $(( ${N_READS} / ${N_THREADS} )) ${READS_TRANSLATED_FILE} "${TMPFILE_PATH}-1-"
for K in $(seq ${MIN_K} ${MAX_K}); do
	SORT_OPTIONS_KMERS=""
	for i in $(seq 1 ${K}); do
		SORT_OPTIONS_KMERS="${SORT_OPTIONS_KMERS} -k ${i},${i}n"
	done
	echo "Collecting ${K}-mers..."
	for FILE in $(find -s ${INPUT_DIR} -name "${TMPFILE_NAME}-1-*"); do
		THREAD_ID=${FILE#${INPUT_DIR}/${TMPFILE_NAME}-1-}
		kmersThread ${K} ${FILE} ${TMPFILE_PATH}-kmers-${K}-${THREAD_ID} &
	done
	wait
	sort --parallel ${N_THREADS} -m -t , ${SORT_OPTIONS_KMERS} ${TMPFILE_PATH}-kmers-${K}-* > ${TMPFILE_PATH}-${K}.txt
	FREQUENT_KMERS_FILE="${INPUT_DIR}/frequent-k${K}.txt"
	echo "Finding frequent ${K}-mers..."
	java ${JAVA_RUNTIME_FLAGS} -classpath "${REVANT_BINARIES}" de.mpi_cbg.revant.apps.CompactKmers ${TMPFILE_PATH}-${K}.txt ${K} ${MIN_FREQUENCY_UNIQUE} -1 ${FREQUENT_KMERS_FILE} 0 null
	echo "Computing $((${K}-1))-mers..."
	K_MINUS_ONE_MERS_FILE="${INPUT_DIR}/kMinusOne-k${K}.txt"
	java ${JAVA_RUNTIME_FLAGS} -classpath "${REVANT_BINARIES}" de.mpi_cbg.revant.apps.GetKMinusOneMers ${ALPHABET_FILE} ${FREQUENT_KMERS_FILE} ${K} ${K_MINUS_ONE_MERS_FILE}
	echo "Disambiguating read ends using contexts of length $((${K}-1))..."
	for FILE in $(find -s ${INPUT_DIR} -name "${TMPFILE_NAME}-$((${K}-1))-*"); do
		THREAD_ID=${FILE#${INPUT_DIR}/${TMPFILE_NAME}-$((${K}-1))-}
		fixThread ${FILE} ${K_MINUS_ONE_MERS_FILE} $((${K}-1)) ${TMPFILE_PATH}-${K}-${THREAD_ID} > ${TMPFILE_PATH}-counts-${K}-${THREAD_ID} &
	done
	wait
	N_FIXED="0"; N_FIXABLE="0"; N_ENDS="0";
	for FILE in $(find -s ${INPUT_DIR} -name "${TMPFILE_NAME}-counts-${K}-*"); do
		N_FIXED=$(( ${N_FIXED} + $(cut -d , -f 1 ${FILE}) ))
		N_FIXABLE=$(( ${N_FIXABLE} + $(cut -d , -f 2 ${FILE}) ))
		N_ENDS=$(( ${N_ENDS} + $(cut -d , -f 3 ${FILE}) ))
	done
	N_ENDS=$((${N_ENDS}*2))
	echo "Disambiguated ${N_FIXED} read ends out of ${N_FIXABLE} fixable (${N_ENDS} total)"
done
rm -f ${READS_DISAMBIGUATED_FILE}
for FILE in $(find -s ${INPUT_DIR} -name "${TMPFILE_NAME}-${MAX_K}-*" ); do
	cat ${FILE} >> ${READS_DISAMBIGUATED_FILE}
done
