#!/bin/bash
# 
# For increasing values of k, the script computes k-mers that: (1) likely occur just once
# in the genome; (2) do not contain any other unique h-mer. All the occurrences of unique
# k-mers in every read are marked, along with an estimate of the number of haplotypes they
# belong to.
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
HAPLOTYPE_COVERAGE="5"  # Of one haplotype
MAX_K="10"  # Stops finding unique k-mers after this length. Should be set using the histogram of recoded lengths.
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
ALPHABET_FILE="${INPUT_DIR}/alphabet-cleaned.txt"
MIN_FREQUENCY_UNIQUE=${HAPLOTYPE_COVERAGE}
MAX_FREQUENCY_UNIQUE=$(( ${HAPLOTYPE_COVERAGE}*${N_HAPLOTYPES} + (${HAPLOTYPE_COVERAGE}-1) ))
# We use open blocks if they match just one character, since they are the only way to
# detect e.g. a transposon that is longer than every read and that occurs just once in
# the genome, or an extremely long satellite that occurs just once in the genome.
UNIQUE_MODE="1"; OPEN_MODE="0"; MULTI_MODE="1"
MAX_HISTOGRAM_COUNT="10000"  # Arbitrary
rm -f ${TMPFILE_PATH}*
rm -f ${INPUT_DIR}/unique-*
rm -f ${INPUT_DIR}/histogram-*

function kmersThread() {
	local LOCAL_K=$1
	local LOCAL_TRANSLATED_READS_FILE=$2
	local LOCAL_K_MINUS_ONE_INTERVALS_FILE=$3
	local LOCAL_KMERS_FILE=$4
	java ${JAVA_RUNTIME_FLAGS} -classpath "${REVANT_BINARIES}" de.mpi_cbg.revant.apps.CollectKmers ${LOCAL_K} ${LOCAL_TRANSLATED_READS_FILE} ${ALPHABET_FILE} ${UNIQUE_MODE} ${OPEN_MODE} ${MULTI_MODE} ${LOCAL_K_MINUS_ONE_INTERVALS_FILE} ${LOCAL_KMERS_FILE}
}

function intervalsThread() {
	local LOCAL_K=$1
	local LOCAL_TRANSLATED_READS_FILE=$2
	local LOCAL_UNIQUE_KMERS_FILE=$3
	local LOCAL_K_MINUS_ONE_INTERVALS_FILE=$4
	local LOCAL_INTERVALS_FILE=$5
	java ${JAVA_RUNTIME_FLAGS} -classpath "${REVANT_BINARIES}" de.mpi_cbg.revant.apps.GetShortestUniqueIntervals ${LOCAL_K} ${LOCAL_TRANSLATED_READS_FILE} ${ALPHABET_FILE} ${UNIQUE_MODE} ${OPEN_MODE} ${MULTI_MODE} ${LOCAL_UNIQUE_KMERS_FILE} ${HAPLOTYPE_COVERAGE} ${LOCAL_K_MINUS_ONE_INTERVALS_FILE} ${LOCAL_INTERVALS_FILE}
}

split -l $(( ${N_READS} / ${N_THREADS} )) ${READS_TRANSLATED_FILE} "${TMPFILE_PATH}-0-"
for K in $(seq 1 ${MAX_K}); do
	SORT_OPTIONS_KMERS=""
	for i in $(seq 1 ${K}); do
		SORT_OPTIONS_KMERS="${SORT_OPTIONS_KMERS} -k ${i},${i}n"
	done
	echo "Collecting ${K}-mers..."
	for FILE in $(find -s ${INPUT_DIR} -name "${TMPFILE_NAME}-0-*" ); do
		THREAD_ID=${FILE#${INPUT_DIR}/${TMPFILE_NAME}-0-}
		if [ ${K} -le 1 ]; then
			PREVIOUS_INTERVALS="null"
		else
			PREVIOUS_INTERVALS="${TMPFILE_PATH}-$((${K}-1))-intervals-${THREAD_ID}"
		fi
		kmersThread ${K} ${FILE} ${PREVIOUS_INTERVALS} ${TMPFILE_PATH}-${K}-kmers-${THREAD_ID} &
	done
	wait
	sort --parallel ${N_THREADS} -m -t , ${SORT_OPTIONS_KMERS} ${TMPFILE_PATH}-${K}-kmers-* > ${TMPFILE_PATH}-${K}.txt
	UNIQUE_KMERS_FILE="${INPUT_DIR}/unique-k${K}.txt"
	OUTPUT_FILE_HISTOGRAM="${INPUT_DIR}/histogram-k${K}.txt"
	echo "Finding unique ${K}-mers..."
	java ${JAVA_RUNTIME_FLAGS} -classpath "${REVANT_BINARIES}" de.mpi_cbg.revant.apps.CompactKmers ${TMPFILE_PATH}-${K}.txt ${K} ${MIN_FREQUENCY_UNIQUE} ${MAX_FREQUENCY_UNIQUE} ${UNIQUE_KMERS_FILE} ${MAX_HISTOGRAM_COUNT} ${OUTPUT_FILE_HISTOGRAM}
	echo "Updating shortest unique intervals file..."
	for FILE in $(find -s ${INPUT_DIR} -name "${TMPFILE_NAME}-0-*" ); do
		THREAD_ID=${FILE#${INPUT_DIR}/${TMPFILE_NAME}-0-}
		if [ ${K} -le 1 ]; then
			PREVIOUS_INTERVALS="null"
		else
			PREVIOUS_INTERVALS="${TMPFILE_PATH}-$((${K}-1))-intervals-${THREAD_ID}"
		fi
		intervalsThread ${K} ${FILE} ${UNIQUE_KMERS_FILE} ${PREVIOUS_INTERVALS} ${TMPFILE_PATH}-${K}-intervals-${THREAD_ID} &
	done
	wait
done
FINAL_INTERVALS_FILE="${INPUT_DIR}/unique-intervals-k1-${MAX_K}.txt"
rm -f ${FINAL_INTERVALS_FILE}
for FILE in $(find -s ${INPUT_DIR} -name "${TMPFILE_PATH}-${MAX_K}-intervals-*" ); do
	cat ${FILE} >> ${FINAL_INTERVALS_FILE}
done
