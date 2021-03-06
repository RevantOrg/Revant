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
N_HAPLOTYPES="1"
HAPLOTYPE_COVERAGE="30"  # Of one haplotype
MAX_K="8"  # Stops looking for unique k-mers after this length. Should be set using the histogram of recoded lengths.
UNIQUE_MODE="1"  # Non-repetitive blocks are allowed in a k-mer, except at the first/last
# position of the k-mer. Usually a good choice.
MULTI_MODE="0"; ONEMER_FILTER="3"
# Good settings for a mostly periodic genome: MULTI_MODE="0"; ONEMER_FILTER="3"
# Good settings for a mostly nonperiodic genome: MULTI_MODE="1"; ONEMER_FILTER="2"
N_THREADS="1"
DELETE_TMP_FILES="0"
# REVANT
JAVA_RUNTIME_FLAGS="-Xms2G -Xmx10G"
# ----------------------------------------------------------------------------------------

set -o pipefail; set -e; set -u
export LC_ALL=C  # To speed up the $sort$ command.
READ_LENGTHS_FILE="${INPUT_DIR}/reads-lengths.txt"
READ_IDS_FILE="${INPUT_DIR}/reads-ids.txt"
N_READS=$(wc -l < ${READ_IDS_FILE})
TMPFILE_NAME="getUniqueSubstrings-tmp"
TMPFILE_PATH="${INPUT_DIR}/${TMPFILE_NAME}"
READS_TRANSLATED_FILE="${INPUT_DIR}/reads-translated-disambiguated.txt"
READS_TRANSLATED_BOUNDARIES="${INPUT_DIR}/reads-translated-boundaries-new.txt"
ALPHABET_FILE="${INPUT_DIR}/alphabet-cleaned.txt"
TANDEMS_FILE="${INPUT_DIR}/tandems.txt"
REPEAT_LENGTHS_FILE="${INPUT_DIR}/repeats-lengths.txt"
N_REPEATS=$(wc -l < ${REPEAT_LENGTHS_FILE})
MIN_FREQUENCY_UNIQUE=$(( ${HAPLOTYPE_COVERAGE} / 2 ))
MAX_FREQUENCY_UNIQUE=$(( ${HAPLOTYPE_COVERAGE}*${N_HAPLOTYPES} + ${HAPLOTYPE_COVERAGE}/2 ))
MAX_HISTOGRAM_COUNT="10000"  # Arbitrary
rm -f ${TMPFILE_PATH}*
rm -f ${INPUT_DIR}/unique-*
rm -f ${INPUT_DIR}/histogram-*

function kmersThread() {
	local LOCAL_K=$1
	local LOCAL_TRANSLATED_READS_FILE=$2
	local LOCAL_BOUNDARIES_FILE=$3
	local LOCAL_READ_LENGTHS_FILE=$4
	local LOCAL_K_MINUS_ONE_INTERVALS_FILE=$5
	local LOCAL_KMERS_FILE=$6
	java ${JAVA_RUNTIME_FLAGS} -classpath "${REVANT_BINARIES}" de.mpi_cbg.revant.apps.CollectKmers ${LOCAL_K} ${LOCAL_TRANSLATED_READS_FILE} ${LOCAL_BOUNDARIES_FILE} ${LOCAL_READ_LENGTHS_FILE} ${ALPHABET_FILE} ${UNIQUE_MODE} ${MULTI_MODE} ${LOCAL_K_MINUS_ONE_INTERVALS_FILE} ${LOCAL_KMERS_FILE}
}

function intervalsThread() {
	local LOCAL_K=$1
	local LOCAL_TRANSLATED_READS_FILE=$2
	local LOCAL_BOUNDARIES_FILE=$3
	local LOCAL_READ_LENGTHS_FILE=$4
	local LOCAL_UNIQUE_KMERS_FILE=$5
	local LOCAL_K_MINUS_ONE_INTERVALS_FILE=$6
	local LOCAL_INTERVALS_FILE=$7
	java ${JAVA_RUNTIME_FLAGS} -classpath "${REVANT_BINARIES}" de.mpi_cbg.revant.apps.GetShortestUniqueIntervals ${LOCAL_K} ${LOCAL_TRANSLATED_READS_FILE} ${LOCAL_BOUNDARIES_FILE} ${LOCAL_READ_LENGTHS_FILE} ${ALPHABET_FILE} ${UNIQUE_MODE} ${MULTI_MODE} ${ONEMER_FILTER} ${LOCAL_UNIQUE_KMERS_FILE} ${HAPLOTYPE_COVERAGE} ${LOCAL_K_MINUS_ONE_INTERVALS_FILE} ${LOCAL_INTERVALS_FILE}
}

FINAL_INTERVALS_FILE="${INPUT_DIR}/unique-intervals-k1-${MAX_K}.txt"
split -l $(( ${N_READS} / ${N_THREADS} )) ${READS_TRANSLATED_FILE} "${TMPFILE_PATH}-0-"
split -l $(( ${N_READS} / ${N_THREADS} )) ${READS_TRANSLATED_BOUNDARIES} "${TMPFILE_PATH}-1-"
split -l $(( ${N_READS} / ${N_THREADS} )) ${READ_LENGTHS_FILE} "${TMPFILE_PATH}-2-"
for K in $(seq 1 ${MAX_K}); do
	SORT_OPTIONS_KMERS=""
	for i in $(seq 1 ${K}); do
		SORT_OPTIONS_KMERS="${SORT_OPTIONS_KMERS} -k ${i},${i}n"
	done
	echo "Collecting ${K}-mers..."
	for FILE in $(find -s ${INPUT_DIR} -name "${TMPFILE_NAME}-0-*"); do
		THREAD_ID=${FILE#${INPUT_DIR}/${TMPFILE_NAME}-0-}
		if [ ${K} -le 1 ]; then
			PREVIOUS_INTERVALS="null"
		else
			PREVIOUS_INTERVALS="${TMPFILE_PATH}-$((${K}-1))-intervals-${THREAD_ID}"
		fi
		kmersThread ${K} ${FILE} ${TMPFILE_PATH}-1-${THREAD_ID} ${TMPFILE_PATH}-2-${THREAD_ID} ${PREVIOUS_INTERVALS} ${TMPFILE_PATH}-${K}-kmers-${THREAD_ID} &
	done
	wait
	sort --parallel=${N_THREADS} -m -t , ${SORT_OPTIONS_KMERS} ${TMPFILE_PATH}-${K}-kmers-* > ${TMPFILE_PATH}-${K}.txt
	if [ ! -s ${TMPFILE_PATH}-${K}.txt ]; then
		MAX_K=$((${K}-1))
		break
	fi
	UNIQUE_KMERS_FILE="${INPUT_DIR}/unique-k${K}.txt"
	OUTPUT_FILE_HISTOGRAM="${INPUT_DIR}/histogram-k${K}.txt"
	echo "Finding unique ${K}-mers..."
	java ${JAVA_RUNTIME_FLAGS} -classpath "${REVANT_BINARIES}" de.mpi_cbg.revant.apps.CompactKmers ${TMPFILE_PATH}-${K}.txt ${K} ${MIN_FREQUENCY_UNIQUE} ${MAX_FREQUENCY_UNIQUE} ${UNIQUE_KMERS_FILE} ${MAX_HISTOGRAM_COUNT} 1 ${OUTPUT_FILE_HISTOGRAM}
	echo "Updating shortest unique intervals file..."
	for FILE in $(find -s ${INPUT_DIR} -name "${TMPFILE_NAME}-0-*"); do
		THREAD_ID=${FILE#${INPUT_DIR}/${TMPFILE_NAME}-0-}
		if [ ${K} -le 1 ]; then
			PREVIOUS_INTERVALS="null"
		else
			PREVIOUS_INTERVALS="${TMPFILE_PATH}-$((${K}-1))-intervals-${THREAD_ID}"
		fi
		intervalsThread ${K} ${FILE} ${TMPFILE_PATH}-1-${THREAD_ID} ${TMPFILE_PATH}-2-${THREAD_ID} ${UNIQUE_KMERS_FILE} ${PREVIOUS_INTERVALS} ${TMPFILE_PATH}-${K}-intervals-${THREAD_ID} &
	done
	wait
done
rm -f ${FINAL_INTERVALS_FILE}
for FILE in $(find -s ${INPUT_DIR} -name "${TMPFILE_NAME}-${MAX_K}-intervals-*" ); do
	cat ${FILE} >> ${FINAL_INTERVALS_FILE}
done

# Collecting tandem intervals
function tandemsThread() {
	local LOCAL_TRANSLATED_READS_FILE=$1
	local LOCAL_BOUNDARIES_FILE=$2
	local LOCAL_READ_LENGTHS_FILE=$3
	local LOCAL_TANDEMS_FILE=$4
	java ${JAVA_RUNTIME_FLAGS} -classpath "${REVANT_BINARIES}" de.mpi_cbg.revant.apps.CollectTandems ${ALPHABET_FILE} ${LOCAL_TRANSLATED_READS_FILE} ${LOCAL_BOUNDARIES_FILE} ${LOCAL_READ_LENGTHS_FILE} ${REPEAT_LENGTHS_FILE} ${N_REPEATS} ${LOCAL_TANDEMS_FILE}
}
echo "Collecting tandems..."
for FILE in $(find -s ${INPUT_DIR} -name "${TMPFILE_NAME}-0-*"); do
	THREAD_ID=${FILE#${INPUT_DIR}/${TMPFILE_NAME}-0-}
	tandemsThread ${FILE} ${TMPFILE_PATH}-1-${THREAD_ID} ${TMPFILE_PATH}-2-${THREAD_ID} ${TMPFILE_PATH}-tandems-${THREAD_ID} &
done
wait
rm -f ${TANDEMS_FILE}
for FILE in $(find -s ${INPUT_DIR} -name "${TMPFILE_NAME}-0-*"); do
	THREAD_ID=${FILE#${INPUT_DIR}/${TMPFILE_NAME}-0-}
	cat ${TMPFILE_PATH}-tandems-${THREAD_ID} >> ${TANDEMS_FILE}
done

# Removing all temp files that are not used downstream
if [ ${DELETE_TMP_FILES} -eq 1 ]; then
	rm -f ${TMPFILE_PATH}*
fi