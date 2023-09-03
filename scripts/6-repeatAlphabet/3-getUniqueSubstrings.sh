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
GENOME_LENGTH=$2
N_HAPLOTYPES=$3
MAX_K=$4  # Stops looking for unique k-mers after this length. Should be set using the
# histogram of recoded lengths.
N_THREADS=$5
DELETE_TMP_FILES=$6
IDENTITY_THRESHOLD=$7
DISTANCE_THRESHOLD=$8
CHARACTER_THRESHOLD=$9
MIN_ALIGNMENT_LENGTH=${10}  # Read-repeat
MAX_KMER_LENGTH_BPS=${11}
UNIQUE_MODE="1"  # Non-repetitive blocks are allowed in a k-mer, except at the first/last
# position of the k-mer. Usually a good choice.
SPANNING_BPS="150"  # Bps before and after a k-mer to consider it observed in a read.
# ------------------------------------ REVANT --------------------------------------------
REVANT_LIBRARIES="${REVANT_BINARIES}/../lib/*"
# ----------------------------------------------------------------------------------------


set -euo pipefail
export LC_ALL=C  # To speed up the $sort$ command.
READ_IDS_FILE="${INPUT_DIR}/reads-ids.txt"
N_READS=$(wc -l < ${READ_IDS_FILE})
READ_LENGTHS_FILE="${INPUT_DIR}/reads-lengths.txt"
AVG_READ_LENGTH=$(paste -sd+ ${READ_LENGTHS_FILE} | bc)
AVG_READ_LENGTH=$(( ${AVG_READ_LENGTH} / ${N_READS} ))
TMPFILE_NAME="getUniqueSubstrings-tmp"
TMPFILE_PATH="${INPUT_DIR}/${TMPFILE_NAME}"
READS_TRANSLATED_FILE="${INPUT_DIR}/reads-translated-disambiguated.txt"
READS_TRANSLATED_BOUNDARIES="${INPUT_DIR}/reads-translated-boundaries-new.txt"
ALPHABET_FILE="${INPUT_DIR}/alphabet-cleaned.txt"
TANDEMS_FILE="${INPUT_DIR}/tandems.txt"
REPEAT_LENGTHS_FILE="${INPUT_DIR}/repeats-lengths.txt"
N_REPEATS=$(wc -l < ${REPEAT_LENGTHS_FILE})
MAX_HISTOGRAM_COUNT="10000"  # Arbitrary
rm -f ${TMPFILE_PATH}*
rm -f ${INPUT_DIR}/unique-*
rm -f ${INPUT_DIR}/histogram-*

function waitAndCheck() {
    local ARRAY_NAME=$1[@]
    
    local PIDS=(${!ARRAY_NAME})
    local LAST_THREAD=$((${#PIDS[@]} - 1))
    N_FAILED="0"
    for THREAD in $(seq 0 ${LAST_THREAD}); do
        wait ${PIDS[${THREAD}]} || N_FAILED=$(( ${N_FAILED} + 1 ))
    done
    if [ ${N_FAILED} -ne 0 ]; then
        exit 1
    fi
}

function enumerateKmersThread() {
	local LOCAL_K=$1
	local LOCAL_TRANSLATED_READS_FILE=$2
	local LOCAL_BOUNDARIES_FILE=$3
	local LOCAL_READ_LENGTHS_FILE=$4
	local LOCAL_K_MINUS_ONE_INTERVALS_FILE=$5
	local LOCAL_KMERS_FILE=$6
	java ${JAVA_RUNTIME_FLAGS} -classpath "${REVANT_BINARIES}" de.mpi_cbg.revant.apps.CollectKmers 0 ${LOCAL_K} ${MAX_KMER_LENGTH_BPS} 2 ${LOCAL_TRANSLATED_READS_FILE} ${LOCAL_BOUNDARIES_FILE} ${LOCAL_READ_LENGTHS_FILE} ${ALPHABET_FILE} ${LOCAL_K_MINUS_ONE_INTERVALS_FILE} null ${LOCAL_KMERS_FILE}
}

function countKmersThread() {
	local LOCAL_K=$1
	local LOCAL_TRANSLATED_READS_FILE=$2
	local LOCAL_BOUNDARIES_FILE=$3
	local LOCAL_READ_LENGTHS_FILE=$4
	local LOCAL_K_MINUS_ONE_INTERVALS_FILE=$5
	local LOCAL_KMERS_FILE_INPUT=$6
    local LOCAL_KMERS_FILE_OUTPUT=$7
	java ${JAVA_RUNTIME_FLAGS} -classpath "${REVANT_BINARIES}" de.mpi_cbg.revant.apps.CollectKmers 1 ${LOCAL_K} ${MAX_KMER_LENGTH_BPS} 2 ${LOCAL_TRANSLATED_READS_FILE} ${LOCAL_BOUNDARIES_FILE} ${LOCAL_READ_LENGTHS_FILE} ${ALPHABET_FILE} ${LOCAL_K_MINUS_ONE_INTERVALS_FILE} ${LOCAL_KMERS_FILE_INPUT} ${LOCAL_KMERS_FILE_OUTPUT}
}

function intervalsThread() {
	local LOCAL_K=$1
	local LOCAL_TRANSLATED_READS_FILE=$2
	local LOCAL_BOUNDARIES_FILE=$3
	local LOCAL_READ_LENGTHS_FILE=$4
	local LOCAL_UNIQUE_KMERS_FILE=$5
	local LOCAL_K_MINUS_ONE_INTERVALS_FILE=$6
	local LOCAL_INTERVALS_FILE=$7
	java ${JAVA_RUNTIME_FLAGS} -classpath "${REVANT_BINARIES}:${REVANT_LIBRARIES}" de.mpi_cbg.revant.apps.GetShortestUniqueIntervals ${LOCAL_K} ${MAX_KMER_LENGTH_BPS} ${LOCAL_TRANSLATED_READS_FILE} ${LOCAL_BOUNDARIES_FILE} ${LOCAL_READ_LENGTHS_FILE} ${ALPHABET_FILE} ${LOCAL_UNIQUE_KMERS_FILE} ${N_READS} ${AVG_READ_LENGTH} ${GENOME_LENGTH} ${N_HAPLOTYPES} ${MIN_ALIGNMENT_LENGTH} ${IDENTITY_THRESHOLD} ${DISTANCE_THRESHOLD} ${CHARACTER_THRESHOLD} ${LOCAL_K_MINUS_ONE_INTERVALS_FILE} ${LOCAL_INTERVALS_FILE}
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
	echo "Enumerating distinct ${K}-mers..."
    PIDS=()
	for FILE in $(find -s ${INPUT_DIR} -maxdepth 1 -name "${TMPFILE_NAME}-0-*"); do
		THREAD_ID=${FILE#${INPUT_DIR}/${TMPFILE_NAME}-0-}
		if [ ${K} -le 1 ]; then
			PREVIOUS_INTERVALS="null"
		else
			PREVIOUS_INTERVALS="${TMPFILE_PATH}-$((${K}-1))-intervals-${THREAD_ID}"
		fi
		enumerateKmersThread ${K} ${FILE} ${TMPFILE_PATH}-1-${THREAD_ID} ${TMPFILE_PATH}-2-${THREAD_ID} ${PREVIOUS_INTERVALS} ${TMPFILE_PATH}-${K}-kmers-${THREAD_ID} &
        PIDS+=($!)
	done
	waitAndCheck PIDS
	sort --parallel=${N_THREADS} -m -u -t , ${SORT_OPTIONS_KMERS} ${TMPFILE_PATH}-${K}-kmers-* > ${TMPFILE_PATH}-${K}-distinct.txt
	if [ ! -s ${TMPFILE_PATH}-${K}-distinct.txt ]; then
		MAX_K=$((${K}-1))
		break
	fi
	echo "Counting ${K}-mer occurrences..."
    PIDS=()
	for FILE in $(find -s ${INPUT_DIR} -maxdepth 1 -name "${TMPFILE_NAME}-0-*"); do
		THREAD_ID=${FILE#${INPUT_DIR}/${TMPFILE_NAME}-0-}
		if [ ${K} -le 1 ]; then
			PREVIOUS_INTERVALS="null"
		else
			PREVIOUS_INTERVALS="${TMPFILE_PATH}-$((${K}-1))-intervals-${THREAD_ID}"
		fi
        rm -f ${TMPFILE_PATH}-${K}-kmers-${THREAD_ID}
		countKmersThread ${K} ${FILE} ${TMPFILE_PATH}-1-${THREAD_ID} ${TMPFILE_PATH}-2-${THREAD_ID} ${PREVIOUS_INTERVALS} ${TMPFILE_PATH}-${K}-distinct.txt ${TMPFILE_PATH}-${K}-kmers-${THREAD_ID} &
        PIDS+=($!)
	done
	waitAndCheck PIDS
    sort --parallel=${N_THREADS} -m -t , ${SORT_OPTIONS_KMERS} ${TMPFILE_PATH}-${K}-kmers-* > ${TMPFILE_PATH}-${K}.txt
    rm -f ${TMPFILE_PATH}-${K}-distinct.txt
	UNIQUE_KMERS_FILE="${INPUT_DIR}/unique-k${K}.txt"
	OUTPUT_FILE_HISTOGRAM="${INPUT_DIR}/histogram-k${K}.txt"
	echo "Finding unique ${K}-mers..."
	java ${JAVA_RUNTIME_FLAGS} -classpath "${REVANT_BINARIES}:${REVANT_LIBRARIES}" de.mpi_cbg.revant.apps.CompactKmers ${TMPFILE_PATH}-${K}.txt ${K} ${GENOME_LENGTH} ${N_HAPLOTYPES} ${N_READS} ${AVG_READ_LENGTH} ${SPANNING_BPS} ${MIN_ALIGNMENT_LENGTH} 1 ${ALPHABET_FILE} 0 ${MAX_HISTOGRAM_COUNT} ${UNIQUE_KMERS_FILE} ${OUTPUT_FILE_HISTOGRAM}
	echo "Updating shortest unique intervals file..."
    PIDS=()
	for FILE in $(find -s ${INPUT_DIR} -maxdepth 1 -name "${TMPFILE_NAME}-0-*"); do
		THREAD_ID=${FILE#${INPUT_DIR}/${TMPFILE_NAME}-0-}
		if [ ${K} -le 1 ]; then
			PREVIOUS_INTERVALS="null"
		else
			PREVIOUS_INTERVALS="${TMPFILE_PATH}-$((${K}-1))-intervals-${THREAD_ID}"
		fi
		intervalsThread ${K} ${FILE} ${TMPFILE_PATH}-1-${THREAD_ID} ${TMPFILE_PATH}-2-${THREAD_ID} ${UNIQUE_KMERS_FILE} ${PREVIOUS_INTERVALS} ${TMPFILE_PATH}-${K}-intervals-${THREAD_ID} &
        PIDS+=($!)
	done
	waitAndCheck PIDS
done
rm -f ${FINAL_INTERVALS_FILE}
for FILE in $(find -s ${INPUT_DIR} -maxdepth 1 -name "${TMPFILE_NAME}-${MAX_K}-intervals-*" ); do
	cat ${FILE} >> ${FINAL_INTERVALS_FILE}
done
INTERVAL_STATS_FILE="${INPUT_DIR}/unique-intervals-k1-${MAX_K}-stats.txt"
rm -f ${INTERVAL_STATS_FILE}
java ${JAVA_RUNTIME_FLAGS} -classpath "${REVANT_BINARIES}" de.mpi_cbg.revant.apps.UniqueIntervalsStats ${FINAL_INTERVALS_FILE} ${READS_TRANSLATED_BOUNDARIES} ${INTERVAL_STATS_FILE}

# Collecting tandem intervals
NONPERIODIC_MODE="0"  # Tightest definition of non-periodic tandem
function tandemsThread() {
	local LOCAL_TRANSLATED_READS_FILE=$1
	local LOCAL_BOUNDARIES_FILE=$2
	local LOCAL_READ_LENGTHS_FILE=$3
	local LOCAL_TANDEMS_FILE=$4
	java ${JAVA_RUNTIME_FLAGS} -classpath "${REVANT_BINARIES}" de.mpi_cbg.revant.apps.CollectTandems 1 1 ${NONPERIODIC_MODE} ${ALPHABET_FILE} ${LOCAL_TRANSLATED_READS_FILE} ${LOCAL_BOUNDARIES_FILE} ${LOCAL_READ_LENGTHS_FILE} ${REPEAT_LENGTHS_FILE} ${N_REPEATS} ${LOCAL_TANDEMS_FILE}
}
echo "Collecting tandems..."
PIDS=()
for FILE in $(find -s ${INPUT_DIR} -maxdepth 1 -name "${TMPFILE_NAME}-0-*"); do
	THREAD_ID=${FILE#${INPUT_DIR}/${TMPFILE_NAME}-0-}
	tandemsThread ${FILE} ${TMPFILE_PATH}-1-${THREAD_ID} ${TMPFILE_PATH}-2-${THREAD_ID} ${TMPFILE_PATH}-tandems-${THREAD_ID} &
    PIDS+=($!)
done
waitAndCheck PIDS
rm -f ${TANDEMS_FILE}
for FILE in $(find -s ${INPUT_DIR} -maxdepth 1 -name "${TMPFILE_NAME}-0-*"); do
	THREAD_ID=${FILE#${INPUT_DIR}/${TMPFILE_NAME}-0-}
	cat ${TMPFILE_PATH}-tandems-${THREAD_ID} >> ${TANDEMS_FILE}
done

# Removing all temp files that are not used downstream
if [ ${DELETE_TMP_FILES} -eq 1 ]; then
	rm -f ${TMPFILE_PATH}*
fi