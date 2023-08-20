#!/bin/bash
# 
# For increasing values of k, the script tries to disambiguate the endpoints of reads
# using a context of length k.
#
# Remark: since collecting k-mers is expensive, one might want to try just a few values of
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
BROKEN_READS=$2  # 1=TRUE
LOW_QUALITY_TYPE=$3  # 1=replacement, 0=insertion.
LOW_QUALITY_LENGTH_TOLERANCE="200"  # bps
MIN_K=$4  # One plus the min length of a context used for disambiguation
MAX_K=$5  # One plus the max length of a context used for disambiguation
N_THREADS=$6
DELETE_TMP_FILES=$7
GENOME_LENGTH=$8
N_HAPLOTYPES=$9
TIGHT_MODE="0"
SPANNING_BPS="150"  # Bps before and after a k-mer to consider it observed in a read.
# ------------------------------------ REVANT --------------------------------------------
REVANT_LIBRARIES="${REVANT_BINARIES}/../lib"
REVANT_LIBRARIES="${REVANT_LIBRARIES}/commons-numbers-gamma-1.1.jar:${REVANT_LIBRARIES}/commons-rng-sampling-1.5.jar:${REVANT_LIBRARIES}/commons-statistics-distribution-1.0.jar"
# ----------------------------------------------------------------------------------------

set -o pipefail; set -e; set -u
export LC_ALL=C  # To speed up the $sort$ command.
TMPFILE_NAME="fixEndBlocks-tmp"
TMPFILE_PATH="${INPUT_DIR}/${TMPFILE_NAME}"
READ_IDS_FILE="${INPUT_DIR}/reads-ids.txt"
N_READS=$(wc -l < ${READ_IDS_FILE})
READ_LENGTHS_FILE="${INPUT_DIR}/reads-lengths.txt"
AVG_READ_LENGTH=$(paste -sd+ ${READ_LENGTHS_FILE} | bc)
AVG_READ_LENGTH=$(( ${AVG_READ_LENGTH} / ${N_READS} ))
READS_TRANSLATED_FILE="${INPUT_DIR}/reads-translated-new.txt"
READS_BOUNDARIES_FILE="${INPUT_DIR}/reads-translated-boundaries-new.txt"
READS_DISAMBIGUATED_FILE="${INPUT_DIR}/reads-translated-disambiguated.txt"
ALPHABET_FILE="${INPUT_DIR}/alphabet-cleaned.txt"
rm -f ${TMPFILE_PATH}*

function kmersThread() {
	local LOCAL_K=$1
	local LOCAL_TRANSLATED_READS_FILE=$2
	local LOCAL_BOUNDARIES_FILE=$3
	local LOCAL_READ_LENGTHS_FILE=$4
	local LOCAL_KMERS_FILE=$5
	java ${JAVA_RUNTIME_FLAGS} -classpath "${REVANT_BINARIES}" de.mpi_cbg.revant.apps.CollectKmers ${LOCAL_K} ${LOCAL_TRANSLATED_READS_FILE} ${LOCAL_BOUNDARIES_FILE} ${LOCAL_READ_LENGTHS_FILE} ${ALPHABET_FILE} null ${LOCAL_KMERS_FILE}
	if [ $? -ne 0 ]; then
		exit
	fi
}

function fixThread() {
	local LOCAL_OLD_TRANSLATED_FILE=$1
	local LOCAL_BOUNDARIES_FILE=$2
	local LOCAL_READ_LENGTHS_FILE=$3
	local LOCAL_KMERS_FILE=$4
	local LOCAL_K=$5
	local LOCAL_NEW_TRANSLATED_FILE=$6
	local LOCAL_STATS_FILE=$7
	java ${JAVA_RUNTIME_FLAGS} -classpath "${REVANT_BINARIES}" de.mpi_cbg.revant.apps.FixEndBlocks ${ALPHABET_FILE} ${LOCAL_OLD_TRANSLATED_FILE} ${LOCAL_BOUNDARIES_FILE} ${LOCAL_READ_LENGTHS_FILE} ${LOCAL_KMERS_FILE} ${LOCAL_K} ${TIGHT_MODE} ${LOCAL_NEW_TRANSLATED_FILE} ${LOCAL_STATS_FILE}
	if [ $? -ne 0 ]; then
		exit
	fi
}

rm -f "${TMPFILE_PATH}-1-*"
split -l $(( ${N_READS} / ${N_THREADS} )) ${READS_TRANSLATED_FILE} "${TMPFILE_PATH}-1-"
split -l $(( ${N_READS} / ${N_THREADS} )) ${READS_BOUNDARIES_FILE} "${TMPFILE_PATH}-z1-"
split -l $(( ${N_READS} / ${N_THREADS} )) ${READ_LENGTHS_FILE} "${TMPFILE_PATH}-z2-"
for K in $(seq ${MIN_K} ${MAX_K}); do
	SORT_OPTIONS_KMERS=""
	for i in $(seq 1 ${K}); do
		SORT_OPTIONS_KMERS="${SORT_OPTIONS_KMERS} -k ${i},${i}n"
	done
	echo "Collecting ${K}-mers..."
	for FILE in $(find -s ${INPUT_DIR} -name "${TMPFILE_NAME}-1-*"); do
		THREAD_ID=${FILE#${INPUT_DIR}/${TMPFILE_NAME}-1-}
		kmersThread ${K} ${FILE} ${TMPFILE_PATH}-z1-${THREAD_ID} ${TMPFILE_PATH}-z2-${THREAD_ID} ${TMPFILE_PATH}-kmers-${K}-${THREAD_ID} &
	done
	wait
	sort --parallel=${N_THREADS} -m -t , ${SORT_OPTIONS_KMERS} ${TMPFILE_PATH}-kmers-${K}-* > ${TMPFILE_PATH}-${K}.txt
	if [ ! -s ${TMPFILE_PATH}-${K}.txt ]; then
		MAX_K=$((${K}-1))
		break
	fi
	FREQUENT_KMERS_FILE="${INPUT_DIR}/frequent-k${K}.txt"
    OUTPUT_FILE_HISTOGRAM="${INPUT_DIR}/histogram-k${K}.txt"
	echo "Finding frequent ${K}-mers..."
	java ${JAVA_RUNTIME_FLAGS} -classpath "${REVANT_BINARIES}:${REVANT_LIBRARIES}" de.mpi_cbg.revant.apps.CompactKmers ${TMPFILE_PATH}-${K}.txt ${K} ${GENOME_LENGTH} ${N_HAPLOTYPES} ${N_READS} ${AVG_READ_LENGTH} ${SPANNING_BPS} 0 ${ALPHABET_FILE} 1 10000 ${FREQUENT_KMERS_FILE} ${OUTPUT_FILE_HISTOGRAM}
	echo "Computing $((${K}-1))-mers..."
	K_MINUS_ONE_MERS_FILE="${INPUT_DIR}/kMinusOne-k${K}.txt"
	java ${JAVA_RUNTIME_FLAGS} -classpath "${REVANT_BINARIES}" de.mpi_cbg.revant.apps.GetKMinusOneMers ${ALPHABET_FILE} ${FREQUENT_KMERS_FILE} ${K} ${K_MINUS_ONE_MERS_FILE}
	echo "Disambiguating read ends using contexts of length $((${K}-1))..."
	for FILE in $(find -s ${INPUT_DIR} -name "${TMPFILE_NAME}-$((${K}-1))-*"); do
		THREAD_ID=${FILE#${INPUT_DIR}/${TMPFILE_NAME}-$((${K}-1))-}
		fixThread ${FILE} ${TMPFILE_PATH}-z1-${THREAD_ID} ${TMPFILE_PATH}-z2-${THREAD_ID} ${K_MINUS_ONE_MERS_FILE} $((${K}-1)) ${TMPFILE_PATH}-${K}-${THREAD_ID} ${TMPFILE_PATH}-counts-${K}-${THREAD_ID} &
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
if [ ${BROKEN_READS} -eq 1 ]; then
	READS_DISAMBIGUATED_FILE_NEW="${INPUT_DIR}/reads-translated-disambiguated-new.txt"
	READS_DISAMBIGUATED_FILE_PRE="${INPUT_DIR}/reads-translated-disambiguated-pre.txt"
	rm -f ${READS_DISAMBIGUATED_FILE_NEW} ${READS_DISAMBIGUATED_FILE_PRE}
	NEW2OLD_FILE="${INPUT_DIR}/broken2unbroken.txt"
	READS_TRANSLATED_BOUNDARIES_NEW="${INPUT_DIR}/reads-translated-boundaries-new.txt"
	java ${JAVA_RUNTIME_FLAGS} -classpath "${REVANT_BINARIES}" de.mpi_cbg.revant.apps.BreakReads3 ${LOW_QUALITY_TYPE} ${LOW_QUALITY_LENGTH_TOLERANCE} ${NEW2OLD_FILE} ${N_READS} ${READ_LENGTHS_FILE} ${ALPHABET_FILE} ${READS_TRANSLATED_FILE} ${READS_DISAMBIGUATED_FILE} ${READS_TRANSLATED_BOUNDARIES_NEW} ${READS_DISAMBIGUATED_FILE_NEW}
	mv ${READS_DISAMBIGUATED_FILE} ${READS_DISAMBIGUATED_FILE_PRE}
	mv ${READS_DISAMBIGUATED_FILE_NEW} ${READS_DISAMBIGUATED_FILE}
fi

# Removing all temp files that are not used downstream
if [ ${DELETE_TMP_FILES} -eq 1 ]; then
	rm -f ${TMPFILE_PATH}*
fi
