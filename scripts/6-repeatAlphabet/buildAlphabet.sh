#!/bin/bash
# 
# 
#
# Remark: the script assumes that environment variable $REVANT_BINARIES$ is already set
# to the directory that contains REVANT's binaries, i.e. to the directory that contains
# the $de/mpi_cbg/revant/...$ subpath.
#
# --------------------------------- CONFIGURATION ----------------------------------------
# This is the only section of the script that needs to be customized.
#
INPUT_DIR=$1
MAX_ALIGNMENT_ERROR="0.2"  # Alignments with more error than this are discarded
SPLIT_IN_PARTS="4"  # For parallelism
# REVANT
JAVA_RUNTIME_FLAGS="-Xms2G -Xmx10G"
# ----------------------------------------------------------------------------------------


READ_LENGTHS_FILE="${INPUT_DIR}/reads-lengths.txt"
READ_IDS_FILE="${INPUT_DIR}/reads-ids.txt"
N_READS=$(wc -l < ${READ_IDS_FILE})
REPEAT_LENGTHS_FILE="${INPUT_DIR}/repeats-lengths.txt"
REPEAT_ISPERIODIC_FILE="${INPUT_DIR}/repeats-isPeriodic.txt"
N_REPEATS=$(wc -l < ${REPEAT_LENGTHS_FILE})
PARTS_PREFIX="${INPUT_DIR}/tmpSplit"
SORT_OPTIONS=""
for i in $(seq 1 9); do  # Should be in sync with the serialization of $Character$.
	SORT_OPTIONS="${SORT_OPTIONS} -k ${i},${i}n"
done
rm -rf ${PARTS_PREFIX}*

echo "Splitting the alignments file..."
ALIGNMENTS_FILE="${INPUT_DIR}/LAshow-reads-repeats.txt"
N_ALIGNMENTS=$(( $(wc -l < ${ALIGNMENTS_FILE}) - 2 ))
LAST_READA_FILE="${INPUT_DIR}/LAshow-lastReadA.txt"
java ${JAVA_RUNTIME_FLAGS} -classpath "${REVANT_BINARIES}" de.mpi_cbg.revant.factorize.SplitAlignments ${N_ALIGNMENTS} ${SPLIT_IN_PARTS} ${ALIGNMENTS_FILE} ${PARTS_PREFIX}-1- ${LAST_READA_FILE}
echo "Alignments filtered and split in ${SPLIT_IN_PARTS} parts"

echo "Collecting character instances..."
function collectionThread() {
	ALIGNMENTS_FILE_ID=$1
	PREFIX_1=$2
	PREFIX_2=$3
	PREFIX_3=$4
	java ${JAVA_RUNTIME_FLAGS} -classpath "${REVANT_BINARIES}" de.mpi_cbg.revant.apps.GetCharacterInstances ${N_READS} ${READ_IDS_FILE} ${READ_LENGTHS_FILE} ${N_REPEATS} ${REPEAT_LENGTHS_FILE} ${REPEAT_ISPERIODIC_FILE} ${PREFIX_1}${ALIGNMENTS_FILE_ID}.txt ${MAX_ALIGNMENT_ERROR} ${PREFIX_2}${ALIGNMENTS_FILE_ID}.txt
	sort --parallel 1 -t , ${SORT_OPTIONS} ${PREFIX_2}${ALIGNMENTS_FILE_ID}.txt | uniq - ${PREFIX_3}${ALIGNMENTS_FILE_ID}.txt
}
if [ -e ${PARTS_PREFIX}-1-${SPLIT_IN_PARTS}.txt ]; then
	TO=${SPLIT_IN_PARTS}
else
	TO=$(( ${SPLIT_IN_PARTS} - 1 ))
fi
for THREAD in $(seq 0 ${TO}); do
	collectionThread ${THREAD} "${PARTS_PREFIX}-1-" "${PARTS_PREFIX}-2-" "${PARTS_PREFIX}-3-" &
done
wait
sort --parallel ${SPLIT_IN_PARTS} -m -n -t , ${SORT_OPTIONS} ${PARTS_PREFIX}-3-*.txt | uniq - ${PARTS_PREFIX}-4.txt
N_INSTANCES=$(wc -l < ${PARTS_PREFIX}-4.txt)
java ${JAVA_RUNTIME_FLAGS} -classpath "${REVANT_BINARIES}" de.mpi_cbg.revant.apps.SplitCharacterInstances ${N_INSTANCES} ${SPLIT_IN_PARTS} ${PARTS_PREFIX}-4.txt ${PARTS_PREFIX}-5-

echo "Compacting character instances..."
function compactionThread() {
	INSTANCES_FILE_ID=$1
	PREFIX_1=$2
	PREFIX_2=$3
	PREFIX_3=$4
	cat ${PREFIX_1}${INSTANCES_FILE_ID}-header.txt ${PREFIX_1}${INSTANCES_FILE_ID}.txt > ${PREFIX_2}${INSTANCES_FILE_ID}.txt
	java ${JAVA_RUNTIME_FLAGS} -classpath "${REVANT_BINARIES}" de.mpi_cbg.revant.apps.CompactCharacterInstances ${PREFIX_2}${INSTANCES_FILE_ID}.txt ${N_REPEATS} ${REPEAT_LENGTHS_FILE} ${PREFIX_3}${INSTANCES_FILE_ID}.txt
}
if [ -e ${PARTS_PREFIX}-5-${SPLIT_IN_PARTS}.txt ]; then
	TO=${SPLIT_IN_PARTS}
else
	TO=$(( ${SPLIT_IN_PARTS} - 1 ))
fi
for THREAD in $(seq 0 ${TO}); do
	compactionThread ${THREAD} "${PARTS_PREFIX}-5-" "${PARTS_PREFIX}-6-" "${PARTS_PREFIX}-7-" &
done
wait
ALPHABET_FILE="${INPUT_DIR}/alphabet.txt"
rm -f ${ALPHABET_FILE}
java ${JAVA_RUNTIME_FLAGS} -classpath "${REVANT_BINARIES}" de.mpi_cbg.revant.apps.MergeAlphabetHeaders ${INPUT_DIR} "tmpSplit-7-" ${ALPHABET_FILE}
for THREAD in $(seq 0 ${TO}); do
	tail -n +2 ${PARTS_PREFIX}-7-${THREAD}.txt >> ${ALPHABET_FILE}
done
ALPHABET_SIZE=$( wc -l < ${ALPHABET_FILE} )
ALPHABET_SIZE=$(( ${ALPHABET_SIZE} - 1 ))

echo "Translating reads..."
function translationThread() {
	ALIGNMENTS_FILE_ID=$1
	LAST_TRANSLATED_READ=$2
	PREFIX_1=$3
	PREFIX_2=$4
	java ${JAVA_RUNTIME_FLAGS} -classpath "${REVANT_BINARIES}" de.mpi_cbg.revant.apps.TranslateReads ${N_READS} ${READ_IDS_FILE} ${READ_LENGTHS_FILE} ${N_REPEATS} ${REPEAT_LENGTHS_FILE} ${REPEAT_ISPERIODIC_FILE} ${PREFIX_1}${ALIGNMENTS_FILE_ID}.txt ${MAX_ALIGNMENT_ERROR} ${ALPHABET_FILE} ${LAST_TRANSLATED_READ} ${PREFIX_2}${ALIGNMENTS_FILE_ID}.txt
}
if [ -e ${PARTS_PREFIX}-1-${SPLIT_IN_PARTS}.txt ]; then
	TO=${SPLIT_IN_PARTS}
else
	TO=$(( ${SPLIT_IN_PARTS} - 1 ))
fi
translationThread 0 -1 "${PARTS_PREFIX}-1-" "${PARTS_PREFIX}-8-" &
for THREAD in $(seq 1 ${TO}); do
	LAST_TRANSLATED_READ=$(tail -n 1 ${PARTS_PREFIX}-1-$(( ${THREAD} - 1 )).txt | awk '{ print $1 }' | tr -d , )
	LAST_TRANSLATED_READ=$(( ${LAST_TRANSLATED_READ} - 1 ))
	translationThread ${THREAD} ${LAST_TRANSLATED_READ} "${PARTS_PREFIX}-1-" "${PARTS_PREFIX}-8-" &
done
wait
READS_TRANSLATED_FILE="${INPUT_DIR}/reads-translated.txt"
rm -f ${READS_TRANSLATED_FILE}
for THREAD in $(seq 0 ${TO}); do
	cat ${PARTS_PREFIX}-8-${THREAD}.txt >> ${READS_TRANSLATED_FILE}
done

echo "Computing character counts..."
COUNTS_FILE="${INPUT_DIR}/alphabet-counts.txt"
HISTOGRAM_FILE="${INPUT_DIR}/alphabet-histogram.txt"
java ${JAVA_RUNTIME_FLAGS} -classpath "${REVANT_BINARIES}" de.mpi_cbg.revant.apps.GetCharacterCounts ${READS_TRANSLATED_FILE} ${ALPHABET_FILE} ${COUNTS_FILE} ${HISTOGRAM_FILE}
