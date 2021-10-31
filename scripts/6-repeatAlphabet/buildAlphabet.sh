#!/bin/bash
# 
# Recodes every read into a sequence of "characters", that are substrings of a repeat in
# a specific orientation. Builds an alphabet of all distinct characters. Discards rare
# characters and recodes the reads again in the cleaned-up alphabet.
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
N_THREADS="4"
MIN_CHARACTER_FREQUENCY="5"  # Should be the coverage of one haplotype
# REVANT
JAVA_RUNTIME_FLAGS="-Xms2G -Xmx10G"
# ----------------------------------------------------------------------------------------


READ_LENGTHS_FILE="${INPUT_DIR}/reads-lengths.txt"
READ_IDS_FILE="${INPUT_DIR}/reads-ids.txt"
N_READS=$(wc -l < ${READ_IDS_FILE})
REPEAT_LENGTHS_FILE="${INPUT_DIR}/repeats-lengths.txt"
REPEAT_ISPERIODIC_FILE="${INPUT_DIR}/repeats-isPeriodic.txt"
N_REPEATS=$(wc -l < ${REPEAT_LENGTHS_FILE})
TMPFILE_NAME="buildAlphabet-tmp"
TMPFILE_PATH="${INPUT_DIR}/${TMPFILE_NAME}"
SORT_OPTIONS=""
for i in $(seq 1 9); do  # Should be in sync with the serialization of $Character$.
	SORT_OPTIONS="${SORT_OPTIONS} -k ${i},${i}n"
done
rm -rf ${TMPFILE_PATH}*

echo "Splitting the alignments file..."
ALIGNMENTS_FILE="${INPUT_DIR}/LAshow-reads-repeats.txt"
N_ALIGNMENTS=$(( $(wc -l < ${ALIGNMENTS_FILE}) - 2 ))
LAST_READA_FILE="${INPUT_DIR}/LAshow-lastReadA.txt"
java ${JAVA_RUNTIME_FLAGS} -classpath "${REVANT_BINARIES}" de.mpi_cbg.revant.factorize.SplitAlignments ${N_ALIGNMENTS} ${N_THREADS} ${ALIGNMENTS_FILE} ${TMPFILE_PATH}-1- ${LAST_READA_FILE}
echo "Alignments filtered and split in ${N_THREADS} parts"

echo "Collecting character instances..."
function collectionThread() {
	ALIGNMENTS_FILE_ID=$1
	PREFIX_1=$2
	PREFIX_2=$3
	PREFIX_3=$4
	java ${JAVA_RUNTIME_FLAGS} -classpath "${REVANT_BINARIES}" de.mpi_cbg.revant.apps.GetCharacterInstances ${N_READS} ${READ_IDS_FILE} ${READ_LENGTHS_FILE} ${N_REPEATS} ${REPEAT_LENGTHS_FILE} ${REPEAT_ISPERIODIC_FILE} ${PREFIX_1}${ALIGNMENTS_FILE_ID}.txt ${MAX_ALIGNMENT_ERROR} ${PREFIX_2}${ALIGNMENTS_FILE_ID}.txt "${PREFIX_2}unique-${ALIGNMENTS_FILE_ID}.txt"
	sort --parallel 1 -t , ${SORT_OPTIONS} ${PREFIX_2}${ALIGNMENTS_FILE_ID}.txt | uniq - ${PREFIX_3}${ALIGNMENTS_FILE_ID}.txt
}
if [ -e ${TMPFILE_PATH}-1-${N_THREADS}.txt ]; then
	TO=${N_THREADS}
else
	TO=$(( ${N_THREADS} - 1 ))
fi
for THREAD in $(seq 0 ${TO}); do
	collectionThread ${THREAD} "${TMPFILE_PATH}-1-" "${TMPFILE_PATH}-2-" "${TMPFILE_PATH}-3-" &
done
wait
sort --parallel ${N_THREADS} -m -t , ${SORT_OPTIONS} ${TMPFILE_PATH}-3-*.txt | uniq - ${TMPFILE_PATH}-4.txt
N_INSTANCES=$(wc -l < ${TMPFILE_PATH}-4.txt)
java ${JAVA_RUNTIME_FLAGS} -classpath "${REVANT_BINARIES}" de.mpi_cbg.revant.apps.SplitCharacterInstances ${N_INSTANCES} ${N_THREADS} ${TMPFILE_PATH}-4.txt ${TMPFILE_PATH}-5-

echo "Compacting character instances..."
function compactionThread() {
	INSTANCES_FILE_ID=$1
	PREFIX_1=$2
	PREFIX_2=$3
	PREFIX_3=$4
	cat ${PREFIX_1}${INSTANCES_FILE_ID}-header.txt ${PREFIX_1}${INSTANCES_FILE_ID}.txt > ${PREFIX_2}${INSTANCES_FILE_ID}.txt
	java ${JAVA_RUNTIME_FLAGS} -classpath "${REVANT_BINARIES}" de.mpi_cbg.revant.apps.CompactCharacterInstances ${PREFIX_2}${INSTANCES_FILE_ID}.txt ${N_REPEATS} ${REPEAT_LENGTHS_FILE} ${PREFIX_3}${INSTANCES_FILE_ID}.txt
}
if [ -e ${TMPFILE_PATH}-5-${N_THREADS}.txt ]; then
	TO=${N_THREADS}
else
	TO=$(( ${N_THREADS} - 1 ))
fi
for THREAD in $(seq 0 ${TO}); do
	compactionThread ${THREAD} "${TMPFILE_PATH}-5-" "${TMPFILE_PATH}-6-" "${TMPFILE_PATH}-7-" &
done
wait
ALPHABET_FILE="${INPUT_DIR}/alphabet.txt"
rm -f ${ALPHABET_FILE}
java ${JAVA_RUNTIME_FLAGS} -classpath "${REVANT_BINARIES}" de.mpi_cbg.revant.apps.MergeAlphabetHeaders ${INPUT_DIR} "${TMPFILE_NAME}-7-" "${TMPFILE_NAME}-2-unique-" ${ALPHABET_FILE}
for THREAD in $(seq 0 ${TO}); do
	tail -n +2 ${TMPFILE_PATH}-7-${THREAD}.txt >> ${ALPHABET_FILE}
done
ALPHABET_SIZE=$( wc -l < ${ALPHABET_FILE} )
ALPHABET_SIZE=$(( ${ALPHABET_SIZE} - 1 ))

echo "Translating reads..."
function translationThread() {
	ALIGNMENTS_FILE_ID=$1
	LAST_TRANSLATED_READ=$2
	PREFIX_1=$3
	PREFIX_2=$4
	PREFIX_3=$5
	PREFIX_4=$6
	java ${JAVA_RUNTIME_FLAGS} -classpath "${REVANT_BINARIES}" de.mpi_cbg.revant.apps.TranslateReads ${N_READS} ${READ_IDS_FILE} ${READ_LENGTHS_FILE} ${N_REPEATS} ${REPEAT_LENGTHS_FILE} ${REPEAT_ISPERIODIC_FILE} ${PREFIX_1}${ALIGNMENTS_FILE_ID}.txt ${MAX_ALIGNMENT_ERROR} ${ALPHABET_FILE} ${LAST_TRANSLATED_READ} ${PREFIX_2}${ALIGNMENTS_FILE_ID}.txt ${PREFIX_3}${ALIGNMENTS_FILE_ID}.txt ${PREFIX_4}${ALIGNMENTS_FILE_ID}.txt
}
if [ -e ${TMPFILE_PATH}-1-${N_THREADS}.txt ]; then
	TO=${N_THREADS}
else
	TO=$(( ${N_THREADS} - 1 ))
fi
translationThread 0 -1 "${TMPFILE_PATH}-1-" "${TMPFILE_PATH}-8-" "${TMPFILE_PATH}-9-" "${TMPFILE_PATH}-9-hist-" &
for THREAD in $(seq 1 ${TO}); do
	LAST_TRANSLATED_READ=$(tail -n 1 ${TMPFILE_PATH}-1-$(( ${THREAD} - 1 )).txt | awk '{ print $1 }' | tr -d , )
	LAST_TRANSLATED_READ=$(( ${LAST_TRANSLATED_READ} - 1 ))
	translationThread ${THREAD} ${LAST_TRANSLATED_READ} "${TMPFILE_PATH}-1-" "${TMPFILE_PATH}-8-" "${TMPFILE_PATH}-9-" "${TMPFILE_PATH}-9-hist-" &
done
wait
READS_TRANSLATED_FILE="${INPUT_DIR}/reads-translated.txt"
rm -f ${READS_TRANSLATED_FILE}
for THREAD in $(seq 0 ${TO}); do
	cat ${TMPFILE_PATH}-8-${THREAD}.txt >> ${READS_TRANSLATED_FILE}
done
READS_TRANSLATED_BOUNDARIES="${INPUT_DIR}/reads-translated-boundaries.txt"
rm -f ${READS_TRANSLATED_BOUNDARIES}
for THREAD in $(seq 0 ${TO}); do
	cat ${TMPFILE_PATH}-9-${THREAD}.txt >> ${READS_TRANSLATED_BOUNDARIES}
done
HISTOGRAM_FILE="${INPUT_DIR}/reads-translated-histogram.txt"
rm -f ${HISTOGRAM_FILE}
PASTE_OPTIONS="awk '{print \$1"
for i in $(seq 2 $((${TO} + 1))); do
	PASTE_OPTIONS="${PASTE_OPTIONS}+\$$i"
done
PASTE_OPTIONS="${PASTE_OPTIONS}}'"
paste ${TMPFILE_PATH}-9-hist-*.txt | eval ${PASTE_OPTIONS} - > ${HISTOGRAM_FILE}

echo "Discarding rare characters..."
COUNTS_FILE="${INPUT_DIR}/alphabet-counts.txt"
HISTOGRAM_FILE="${INPUT_DIR}/alphabet-histogram.txt"
java ${JAVA_RUNTIME_FLAGS} -classpath "${REVANT_BINARIES}" de.mpi_cbg.revant.apps.GetCharacterCounts ${READS_TRANSLATED_FILE} ${ALPHABET_FILE} ${COUNTS_FILE} ${HISTOGRAM_FILE}
function cleaningThread() {
	TRANSLATED_CHARACTERS=$1
	TRANSLATED_BOUNDARIES=$2
	PREFIX_1=$3
	PREFIX_2=$4
	ID=$5
	java ${JAVA_RUNTIME_FLAGS} -classpath "${REVANT_BINARIES}" de.mpi_cbg.revant.apps.CleanTranslatedReads1 ${ALPHABET_FILE} ${COUNTS_FILE} ${N_READS} ${READ_IDS_FILE} ${READ_LENGTHS_FILE} ${TRANSLATED_CHARACTERS} ${TRANSLATED_BOUNDARIES} ${MIN_CHARACTER_FREQUENCY} ${PREFIX_1}${ID}.txt > ${PREFIX_1}unique-${ID}.txt
	sort --parallel 1 -t , ${SORT_OPTIONS} ${PREFIX_1}${ID}.txt | uniq - ${PREFIX_2}${ID}.txt
}
split -l $(( ${N_READS} / ${N_THREADS} )) ${READS_TRANSLATED_FILE} "${TMPFILE_PATH}-10-"
split -l $(( ${N_READS} / ${N_THREADS} )) ${READS_TRANSLATED_BOUNDARIES} "${TMPFILE_PATH}-11-"
for FILE in $(find ${INPUT_DIR} -name "${TMPFILE_NAME}-10-*" ); do
	THREAD_ID=${FILE#${INPUT_DIR}/${TMPFILE_NAME}-10-}
	cleaningThread "${INPUT_DIR}/${TMPFILE_NAME}-10-${THREAD_ID}" "${INPUT_DIR}/${TMPFILE_NAME}-11-${THREAD_ID}" "${TMPFILE_PATH}-12-" "${TMPFILE_PATH}-13-" ${THREAD_ID} &
done
wait
sort --parallel ${N_THREADS} -m -t , ${SORT_OPTIONS} ${TMPFILE_PATH}-13-*.txt | uniq - ${TMPFILE_PATH}-13.txt
cat ${TMPFILE_PATH}-12-unique-*.txt | sort -n -r > ${TMPFILE_PATH}-12-unique.txt
ALPHABET_FILE_CLEANED="${INPUT_DIR}/alphabet-cleaned.txt"
rm -f ${ALPHABET_FILE_CLEANED}
OLD2NEW_FILE="${INPUT_DIR}/alphabet-old2new.txt"
rm -f ${OLD2NEW_FILE}
java ${JAVA_RUNTIME_FLAGS} -classpath "${REVANT_BINARIES}" de.mpi_cbg.revant.apps.CleanTranslatedReads2 ${ALPHABET_FILE} ${COUNTS_FILE} $(wc -l < ${TMPFILE_PATH}-13.txt) ${TMPFILE_PATH}-13.txt ${MIN_CHARACTER_FREQUENCY} ${TMPFILE_PATH}-12-unique.txt ${ALPHABET_FILE_CLEANED} ${OLD2NEW_FILE}
function cleaningThread2() {
	TRANSLATED_CHARACTERS_OLD=$1
	TRANSLATED_BOUNDARIES_OLD=$2
	TRANSLATED_CHARACTERS_NEW=$3
	TRANSLATED_BOUNDARIES_NEW=$4
	HISTOGRAM_FILE=$5
	java ${JAVA_RUNTIME_FLAGS} -classpath "${REVANT_BINARIES}" de.mpi_cbg.revant.apps.CleanTranslatedReads3 ${N_READS} ${READ_IDS_FILE} ${READ_LENGTHS_FILE} ${ALPHABET_FILE} ${COUNTS_FILE} ${TRANSLATED_CHARACTERS_OLD} ${TRANSLATED_BOUNDARIES_OLD} ${MIN_CHARACTER_FREQUENCY} ${ALPHABET_FILE_CLEANED} ${OLD2NEW_FILE} ${TRANSLATED_CHARACTERS_NEW} ${TRANSLATED_BOUNDARIES_NEW} ${HISTOGRAM_FILE}
}
for FILE in $(find ${INPUT_DIR} -name "${TMPFILE_NAME}-10-*" ); do
	THREAD_ID=${FILE#${INPUT_DIR}/${TMPFILE_NAME}-10-}
	cleaningThread2 "${INPUT_DIR}/${TMPFILE_NAME}-10-${THREAD_ID}" "${INPUT_DIR}/${TMPFILE_NAME}-11-${THREAD_ID}" "${INPUT_DIR}/${TMPFILE_NAME}-14-${THREAD_ID}" "${INPUT_DIR}/${TMPFILE_NAME}-15-${THREAD_ID}" "${INPUT_DIR}/${TMPFILE_NAME}-16-${THREAD_ID}" &
done
wait
READS_TRANSLATED_FILE_NEW="${INPUT_DIR}/reads-translated-new.txt"
READS_TRANSLATED_BOUNDARIES_NEW="${INPUT_DIR}/reads-translated-boundaries-new.txt"
HISTOGRAM_FILE="${INPUT_DIR}/reads-translated-histogram-new.txt"
rm -f ${READS_TRANSLATED_FILE_NEW} ${READS_TRANSLATED_BOUNDARIES_NEW}
for FILE in $(find -s ${INPUT_DIR} -name "${TMPFILE_NAME}-14-*" ); do
	THREAD_ID=${FILE#${INPUT_DIR}/${TMPFILE_NAME}-14-}
	cat ${INPUT_DIR}/${TMPFILE_NAME}-14-${THREAD_ID} >> ${READS_TRANSLATED_FILE_NEW}
	cat ${INPUT_DIR}/${TMPFILE_NAME}-15-${THREAD_ID} >> ${READS_TRANSLATED_BOUNDARIES_NEW}
done
PASTE_OPTIONS="awk '{print \$1"
for i in $(seq 2 ${N_THREADS}); do
	PASTE_OPTIONS="${PASTE_OPTIONS}+\$$i"
done
PASTE_OPTIONS="${PASTE_OPTIONS}}'"
paste ${INPUT_DIR}/${TMPFILE_NAME}-16-* | eval ${PASTE_OPTIONS} - > ${HISTOGRAM_FILE}
