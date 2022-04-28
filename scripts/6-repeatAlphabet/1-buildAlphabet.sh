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
BROKEN_READS=$2  # 1=TRUE
KEEP_PERIODIC="1"  # 1=do not remove rare characters if they are periodic
MAX_ALIGNMENT_ERROR="0.3"  # Repeat-read alignments with error > this are discarded
MIN_ALIGNMENT_LENGTH="500"  # Repeat-read alignments with length < this are discarded
HAPLOTYPE_COVERAGE="30"  # Of one haplotype
MAX_SPACER_LENGTH="400"  # 0=assume that the endpoints of periodic repeats are accurate
N_THREADS="4"
DELETE_TMP_FILES="1"
# REVANT
JAVA_RUNTIME_FLAGS="-Xms2G -Xmx10G"
# ----------------------------------------------------------------------------------------

set -o pipefail; set -e; set -u
export LC_ALL=C  # To speed up the $sort$ command.
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
MIN_CHARACTER_FREQUENCY=$(( ${HAPLOTYPE_COVERAGE} / 2 ))
rm -f ${TMPFILE_PATH}*

echo "Splitting the alignments file..."
if [ ${BROKEN_READS} -eq 1 ]; then
	# Reusing the chunks of the read-repeat alignments file that are already there (we
	# assume that they all have the header).
	for FILE in $(ls ${INPUT_DIR}/breakReads-tmp-4-*.txt ); do
		ID=$(basename ${FILE} .txt)
		ID=${ID#breakReads-tmp-4-}
		mv ${INPUT_DIR}/breakReads-tmp-4-${ID}.txt ${TMPFILE_PATH}-1-${ID}.txt
	done
else
	ALIGNMENTS_FILE="${INPUT_DIR}/LAshow-reads-repeats.txt"
	N_ALIGNMENTS=$(( $(wc -l < ${ALIGNMENTS_FILE}) - 2 ))
	LAST_READA_FILE="${INPUT_DIR}/LAshow-lastReadA.txt"
	java ${JAVA_RUNTIME_FLAGS} -classpath "${REVANT_BINARIES}" de.mpi_cbg.revant.factorize.SplitAlignments ${N_ALIGNMENTS} ${N_THREADS} ${ALIGNMENTS_FILE} ${TMPFILE_PATH}-1- ${LAST_READA_FILE}
fi
echo "Alignments filtered and split in ${N_THREADS} parts"

echo "Collecting character instances..."
function collectionThread() {
	local ALIGNMENTS_FILE_ID=$1
	local PREFIX_1=$2
	local PREFIX_2=$3
	local PREFIX_3=$4
	java ${JAVA_RUNTIME_FLAGS} -classpath "${REVANT_BINARIES}" de.mpi_cbg.revant.apps.GetCharacterInstances ${N_READS} ${READ_IDS_FILE} ${READ_LENGTHS_FILE} ${N_REPEATS} ${REPEAT_LENGTHS_FILE} ${REPEAT_ISPERIODIC_FILE} ${PREFIX_1}${ALIGNMENTS_FILE_ID}.txt ${MAX_ALIGNMENT_ERROR} ${MIN_ALIGNMENT_LENGTH} ${PREFIX_2}${ALIGNMENTS_FILE_ID}.txt "${PREFIX_2}unique-${ALIGNMENTS_FILE_ID}.txt"
	sort --parallel=1 -t , ${SORT_OPTIONS} ${PREFIX_2}${ALIGNMENTS_FILE_ID}.txt | uniq - ${PREFIX_3}${ALIGNMENTS_FILE_ID}.txt
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
sort --parallel=${N_THREADS} -m -t , ${SORT_OPTIONS} ${TMPFILE_PATH}-3-*.txt | uniq - ${TMPFILE_PATH}-4.txt
N_INSTANCES=$(wc -l < ${TMPFILE_PATH}-4.txt)
java ${JAVA_RUNTIME_FLAGS} -classpath "${REVANT_BINARIES}" de.mpi_cbg.revant.apps.SplitCharacterInstances ${N_INSTANCES} ${N_THREADS} ${TMPFILE_PATH}-4.txt ${TMPFILE_PATH}-5-

echo "Compacting character instances..."
function compactionThread() {
	local INSTANCES_FILE_ID=$1
	local PREFIX_1=$2
	local PREFIX_2=$3
	local PREFIX_3=$4
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
	local ALIGNMENTS_FILE_ID=$1
	local LAST_TRANSLATED_READ=$2
	local PREFIX_1=$3
	local PREFIX_2=$4
	local PREFIX_3=$5
	local PREFIX_4=$6
	local PREFIX_5=$7
	local PREFIX_6=$8
	local IS_LAST_CHUNK=$9
	java ${JAVA_RUNTIME_FLAGS} -classpath "${REVANT_BINARIES}" de.mpi_cbg.revant.apps.TranslateReads ${N_READS} ${READ_IDS_FILE} ${READ_LENGTHS_FILE} ${N_REPEATS} ${REPEAT_LENGTHS_FILE} ${REPEAT_ISPERIODIC_FILE} ${PREFIX_1}${ALIGNMENTS_FILE_ID}.txt ${MAX_ALIGNMENT_ERROR} ${MIN_ALIGNMENT_LENGTH} ${ALPHABET_FILE} ${LAST_TRANSLATED_READ} ${PREFIX_2}${ALIGNMENTS_FILE_ID}.txt ${PREFIX_3}${ALIGNMENTS_FILE_ID}.txt ${PREFIX_4}${ALIGNMENTS_FILE_ID}.txt ${PREFIX_5}${ALIGNMENTS_FILE_ID}.txt ${PREFIX_6}${ALIGNMENTS_FILE_ID}.txt ${IS_LAST_CHUNK}
}
if [ -e ${TMPFILE_PATH}-1-${N_THREADS}.txt ]; then
	TO=${N_THREADS}
else
	TO=$(( ${N_THREADS} - 1 ))
fi
if [ ${TO} -ge 1 ]; then
	translationThread 0 -1 "${TMPFILE_PATH}-1-" "${TMPFILE_PATH}-8-" "${TMPFILE_PATH}-9-" "${TMPFILE_PATH}-9-hist-" "${TMPFILE_PATH}-10-" "${TMPFILE_PATH}-11-" 0 &
	if [ ${TO} -ge 2 ]; then
		for THREAD in $(seq 1 $((${TO} - 1))); do
			LAST_TRANSLATED_READ=$(tail -n 1 ${TMPFILE_PATH}-1-$(( ${THREAD} - 1 )).txt | awk '{ print $1 }' | tr -d , )
			LAST_TRANSLATED_READ=$(( ${LAST_TRANSLATED_READ} - 1 ))
			translationThread ${THREAD} ${LAST_TRANSLATED_READ} "${TMPFILE_PATH}-1-" "${TMPFILE_PATH}-8-" "${TMPFILE_PATH}-9-" "${TMPFILE_PATH}-9-hist-" "${TMPFILE_PATH}-10-" "${TMPFILE_PATH}-11-" 0 &
		done
	fi
	LAST_TRANSLATED_READ=$(tail -n 1 ${TMPFILE_PATH}-1-$((${TO} - 1)).txt | awk '{ print $1 }' | tr -d , )
	LAST_TRANSLATED_READ=$(( ${LAST_TRANSLATED_READ} - 1 ))
	translationThread ${TO} ${LAST_TRANSLATED_READ} "${TMPFILE_PATH}-1-" "${TMPFILE_PATH}-8-" "${TMPFILE_PATH}-9-" "${TMPFILE_PATH}-9-hist-" "${TMPFILE_PATH}-10-" "${TMPFILE_PATH}-11-" 1 &
else
	translationThread 0 -1 "${TMPFILE_PATH}-1-" "${TMPFILE_PATH}-8-" "${TMPFILE_PATH}-9-" "${TMPFILE_PATH}-9-hist-" "${TMPFILE_PATH}-10-" "${TMPFILE_PATH}-11-" 1 &
fi
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
if [ $((${TO} + 1)) -ge 2 ]; then
	for i in $(seq 2 $((${TO} + 1))); do
		PASTE_OPTIONS="${PASTE_OPTIONS}+\$$i"
	done
fi
PASTE_OPTIONS="${PASTE_OPTIONS}}'"
paste ${TMPFILE_PATH}-9-hist-*.txt | eval ${PASTE_OPTIONS} - > ${HISTOGRAM_FILE}
FULLY_UNIQUE_FILE="${INPUT_DIR}/reads-fullyUnique.txt"
rm -f ${FULLY_UNIQUE_FILE}
for THREAD in $(seq 0 ${TO}); do
	cat ${TMPFILE_PATH}-10-${THREAD}.txt >> ${FULLY_UNIQUE_FILE}
done
FULLY_CONTAINED_FILE="${INPUT_DIR}/reads-fullyContained.txt"
rm -f ${FULLY_CONTAINED_FILE}
for THREAD in $(seq 0 ${TO}); do
	cat ${TMPFILE_PATH}-11-${THREAD}.txt >> ${FULLY_CONTAINED_FILE}
done




echo "Fixing periodic endpoints..."
if [ ${MAX_SPACER_LENGTH} -ne 0 ]; then
	READ_READ_ALIGNMENTS_FILE="${INPUT_DIR}/LAshow-reads-reads.txt"
	echo "Splitting the alignments file..."
	LAST_READA_FILE="${INPUT_DIR}/LAshow-reads-reads-lastReadA.txt"
	if [ ${BROKEN_READS} -eq 1 ]; then
		# Reusing the chunks of the read-read alignments file that are already there (we
		# assume that they all have the header).
		for FILE in $(ls ${INPUT_DIR}/breakReads-tmp-2-*.txt ); do
			ID=$(basename ${FILE} .txt)
			ID=${ID#breakReads-tmp-2-}
			mv ${INPUT_DIR}/breakReads-tmp-2-${ID}.txt ${TMPFILE_PATH}-spacers-1-${ID}.txt
		done
	else
		N_ALIGNMENTS=$(( $(wc -l < ${READ_READ_ALIGNMENTS_FILE}) - 2 ))
		java ${JAVA_RUNTIME_FLAGS} -classpath "${REVANT_BINARIES}" de.mpi_cbg.revant.factorize.SplitAlignments ${N_ALIGNMENTS} ${N_THREADS} ${READ_READ_ALIGNMENTS_FILE} ${TMPFILE_PATH}-spacers-1- ${LAST_READA_FILE}
	fi
	echo "Read-read alignments filtered and split in ${N_THREADS} parts"
	echo "Collecting spacers and assigning breakpoints to them..."
	java ${JAVA_RUNTIME_FLAGS} -classpath "${REVANT_BINARIES}" de.mpi_cbg.revant.apps.FixPeriodicEndpoints1 ${MAX_SPACER_LENGTH} ${N_READS} ${READ_IDS_FILE} ${READ_LENGTHS_FILE} ${READS_TRANSLATED_FILE} ${READS_TRANSLATED_BOUNDARIES} ${READ_READ_ALIGNMENTS_FILE} ${N_THREADS} ${LAST_READA_FILE} ${TMPFILE_PATH}-spacers-2-
	if [ -e ${TMPFILE_PATH}-spacers-2-${N_THREADS}.txt ]; then
		TO=${N_THREADS}
	else
		TO=$(( ${N_THREADS} - 1 ))
	fi
	SPACERS_FILE="${TMPFILE_PATH}-spacers.txt"
	rm -f ${SPACERS_FILE}
	for i in $(seq 0 ${TO}); do
		cat ${TMPFILE_PATH}-spacers-2-${i}.txt >> ${SPACERS_FILE}
	done
	N_SPACERS=$(wc -l < ${SPACERS_FILE})
	echo "Collecting instances of spacer-induced characters..."
	function collectionThread_spacers() {
		local SPACERS_FILE_ID=$1
		local N_SPACERS = $(wc -l < ${SPACERS_FILE_ID})
		java ${JAVA_RUNTIME_FLAGS} -classpath "${REVANT_BINARIES}" de.mpi_cbg.revant.apps.FixPeriodicEndpoints2 ${MAX_SPACER_LENGTH} ${TMPFILE_PATH}-spacers-2-${SPACERS_FILE_ID}.txt ${N_SPACERS} ${ALPHABET_FILE} ${TMPFILE_PATH}-spacers-2-ids-${THREAD}.txt ${TMPFILE_PATH}-spacers-2-lengths-${THREAD}.txt ${TMPFILE_PATH}-8-${THREAD}.txt ${TMPFILE_PATH}-9-${THREAD}.txt ${TMPFILE_PATH}-spacers-3-${SPACERS_FILE_ID}.txt ${TMPFILE_PATH}-spacers-3-unique-${SPACERS_FILE_ID}.txt
		sort --parallel=1 -t , ${SORT_OPTIONS} ${TMPFILE_PATH}-spacers-3-${SPACERS_FILE_ID}.txt | uniq - ${TMPFILE_PATH}-spacers-4-${SPACERS_FILE_ID}.txt
	}
	for THREAD in $(seq 0 ${TO}); do
		collectionThread_spacers ${THREAD} &
	done
	wait
	sort --parallel=${N_THREADS} -m -t , ${SORT_OPTIONS} ${TMPFILE_PATH}-spacers-4-*.txt | uniq - ${TMPFILE_PATH}-spacers-5.txt
	N_INSTANCES=$(wc -l < ${TMPFILE_PATH}-spacers-5.txt)
	java ${JAVA_RUNTIME_FLAGS} -classpath "${REVANT_BINARIES}" de.mpi_cbg.revant.apps.SplitCharacterInstances ${N_INSTANCES} ${N_THREADS} ${TMPFILE_PATH}-spacers-5.txt ${TMPFILE_PATH}-spacers-6-
	echo "Compacting character instances..."
	for THREAD in $(seq 0 ${TO}); do
		compactionThread ${THREAD} "${TMPFILE_PATH}-spacers-6-" "${TMPFILE_PATH}-spacers-7-" "${TMPFILE_PATH}-spacers-8-" &
	done
	wait
	ALPHABET_FILE_SPACERS="${INPUT_DIR}/alphabet-spacers.txt"
	rm -f ${ALPHABET_FILE_SPACERS}
	java ${JAVA_RUNTIME_FLAGS} -classpath "${REVANT_BINARIES}" de.mpi_cbg.revant.apps.MergeAlphabetHeaders ${INPUT_DIR} "${TMPFILE_NAME}-spacers-8-" "${TMPFILE_NAME}-spacers-3-unique-" ${ALPHABET_FILE_SPACERS}
	for THREAD in $(seq 0 ${TO}); do
		tail -n +2 ${TMPFILE_PATH}-spacers-8-${THREAD}.txt >> ${ALPHABET_FILE_SPACERS}
	done
	ALPHABET_SIZE_SPACERS=$( wc -l < ${ALPHABET_FILE_SPACERS} )
	ALPHABET_SIZE_SPACERS=$(( ${ALPHABET_SIZE_SPACERS} - 1 ))
	echo "Translating reads into the spacers-induced alphabet..."
	function translationThread_spacers() {
		local THREAD_ID=$1
		local READ2CHARACTERS_FILE_NEW=$2
		local READ2BOUNDARIES_FILE_NEW=$3
		java ${JAVA_RUNTIME_FLAGS} -classpath "${REVANT_BINARIES}" de.mpi_cbg.revant.apps.FixPeriodicEndpoints3 ${MAX_SPACER_LENGTH}  --------> make it work on chunk    ${SPACERS_FILE} ${N_SPACERS} ${ALPHABET_FILE} ${ALPHABET_FILE_SPACERS} ${READ_IDS_FILE} ${READ_LENGTHS_FILE} ${TMPFILE_PATH}-8-${THREAD_ID}.txt ${TMPFILE_PATH}-9-${THREAD_ID}.txt ${READ2CHARACTERS_FILE_NEW} ${READ2BOUNDARIES_FILE_NEW}	
	}
	for THREAD in $(seq 0 ${TO}); do
		translationThread_spacers ${THREAD} "${TMPFILE_PATH}-spacers-9-" "${TMPFILE_PATH}-spacers-10-" &
	done
	wait
	rm -f ${READS_TRANSLATED_FILE} ${READS_TRANSLATED_BOUNDARIES}
	for THREAD in $(seq 0 ${TO}); do
		cat ${TMPFILE_PATH}-spacers-9-${THREAD}.txt >> ${READS_TRANSLATED_FILE}
		cat ${TMPFILE_PATH}-spacers-10-${THREAD}.txt >> ${READS_TRANSLATED_BOUNDARIES}
	done
fi
echo "Periodic endpoints fixed"




echo "Discarding rare characters..."
COUNTS_FILE="${INPUT_DIR}/alphabet-counts.txt"
HISTOGRAM_FILE="${INPUT_DIR}/alphabet-histogram.txt"
java ${JAVA_RUNTIME_FLAGS} -classpath "${REVANT_BINARIES}" de.mpi_cbg.revant.apps.GetCharacterCounts ${READS_TRANSLATED_FILE} ${ALPHABET_FILE} ${COUNTS_FILE} ${HISTOGRAM_FILE}
function cleaningThread1() {
	local TRANSLATED_CHARACTERS=$1
	local TRANSLATED_BOUNDARIES=$2
	local PREFIX_1=$3
	local PREFIX_2=$4
	local ID=$5
	java ${JAVA_RUNTIME_FLAGS} -classpath "${REVANT_BINARIES}" de.mpi_cbg.revant.apps.CleanTranslatedReads1 ${ALPHABET_FILE} ${COUNTS_FILE} ${N_READS} ${READ_IDS_FILE} ${READ_LENGTHS_FILE} ${TRANSLATED_CHARACTERS} ${TRANSLATED_BOUNDARIES} ${MIN_CHARACTER_FREQUENCY} ${KEEP_PERIODIC} ${PREFIX_1}${ID}.txt > ${PREFIX_1}unique-${ID}.txt
	sort --parallel=1 -t , ${SORT_OPTIONS} ${PREFIX_1}${ID}.txt | uniq - ${PREFIX_2}${ID}.txt
}
split -l $(( ${N_READS} / ${N_THREADS} )) ${READS_TRANSLATED_FILE} "${TMPFILE_PATH}-12-"
split -l $(( ${N_READS} / ${N_THREADS} )) ${READS_TRANSLATED_BOUNDARIES} "${TMPFILE_PATH}-13-"
for FILE in $(find ${INPUT_DIR} -name "${TMPFILE_NAME}-12-*" ); do
	THREAD_ID=${FILE#${INPUT_DIR}/${TMPFILE_NAME}-12-}
	cleaningThread1 "${INPUT_DIR}/${TMPFILE_NAME}-12-${THREAD_ID}" "${INPUT_DIR}/${TMPFILE_NAME}-13-${THREAD_ID}" "${TMPFILE_PATH}-14-" "${TMPFILE_PATH}-15-" ${THREAD_ID} &
done
wait
sort --parallel=${N_THREADS} -m -t , ${SORT_OPTIONS} ${TMPFILE_PATH}-15-*.txt | uniq - ${TMPFILE_PATH}-15.txt
cat ${TMPFILE_PATH}-14-unique-*.txt | sort --parallel=${N_THREADS} -n -r > ${TMPFILE_PATH}-14-unique.txt
ALPHABET_FILE_CLEANED="${INPUT_DIR}/alphabet-cleaned.txt"
rm -f ${ALPHABET_FILE_CLEANED}
OLD2NEW_FILE="${INPUT_DIR}/alphabet-old2new.txt"
rm -f ${OLD2NEW_FILE}
java ${JAVA_RUNTIME_FLAGS} -classpath "${REVANT_BINARIES}" de.mpi_cbg.revant.apps.CleanTranslatedReads2 ${ALPHABET_FILE} ${COUNTS_FILE} $(wc -l < ${TMPFILE_PATH}-15.txt) ${TMPFILE_PATH}-15.txt ${MIN_CHARACTER_FREQUENCY} ${KEEP_PERIODIC} ${TMPFILE_PATH}-14-unique.txt ${ALPHABET_FILE_CLEANED} ${OLD2NEW_FILE}
function cleaningThread3() {
	local TRANSLATED_CHARACTERS_OLD=$1
	local TRANSLATED_BOUNDARIES_OLD=$2
	local TRANSLATED_CHARACTERS_NEW=$3
	local TRANSLATED_BOUNDARIES_NEW=$4
	local HISTOGRAM_FILE=$5
	local FULLY_UNIQUE_FILE_NEW=$6
	local LAST_TRANSLATED_READ=$7
	java ${JAVA_RUNTIME_FLAGS} -classpath "${REVANT_BINARIES}" de.mpi_cbg.revant.apps.CleanTranslatedReads3 ${N_READS} ${READ_IDS_FILE} ${READ_LENGTHS_FILE} ${ALPHABET_FILE} ${COUNTS_FILE} ${TRANSLATED_CHARACTERS_OLD} ${TRANSLATED_BOUNDARIES_OLD} ${MIN_CHARACTER_FREQUENCY} ${KEEP_PERIODIC} ${ALPHABET_FILE_CLEANED} ${OLD2NEW_FILE} ${TRANSLATED_CHARACTERS_NEW} ${TRANSLATED_BOUNDARIES_NEW} ${HISTOGRAM_FILE} ${FULLY_UNIQUE_FILE_NEW} ${LAST_TRANSLATED_READ}
}
LAST_TRANSLATED_READ="-1"
for FILE in $(find -s ${INPUT_DIR} -name "${TMPFILE_NAME}-12-*" ); do
	THREAD_ID=${FILE#${INPUT_DIR}/${TMPFILE_NAME}-12-}
	cleaningThread3 "${INPUT_DIR}/${TMPFILE_NAME}-12-${THREAD_ID}" "${INPUT_DIR}/${TMPFILE_NAME}-13-${THREAD_ID}" "${INPUT_DIR}/${TMPFILE_NAME}-16-${THREAD_ID}" "${INPUT_DIR}/${TMPFILE_NAME}-17-${THREAD_ID}" "${INPUT_DIR}/${TMPFILE_NAME}-18-${THREAD_ID}" "${INPUT_DIR}/${TMPFILE_NAME}-19-${THREAD_ID}" ${LAST_TRANSLATED_READ} &
	LAST_TRANSLATED_READ=$(( ${LAST_TRANSLATED_READ} + $(wc -l < "${INPUT_DIR}/${TMPFILE_NAME}-12-${THREAD_ID}") ))
done
wait
READS_TRANSLATED_FILE_NEW="${INPUT_DIR}/reads-translated-new.txt"
READS_TRANSLATED_BOUNDARIES_NEW="${INPUT_DIR}/reads-translated-boundaries-new.txt"
HISTOGRAM_FILE="${INPUT_DIR}/reads-translated-histogram-new.txt"
FULLY_UNIQUE_FILE_NEW="${INPUT_DIR}/reads-fullyUnique-new.txt"
FULLY_CONTAINED_FILE_NEW="${INPUT_DIR}/reads-fullyContained-new.txt"
rm -f ${READS_TRANSLATED_FILE_NEW} ${READS_TRANSLATED_BOUNDARIES_NEW} ${FULLY_UNIQUE_FILE_NEW} ${FULLY_CONTAINED_FILE_NEW}
for FILE in $(find -s ${INPUT_DIR} -name "${TMPFILE_NAME}-16-*" ); do
	THREAD_ID=${FILE#${INPUT_DIR}/${TMPFILE_NAME}-16-}
	cat ${INPUT_DIR}/${TMPFILE_NAME}-16-${THREAD_ID} >> ${READS_TRANSLATED_FILE_NEW}
	cat ${INPUT_DIR}/${TMPFILE_NAME}-17-${THREAD_ID} >> ${READS_TRANSLATED_BOUNDARIES_NEW}
	cat ${INPUT_DIR}/${TMPFILE_NAME}-19-${THREAD_ID} >> ${FULLY_UNIQUE_FILE_NEW}
done
PASTE_OPTIONS="awk '{print \$1"
if [ ${N_THREADS} -ge 2 ]; then
	for i in $(seq 2 ${N_THREADS}); do
		PASTE_OPTIONS="${PASTE_OPTIONS}+\$$i"
	done
fi
PASTE_OPTIONS="${PASTE_OPTIONS}}'"
paste ${INPUT_DIR}/${TMPFILE_NAME}-18-* | eval ${PASTE_OPTIONS} - > ${HISTOGRAM_FILE}
java ${JAVA_RUNTIME_FLAGS} -classpath "${REVANT_BINARIES}" de.mpi_cbg.revant.apps.Setminus ${FULLY_CONTAINED_FILE} ${FULLY_UNIQUE_FILE_NEW} > ${FULLY_CONTAINED_FILE_NEW}

# Removing all temp files that are not used downstream
if [ ${DELETE_TMP_FILES} -eq 1 ]; then
	rm -f ${TMPFILE_PATH}*
fi
