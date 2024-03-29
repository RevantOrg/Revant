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
# ----------------------------------- CONFIGURATION --------------------------------------
# This is the only section of the script that needs to be customized.
#
INPUT_DIR=$1
BROKEN_READS=$2  # 1=TRUE
MAX_ALIGNMENT_ERROR=$3  # Repeat-read alignments with error > this are discarded
MIN_ALIGNMENT_LENGTH=$4  # Repeat-read alignments with length < this are discarded
N_HAPLOTYPES=$5
N_THREADS=$6
DELETE_TMP_FILES=$7
MAX_SPACER_LENGTH=$8  # 0=assume that the endpoints of periodic repeats are accurate
WOBBLE_LENGTH=$9  # 0=do not wobble
# Good settings for mostly periodic genome: MAX_SPACER_LENGTH="10000"; WOBBLE_LENGTH="100"
# If spacers are real in the read translations (i.e. if they do not come from
# idiosyncrasies of the aligner) the translations are automatically untouched.
# WOBBLE_LENGTH values have to be chosen experimentally. Smaller values might strike a
# better balance between increasing frequency and not making truly unique k-mers too
# frequent as to be classified as repeats.
FIX_TANDEM_SPACERS=${10}  # 0=assume that non-repetitive blocks near tandems are real.
CONCATENATE_BLOCKS=${11}  # 0=do not try to merge adjacent blocks from the same repeat.
AVG_READ_LENGTH=${12}
GENOME_LENGTH=${13}  # Of one haplotype
KEEP_PERIODIC="1"  # 1=do not remove rare characters if they are periodic. Usually good.
SPANNING_BPS="150"  # Bps before and after a character to consider it observed in a read.
# ---------------------------------- TANDEM SPACERS --------------------------------------
TANDEM_SPACERS_ITERATIONS="1"  # >=1
NONREPETITIVE_BLOCKS_MODE="2"  # Should be probably set to 1 in production, 2 is the most aggressive.
CONCATENATE_THRESHOLD="200"
LONG_SPACER_LENGTH=$(( ${AVG_READ_LENGTH} / 3 ))  # Arbitrary
# ------------------------------------ REVANT --------------------------------------------
REVANT_LIBRARIES="${REVANT_BINARIES}/../lib/*"
# ----------------------------------------------------------------------------------------

set -euo pipefail
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
rm -f ${TMPFILE_PATH}*

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
	sort --parallel=1 -t , -u ${SORT_OPTIONS} ${PREFIX_2}${ALIGNMENTS_FILE_ID}.txt > ${PREFIX_3}${ALIGNMENTS_FILE_ID}.txt
}
if [ -e ${TMPFILE_PATH}-1-${N_THREADS}.txt ]; then
	TO=${N_THREADS}
else
	TO=$(( ${N_THREADS} - 1 ))
fi
PIDS=()
for THREAD in $(seq 0 ${TO}); do
	collectionThread ${THREAD} "${TMPFILE_PATH}-1-" "${TMPFILE_PATH}-2-" "${TMPFILE_PATH}-3-" &
    PIDS+=($!)
done
waitAndCheck PIDS
sort --parallel=${N_THREADS} -m -t , -u ${SORT_OPTIONS} ${TMPFILE_PATH}-3-*.txt > ${TMPFILE_PATH}-4.txt
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
PIDS=()
for THREAD in $(seq 0 ${TO}); do
	compactionThread ${THREAD} "${TMPFILE_PATH}-5-" "${TMPFILE_PATH}-6-" "${TMPFILE_PATH}-7-" &
    PIDS+=($!)
done
waitAndCheck PIDS
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
PIDS=()
if [ ${TO} -ge 1 ]; then
	translationThread 0 -1 "${TMPFILE_PATH}-1-" "${TMPFILE_PATH}-8-" "${TMPFILE_PATH}-9-" "${TMPFILE_PATH}-9-hist-" "${TMPFILE_PATH}-10-" "${TMPFILE_PATH}-11-" 0 &
    PIDS+=($!)
	if [ ${TO} -ge 2 ]; then
		for THREAD in $(seq 1 $((${TO} - 1))); do
			LAST_TRANSLATED_READ=$(tail -n 1 ${TMPFILE_PATH}-1-$(( ${THREAD} - 1 )).txt | awk '{ print $1 }' | tr -d , )
			LAST_TRANSLATED_READ=$(( ${LAST_TRANSLATED_READ} - 1 ))
			translationThread ${THREAD} ${LAST_TRANSLATED_READ} "${TMPFILE_PATH}-1-" "${TMPFILE_PATH}-8-" "${TMPFILE_PATH}-9-" "${TMPFILE_PATH}-9-hist-" "${TMPFILE_PATH}-10-" "${TMPFILE_PATH}-11-" 0 &
            PIDS+=($!)
		done
	fi
	LAST_TRANSLATED_READ=$(tail -n 1 ${TMPFILE_PATH}-1-$((${TO} - 1)).txt | awk '{ print $1 }' | tr -d , )
	LAST_TRANSLATED_READ=$(( ${LAST_TRANSLATED_READ} - 1 ))
	translationThread ${TO} ${LAST_TRANSLATED_READ} "${TMPFILE_PATH}-1-" "${TMPFILE_PATH}-8-" "${TMPFILE_PATH}-9-" "${TMPFILE_PATH}-9-hist-" "${TMPFILE_PATH}-10-" "${TMPFILE_PATH}-11-" 1 &
    PIDS+=($!)
else
	translationThread 0 -1 "${TMPFILE_PATH}-1-" "${TMPFILE_PATH}-8-" "${TMPFILE_PATH}-9-" "${TMPFILE_PATH}-9-hist-" "${TMPFILE_PATH}-10-" "${TMPFILE_PATH}-11-" 1 &
    PIDS+=($!)
fi
waitAndCheck PIDS
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




echo "Trying to fix tandem spacers if needed..."
if [ ${TANDEM_SPACERS_ITERATIONS} -gt 0 ]; then
	java ${JAVA_RUNTIME_FLAGS} -classpath "${REVANT_BINARIES}" de.mpi_cbg.revant.apps.SplitSpacers ${LAST_READA_FILE} ${N_THREADS} null ${READ_IDS_FILE} ${READ_LENGTHS_FILE} ${TMPFILE_PATH}-stash-
fi
function tandemsThread() {
	local LOCAL_TRANSLATED_READS_FILE=$1
	local LOCAL_BOUNDARIES_FILE=$2
	local LOCAL_READ_LENGTHS_FILE=$3
	local LOCAL_TANDEMS_FILE=$4
	java ${JAVA_RUNTIME_FLAGS} -classpath "${REVANT_BINARIES}" de.mpi_cbg.revant.apps.CollectTandems 0 1 2 ${ALPHABET_FILE} ${LOCAL_TRANSLATED_READS_FILE} ${LOCAL_BOUNDARIES_FILE} ${LOCAL_READ_LENGTHS_FILE} ${REPEAT_LENGTHS_FILE} ${N_REPEATS} ${LOCAL_TANDEMS_FILE}
}
ITER="1";
while [ ${ITER} -le ${TANDEM_SPACERS_ITERATIONS} ]; do
	rm -f ${TMPFILE_PATH}-tspacers-*
	for FILE in $(find -s ${INPUT_DIR} -maxdepth 1 -name "${TMPFILE_NAME}-stash-*"); do
		SUFFIX=${FILE#${TMPFILE_PATH}-stash-}
		cp ${FILE} ${TMPFILE_PATH}-tspacers-1-${SUFFIX}
	done
	echo "Computing tandem track..."
    PIDS=()
	for FILE in $(find -s ${INPUT_DIR} -maxdepth 1 -name "${TMPFILE_NAME}-8-*.txt"); do
		THREAD_ID=${FILE#${TMPFILE_PATH}-8-}
		THREAD_ID=${THREAD_ID%.txt}
		tandemsThread ${FILE} ${TMPFILE_PATH}-9-${THREAD_ID}.txt ${TMPFILE_PATH}-tspacers-1-lengths-${THREAD_ID}.txt ${TMPFILE_PATH}-tspacers-2-${THREAD_ID}.txt &
        PIDS+=($!)
	done
	waitAndCheck PIDS
	TANDEMS_FILE="${TMPFILE_PATH}-tspacers-3.txt"
	rm -f ${TANDEMS_FILE}
	for FILE in $(find -s ${INPUT_DIR} -maxdepth 1 -name "${TMPFILE_NAME}-tspacers-2-*.txt"); do
		cat ${FILE} >> ${TANDEMS_FILE}
	done
	echo "Collecting tandem spacers and assigning solutions to them..."
	N_FULLY_UNIQUE=$(wc -l < ${FULLY_UNIQUE_FILE})
	N_FULLY_CONTAINED=$(wc -l < ${FULLY_CONTAINED_FILE})
	READ_READ_ALIGNMENTS_FILE="${INPUT_DIR}/LAshow-reads-reads.txt"
	SPACERS_ARE_CORRECT=$(java ${JAVA_RUNTIME_FLAGS} -classpath "${REVANT_BINARIES}" de.mpi_cbg.revant.apps.FixTandemSpacers1 ${N_READS} ${READ_IDS_FILE} ${READ_LENGTHS_FILE} ${ALPHABET_FILE} ${READS_TRANSLATED_FILE} ${READS_TRANSLATED_BOUNDARIES} ${FULLY_UNIQUE_FILE} ${N_FULLY_UNIQUE} ${FULLY_CONTAINED_FILE} ${N_FULLY_CONTAINED} ${READ_READ_ALIGNMENTS_FILE} ${TANDEMS_FILE} ${NONREPETITIVE_BLOCKS_MODE} ${LONG_SPACER_LENGTH} ${TMPFILE_PATH}-tspacers-4.txt)
	if [ ${SPACERS_ARE_CORRECT} -ne 0 ]; then
		break
	else
		java ${JAVA_RUNTIME_FLAGS} -classpath "${REVANT_BINARIES}" de.mpi_cbg.revant.apps.SplitSpacers ${LAST_READA_FILE} ${N_THREADS} ${TMPFILE_PATH}-tspacers-4.txt ${READ_IDS_FILE} ${READ_LENGTHS_FILE} ${TMPFILE_PATH}-tspacers-5-
		echo "Collecting instances of characters induced by solved tandem spacers..."
		function collectionThread_tspacers() {
			local SPACERS_FILE_ID=$1
			local LOCAL_N_SPACERS=$(wc -l < ${TMPFILE_PATH}-tspacers-5-${SPACERS_FILE_ID}.txt)
			java ${JAVA_RUNTIME_FLAGS} -classpath "${REVANT_BINARIES}" de.mpi_cbg.revant.apps.FixTandemSpacers2 ${TMPFILE_PATH}-tspacers-5-${SPACERS_FILE_ID}.txt ${LOCAL_N_SPACERS} ${ALPHABET_FILE} ${TMPFILE_PATH}-tspacers-5-ids-${SPACERS_FILE_ID}.txt ${TMPFILE_PATH}-tspacers-5-lengths-${SPACERS_FILE_ID}.txt ${TMPFILE_PATH}-8-${SPACERS_FILE_ID}.txt ${TMPFILE_PATH}-9-${SPACERS_FILE_ID}.txt ${TMPFILE_PATH}-tspacers-2-${SPACERS_FILE_ID}.txt ${REPEAT_LENGTHS_FILE} ${N_REPEATS} ${NONREPETITIVE_BLOCKS_MODE} ${TMPFILE_PATH}-tspacers-6-${SPACERS_FILE_ID}.txt
			sort --parallel=1 -t , -u ${SORT_OPTIONS} ${TMPFILE_PATH}-tspacers-6-${SPACERS_FILE_ID}.txt > ${TMPFILE_PATH}-tspacers-7-${SPACERS_FILE_ID}.txt
		}
		if [ -e ${TMPFILE_PATH}-tspacers-5-${N_THREADS}.txt ]; then
			TO=${N_THREADS}
		else
			TO=$(( ${N_THREADS} - 1 ))
		fi
        PIDS=()
		for THREAD in $(seq 0 ${TO}); do
			collectionThread_tspacers ${THREAD} &
            PIDS+=($!)
		done
		waitAndCheck PIDS
		sort --parallel=${N_THREADS} -m -t , -u ${SORT_OPTIONS} ${TMPFILE_PATH}-tspacers-7-*.txt > ${TMPFILE_PATH}-tspacers-8.txt
		N_INSTANCES=$(wc -l < ${TMPFILE_PATH}-tspacers-8.txt)
		java ${JAVA_RUNTIME_FLAGS} -classpath "${REVANT_BINARIES}" de.mpi_cbg.revant.apps.SplitCharacterInstances ${N_INSTANCES} ${N_THREADS} ${TMPFILE_PATH}-tspacers-8.txt ${TMPFILE_PATH}-tspacers-9-
		echo "Compacting character instances..."
        PIDS=()
		for THREAD in $(seq 0 ${TO}); do
			compactionThread ${THREAD} "${TMPFILE_PATH}-tspacers-9-" "${TMPFILE_PATH}-tspacers-10-" "${TMPFILE_PATH}-tspacers-11-" &
            PIDS+=($!)
		done
		waitAndCheck PIDS
		head -n 1 ${ALPHABET_FILE} | cut -d , -f 4 > "${TMPFILE_PATH}-tspacers-unique.txt"
		ALPHABET_FILE_SPACERS="${INPUT_DIR}/alphabet-tspacers.txt"
		rm -f ${ALPHABET_FILE_SPACERS}
		java ${JAVA_RUNTIME_FLAGS} -classpath "${REVANT_BINARIES}" de.mpi_cbg.revant.apps.MergeAlphabetHeaders ${INPUT_DIR} "${TMPFILE_NAME}-tspacers-11-" "${TMPFILE_PATH}-tspacers-unique" ${ALPHABET_FILE_SPACERS}
		for THREAD in $(seq 0 ${TO}); do
			tail -n +2 ${TMPFILE_PATH}-tspacers-11-${THREAD}.txt >> ${ALPHABET_FILE_SPACERS}
		done
		ALPHABET_SIZE_SPACERS=$( wc -l < ${ALPHABET_FILE_SPACERS} )
		ALPHABET_SIZE_SPACERS=$(( ${ALPHABET_SIZE_SPACERS} - 1 ))
		echo "Translating reads into the alphabet induced by solved tandem spacers..."
		function translationThread_tspacers() {
			local THREAD_ID=$1
			local READ2CHARACTERS_FILE_NEW=$2
			local READ2BOUNDARIES_FILE_NEW=$3
			local LOCAL_SPACERS_FILE="${TMPFILE_PATH}-tspacers-5-${THREAD_ID}.txt"
			local LOCAL_N_SPACERS=$(wc -l < ${LOCAL_SPACERS_FILE})
			java ${JAVA_RUNTIME_FLAGS} -classpath "${REVANT_BINARIES}" de.mpi_cbg.revant.apps.FixTandemSpacers3 ${LOCAL_SPACERS_FILE} ${LOCAL_N_SPACERS} ${ALPHABET_FILE} ${ALPHABET_FILE_SPACERS} ${TMPFILE_PATH}-tspacers-5-ids-${THREAD_ID}.txt ${TMPFILE_PATH}-tspacers-5-lengths-${THREAD_ID}.txt ${TMPFILE_PATH}-8-${THREAD_ID}.txt ${TMPFILE_PATH}-9-${THREAD_ID}.txt ${READ2CHARACTERS_FILE_NEW}${THREAD_ID}.txt ${READ2BOUNDARIES_FILE_NEW}${THREAD_ID}.txt ${TMPFILE_PATH}-tspacers-2-${THREAD_ID}.txt ${REPEAT_LENGTHS_FILE} ${N_REPEATS} ${NONREPETITIVE_BLOCKS_MODE}
		}
        PIDS=()
		for THREAD in $(seq 0 ${TO}); do
			translationThread_tspacers ${THREAD} "${TMPFILE_PATH}-tspacers-12-" "${TMPFILE_PATH}-tspacers-13-" &
            PIDS+=($!)
		done
		waitAndCheck PIDS
		mv ${READS_TRANSLATED_FILE} ${READS_TRANSLATED_FILE}-preTspacers
		mv ${READS_TRANSLATED_BOUNDARIES} ${READS_TRANSLATED_BOUNDARIES}-preTspacers
		mv ${ALPHABET_FILE} ${ALPHABET_FILE}-preTspacers
        mv ${FULLY_UNIQUE_FILE} ${FULLY_UNIQUE_FILE}-preTspacers
		for THREAD in $(seq 0 ${TO}); do
			cat ${TMPFILE_PATH}-tspacers-12-${THREAD}.txt >> ${READS_TRANSLATED_FILE}
			cat ${TMPFILE_PATH}-tspacers-13-${THREAD}.txt >> ${READS_TRANSLATED_BOUNDARIES}
		done
        touch ${FULLY_UNIQUE_FILE}
        for i in $(sed -n '/^$/=' ${READS_TRANSLATED_FILE}); do
            echo $(( $i - 1 )) >> ${FULLY_UNIQUE_FILE}
        done
		mv ${ALPHABET_FILE_SPACERS} ${ALPHABET_FILE}
		TANDEM_SPACERS_FIXED="1"
		echo "Tandem spacers fixed"
		
		#exit
		
		rm -f ${TMPFILE_PATH}-concatenate-*
		for FILE in $(find -s ${INPUT_DIR} -maxdepth 1 -name "${TMPFILE_NAME}-stash-*"); do
			SUFFIX=${FILE#${TMPFILE_PATH}-stash-}
			cp ${FILE} ${TMPFILE_PATH}-concatenate-1-${SUFFIX}
		done
		echo "Concatenating adjacent blocks from the same non-periodic repeat"
		java ${JAVA_RUNTIME_FLAGS} -classpath "${REVANT_BINARIES}" de.mpi_cbg.revant.apps.SplitTranslations ${READ_IDS_FILE} ${READS_TRANSLATED_FILE} ${READS_TRANSLATED_BOUNDARIES} ${LAST_READA_FILE} ${TMPFILE_PATH}-concatenate-2- ${TMPFILE_PATH}-concatenate-3-
		function concatenateThread() {
			local THREAD_ID=$1
			java ${JAVA_RUNTIME_FLAGS} -classpath "${REVANT_BINARIES}" de.mpi_cbg.revant.apps.ConcatenateBlocks1 ${ALPHABET_FILE} ${TMPFILE_PATH}-concatenate-1-lengths-${THREAD_ID}.txt ${TMPFILE_PATH}-concatenate-2-${THREAD_ID}.txt ${TMPFILE_PATH}-concatenate-3-${THREAD_ID}.txt ${REPEAT_LENGTHS_FILE} ${N_REPEATS} ${CONCATENATE_THRESHOLD} ${TMPFILE_PATH}-concatenate-4-${THREAD_ID}.txt
			sort --parallel=1 -t , -u ${SORT_OPTIONS} ${TMPFILE_PATH}-concatenate-4-${THREAD_ID}.txt > ${TMPFILE_PATH}-concatenate-5-${THREAD_ID}.txt
			rm -f ${TMPFILE_PATH}-concatenate-4-${THREAD_ID}.txt
		}
		echo "Collecting new characters from concatenation..."
		if [ -e ${TMPFILE_PATH}-concatenate-1-lengths-${N_THREADS}.txt ]; then
			TO=${N_THREADS}
		else
			TO=$(( ${N_THREADS} - 1 ))
		fi
        PIDS=()
		for THREAD in $(seq 0 ${TO}); do
			concatenateThread ${THREAD} &
            PIDS+=($!)
		done
		waitAndCheck PIDS
		sort --parallel=${N_THREADS} -m -t , -u ${SORT_OPTIONS} ${TMPFILE_PATH}-concatenate-5-*.txt > ${TMPFILE_PATH}-concatenate-6.txt
		N_INSTANCES=$(wc -l < ${TMPFILE_PATH}-concatenate-6.txt)
		java ${JAVA_RUNTIME_FLAGS} -classpath "${REVANT_BINARIES}" de.mpi_cbg.revant.apps.SplitCharacterInstances ${N_INSTANCES} ${N_THREADS} ${TMPFILE_PATH}-concatenate-6.txt ${TMPFILE_PATH}-concatenate-7-
		echo "Compacting character instances..."
        PIDS=()
		for THREAD in $(seq 0 ${TO}); do
			compactionThread ${THREAD} "${TMPFILE_PATH}-concatenate-7-" "${TMPFILE_PATH}-concatenate-8-" "${TMPFILE_PATH}-concatenate-9-" &
            PIDS+=($!)
		done
		waitAndCheck PIDS
		head -n 1 ${ALPHABET_FILE} | cut -d , -f 4 > "${TMPFILE_PATH}-concatenate-unique.txt"
		ALPHABET_FILE_CONCATENATE="${INPUT_DIR}/alphabet-concatenate.txt"
		rm -f ${ALPHABET_FILE_CONCATENATE}
		java ${JAVA_RUNTIME_FLAGS} -classpath "${REVANT_BINARIES}" de.mpi_cbg.revant.apps.MergeAlphabetHeaders ${INPUT_DIR} "${TMPFILE_NAME}-concatenate-9-" "${TMPFILE_PATH}-concatenate-unique" ${ALPHABET_FILE_CONCATENATE}
		for THREAD in $(seq 0 ${TO}); do
			tail -n +2 ${TMPFILE_PATH}-concatenate-9-${THREAD}.txt >> ${ALPHABET_FILE_CONCATENATE}
		done
		ALPHABET_SIZE_CONCATENATE=$( wc -l < ${ALPHABET_FILE_CONCATENATE} )
		ALPHABET_SIZE_CONCATENATE=$(( ${ALPHABET_SIZE_CONCATENATE} - 1 ))
		echo "Translating reads into the alphabet induced by concatenation..."
		function translationThread_concatenate() {
			local THREAD_ID=$1
			java ${JAVA_RUNTIME_FLAGS} -classpath "${REVANT_BINARIES}" de.mpi_cbg.revant.apps.ConcatenateBlocks2 ${ALPHABET_FILE} ${ALPHABET_FILE_CONCATENATE} ${TMPFILE_PATH}-concatenate-1-lengths-${THREAD_ID}.txt ${REPEAT_LENGTHS_FILE} ${N_REPEATS} ${TMPFILE_PATH}-concatenate-2-${THREAD_ID}.txt ${TMPFILE_PATH}-concatenate-3-${THREAD_ID}.txt ${CONCATENATE_THRESHOLD} ${TMPFILE_PATH}-concatenate-10-${THREAD_ID}.txt ${TMPFILE_PATH}-concatenate-11-${THREAD_ID}.txt ${TMPFILE_PATH}-concatenate-12-${THREAD_ID}.txt
		}
        PIDS=()
		for THREAD in $(seq 0 ${TO}); do
			translationThread_concatenate ${THREAD} &
            PIDS+=($!)
		done
		waitAndCheck PIDS
		mv ${READS_TRANSLATED_FILE} ${READS_TRANSLATED_FILE}-preConcatenation
		mv ${READS_TRANSLATED_BOUNDARIES} ${READS_TRANSLATED_BOUNDARIES}-preConcatenation
		mv ${FULLY_CONTAINED_FILE} ${FULLY_CONTAINED_FILE}-preConcatenation
		mv ${ALPHABET_FILE} ${ALPHABET_FILE}-preConcatenation
		for THREAD in $(seq 0 ${TO}); do
			cat ${TMPFILE_PATH}-concatenate-10-${THREAD}.txt >> ${READS_TRANSLATED_FILE}
			cat ${TMPFILE_PATH}-concatenate-11-${THREAD}.txt >> ${READS_TRANSLATED_BOUNDARIES}
			cat ${TMPFILE_PATH}-concatenate-12-${THREAD}.txt >> ${FULLY_CONTAINED_FILE}
		done
		mv ${ALPHABET_FILE_CONCATENATE} ${ALPHABET_FILE}
		echo "Concatenation completed"
		
		# Next iteration
		for THREAD in $(seq 0 ${TO}); do
			mv ${TMPFILE_PATH}-concatenate-10-${THREAD}.txt ${TMPFILE_PATH}-8-${THREAD}.txt
			mv ${TMPFILE_PATH}-concatenate-11-${THREAD}.txt ${TMPFILE_PATH}-9-${THREAD}.txt
		done
		ITER=$(( ${ITER} + 1 ))
	fi
done
echo "Tandem spacers fixed"
if [ ${WOBBLE_LENGTH} -ne 0 ]; then
	echo "Wobbling..."
	java ${JAVA_RUNTIME_FLAGS} -classpath "${REVANT_BINARIES}" de.mpi_cbg.revant.apps.SplitTranslations ${READ_IDS_FILE} ${READS_TRANSLATED_FILE} ${READS_TRANSLATED_BOUNDARIES} ${LAST_READA_FILE} ${TMPFILE_PATH}-wobble-1- ${TMPFILE_PATH}-wobble-2-
	java ${JAVA_RUNTIME_FLAGS} -classpath "${REVANT_BINARIES}" de.mpi_cbg.revant.apps.SplitSpacers ${LAST_READA_FILE} ${N_THREADS} null ${READ_IDS_FILE} ${READ_LENGTHS_FILE} ${TMPFILE_PATH}-wobble-3-
	if [ -e ${TMPFILE_PATH}-wobble-1-${N_THREADS}.txt ]; then
		TO=${N_THREADS}
	else
		TO=$(( ${N_THREADS} - 1 ))
	fi
	echo "Computing tandem track..."
    PIDS=()
	for THREAD_ID in $(seq 0 ${TO}); do		 
		tandemsThread ${TMPFILE_PATH}-wobble-1-${THREAD_ID}.txt ${TMPFILE_PATH}-wobble-2-${THREAD_ID}.txt ${TMPFILE_PATH}-wobble-3-lengths-${THREAD_ID}.txt ${TMPFILE_PATH}-wobble-4-${THREAD_ID}.txt &
        PIDS+=($!)
	done
	waitAndCheck PIDS
	function wobbleThreadCreate() {
		local WOBBLE_FILE_ID=$1
		java ${JAVA_RUNTIME_FLAGS} -classpath "${REVANT_BINARIES}" de.mpi_cbg.revant.apps.WobbleLongPeriodCreateAlphabet1 ${ALPHABET_FILE} ${TMPFILE_PATH}-wobble-1-${WOBBLE_FILE_ID}.txt ${TMPFILE_PATH}-wobble-4-${WOBBLE_FILE_ID}.txt ${TMPFILE_PATH}-wobble-3-lengths-${WOBBLE_FILE_ID}.txt ${TMPFILE_PATH}-wobble-5-flags-${WOBBLE_FILE_ID}.txt
	}
    PIDS=()
	for THREAD in $(seq 0 ${TO}); do
		wobbleThreadCreate ${THREAD} &
        PIDS+=($!)
	done
	waitAndCheck PIDS
	WOBBLE_ALPHABET="${INPUT_DIR}/alphabet-wobble.txt"
	WOBBLE_OLD2NEW="${INPUT_DIR}/alphabet-wobble-old2new.txt"
	java ${JAVA_RUNTIME_FLAGS} -classpath "${REVANT_BINARIES}" de.mpi_cbg.revant.apps.WobbleCreateAlphabet2 ${ALPHABET_FILE} ${WOBBLE_LENGTH} ${MIN_ALIGNMENT_LENGTH} ${REPEAT_LENGTHS_FILE} ${N_REPEATS} ${TMPFILE_PATH}-wobble-5-flags ${TO} 0 ${WOBBLE_ALPHABET} ${WOBBLE_OLD2NEW}
	function wobbleThread() {
		local WOBBLE_FILE_ID=$1
		java ${JAVA_RUNTIME_FLAGS} -classpath "${REVANT_BINARIES}" de.mpi_cbg.revant.apps.WobbleLongPeriod ${TMPFILE_PATH}-wobble-1-${WOBBLE_FILE_ID}.txt ${WOBBLE_LENGTH} ${ALPHABET_FILE} ${WOBBLE_ALPHABET} ${WOBBLE_OLD2NEW} ${REPEAT_LENGTHS_FILE} ${N_REPEATS} ${TMPFILE_PATH}-wobble-4-${WOBBLE_FILE_ID}.txt ${TMPFILE_PATH}-wobble-6-${WOBBLE_FILE_ID}.txt
	}
    PIDS=()
	for THREAD in $(seq 0 ${TO}); do
		wobbleThread ${THREAD} &
        PIDS+=($!)
	done
	waitAndCheck PIDS
	mv ${READS_TRANSLATED_FILE} ${READS_TRANSLATED_FILE}-prewobble
	for THREAD in $(seq 0 ${TO}); do
		cat ${TMPFILE_PATH}-wobble-6-${THREAD}.txt >> ${READS_TRANSLATED_FILE}
	done
	mv ${ALPHABET_FILE} ${ALPHABET_FILE}-prewobble
	mv ${WOBBLE_ALPHABET} ${ALPHABET_FILE}
	echo "Wobbling completed"
fi












PERIODIC_ENDPOINTS_FIXED="0"
if [ ${MAX_SPACER_LENGTH} -ne 0 ]; then
	echo "Trying to fix periodic endpoints if needed..."
	READ_READ_ALIGNMENTS_FILE="${INPUT_DIR}/LAshow-reads-reads.txt"
	echo "Collecting spacers and assigning breakpoints to them..."
	N_FULLY_UNIQUE=$(wc -l < ${FULLY_UNIQUE_FILE})
	N_FULLY_CONTAINED=$(wc -l < ${FULLY_CONTAINED_FILE})
	SPACERS_ARE_CORRECT=$(java ${JAVA_RUNTIME_FLAGS} -classpath "${REVANT_BINARIES}" de.mpi_cbg.revant.apps.FixPeriodicEndpoints1 ${MAX_SPACER_LENGTH} ${MIN_ALIGNMENT_LENGTH} ${N_READS} ${READ_IDS_FILE} ${READ_LENGTHS_FILE} ${ALPHABET_FILE} ${READS_TRANSLATED_FILE} ${READS_TRANSLATED_BOUNDARIES} ${FULLY_UNIQUE_FILE} ${N_FULLY_UNIQUE} ${FULLY_CONTAINED_FILE} ${N_FULLY_CONTAINED} ${READ_READ_ALIGNMENTS_FILE} ${AVG_READ_LENGTH} ${GENOME_LENGTH} ${TMPFILE_PATH}-spacers.txt)
	if [ ${SPACERS_ARE_CORRECT} -eq 1 ]; then
		PERIODIC_ENDPOINTS_FIXED="0"
	else
		echo "Splitting the alignments file..."
		LAST_READA_FILE="${INPUT_DIR}/LAshow-reads-reads-lastReadA.txt"
		if [ ${BROKEN_READS} -eq 1 ]; then
			# Reusing the chunks of the read-read alignments file that are already there
			# (we assume that they all have the header).
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
		# Splitting read translations based on read-read alignments. This is necessary,
		# since the previous translated file chunks were split by balancing read-repeat
		# alignments, rather than read-read alignments.
		java ${JAVA_RUNTIME_FLAGS} -classpath "${REVANT_BINARIES}" de.mpi_cbg.revant.apps.SplitTranslations ${READ_IDS_FILE} ${READS_TRANSLATED_FILE} ${READS_TRANSLATED_BOUNDARIES} ${LAST_READA_FILE} ${TMPFILE_PATH}-spacers-1-1- ${TMPFILE_PATH}-spacers-1-2-
		java ${JAVA_RUNTIME_FLAGS} -classpath "${REVANT_BINARIES}" de.mpi_cbg.revant.apps.SplitSpacers ${LAST_READA_FILE} ${N_THREADS} ${TMPFILE_PATH}-spacers.txt ${READ_IDS_FILE} ${READ_LENGTHS_FILE} ${TMPFILE_PATH}-spacers-2-
		echo "Collecting instances of spacer-induced characters..."
		function collectionThread_spacers() {
			local SPACERS_FILE_ID=$1
			local LOCAL_N_SPACERS=$(wc -l < ${TMPFILE_PATH}-spacers-2-${SPACERS_FILE_ID}.txt)
			java ${JAVA_RUNTIME_FLAGS} -classpath "${REVANT_BINARIES}" de.mpi_cbg.revant.apps.FixPeriodicEndpoints2 ${MAX_SPACER_LENGTH} ${TMPFILE_PATH}-spacers-2-${SPACERS_FILE_ID}.txt ${LOCAL_N_SPACERS} ${ALPHABET_FILE} ${TMPFILE_PATH}-spacers-2-ids-${THREAD}.txt ${TMPFILE_PATH}-spacers-2-lengths-${THREAD}.txt ${TMPFILE_PATH}-spacers-1-1-${THREAD}.txt ${TMPFILE_PATH}-spacers-1-2-${THREAD}.txt ${TMPFILE_PATH}-spacers-3-${SPACERS_FILE_ID}.txt ${TMPFILE_PATH}-spacers-3-unique-${SPACERS_FILE_ID}.txt
			sort --parallel=1 -t , -u ${SORT_OPTIONS} ${TMPFILE_PATH}-spacers-3-${SPACERS_FILE_ID}.txt > ${TMPFILE_PATH}-spacers-4-${SPACERS_FILE_ID}.txt
		}
		if [ -e ${TMPFILE_PATH}-spacers-2-${N_THREADS}.txt ]; then
			TO=${N_THREADS}
		else
			TO=$(( ${N_THREADS} - 1 ))
		fi
        PIDS=()
		for THREAD in $(seq 0 ${TO}); do
			collectionThread_spacers ${THREAD} &
            PIDS+=($!)
		done
		waitAndCheck PIDS
		sort --parallel=${N_THREADS} -m -t , -u ${SORT_OPTIONS} ${TMPFILE_PATH}-spacers-4-*.txt > ${TMPFILE_PATH}-spacers-5.txt
		N_INSTANCES=$(wc -l < ${TMPFILE_PATH}-spacers-5.txt)
		java ${JAVA_RUNTIME_FLAGS} -classpath "${REVANT_BINARIES}" de.mpi_cbg.revant.apps.SplitCharacterInstances ${N_INSTANCES} ${N_THREADS} ${TMPFILE_PATH}-spacers-5.txt ${TMPFILE_PATH}-spacers-6-
		echo "Compacting character instances..."
        PIDS=()
		for THREAD in $(seq 0 ${TO}); do
			compactionThread ${THREAD} "${TMPFILE_PATH}-spacers-6-" "${TMPFILE_PATH}-spacers-7-" "${TMPFILE_PATH}-spacers-8-" &
            PIDS+=($!)
		done
		waitAndCheck PIDS
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
			local FULLYCONTAINED_FILE_NEW=$4
			local LOCAL_SPACERS_FILE="${TMPFILE_PATH}-spacers-2-${THREAD_ID}.txt"
			local LOCAL_N_SPACERS=$(wc -l < ${LOCAL_SPACERS_FILE})
			java ${JAVA_RUNTIME_FLAGS} -classpath "${REVANT_BINARIES}" de.mpi_cbg.revant.apps.FixPeriodicEndpoints3 ${MAX_SPACER_LENGTH} ${LOCAL_SPACERS_FILE} ${LOCAL_N_SPACERS} ${ALPHABET_FILE} ${ALPHABET_FILE_SPACERS} ${TMPFILE_PATH}-spacers-2-ids-${THREAD_ID}.txt ${TMPFILE_PATH}-spacers-2-lengths-${THREAD_ID}.txt ${TMPFILE_PATH}-spacers-1-1-${THREAD_ID}.txt ${TMPFILE_PATH}-spacers-1-2-${THREAD_ID}.txt ${READ2CHARACTERS_FILE_NEW}${THREAD_ID}.txt ${READ2BOUNDARIES_FILE_NEW}${THREAD_ID}.txt ${FULLYCONTAINED_FILE_NEW}${THREAD_ID}.txt
		}
        PIDS=()
		for THREAD in $(seq 0 ${TO}); do
			translationThread_spacers ${THREAD} "${TMPFILE_PATH}-spacers-9-" "${TMPFILE_PATH}-spacers-10-" "${TMPFILE_PATH}-spacers-11-" &
            PIDS+=($!)
		done
		waitAndCheck PIDS
		mv ${READS_TRANSLATED_FILE} ${READS_TRANSLATED_FILE}-prespacers
		mv ${READS_TRANSLATED_BOUNDARIES} ${READS_TRANSLATED_BOUNDARIES}-prespacers
		mv ${ALPHABET_FILE} ${ALPHABET_FILE}-prespacers
		cp ${FULLY_CONTAINED_FILE} ${TMPFILE_PATH}-spacers-12.txt
		mv ${FULLY_CONTAINED_FILE} ${FULLY_CONTAINED_FILE}-prespacers
		for THREAD in $(seq 0 ${TO}); do
			cat ${TMPFILE_PATH}-spacers-9-${THREAD}.txt >> ${READS_TRANSLATED_FILE}
			cat ${TMPFILE_PATH}-spacers-10-${THREAD}.txt >> ${READS_TRANSLATED_BOUNDARIES}
			cat ${TMPFILE_PATH}-spacers-11-${THREAD}.txt >> ${TMPFILE_PATH}-spacers-12.txt
		done
		sort --parallel=${N_THREADS} -n ${TMPFILE_PATH}-spacers-12.txt > ${FULLY_CONTAINED_FILE}
		mv ${ALPHABET_FILE_SPACERS} ${ALPHABET_FILE}
		PERIODIC_ENDPOINTS_FIXED="1"
		echo "Periodic endpoints fixed"
		if [ ${WOBBLE_LENGTH} -ne 0 ]; then
			echo "Wobbling..."
			if [ ${MAX_SPACER_LENGTH} -ne 0 ]; then
				WOBBLE_PREFIX="${TMPFILE_PATH}-spacers-9"
			else
				# Splitting read translations based on read-read alignments. This is
				# necessary, since the previous translated file chunks were split by
				# balancing read-repeat alignments, rather than read-read alignments.
				java ${JAVA_RUNTIME_FLAGS} -classpath "${REVANT_BINARIES}" de.mpi_cbg.revant.apps.SplitTranslations ${READ_IDS_FILE} ${READS_TRANSLATED_FILE} ${READS_TRANSLATED_BOUNDARIES} ${LAST_READA_FILE} ${TMPFILE_PATH}-wobble-1- ${TMPFILE_PATH}-wobble-2-
				WOBBLE_PREFIX="${TMPFILE_PATH}-wobble-1"
			fi	
			function wobbleThreadCreate() {
				local WOBBLE_FILE_ID=$1
				java ${JAVA_RUNTIME_FLAGS} -classpath "${REVANT_BINARIES}" de.mpi_cbg.revant.apps.WobbleCreateAlphabet1 ${ALPHABET_FILE} ${WOBBLE_PREFIX}-${WOBBLE_FILE_ID}.txt ${WOBBLE_PREFIX}-flags-${WOBBLE_FILE_ID}.txt
			}
			if [ -e ${WOBBLE_PREFIX}-${N_THREADS}.txt ]; then
				TO=${N_THREADS}
			else
				TO=$(( ${N_THREADS} - 1 ))
			fi
            PIDS=()
			for THREAD in $(seq 0 ${TO}); do
				wobbleThreadCreate ${THREAD} &
                PIDS+=($!)
			done
			waitAndCheck PIDS
			WOBBLE_ALPHABET="${INPUT_DIR}/alphabet-wobble.txt"
			WOBBLE_OLD2NEW="${INPUT_DIR}/alphabet-wobble-old2new.txt"
			java ${JAVA_RUNTIME_FLAGS} -classpath "${REVANT_BINARIES}" de.mpi_cbg.revant.apps.WobbleCreateAlphabet2 ${ALPHABET_FILE} ${WOBBLE_LENGTH} ${MIN_ALIGNMENT_LENGTH} ${REPEAT_LENGTHS_FILE} ${N_REPEATS} ${WOBBLE_PREFIX}-flags ${TO} 1 ${WOBBLE_ALPHABET} ${WOBBLE_OLD2NEW}
			function wobbleThread() {
				local WOBBLE_FILE_ID=$1
				java ${JAVA_RUNTIME_FLAGS} -classpath "${REVANT_BINARIES}" de.mpi_cbg.revant.apps.Wobble ${WOBBLE_PREFIX}-${WOBBLE_FILE_ID}.txt ${WOBBLE_LENGTH} ${ALPHABET_FILE} ${WOBBLE_ALPHABET} ${WOBBLE_OLD2NEW} ${REPEAT_LENGTHS_FILE} ${N_REPEATS} ${TMPFILE_PATH}-wobble-3-${WOBBLE_FILE_ID}.txt
			}
            PIDS=()
			for THREAD in $(seq 0 ${TO}); do
				wobbleThread ${THREAD} &
                PIDS+=($!)
			done
			waitAndCheck PIDS
			mv ${READS_TRANSLATED_FILE} ${READS_TRANSLATED_FILE}-prewobble
			for THREAD in $(seq 0 ${TO}); do
				cat ${TMPFILE_PATH}-wobble-3-${THREAD}.txt >> ${READS_TRANSLATED_FILE}
			done
			mv ${ALPHABET_FILE} ${ALPHABET_FILE}-prewobble
			mv ${WOBBLE_ALPHABET} ${ALPHABET_FILE}
			echo "Wobbling completed"
		fi
	fi
fi


echo "Discarding rare characters..."
COUNTS_FILE="${INPUT_DIR}/alphabet-counts.txt"
HISTOGRAM_FILE="${INPUT_DIR}/alphabet-histogram.txt"
java ${JAVA_RUNTIME_FLAGS} -classpath "${REVANT_BINARIES}" de.mpi_cbg.revant.apps.GetCharacterCounts ${READS_TRANSLATED_FILE} ${READS_TRANSLATED_BOUNDARIES} ${READ_LENGTHS_FILE} ${ALPHABET_FILE} ${COUNTS_FILE} ${HISTOGRAM_FILE}
function cleaningThread1() {
	local TRANSLATED_CHARACTERS=$1
	local TRANSLATED_BOUNDARIES=$2
	local PREFIX_1=$3
	local PREFIX_2=$4
	local ID=$5
	java ${JAVA_RUNTIME_FLAGS} -classpath "${REVANT_BINARIES}:${REVANT_LIBRARIES}" de.mpi_cbg.revant.apps.CleanTranslatedReads1 ${ALPHABET_FILE} ${COUNTS_FILE} ${N_READS} ${READ_IDS_FILE} ${READ_LENGTHS_FILE} ${AVG_READ_LENGTH} ${SPANNING_BPS} ${MIN_ALIGNMENT_LENGTH} ${GENOME_LENGTH} ${N_HAPLOTYPES} ${TRANSLATED_CHARACTERS} ${TRANSLATED_BOUNDARIES} ${KEEP_PERIODIC} ${PREFIX_1}${ID}.txt > ${PREFIX_1}unique-${ID}.txt
	sort --parallel=1 -t , -u ${SORT_OPTIONS} ${PREFIX_1}${ID}.txt > ${PREFIX_2}${ID}.txt
}
split -l $(( ${N_READS} / ${N_THREADS} )) ${READS_TRANSLATED_FILE} "${TMPFILE_PATH}-12-"
split -l $(( ${N_READS} / ${N_THREADS} )) ${READS_TRANSLATED_BOUNDARIES} "${TMPFILE_PATH}-13-"
PIDS=()
for FILE in $( find ${INPUT_DIR} -maxdepth 1 -name "${TMPFILE_NAME}-12-*" ); do
	THREAD_ID=${FILE#${INPUT_DIR}/${TMPFILE_NAME}-12-}
	cleaningThread1 "${INPUT_DIR}/${TMPFILE_NAME}-12-${THREAD_ID}" "${INPUT_DIR}/${TMPFILE_NAME}-13-${THREAD_ID}" "${TMPFILE_PATH}-14-" "${TMPFILE_PATH}-15-" ${THREAD_ID} &
    PIDS+=($!)
done
waitAndCheck PIDS
sort --parallel=${N_THREADS} -m -t , -u ${SORT_OPTIONS} ${TMPFILE_PATH}-15-*.txt > ${TMPFILE_PATH}-15.txt
cat ${TMPFILE_PATH}-14-unique-*.txt | sort --parallel=${N_THREADS} -n -r > ${TMPFILE_PATH}-14-unique.txt
ALPHABET_FILE_CLEANED="${INPUT_DIR}/alphabet-cleaned.txt"
rm -f ${ALPHABET_FILE_CLEANED}
OLD2NEW_FILE="${INPUT_DIR}/alphabet-old2new.txt"
rm -f ${OLD2NEW_FILE}
java ${JAVA_RUNTIME_FLAGS} -classpath "${REVANT_BINARIES}:${REVANT_LIBRARIES}" de.mpi_cbg.revant.apps.CleanTranslatedReads2 ${ALPHABET_FILE} ${COUNTS_FILE} ${N_READS} ${AVG_READ_LENGTH} ${SPANNING_BPS} ${MIN_ALIGNMENT_LENGTH} ${GENOME_LENGTH} ${N_HAPLOTYPES} $(wc -l < ${TMPFILE_PATH}-15.txt) ${TMPFILE_PATH}-15.txt ${KEEP_PERIODIC} ${TMPFILE_PATH}-14-unique.txt ${ALPHABET_FILE_CLEANED} ${OLD2NEW_FILE}
function cleaningThread3() {
	local TRANSLATED_CHARACTERS_OLD=$1
	local TRANSLATED_BOUNDARIES_OLD=$2
	local TRANSLATED_CHARACTERS_NEW=$3
	local TRANSLATED_BOUNDARIES_NEW=$4
	local HISTOGRAM_FILE=$5
	local FULLY_UNIQUE_FILE_NEW=$6
	local LAST_TRANSLATED_READ=$7
	java ${JAVA_RUNTIME_FLAGS} -classpath "${REVANT_BINARIES}:${REVANT_LIBRARIES}" de.mpi_cbg.revant.apps.CleanTranslatedReads3 ${N_READS} ${READ_IDS_FILE} ${READ_LENGTHS_FILE} ${ALPHABET_FILE} ${COUNTS_FILE} ${TRANSLATED_CHARACTERS_OLD} ${TRANSLATED_BOUNDARIES_OLD} ${AVG_READ_LENGTH} ${SPANNING_BPS} ${MIN_ALIGNMENT_LENGTH} ${GENOME_LENGTH} ${N_HAPLOTYPES} ${KEEP_PERIODIC} ${ALPHABET_FILE_CLEANED} ${OLD2NEW_FILE} ${TRANSLATED_CHARACTERS_NEW} ${TRANSLATED_BOUNDARIES_NEW} ${HISTOGRAM_FILE} ${FULLY_UNIQUE_FILE_NEW} ${LAST_TRANSLATED_READ}
}
LAST_TRANSLATED_READ="-1"
PIDS=()
for FILE in $(find -s ${INPUT_DIR} -maxdepth 1 -name "${TMPFILE_NAME}-12-*" ); do
	THREAD_ID=${FILE#${INPUT_DIR}/${TMPFILE_NAME}-12-}
	cleaningThread3 "${INPUT_DIR}/${TMPFILE_NAME}-12-${THREAD_ID}" "${INPUT_DIR}/${TMPFILE_NAME}-13-${THREAD_ID}" "${INPUT_DIR}/${TMPFILE_NAME}-16-${THREAD_ID}" "${INPUT_DIR}/${TMPFILE_NAME}-17-${THREAD_ID}" "${INPUT_DIR}/${TMPFILE_NAME}-18-${THREAD_ID}" "${INPUT_DIR}/${TMPFILE_NAME}-19-${THREAD_ID}" ${LAST_TRANSLATED_READ} &
    PIDS+=($!)
	LAST_TRANSLATED_READ=$(( ${LAST_TRANSLATED_READ} + $(wc -l < "${INPUT_DIR}/${TMPFILE_NAME}-12-${THREAD_ID}") ))
done
waitAndCheck PIDS
READS_TRANSLATED_FILE_NEW="${INPUT_DIR}/reads-translated-new.txt"
READS_TRANSLATED_BOUNDARIES_NEW="${INPUT_DIR}/reads-translated-boundaries-new.txt"
HISTOGRAM_FILE="${INPUT_DIR}/reads-translated-histogram-new.txt"
FULLY_UNIQUE_FILE_NEW="${INPUT_DIR}/reads-fullyUnique-new.txt"
FULLY_CONTAINED_FILE_NEW="${INPUT_DIR}/reads-fullyContained-new.txt"
rm -f ${READS_TRANSLATED_FILE_NEW} ${READS_TRANSLATED_BOUNDARIES_NEW} ${FULLY_UNIQUE_FILE_NEW} ${FULLY_CONTAINED_FILE_NEW}
for FILE in $(find -s ${INPUT_DIR} -maxdepth 1 -name "${TMPFILE_NAME}-16-*" ); do
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
	SPACERS_NUMBER=$( find ${INPUT_DIR} -maxdepth 1 -name "${TMPFILE_NAME}-spacers-1-*.txt" | wc -l )
	if [ ${SPACERS_NUMBER} -gt 0 ]; then
		rm -rf ${INPUT_DIR}/preserve
		mkdir ${INPUT_DIR}/preserve
		mv ${TMPFILE_PATH}-spacers-1-*.txt ${INPUT_DIR}/preserve
	fi
	rm -f ${TMPFILE_PATH}*
	if [ ${SPACERS_NUMBER} -gt 0 ]; then
		mv ${INPUT_DIR}/preserve/*.txt ${INPUT_DIR}
		rm -rf ${INPUT_DIR}/preserve
	fi
fi

# Return value
echo ${PERIODIC_ENDPOINTS_FIXED} > ${TMPFILE_PATH}-return.txt
