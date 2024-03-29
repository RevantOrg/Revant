#!/bin/bash
# 
# This step is needed only for CLR: if the data is not CLR, start from Step 1 instead.
# Breaks reads at long low-quality intervals before feeding them to the following steps
# of the pipeline.
#
# Remark: the script assumes that environment variable $REVANT_BINARIES$ is already set
# to the directory that contains REVANT's binaries, i.e. to the directory that contains
# the $de/mpi_cbg/revant/...$ subpath.
#
# --------------------------------- CONFIGURATION ----------------------------------------
# This is the only section of the script that needs to be customized.
#
INPUT_DIR=$1
LOW_QUALITY_LENGTH=$2  # >=100
N_THREADS=$3
DELETE_TMP_FILES=$4
QUALITY_TRACK_FILE="${INPUT_DIR}/reads-phred.dbdump"
# ----------------------------------------------------------------------------------------

set -o pipefail; set -e; set -u
QUALITY_THRESHOLDS_FILE="${INPUT_DIR}/qualityThresholds.txt"
READ_LENGTHS_FILE="${INPUT_DIR}/reads-lengths.txt"
N_READS=$(wc -l < ${READ_LENGTHS_FILE})
TMPFILE_NAME="breakReads-tmp"
TMPFILE_PATH="${INPUT_DIR}/${TMPFILE_NAME}"
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


echo "Breaking reads..."
OLD2NEW_FILE="${INPUT_DIR}/unbroken2broken.txt"
NEW2OLD_FILE="${INPUT_DIR}/broken2unbroken.txt"
READ_LENGTHS_FILE_NEW="${INPUT_DIR}/reads-lengths-broken.txt"
rm -f ${OLD2NEW_FILE} ${NEW2OLD_FILE} ${READ_LENGTHS_FILE_NEW}
java ${JAVA_RUNTIME_FLAGS} -classpath "${REVANT_BINARIES}" de.mpi_cbg.revant.apps.BreakReads1 ${N_READS} ${READ_LENGTHS_FILE} ${QUALITY_TRACK_FILE} ${QUALITY_THRESHOLDS_FILE} ${LOW_QUALITY_LENGTH} ${OLD2NEW_FILE} ${NEW2OLD_FILE} ${READ_LENGTHS_FILE_NEW} > ${TMPFILE_PATH}-nReads
N_READS_BROKEN=$(cat ${TMPFILE_PATH}-nReads)

echo "Translating read-read alignments..."
function alignmentsThread1() {
	local ALIGNMENTS_FILE_ID=$1
	local WRITE_HEADER=$2
	java ${JAVA_RUNTIME_FLAGS} -classpath "${REVANT_BINARIES}" de.mpi_cbg.revant.apps.BreakReads2 ${N_READS} ${READ_LENGTHS_FILE} ${OLD2NEW_FILE} 1 ${WRITE_HEADER} ${TMPFILE_PATH}-1-${ALIGNMENTS_FILE_ID}.txt ${TMPFILE_PATH}-2-${ALIGNMENTS_FILE_ID}.txt
}
ALIGNMENTS_FILE="${INPUT_DIR}/LAshow-reads-reads.txt"
N_ALIGNMENTS=$(( $(wc -l < ${ALIGNMENTS_FILE}) - 2 ))
LAST_READA_FILE="${INPUT_DIR}/LAshow-reads-reads-lastReadA.txt"
java ${JAVA_RUNTIME_FLAGS} -classpath "${REVANT_BINARIES}" de.mpi_cbg.revant.factorize.SplitAlignments ${N_ALIGNMENTS} ${N_THREADS} ${ALIGNMENTS_FILE} ${TMPFILE_PATH}-1- ${LAST_READA_FILE}
echo "Read-read alignments split in ${N_THREADS} parts"
if [ -e ${TMPFILE_PATH}-1-${N_THREADS}.txt ]; then
	TO=${N_THREADS}
else
	TO=$(( ${N_THREADS} - 1 ))
fi
PIDS=()
alignmentsThread1 0 1 &
PIDS+=($!)
for THREAD in $(seq 1 ${TO}); do
	alignmentsThread1 ${THREAD} 1 &  # We always write the header for future scripts.
    PIDS+=($!)
done
waitAndCheck PIDS
echo "Read-read alignments translated successfully"
# We do not need to concatenate the chunks in a single file, since they will be used
# directly as chunks by the following scripts.
# Remark: if we really want to concatenate, we should remove headers from all chunks,
# except the first one, in the call above.
#ALIGNMENTS_FILE_BROKEN="${INPUT_DIR}/LAshow-reads-reads-broken.txt"
#rm -f ${ALIGNMENTS_FILE_BROKEN}
#for THREAD in $(seq 0 ${TO}); do
#	cat ${TMPFILE_PATH}-2-${THREAD}.txt >> ${ALIGNMENTS_FILE_BROKEN}
#done

echo "Translating read-repeat alignments..."
function alignmentsThread2() {
	local ALIGNMENTS_FILE_ID=$1
	local WRITE_HEADER=$2
	java ${JAVA_RUNTIME_FLAGS} -classpath "${REVANT_BINARIES}" de.mpi_cbg.revant.apps.BreakReads2 ${N_READS} ${READ_LENGTHS_FILE} ${OLD2NEW_FILE} 0 ${WRITE_HEADER} ${TMPFILE_PATH}-3-${ALIGNMENTS_FILE_ID}.txt ${TMPFILE_PATH}-4-${ALIGNMENTS_FILE_ID}.txt
}
ALIGNMENTS_FILE="${INPUT_DIR}/LAshow-reads-repeats.txt"
N_ALIGNMENTS=$(( $(wc -l < ${ALIGNMENTS_FILE}) - 2 ))
LAST_READA_FILE="${INPUT_DIR}/LAshow-reads-repeats-lastReadA.txt"
java ${JAVA_RUNTIME_FLAGS} -classpath "${REVANT_BINARIES}" de.mpi_cbg.revant.factorize.SplitAlignments ${N_ALIGNMENTS} ${N_THREADS} ${ALIGNMENTS_FILE} ${TMPFILE_PATH}-3- ${LAST_READA_FILE}
echo "Read-repeat alignments split in ${N_THREADS} parts"
if [ -e ${TMPFILE_PATH}-3-${N_THREADS}.txt ]; then
	TO=${N_THREADS}
else
	TO=$(( ${N_THREADS} - 1 ))
fi
PIDS=()
alignmentsThread2 0 1 &
PIDS+=($!)
for THREAD in $(seq 1 ${TO}); do
	alignmentsThread2 ${THREAD} 1 &  # We always write the header for future scripts.
    PIDS+=($!)
done
waitAndCheck PIDS
echo "Read-repeat alignments translated successfully"
# We do not need to concatenate the chunks in a single file, since they will be used
# directly as chunks by the following scripts.
# Remark: if we really want to concatenate, we should remove headers from all chunks,
# except the first one, in the call above.
#ALIGNMENTS_FILE_BROKEN="${INPUT_DIR}/LAshow-reads-repeats-broken.txt"
#rm -f ${ALIGNMENTS_FILE_BROKEN}
#for THREAD in $(seq 0 ${TO}); do
#	cat ${TMPFILE_PATH}-4-${THREAD}.txt >> ${ALIGNMENTS_FILE_BROKEN}
#done

# Renaming files for the following steps of the pipeline
mv "${INPUT_DIR}/LAshow-reads-reads.txt" "${INPUT_DIR}/LAshow-reads-reads-unbroken.txt"
#mv "${INPUT_DIR}/LAshow-reads-reads-broken.txt" "${INPUT_DIR}/LAshow-reads-reads.txt"
mv "${INPUT_DIR}/LAshow-reads-repeats.txt" "${INPUT_DIR}/LAshow-reads-repeats-unbroken.txt"
#mv "${INPUT_DIR}/LAshow-reads-repeats-broken.txt" "${INPUT_DIR}/LAshow-reads-repeats.txt"
mv "${INPUT_DIR}/reads-lengths.txt" "${INPUT_DIR}/reads-lengths-unbroken.txt"
mv "${INPUT_DIR}/reads-lengths-broken.txt" "${INPUT_DIR}/reads-lengths.txt"
mv "${INPUT_DIR}/reads-ids.txt" "${INPUT_DIR}/reads-ids-unbroken.txt"
mv "${INPUT_DIR}/LAshow-reads-reads-lastReadA.txt" "${INPUT_DIR}/LAshow-reads-reads-lastReadA-unbroken.txt"
mv "${INPUT_DIR}/LAshow-reads-repeats-lastReadA.txt" "${INPUT_DIR}/LAshow-reads-repeats-lastReadA-unbroken.txt"
seq 0 $(( ${N_READS_BROKEN} - 1 )) > "${INPUT_DIR}/reads-ids.txt"

# Removing all temp files that are not used downstream
if [ ${DELETE_TMP_FILES} -eq 1 ]; then
	rm -f ${TMPFILE_PATH}-3-*
fi