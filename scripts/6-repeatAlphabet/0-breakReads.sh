#!/bin/bash
# 
# This step is needed only for CLR. Breaks reads at long low-quality intervals before
# feeding them to the following steps of the pipeline.
#
# Remark: the script assumes that environment variable $REVANT_BINARIES$ is already set
# to the directory that contains REVANT's binaries, i.e. to the directory that contains
# the $de/mpi_cbg/revant/...$ subpath.
#
# --------------------------------- CONFIGURATION ----------------------------------------
# This is the only section of the script that needs to be customized.
#
INPUT_DIR=$1
QUALITY_TRACK_FILE="${INPUT_DIR}/reads-phred.dbdump"
LOW_QUALITY_LENGTH="100"  # >=100
N_THREADS="4"
# REVANT
JAVA_RUNTIME_FLAGS="-Xms2G -Xmx10G"
# ----------------------------------------------------------------------------------------


QUALITY_THRESHOLDS_FILE="${INPUT_DIR}/qualityThresholds.txt"
READ_LENGTHS_FILE="${INPUT_DIR}/reads-lengths.txt"
N_READS=$(wc -l < ${READ_LENGTHS_FILE})
TMPFILE_NAME="breakReads-tmp"
TMPFILE_PATH="${INPUT_DIR}/${TMPFILE_NAME}"
rm -f ${TMPFILE_PATH}*

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
alignmentsThread1 0 1 &
for THREAD in $(seq 1 ${TO}); do
	alignmentsThread1 ${THREAD} 0 &
done
wait
echo "Read-read alignments translated successfully"
ALIGNMENTS_FILE_BROKEN="${INPUT_DIR}/LAshow-reads-reads-broken.txt"
rm -f ${ALIGNMENTS_FILE_BROKEN}
for THREAD in $(seq 0 ${TO}); do
	cat ${TMPFILE_PATH}-2-${THREAD}.txt >> ${ALIGNMENTS_FILE_BROKEN}
done

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
alignmentsThread2 0 1 &
for THREAD in $(seq 1 ${TO}); do
	alignmentsThread2 ${THREAD} 0 &
done
wait
echo "Read-repeat alignments translated successfully"
ALIGNMENTS_FILE_BROKEN="${INPUT_DIR}/LAshow-reads-repeats-broken.txt"
rm -f ${ALIGNMENTS_FILE_BROKEN}
for THREAD in $(seq 0 ${TO}); do
	cat ${TMPFILE_PATH}-4-${THREAD}.txt >> ${ALIGNMENTS_FILE_BROKEN}
done

# Renaming files for the following steps of the pipeline
mv "${INPUT_DIR}/LAshow-reads-reads.txt" "${INPUT_DIR}/LAshow-reads-reads-unbroken.txt"
mv "${INPUT_DIR}/LAshow-reads-reads-broken.txt" "${INPUT_DIR}/LAshow-reads-reads.txt"
mv "${INPUT_DIR}/LAshow-reads-repeats.txt" "${INPUT_DIR}/LAshow-reads-repeats-unbroken.txt"
mv "${INPUT_DIR}/LAshow-reads-repeats-broken.txt" "${INPUT_DIR}/LAshow-reads-repeats.txt"
mv "${INPUT_DIR}/reads-lengths.txt" "${INPUT_DIR}/reads-lengths-unbroken.txt"
mv "${INPUT_DIR}/reads-lengths-broken.txt" "${INPUT_DIR}/reads-lengths.txt"
mv "${INPUT_DIR}/reads-ids.txt" "${INPUT_DIR}/reads-ids-unbroken.txt"
mv "${INPUT_DIR}/LAshow-reads-reads-lastReadA.txt" "${INPUT_DIR}/LAshow-reads-reads-lastReadA-unbroken.txt"
mv "${INPUT_DIR}/LAshow-reads-repeats-lastReadA.txt" "${INPUT_DIR}/LAshow-reads-repeats-lastReadA-unbroken.txt"
seq 0 $(( ${N_READS_BROKEN} - 1 )) > "${INPUT_DIR}/reads-ids.txt"