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
PARTS_PREFIX="${INPUT_DIR}/LAshow-"


# ------------------------------------ A THREAD ------------------------------------------
function run() {
	ALIGNMENTS_FILE_ID=$1
	java ${JAVA_RUNTIME_FLAGS} -classpath "${REVANT_BINARIES}" de.mpi_cbg.revant.apps.CollectCharacterInstances ${N_READS} ${READ_IDS_FILE} ${READ_LENGTHS_FILE} ${N_REPEATS} ${REPEAT_LENGTHS_FILE} ${REPEAT_ISPERIODIC_FILE} ${PARTS_PREFIX}${ALIGNMENTS_FILE_ID}.txt ${MAX_ALIGNMENT_ERROR} ${PARTS_PREFIX}${ALIGNMENTS_FILE_ID}-instances.txt
	sort -n -t , -o ${PARTS_PREFIX}${ALIGNMENTS_FILE_ID}-instances-sorted.txt ${PARTS_PREFIX}${ALIGNMENTS_FILE_ID}-instances.txt
}



# ---------------------------------- MAIN PROGRAM ----------------------------------------

echo "Splitting the alignments file..."
ALIGNMENTS_FILE="${INPUT_DIR}/LAshow-reads-repeats.txt"
N_ALIGNMENTS=$(( $(wc -l < ${ALIGNMENTS_FILE}) - 2 ))
LAST_READA_FILE="${INPUT_DIR}/LAshow-lastReadA.txt"
java ${JAVA_RUNTIME_FLAGS} -classpath "${REVANT_BINARIES}" de.mpi_cbg.revant.factorize.SplitAlignments ${N_ALIGNMENTS} ${SPLIT_IN_PARTS} ${ALIGNMENTS_FILE} ${PARTS_PREFIX} ${LAST_READA_FILE}
echo "Alignments filtered and split in ${SPLIT_IN_PARTS} parts"

echo "Collecting character instances..."
for THREAD in $(seq 0 $(( ${SPLIT_IN_PARTS} - 1 ))); do
	run ${THREAD} &
done
wait
sort -m -n -t , -o ${PARTS_PREFIX}instances-sorted.txt ${PARTS_PREFIX}*-instances-sorted.txt
