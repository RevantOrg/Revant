#!/bin/bash
#
# Given the directory of a project, the script builds the graph of all intervals and
# partitions it into connected components. Every connected component can then be
# processed in parallel.
#
# Remark: the script assumes that environment variable $REVANT_BINARIES$ is already set
# to the directory that contains REVANT's binaries, i.e. to the directory that contains
# the $de/mpi_cbg/revant/...$ subpath.
#
# --------------------------------- CONFIGURATION ----------------------------------------
# This is the only section of the script that needs to be customized.
#
PROJECT_DIR=$1
MIN_ALIGNMENT_LENGTH="1500"
MIN_OCCURRENCES_IN_GENOME="3"
FIX_DANGLING_EDGES="1"
FIX_UNASSIGNED_INTERVALS="1"
ADD_SAME_READ_EDGES="0"
GENOME_COVERAGE="1"
DISCARD_CONTAINED_IN_SHORT_PERIOD="0"
DISCARD_CONTAINED_IN_LONG_PERIOD="0"
DISCARD_PERIOD_IN_PERIOD="0"
DISCARD_DENSE_IN_DENSE="0"
# REVANT
JAVA_RUNTIME_FLAGS="-Xms10G -Xmx16G"
# ----------------------------------------------------------------------------------------




N_READS=$(wc -l < ${PROJECT_DIR}/reads-lengths.txt)
TAGS_DIR="${OUTPUT_DIR}/finalOutput"
TAGS_FILE="${TAGS_DIR}/tags-root.txt"
OUTPUT_DIR="${PROJECT_DIR}/step1"

rm -rf ${OUTPUT_DIR}
mkdir ${OUTPUT_DIR}
mkdir ${TAGS_DIR}
java ${JAVA_RUNTIME_FLAGS} -classpath "${REVANT_BINARIES}" de.mpi_cbg.revant.intervalgraph.IntervalGraphStep1 ${PROJECT_DIR} ${N_READS} ${GENOME_COVERAGE} ${MIN_OCCURRENCES_IN_GENOME} ${MIN_ALIGNMENT_LENGTH} ${DISCARD_CONTAINED_IN_SHORT_PERIOD} ${DISCARD_CONTAINED_IN_LONG_PERIOD} ${DISCARD_PERIOD_IN_PERIOD} ${DISCARD_DENSE_IN_DENSE} ${FIX_DANGLING_EDGES} ${FIX_UNASSIGNED_INTERVALS} ${ADD_SAME_READ_EDGES}
TMP_FILE="${TAGS_DIR}/tags-root.tmp"
sort -t , -k 1,1n -k 2,2n -k 3,3n -o ${TMP_FILE} ${TAGS_FILE}
rm -f ${TAGS_FILE}
mv ${TMP_FILE} ${TAGS_FILE}
