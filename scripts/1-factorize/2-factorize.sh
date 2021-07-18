#!/bin/bash
#
# Given the directory of a project, the script performs read factorization on a chunk of
# input alignments. Multiple instances of this script can be executed in parallel on
# different chunks.
#
# Remark: the script assumes that environment variable $REVANT_BINARIES$ is already set
# to the directory that contains REVANT's binaries, i.e. to the directory that contains
# the $de/mpi_cbg/revant/...$ subpath.
#
# Remark: the script assumes that environment variable $REVANT_BIOLOGY$ is already set
# to the directory that contains REVANT's biology constants.
#
# --------------------------------- CONFIGURATION ----------------------------------------
# This is the only section of the script that needs to be customized.
#
PROJECT_DIR=$1
CHUNK=$1
MIN_ALIGNMENT_LENGTH="1500"
MAX_READ_LENGTH="261338"   # Median=28005
GENOME_COVERAGE="1"  # Use 1 even when smaller than one.
MIN_OCCURRENCES_IN_GENOME="3"  # For calling something a repeat
INCLUDES_SELFALIGNMENTS="0"  # 1=the alignments file contains alignments of a read to itself.
MIN_FACTOR_LENGTH="250"
MIN_LONG_INSERTION_LENGTH="1000"  # Min length of a long random insertion in CCS reads.
# REVANT
JAVA_RUNTIME_FLAGS="-Xms2G -Xmx10G"
# ----------------------------------------------------------------------------------------




OUTPUT_ALL_PERIODIC="1"
READ_LENGTHS_FILE="${PROJECT_DIR}/reads-lengths.txt"
N_READS=$(wc -l < ${READ_LENGTHS_FILE})
READ_IDS_FILE="${PROJECT_DIR}/reads-ids.txt"
ALIGNMENTS_FILE="${PROJECT_DIR}/LAshow-${CHUNK}.txt"
N_ROWS=$(wc -l < ${ALIGNMENTS_FILE})
N_ALIGNMENTS=$(( ${N_ROWS} - 2 ))
QUALITY_THRESHOLDS_FILE="${PROJECT_DIR}/qualityThresholds.txt"
FACTORS_FILE="${PROJECT_DIR}/factors-${CHUNK}.txt"
FACTORS_INTERVALS_DENSE_FILE="${PROJECT_DIR}/factors-intervals-dense-${CHUNK}.txt"
FACTORS_INTERVALS_PERIODIC_FILE="${PROJECT_DIR}/factors-intervals-periodic-${CHUNK}.txt"
FACTORS_INTERVALS_ALIGNMENT_FILE="${PROJECT_DIR}/factors-intervals-alignment-${CHUNK}.txt"
INTERVALS_DENSE_FILE="${PROJECT_DIR}/intervals-dense-${CHUNK}.txt"
INTERVALS_PERIODIC_FILE="${PROJECT_DIR}/intervals-periodic-${CHUNK}.txt"
INTERVALS_ALIGNMENTS_FILE="${PROJECT_DIR}/intervals-alignments-${CHUNK}.txt"
INTERVALS_CONNECTION_FILE="${PROJECT_DIR}/connection-${CHUNK}.txt"
QUALITIES_FILE="${PROJECT_DIR}/reads-phred.damar"
if [ ! -e ${QUALITIES_FILE} ]; then
	QUALITIES_FILE="${PROJECT_DIR}/reads-phred.dbdump"
fi

java ${JAVA_RUNTIME_FLAGS} -classpath "${REVANT_BINARIES}" de.mpi_cbg.revant.factorize.Factorize ${N_READS} ${MAX_READ_LENGTH} ${READ_LENGTHS_FILE} ${READ_IDS_FILE} ${N_ALIGNMENTS} ${ALIGNMENTS_FILE} ${QUALITY_THRESHOLDS_FILE} ${QUALITIES_FILE} ${GENOME_COVERAGE} ${MIN_OCCURRENCES_IN_GENOME} ${INCLUDES_SELFALIGNMENTS} ${MIN_ALIGNMENT_LENGTH} ${MIN_FACTOR_LENGTH} ${MIN_LONG_INSERTION_LENGTH} ${FACTORS_FILE} ${FACTORS_INTERVALS_DENSE_FILE} ${FACTORS_INTERVALS_PERIODIC_FILE} ${FACTORS_INTERVALS_ALIGNMENT_FILE} ${INTERVALS_DENSE_FILE} ${INTERVALS_PERIODIC_FILE} ${INTERVALS_ALIGNMENTS_FILE} ${OUTPUT_ALL_PERIODIC} ${INTERVALS_CONNECTION_FILE} "${REVANT_BIOLOGY}"
java ${JAVA_RUNTIME_FLAGS} -classpath "${REVANT_BINARIES}" de.mpi_cbg.revant.factorize.CheckFactorization ${ALIGNMENTS_FILE} ${INTERVALS_CONNECTION_FILE} ${INTERVALS_ALIGNMENTS_FILE} ${INTERVALS_DENSE_FILE} ${INTERVALS_PERIODIC_FILE}
