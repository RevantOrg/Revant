#!/bin/bash
#
# Given the directory of a project, the script draws a graphical representation of every
# read in a given chunk. One image file is created in the current directory per read.
#
# Remark: the script assumes that environment variable $REVANT_BINARIES$ is already set
# to the directory that contains REVANT's binaries, i.e. to the directory that contains
# the $de/mpi_cbg/revant/...$ subpath.
#
# --------------------------------- CONFIGURATION ----------------------------------------
# This is the only section of the script that needs to be customized.
#
PROJECT_DIR=$1
CHUNK=$2
MIN_ALIGNMENT_LENGTH="1500"
MAX_ERROR_RATE="0.30"
MAX_READ_LENGTH="261338"
GENOME_COVERAGE="1"
MIN_OCCURRENCES_IN_GENOME="3"
INCLUDES_SELFALIGNMENTS="0"
MIN_FACTOR_LENGTH="250"
MIN_LONG_INSERTION_LENGTH="1000"
STEP="0"  # Draw the information produced: 0=after factorization; {3,4}=after Step{3,4} of the Interval Graph pipeline.
# Ground-truth repeats
TRANSPARENT="0"  # 0/1: 1=does not fill the ground-truth repeat marks.
REFERENCE_ALIGNMENTS="null"
REFERENCE_NAMES="null"
REFERENCE_LENGTHS="null"
# Another family of ground-truth repeats (drawn in a different color).
REFERENCE_ALIGNMENTS_PRIME="null"
REFERENCE_NAMES_PRIME="null"
REFERENCE_LENGTHS_PRIME="null"
# Souce files for plotting additional things
AUXILIARY_COVERAGE_FILE="null"  # Another coverage histogram
DAZZLER_TRACK_LOWCOMPLEXITY_FILE="null"
DAZZLER_TRACK_TANDEM_FILE="null"
# REVANT
JAVA_RUNTIME_FLAGS="-Xms2G -Xmx10G"
# ----------------------------------------------------------------------------------------




OUTPUT_ALL_PERIODIC="1"
N_READS=$(wc -l < "${PROJECT_DIR}/reads-ids.txt")
N_ALIGNMENTS=$(wc -l < "${PROJECT_DIR}/LAshow-${CHUNK}.txt")
N_ALIGNMENTS=$(( ${N_ALIGNMENTS} - 2 ))  # Removing LAshow header
READ_LENGTHS_FILE="${PROJECT_DIR}/reads-lengths.txt"
READ_IDS_FILE="${PROJECT_DIR}/reads-ids.txt"
ALIGNMENTS_FILE="${PROJECT_DIR}/LAshow-${CHUNK}.txt"
QUALITY_THRESHOLDS_FILE="${PROJECT_DIR}/qualityThresholds.txt"
QUALITIES_FILE="${PROJECT_DIR}/reads-phred.damar"
if [ ! -e ${QUALITIES_FILE} ]; then
	QUALITIES_FILE="${PROJECT_DIR}/reads-phred.dbdump"
fi
INTERVALS_DENSE_FILE="${PROJECT_DIR}/intervals-dense-${CHUNK}.txt"
INTERVALS_PERIODIC_FILE="${PROJECT_DIR}/intervals-periodic-${CHUNK}.txt"
INTERVALS_ALIGNMENTS_FILE="${PROJECT_DIR}/intervals-alignments-${CHUNK}.txt"
INTERVALS_CONNECTION_FILE="${PROJECT_DIR}/connection-${CHUNK}.txt"
if [ ${STEP} -eq 0 ]; then
	KERNEL_TAGS_FILE="null"
	KERNEL2DESCRIPTOR_DIR="null"
elif [ ${STEP} -eq 3 ]; then
	KERNEL_TAGS_FILE="${PROJECT_DIR}/step1/finalOutput/allTags.txt"
	KERNEL2DESCRIPTOR_DIR="null"
elif [ ${STEP} -eq 4 ]; then
	KERNEL_TAGS_FILE="${PROJECT_DIR}/step1/finalOutput/step4/allTags.txt"
	KERNEL2DESCRIPTOR_DIR="${PROJECT_DIR}/step1/finalOutput/step4"
elif
	
java ${JAVA_RUNTIME_FLAGS} -classpath "${REVANT_BINARIES}" de.mpi_cbg.revant.graphics.PrintReadAlignments $N_READS $MAX_READ_LENGTH $READ_LENGTHS_FILE ${READ_IDS_FILE} $N_ALIGNMENTS $ALIGNMENTS_FILE $QUALITY_THRESHOLDS_FILE $QUALITIES_FILE $GENOME_COVERAGE $MIN_OCCURRENCES_IN_GENOME $INCLUDES_SELFALIGNMENTS $MIN_ALIGNMENT_LENGTH $MIN_FACTOR_LENGTH $MIN_LONG_INSERTION_LENGTH null null null null $INTERVALS_DENSE_FILE $INTERVALS_PERIODIC_FILE $INTERVALS_ALIGNMENTS_FILE $OUTPUT_ALL_PERIODIC $INTERVALS_CONNECTION_FILE ${REFERENCE_ALIGNMENTS} ${REFERENCE_LENGTHS} ${REFERENCE_NAMES} ${REFERENCE_ALIGNMENTS_PRIME} ${REFERENCE_LENGTHS_PRIME} ${REFERENCE_NAMES_PRIME} ${AUXILIARY_COVERAGE_FILE} ${DAZZLER_TRACK_LOWCOMPLEXITY_FILE} ${DAZZLER_TRACK_TANDEM_FILE} ${KERNEL_TAGS_FILE} ${MAX_ERROR_RATE} ${TRANSPARENT} ${KERNEL2DESCRIPTOR_DIR}
