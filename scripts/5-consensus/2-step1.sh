#!/bin/bash
#
# Checks for consistency all consensus alignments built by the previous script. See
# $ConsensusStep1.java$ for details.
#
# Remark: the script assumes that environment variable $REVANT_BINARIES$ is already set
# to the directory that contains REVANT's binaries, i.e. to the directory that contains
# the $de/mpi_cbg/revant/...$ subpath.
#
# --------------------------------- CONFIGURATION ----------------------------------------
# This is the only part of the script that needs to be customized.
#
PROJECT_DIR=$1
CONSENSI_REPEAT_TRACK=$2
CONSENSI_REPEAT_TRACK_FORMAT="0"  # 0=DAZZLER, 1=DAMAR.
CONSENSUS_SUFFIX=".filt.daccord.consensus"  # Excluding ".fasta"
FRAGMENT_CONSENSUS_ALIGNMENTS_FORMAT="1"  # Consensi, fragments in distinct DBs.
FRAGMENT_DOMINATOR_ALIGNMENTS_FORMAT="0"  # Consensi, fragments in the same DB.
MAX_READ_LENGTH="200000"  # Max length of a consensus
FRAGMENT_CONSENSUS_ERROR="0.3"
CONSENSUS_CONSENSUS_ERROR="0.15"
MIN_ALIGNMENT_LENGTH_CONSENSUS_CONSENSUS="500"
MIN_ALIGNMENT_LENGTH_REPEAT_INFERENCE="1500"
BLOCK_CONSENSUS_ANALYSIS="0"
# REVANT
JAVA_RUNTIME_FLAGS="-Xms2G -Xmx10G"
# ----------------------------------------------------------------------------------------




STEP1_DIR="${PROJECT_DIR}/step1"
QUALITY_THRESHOLDS_FILE="${PROJECT_DIR}/qualityThresholds.txt"
CONSENSI_STRINGS_DIR="${STEP1_DIR}/finalOutput/step4/step5/fragments-strings-alignments/consensi/consensi-strings"
CONSENSI_LIST="${CONSENSI_STRINGS_DIR}/list-consensi.txt"
CONSENSUS_CONSENSUS_ALIGNMENTS="${CONSENSI_STRINGS_DIR}/LAshow-consensi-consensi.txt"
READ_LENGTHS_FILE="${CONSENSI_STRINGS_DIR}/reads-lengths.txt"
READ_IDS_FILE="${CONSENSI_STRINGS_DIR}/reads-ids.txt"
CONSENSUS_CONSENSUS_ALIGNMENTS_FILTERED="${CONSENSI_STRINGS_DIR}/LAshow-consensi-consensi-filtered.txt"

# Filtering consensus-consensus alignments
rm -rf ${READ_LENGTHS_FILE}
for INPUT_FILE in $(cat ${CONSENSI_LIST}); do
	BASE_NAME=$(basename ${INPUT_FILE} .fasta)
	java ${JAVA_RUNTIME_FLAGS} -classpath "${REVANT_BINARIES}" de.mpi_cbg.revant.util.Fasta2Lengths ${CONSENSI_STRINGS_DIR}/${INPUT_FILE} > ${CONSENSI_STRINGS_DIR}/${BASE_NAME}-lengths.txt
	cat ${CONSENSI_STRINGS_DIR}/${BASE_NAME}-lengths.txt >> ${READ_LENGTHS_FILE}
done
N_ALIGNMENTS=$(wc -l < ${CONSENSUS_CONSENSUS_ALIGNMENTS})  # Just total n. of rows
N_READS=$(wc -l < ${READ_LENGTHS_FILE})
seq 0 $(( ${N_READS} - 1 )) > ${READ_IDS_FILE}
java ${JAVA_RUNTIME_FLAGS} -classpath "${REVANT_BINARIES}" de.mpi_cbg.revant.factorize.FilterAlignments ${CONSENSUS_CONSENSUS_ERROR} ${MIN_ALIGNMENT_LENGTH_CONSENSUS_CONSENSUS} ${CONSENSUS_CONSENSUS_ALIGNMENTS} ${CONSENSUS_CONSENSUS_ALIGNMENTS_FILTERED} ${N_READS} ${MAX_READ_LENGTH} ${READ_LENGTHS_FILE} ${READ_IDS_FILE} null null 0 0 $(( ${N_READS} - 1 )) ${N_ALIGNMENTS} ${QUALITY_THRESHOLDS_FILE} null 

# Computing statistics
N_READS=$(wc -l < ${PROJECT_DIR}/reads-lengths.txt)
java ${JAVA_RUNTIME_FLAGS} -classpath "${REVANT_BINARIES}" de.mpi_cbg.revant.consensus.ConsensusStep1 ${STEP1_DIR} ${FRAGMENT_CONSENSUS_ERROR} ${MIN_ALIGNMENT_LENGTH_CONSENSUS_CONSENSUS} ${CONSENSUS_SUFFIX} ${N_READS} ${MIN_ALIGNMENT_LENGTH_REPEAT_INFERENCE} ${BLOCK_CONSENSUS_ANALYSIS} ${CONSENSUS_CONSENSUS_ALIGNMENTS_FILTERED} ${FRAGMENT_DOMINATOR_ALIGNMENTS_FORMAT} ${FRAGMENT_CONSENSUS_ALIGNMENTS_FORMAT} ${CONSENSI_REPEAT_TRACK} ${CONSENSI_REPEAT_TRACK_FORMAT}
