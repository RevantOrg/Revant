#!/bin/bash
# 
# Main script that controls the entire alignments filtering pipeline.
#
# Remark: the script assumes that environment variable $REVANT_BINARIES$ is already set
# to the directory that contains REVANT's binaries, i.e. to the directory that contains
# the $de/mpi_cbg/revant/...$ subpath.
#
# --------------------------------- CONFIGURATION ----------------------------------------
# This is the only section of the script that needs to be customized.
#
INPUT_DIR=$1
# ---------------------------- Properties of the genome ----------------------------------
N_HAPLOTYPES="1"
HAPLOTYPE_COVERAGE="30"  # Of one haplotype
# ----------------------------- Properties of the reads ----------------------------------
MAX_ALIGNMENT_ERROR="0.3"
MIN_ALIGNMENT_LENGTH_READ_REPEAT="500"
MIN_ALIGNMENT_LENGTH_READ_READ="500"
# ----------------------------- Properties of CLR reads ----------------------------------
IS_CLR="0"
LOW_QUALITY_LENGTH="100"  # >=100. Break reads at low-quality regions >= this length.
LOW_QUALITY_TYPE="1"  # 1=replacement, 0=insertion.
# ---------------------------- Properties of the aligner ---------------------------------
MAX_SPACER_LENGTH="10000"  # 0=assume that the endpoints of periodic repeats are accurate
WOBBLE_LENGTH="100"  # 0=do not wobble.
# Good settings for mostly periodic genome: MAX_SPACER_LENGTH="10000"; WOBBLE_LENGTH="100"
# ------------------------- Properties of genome addresses -------------------------------
IDENTITY_THRESHOLD="100"  # For ambiguous characters in first/last block. >=WOBBLE_LENGTH.
DISTANCE_THRESHOLD=$(( ${IDENTITY_THRESHOLD} * 4 ))  # >=IDENTITY_THRESHOLD
MULTI_MODE="0"; ONEMER_FILTER="3"
# Good settings for a mostly periodic genome: MULTI_MODE="0"; ONEMER_FILTER="3"
# Good settings for a mostly nonperiodic genome: MULTI_MODE="1"; ONEMER_FILTER="2"
MAX_K_UNIQUE="8"  # Use k-mers up to this length as unique addresses
MIN_INTERSECTION_NONREPETITIVE="100000"  # Non-repetitive regions shorter than this n. of
# bps are not considered trustworthy addresses on the genome.
# Good settings for a mostly periodic genome: MIN_INTERSECTION_NONREPETITIVE="100000"
# Good settings for a mostly nonperiodic genome: MIN_INTERSECTION_NONREPETITIVE="500"
# ------------------------ Properties of alignment filters -------------------------------
ALIGNMENT_FILTERING_MODE="0"  # 0=loose, 1=tight, 2=tight with matching characters.
# ----------------------------------- Resources ------------------------------------------
N_THREADS="4"
JAVA_RUNTIME_FLAGS="-Xms2G -Xmx10G"
DELETE_TMP_FILES="1"
# ----------------------------------------------------------------------------------------

set -o pipefail; set -e; set -u
export JAVA_RUNTIME_FLAGS
if [ ${IS_CLR} -eq 1 ]; then
	BROKEN_READS="1"
	./0-breakReads.sh ${INPUT_DIR} ${LOW_QUALITY_LENGTH} ${N_THREADS} ${DELETE_TMP_FILES}
else
	BROKEN_READS="0"
fi
./1-buildAlphabet.sh ${INPUT_DIR} ${BROKEN_READS} ${MAX_ALIGNMENT_ERROR} ${MIN_ALIGNMENT_LENGTH_READ_REPEAT} ${N_HAPLOTYPES} ${HAPLOTYPE_COVERAGE} ${N_THREADS} ${DELETE_TMP_FILES} ${MAX_SPACER_LENGTH} ${WOBBLE_LENGTH}
MIN_K_FOR_DISAMBIGUATION="2"; MAX_K_FOR_DISAMBIGUATION="4"
./2-fixEndBlocks.sh ${INPUT_DIR} ${BROKEN_READS} ${HAPLOTYPE_COVERAGE} ${LOW_QUALITY_TYPE} ${MIN_K_FOR_DISAMBIGUATION} ${MAX_K_FOR_DISAMBIGUATION} ${N_THREADS} ${DELETE_TMP_FILES}
./3-getUniqueSubstrings.sh ${INPUT_DIR} ${N_HAPLOTYPES} ${HAPLOTYPE_COVERAGE} ${MAX_K_UNIQUE} ${ONEMER_FILTER} ${N_THREADS} ${DELETE_TMP_FILES} ${IDENTITY_THRESHOLD} ${DISTANCE_THRESHOLD}
./4-filterAlignments.sh ${INPUT_DIR} ${BROKEN_READS} ${MAX_SPACER_LENGTH} ${MIN_ALIGNMENT_LENGTH_READ_READ} ${MIN_ALIGNMENT_LENGTH_READ_REPEAT} ${MAX_K_UNIQUE} ${ALIGNMENT_FILTERING_MODE} ${MIN_INTERSECTION_NONREPETITIVE} ${N_THREADS} ${DELETE_TMP_FILES}
READ_LENGTHS_FILE="${INPUT_DIR}/reads-lengths.txt"
N_READS=$(wc -l < ${READ_LENGTHS_FILE})
GRAPH_DIR="${INPUT_DIR}/components-mode-${ALIGNMENT_FILTERING_MODE}"
rm -rf ${GRAPH_DIR}; mkdir ${GRAPH_DIR}
java ${JAVA_RUNTIME_FLAGS} -classpath "${REVANT_BINARIES}" de.mpi_cbg.revant.apps.BuildAssemblyGraph ${INPUT_DIR} ${N_READS} ${ALIGNMENT_FILTERING_MODE} 2 ${MAX_ALIGNMENT_ERROR} ${GRAPH_DIR}
