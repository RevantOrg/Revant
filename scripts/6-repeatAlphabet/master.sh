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
GENOME_LENGTH=$(wc -c < ${INPUT_DIR}/genome.fasta)
N_HAPLOTYPES="1"
# ----------------------------- Properties of the reads ----------------------------------
MAX_ALIGNMENT_ERROR="0.3"
MIN_ALIGNMENT_LENGTH_READ_REPEAT="500"
MIN_ALIGNMENT_LENGTH_READ_READ="500"
AVG_READ_LENGTH="30000"
READ_READ_ALIGNMENTS_FILE="${INPUT_DIR}/LAshow-reads-reads.txt"
# ----------------------------- Properties of CLR reads ----------------------------------
IS_CLR="0"
LOW_QUALITY_LENGTH="100"  # >=100. Break reads at low-quality regions >= this length.
LOW_QUALITY_TYPE="1"  # 1=replacement, 0=insertion.
# ---------------------------- Properties of the aligner ---------------------------------
MAX_SPACER_LENGTH="10000"  # 0=assume that the endpoints of periodic repeats are accurate
WOBBLE_LENGTH="100"  # 0=do not wobble.
# Good settings for mostly periodic genome: MAX_SPACER_LENGTH="10000"; WOBBLE_LENGTH="100"
FIX_TANDEM_SPACERS="1"
CONCATENATE_BLOCKS="1"
# ------------------------- Properties of genome addresses -------------------------------
IDENTITY_THRESHOLD="100"  # For ambiguous characters in first/last block. >=WOBBLE_LENGTH.
DISTANCE_THRESHOLD=$(( ${IDENTITY_THRESHOLD} * 4 ))  # >=IDENTITY_THRESHOLD
CHARACTER_THRESHOLD="0.9"  # For ambiguous characters in first/last block. Arbitrary.
MAX_K_UNIQUE="6"  # Use k-mers up to this length as unique addresses
MIN_INTERSECTION_NONREPETITIVE="100000"  # Non-repetitive regions shorter than this n. of
# bps are not considered trustworthy addresses on the genome.
# Good settings for a mostly periodic genome: MIN_INTERSECTION_NONREPETITIVE="100000"
# Good settings for a mostly nonperiodic genome: MIN_INTERSECTION_NONREPETITIVE="500"
# ------------------------ Properties of alignment filters -------------------------------
ALIGNMENT_FILTERING_MODE="0"  # 0=loose, 1=tight, 2=tight with matching characters.
# ------------------------ Properties of the assembly graph ------------------------------
SIMPLIFY_ASSEMBLY_GRAPH="0"
# ----------------------------------- Resources ------------------------------------------
N_THREADS="1"
JAVA_RUNTIME_FLAGS="-Xms2G -Xmx10G"
DELETE_TMP_FILES="0"
# ----------------------------------------------------------------------------------------

set -euo pipefail
export JAVA_RUNTIME_FLAGS
MAX_KMER_LENGTH_BPS=$(java ${JAVA_RUNTIME_FLAGS} -classpath "${REVANT_BINARIES}" de.mpi_cbg.revant.apps.GetAlignmentLengthThreshold ${READ_READ_ALIGNMENTS_FILE} 0 ${AVG_READ_LENGTH} ${INPUT_DIR}/histogram-alignmentLength-readRead.txt)
echo "Histogram of suffix-prefix read-read alignment lengths:"
cat ${INPUT_DIR}/histogram-alignmentLength-readRead.txt
if [ ${IS_CLR} -eq 1 ]; then
	BROKEN_READS="1"
	./0-breakReads.sh ${INPUT_DIR} ${LOW_QUALITY_LENGTH} ${N_THREADS} ${DELETE_TMP_FILES}
else
	BROKEN_READS="0"
fi
./1-buildAlphabet.sh ${INPUT_DIR} ${BROKEN_READS} ${MAX_ALIGNMENT_ERROR} ${MIN_ALIGNMENT_LENGTH_READ_REPEAT} ${N_HAPLOTYPES} ${N_THREADS} ${DELETE_TMP_FILES} ${MAX_SPACER_LENGTH} ${WOBBLE_LENGTH} ${FIX_TANDEM_SPACERS} ${CONCATENATE_BLOCKS} ${AVG_READ_LENGTH} ${GENOME_LENGTH}
PERIODIC_ENDPOINTS_FIXED=$(cat ${INPUT_DIR}/buildAlphabet-tmp-return.txt)
MIN_K_FOR_DISAMBIGUATION="2"; MAX_K_FOR_DISAMBIGUATION="4"
./2-fixEndBlocks.sh ${INPUT_DIR} ${BROKEN_READS} ${LOW_QUALITY_TYPE} ${MIN_K_FOR_DISAMBIGUATION} ${MAX_K_FOR_DISAMBIGUATION} ${N_THREADS} ${DELETE_TMP_FILES} ${GENOME_LENGTH} ${N_HAPLOTYPES} ${MAX_KMER_LENGTH_BPS}
./3-getUniqueSubstrings.sh ${INPUT_DIR} ${GENOME_LENGTH} ${N_HAPLOTYPES} ${MAX_K_UNIQUE} ${N_THREADS} ${DELETE_TMP_FILES} ${IDENTITY_THRESHOLD} ${DISTANCE_THRESHOLD} ${CHARACTER_THRESHOLD} ${MIN_ALIGNMENT_LENGTH_READ_REPEAT} ${MAX_KMER_LENGTH_BPS}
./4-filterAlignments.sh ${INPUT_DIR} ${BROKEN_READS} ${PERIODIC_ENDPOINTS_FIXED} ${MIN_ALIGNMENT_LENGTH_READ_READ} ${MIN_ALIGNMENT_LENGTH_READ_REPEAT} ${MAX_K_UNIQUE} ${ALIGNMENT_FILTERING_MODE} ${MIN_INTERSECTION_NONREPETITIVE} ${N_THREADS} ${DELETE_TMP_FILES}
READ_LENGTHS_FILE="${INPUT_DIR}/reads-lengths.txt"
N_READS=$(wc -l < ${READ_LENGTHS_FILE})
GRAPH_DIR="${INPUT_DIR}/components-mode-${ALIGNMENT_FILTERING_MODE}"
rm -rf ${GRAPH_DIR}; mkdir ${GRAPH_DIR}
java ${JAVA_RUNTIME_FLAGS} -classpath "${REVANT_BINARIES}" de.mpi_cbg.revant.apps.BuildAssemblyGraph ${INPUT_DIR} ${N_READS} ${ALIGNMENT_FILTERING_MODE} 2 ${MAX_ALIGNMENT_ERROR} ${AVG_READ_LENGTH} ${SIMPLIFY_ASSEMBLY_GRAPH} ${GRAPH_DIR}
