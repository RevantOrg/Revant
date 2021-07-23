#!/bin/bash
# 
# Given the directory of a project, the script reads its $input$ subfolder and converts
# its files in the format expected by the factorization stage. The result of the script
# is a set of files in the main directory of the project. The script can work on just a
# compact range of all the reads in the dataset.
#
# Remark: the script assumes that environment variable $REVANT_BINARIES$ is already set
# to the directory that contains REVANT's binaries, i.e. to the directory that contains
# the $de/mpi_cbg/revant/...$ subpath.
#
# Remark: every input file is a human-readable text file. $DAZZLER_*$ track files are
# text files produced by DAZZLER's $DBdump -r$ tool. $DAMAR_*$ are text files produced by
# DAMAR tools. Repeat tracks are text files as well.
#
# --------------------------------- CONFIGURATION ----------------------------------------
# This is the only section of the script that needs to be customized.
#
PROJECT_DIR=$1  # Contains the $input$ folder
FIRST_READ_IN_RANGE="971409"  # Start of read range, zero-based.
LAST_READ_IN_RANGE="1978662"  # End of read range, zero-based.
MAX_ALIGNMENT_ERROR="0.3"  # Alignments with more error than this are discarded
MIN_ALIGNMENT_LENGTH="1500"  # Alignments shorter than this are discarded
MAX_READ_LENGTH="200000"  # Coarse upper bound
SPLIT_IN_PARTS="4"  # For factorizing in parallel
# Track files
INPUT_DIR="${PROJECT_DIR}/input"
DAMAR_QUALITY_TRACK="${INPUT_DIR}/DAMAR_Qtrack.txt"  # Use "null" to discard it.
DAZZLER_QUALITY_TRACK="null"  # Use "null" to discard it.
DAZZLER_DUST_TRACK="null"  # Use "null" to discard it.
DAZZLER_TANDEM_TRACK="null"  # Use "null" to discard it.
REPEAT_TRACK_DIR="${INPUT_DIR}"  # Use "null" to discard it.
REPEAT_TRACK_PREFIX="repmod"  # Use "null" to discard it.
REPEAT_TRACK_FORMAT="0"  # 0=DAZZLER, 1=DAMAR format
VERBOSE="0"  # 0/1
# REVANT
JAVA_RUNTIME_FLAGS="-Xms2G -Xmx10G"
# ----------------------------------------------------------------------------------------




ALIGNMENTS_FILE="${INPUT_DIR}/LAshow.txt"
QUALITY_THRESHOLDS_FILE="${INPUT_DIR}/qualityThresholds.txt"
DBDUMP_READ_LENGTHS="${INPUT_DIR}/DBdump-lengths.txt"
READ_LENGTHS_FILE="${PROJECT_DIR}/reads-lengths.txt"
READ_IDS_FILE="${PROJECT_DIR}/reads-ids.txt"
FILTERED_ALIGNMENTS_FILE="${PROJECT_DIR}/LAshow.txt"
OUTPUT_QUALITY_THRESHOLDS_FILE="${PROJECT_DIR}/qualityThresholds.txt"
MAX_QUALITY_FOR_HISTOGRAM="51"  # Just for display purposes
MIN_QUALITY_FOR_MERGE="45"  # We try to replace (with a good value) all quality values this much or more that are in tandem but not in dust.

echo "Building the read lengths file..."
java ${JAVA_RUNTIME_FLAGS} -classpath "${REVANT_BINARIES}" de.mpi_cbg.revant.factorize.DBdump2ReadIDs ${DBDUMP_READ_LENGTHS} > ${READ_IDS_FILE}
java ${JAVA_RUNTIME_FLAGS} -classpath "${REVANT_BINARIES}" de.mpi_cbg.revant.factorize.DBdump2ReadLengths ${DBDUMP_READ_LENGTHS} > ${READ_LENGTHS_FILE}
N_READS=$(wc -l < ${READ_LENGTHS_FILE})
echo "Read lengths file built"

echo "Building the quality track file..."
PHRED_FILE="null"
if [ ${DAMAR_QUALITY_TRACK} != "null" ]; then  # Priority to DAMAR
	PHRED_FILE="${PROJECT_DIR}/reads-phred.damar"
	java ${JAVA_RUNTIME_FLAGS} -classpath "${REVANT_BINARIES}" de.mpi_cbg.revant.factorize.Damar2Phred ${DAMAR_QUALITY_TRACK} ${READ_LENGTHS_FILE} ${READ_IDS_FILE} ${N_READS} ${QUALITY_THRESHOLDS_FILE} > ${PHRED_FILE}
	echo "Quality histogram:"
	java ${JAVA_RUNTIME_FLAGS} -classpath "${REVANT_BINARIES}" de.mpi_cbg.revant.factorize.QualityHistogram ${PHRED_FILE} ${MAX_QUALITY_FOR_HISTOGRAM}
elif [ ${DAZZLER_QUALITY_TRACK} != "null" ]; then
	PHRED_FILE="${PROJECT_DIR}/reads-phred.dbdump"
	java ${JAVA_RUNTIME_FLAGS} -classpath "${REVANT_BINARIES}" de.mpi_cbg.revant.factorize.DBdump2Phred ${DAZZLER_QUALITY_TRACK} > ${PHRED_FILE}
	echo "Quality histogram:"
	java ${JAVA_RUNTIME_FLAGS} -classpath "${REVANT_BINARIES}" de.mpi_cbg.revant.factorize.QualityHistogram ${PHRED_FILE} ${MAX_QUALITY_FOR_HISTOGRAM}
else
	echo "Input for quality track construction not found"
	#exit 1
fi
DUST_FILE="${PROJECT_DIR}/reads-dust.txt"
if [ ${DAZZLER_DUST_TRACK} != "null" ]; then
	java ${JAVA_RUNTIME_FLAGS} -classpath "${REVANT_BINARIES}" de.mpi_cbg.revant.factorize.DBdump2Track ${DAZZLER_DUST_TRACK} > ${DUST_FILE}
fi
TANDEM_FILE="${PROJECT_DIR}/reads-tandem.txt"
if [ ${DAZZLER_TANDEM_TRACK} != "null" ]; then
	java ${JAVA_RUNTIME_FLAGS} -classpath "${REVANT_BINARIES}" de.mpi_cbg.revant.factorize.DBdump2Track ${DAZZLER_TANDEM_TRACK} > ${TANDEM_FILE}
fi
if [ ${DAZZLER_DUST_TRACK} != "null" -a ${DAZZLER_TANDEM_TRACK} != "null" ]; then
	TMP_FILE="${PROJECT_DIR}/tmp.txt"
	java ${JAVA_RUNTIME_FLAGS} -classpath "${REVANT_BINARIES}" de.mpi_cbg.revant.factorize.AddTracks2Phred ${MIN_QUALITY_FOR_MERGE} ${PHRED_FILE} ${QUALITY_THRESHOLDS_FILE} ${DUST_FILE} ${TANDEM_FILE} ${TMP_FILE}
	mv ${PHRED_FILE} "${PHRED_FILE}.beforeMerge"
	mv ${TMP_FILE} ${PHRED_FILE}
	echo "Quality histogram after merge:"
	java ${JAVA_RUNTIME_FLAGS} -classpath "${REVANT_BINARIES}" de.mpi_cbg.revant.factorize.QualityHistogram ${PHRED_FILE} ${MAX_QUALITY_FOR_HISTOGRAM}
fi
cp ${QUALITY_THRESHOLDS_FILE} ${OUTPUT_QUALITY_THRESHOLDS_FILE}
echo "Quality track file built: ${PHRED_FILE}"

echo "Filtering and splitting alignments..."
N_ALIGNMENTS=$(( $(wc -l < ${ALIGNMENTS_FILE}) - 2 ))
java ${JAVA_RUNTIME_FLAGS} -classpath "${REVANT_BINARIES}" de.mpi_cbg.revant.factorize.FilterAlignments ${MAX_ALIGNMENT_ERROR} ${MIN_ALIGNMENT_LENGTH} ${ALIGNMENTS_FILE} ${FILTERED_ALIGNMENTS_FILE} ${N_READS} ${MAX_READ_LENGTH} ${READ_LENGTHS_FILE} ${READ_IDS_FILE} ${REPEAT_TRACK_DIR} ${REPEAT_TRACK_PREFIX} ${REPEAT_TRACK_FORMAT} ${FIRST_READ_IN_RANGE} ${LAST_READ_IN_RANGE} ${N_ALIGNMENTS} ${QUALITY_THRESHOLDS_FILE} ${PHRED_FILE} ${VERBOSE}
EXIT_STATUS=$?
if [ ${EXIT_STATUS} -ne 0 ]; then
    echo "The following command returned error ${EXIT_STATUS}:"
	echo "java ${JAVA_RUNTIME_FLAGS} -classpath "${REVANT_BINARIES}" de.mpi_cbg.revant.factorize.FilterAlignments ${MAX_ALIGNMENT_ERROR} ${MIN_ALIGNMENT_LENGTH} ${ALIGNMENTS_FILE} ${FILTERED_ALIGNMENTS_FILE} ${N_READS} ${MAX_READ_LENGTH} ${READ_LENGTHS_FILE} ${READ_IDS_FILE} ${REPEAT_TRACK_DIR} ${REPEAT_TRACK_PREFIX} ${FIRST_READ_IN_RANGE} ${LAST_READ_IN_RANGE} ${N_ALIGNMENTS}"
	exit ${EXIT_STATUS}
fi
N_ALIGNMENTS=$(( $(wc -l < ${FILTERED_ALIGNMENTS_FILE}) - 2 ))
PARTS_PREFIX="${PROJECT_DIR}/LAshow-"
LAST_READA_FILE="${PROJECT_DIR}/LAshow-lastReadA.txt"
java ${JAVA_RUNTIME_FLAGS} -classpath "${REVANT_BINARIES}" de.mpi_cbg.revant.factorize.SplitAlignments ${N_ALIGNMENTS} ${SPLIT_IN_PARTS} ${FILTERED_ALIGNMENTS_FILE} ${PARTS_PREFIX} ${LAST_READA_FILE}
echo "Alignments filtered and split in ${SPLIT_IN_PARTS} parts"
