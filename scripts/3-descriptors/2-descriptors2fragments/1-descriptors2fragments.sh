#!/bin/bash
#
# Assume that we have a list of descriptor files, and that we want to convert every
# interval in every descriptor file into strings. The script builds a list of all the
# substrings of reads that are needed for the conversion.
#
# Remark: at the end of the descriptors stage, there could be more reference files
# than fragment files. This is not necessarily an error.
#
# Remark: the script assumes that environment variable $REVANT_BINARIES$ is already set
# to the directory that contains REVANT's binaries, i.e. to the directory that contains
# the $de/mpi_cbg/revant/...$ subpath.
#
# --------------------------------- CONFIGURATION ----------------------------------------
# This is the only part of the script that needs to be customized.
#
PROJECT_DIR=$1
MIN_ALIGNMENT_LENGTH="1500"  # In repeat inference
MIN_N_FRAGMENTS="10"
# REVANT
JAVA_RUNTIME_FLAGS="-Xms2G -Xmx10G"
# ----------------------------------------------------------------------------------------




STEP1_DIR="${PROJECT_DIR}/step1"
FRAGMENTS_DIR="${STEP1_DIR}/finalOutput/step4/step5/fragments"
if [ ! -d ${FRAGMENTS_DIR} ]; then
	mkdir ${FRAGMENTS_DIR}
fi
java ${JAVA_RUNTIME_FLAGS} -classpath "${REVANT_BINARIES}" de.mpi_cbg.revant.fragments.BasinDescriptor2Fragments ${STEP1_DIR} ${MIN_N_FRAGMENTS} ${MIN_ALIGNMENT_LENGTH} 0 1
