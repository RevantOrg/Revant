#!/bin/bash
#
# Given a set of string labels of nodes of bidirected graphs, and the set of all sequences
# of reads that appear in some descriptor file, the script builds the string of every 
# descriptor reference (which is a path in some bidirected graph).
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
MIN_ALIGNMENT_LENGTH="1500"
MIN_N_FRAGMENTS="10"
# REVANT
JAVA_RUNTIME_FLAGS="-Xms2G -Xmx10G"
# ----------------------------------------------------------------------------------------




STEP1_DIR="${PROJECT_DIR}/step1"
REFERENCE_STRINGS_DIR="${STEP1_DIR}/finalOutput/step4/step5/references-strings"
if [ ! -d ${REFERENCE_STRINGS_DIR} ]; then
	mkdir ${REFERENCE_STRINGS_DIR}
fi
java ${JAVA_RUNTIME_FLAGS} -classpath "${REVANT_BINARIES}" de.mpi_cbg.revant.fragments.BasinDescriptor2Reference ${STEP1_DIR} ${MIN_N_FRAGMENTS} ${MIN_ALIGNMENT_LENGTH}