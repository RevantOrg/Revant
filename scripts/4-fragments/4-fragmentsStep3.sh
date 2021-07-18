#!/bin/bash
#
# Assume that we are given the set of all alignments of fragments to the new references 
# (computed by the previous script), and assume that we want to use such star alignments
# to compute the consensus sequence of every repeat. This script keeps just some of those
# alignments, selected in such a way that the same character from a fragment votes just
# once for the sequence of the consensus. The script keeps the original alignment files
# intact, and it outputs a text file with a zero or a one at line $i$ if the $i$-th
# alignment should be kept.
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
LASHOW_FORMAT="0"  # Pippel's format
# REVANT
JAVA_RUNTIME_FLAGS="-Xms2G -Xmx10G"
# ----------------------------------------------------------------------------------------




STEP1_DIR="${PROJECT_DIR}/step1"
FILTER_DIR="${STEP1_DIR}/finalOutput/step4/step5/fragments-strings-alignments/fragments-strings-alignments-new"

if [ ! -d ${FILTER_DIR}/stats ]; then
	mkdir ${FILTER_DIR}/stats
fi
java ${JAVA_RUNTIME_FLAGS} -classpath "${REVANT_BINARIES}" de.mpi_cbg.revant.fragments.FragmentsStep3 ${STEP1_DIR} ${MIN_ALIGNMENT_LENGTH} 1 0 ${LASHOW_FORMAT}
