#!/bin/bash
#
STEP1_DIR=$1
FILTER_DIR="${STEP1_DIR}/finalOutput/step4/step5/fragments-strings-alignments/fragments-strings-alignments-new"
MIN_ALIGNMENT_LENGTH="1500"
LASHOW_FORMAT="0"  # Pippel's format
CODE_DIR="/Users/ramseysnow/Dropbox/dropbox private/collaborators/myers/NEW/sw"

if [ ! -d ${FILTER_DIR}/stats ]; then
	mkdir ${FILTER_DIR}/stats
fi
cd "${CODE_DIR}"
java FragmentsStep3 ${STEP1_DIR} ${MIN_ALIGNMENT_LENGTH} 1 0 ${LASHOW_FORMAT}
cd -
