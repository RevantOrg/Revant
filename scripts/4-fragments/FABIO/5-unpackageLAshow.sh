#!/bin/bash
#
# Copies the LAshow files created by <pippel-printLAshow.sh> to some directory.
#
FROM_DIR=$1
STEP1_DIR=$2
TO_DIR="${STEP1_DIR}/finalOutput/step4/step5/fragments-strings-alignments/fragments-strings-alignments-new"
TO_PREFIX="test-basin-fragments"
LASHOW_FILE_1="LAshow-fragments-reference.txt"
LASHOW_FILE_2="LAshow-reference-fragments.txt"

cd ${FROM_DIR}
for MODULE in $(ls -d *_*_*/); do
	MODULE=${MODULE%?}
	DEST_DIR="${TO_DIR}/${TO_PREFIX}-${MODULE//_/-}/"
	if [ ! -d ${DEST_DIR} ]; then
		mkdir ${DEST_DIR}
	fi
	cp "${MODULE}/${LASHOW_FILE_1}" ${DEST_DIR}
	cp "${MODULE}/${LASHOW_FILE_2}" ${DEST_DIR}
done
cd -