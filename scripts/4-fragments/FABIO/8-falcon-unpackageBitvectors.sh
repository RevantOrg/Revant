#!/bin/bash
#
# Moves the packaged files built by <pippel-printLAshow.sh> to their destination 
# directories.
#
INPUT_DIR=$1
PACKAGE_DIR="${INPUT_DIR}/packaged"

cd ${INPUT_DIR}
for MODULE in $(ls -d *_*_*/); do
	cp ${PACKAGE_DIR}/${MODULE}/LAshow-fragments-reference-keep.txt ${INPUT_DIR}/${MODULE}/
	cp ${PACKAGE_DIR}/${MODULE}/LAshow-reference-fragments-keep.txt ${INPUT_DIR}/${MODULE}/
done
