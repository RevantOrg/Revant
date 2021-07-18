#!/bin/bash
#
# Packages for transfer just the file built by <pippel-printLAshow.sh>.
#
INPUT_DIR=$1
OUTPUT_DIR="${INPUT_DIR}/packaged"

rm -rf ${OUTPUT_DIR}
mkdir ${OUTPUT_DIR}
cd ${INPUT_DIR}
for MODULE in $(ls -d *_*_*/); do
	PACKAGE_DIR="${OUTPUT_DIR}/${MODULE}"
	mkdir ${PACKAGE_DIR}
	cd ${MODULE}
	cp LAshow-fragments-reference.txt ${PACKAGE_DIR}
	cp LAshow-reference-fragments.txt ${PACKAGE_DIR}
	cd ..
done
