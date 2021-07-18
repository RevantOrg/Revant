#!/bin/bash
#
INPUT_DIR=$1
OUTPUT_DIR="${INPUT_DIR}/packaged"

rm -rf ${OUTPUT_DIR}
mkdir ${OUTPUT_DIR}
cd ${INPUT_DIR}
for MODULE in $(ls -d *_*_*/); do
	PACKAGE_DIR="${OUTPUT_DIR}/${MODULE//_/-}"
	mkdir ${PACKAGE_DIR}
	cd ${MODULE}
	cp LAshow-fragments-reference.txt ${PACKAGE_DIR}
	cp LAshow-reference-fragments.txt ${PACKAGE_DIR}
	cp LAshow-fragments-reference-keep.txt ${PACKAGE_DIR}
	cp LAshow-reference-fragments-keep.txt ${PACKAGE_DIR}
	cd ..
done
