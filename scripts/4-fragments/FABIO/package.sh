#!/bin/bash
#
INPUT_DIR=$1
OUTPUT_DIR="${INPUT_DIR}/packaged"

rm -rf ${OUTPUT_DIR}
mkdir ${OUTPUT_DIR}
cd ${INPUT_DIR}
for MODULE in $(ls -d test-basin-fragments-*-*-*/); do
	PACKAGE_DIR="${OUTPUT_DIR}/${MODULE}"
	mkdir ${PACKAGE_DIR}
	cd ${MODULE}
	cp LAshow.txt ${PACKAGE_DIR}
	cp LAshow-fragments-reference.txt ${PACKAGE_DIR}
	cp output-DBdump.txt ${PACKAGE_DIR}
	cd ..
done
