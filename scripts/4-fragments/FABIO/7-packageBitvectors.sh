#!/bin/bash
#
# Packages for transfer just the bitvectors built by $FragmentsStep3.java$.
#
STEP1_DIR=$1
OUTPUT_DIR=$2
INPUT_DIR="${STEP1_DIR}/finalOutput/step4/step5/fragments-strings-alignments/fragments-strings-alignments-new"

rm -rf ${OUTPUT_DIR}
mkdir ${OUTPUT_DIR}
cd ${INPUT_DIR}
for MODULE in $(ls -d test-basin-fragments-*/); do
	NEW_MODULE=${MODULE#"test-basin-fragments-"}
	NEW_MODULE=${NEW_MODULE//-/_}
	PACKAGE_DIR="${OUTPUT_DIR}/${NEW_MODULE}"
	mkdir ${PACKAGE_DIR}
	cd ${MODULE}
	cp LAshow-fragments-reference-keep.txt ${PACKAGE_DIR}
	cp LAshow-reference-fragments-keep.txt ${PACKAGE_DIR}
	cd ..
done
