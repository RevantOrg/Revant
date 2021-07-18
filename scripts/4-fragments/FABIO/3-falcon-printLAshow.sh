#!/bin/bash
#
# Prints the LAshow files from Pippel's fragment-reference alignments. Such files will be
# the input of FragmentsStep3. 
#
INPUT_DIR=$1
PATH="/projects/dazzler/pippel/prog/dazzlerGIT/TRACE_XOVR_125/bin:$PATH"
export PATH

cd ${INPUT_DIR}
for MODULE in $(ls -d *_*_*/); do
	MODULE=${MODULE%?}
	cd ${MODULE}
	DB_FILE="${MODULE}_Z.db"
	ALIGNMENTS_FILE="${MODULE}_Z.2.${MODULE}_Z.1.las"
	LAshow ${DB_FILE} ${DB_FILE} ${ALIGNMENTS_FILE} > "LAshow-fragments-reference.txt"
	ALIGNMENTS_FILE="${MODULE}_Z.1.${MODULE}_Z.2.las"
	LAshow ${DB_FILE} ${DB_FILE} ${ALIGNMENTS_FILE} > "LAshow-reference-fragments.txt"
	cd ..
done
