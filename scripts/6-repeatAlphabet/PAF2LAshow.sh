#!/bin/bash
#
INPUT_DIR=$1
MIN_ALIGNMENT_LENGTH_READ_READ="500"
MIN_ALIGNMENT_LENGTH_READ_REPEAT="500"
MAX_ERROR_RATE_READ_READ="0.02"
MAX_ERROR_RATE_READ_REPEAT="0.3"
N_THREADS="4"

REVANT_BIN="/Users/ramseysnow/git/Revant/bin"
LENGTHS_FILE="${INPUT_DIR}/reads-lengths-pre.txt"
CONVERSION_FILE="${INPUT_DIR}/conversion.txt"
TMP_1="${INPUT_DIR}/tmp1.txt"
TMP_2="${INPUT_DIR}/tmp2.txt"

cut -f 1,2 ${INPUT_DIR}/reads-reads.paf | sort --parallel ${N_THREADS} | uniq > ${TMP_1}
cut -f 6,7 ${INPUT_DIR}/reads-reads.paf | sort --parallel ${N_THREADS} | uniq > ${TMP_2}
sort --parallel ${N_THREADS} -m ${TMP_1} ${TMP_2} | uniq > ${LENGTHS_FILE}
rm -f ${TMP_1} ${TMP_2}
 
java -cp ${REVANT_BIN} de.mpi_cbg.revant.factorize.PAF2LAshowReadRead ${INPUT_DIR}/reads-reads.paf ${LENGTHS_FILE} ${INPUT_DIR}/LAshow-reads-reads.txt ${INPUT_DIR}/reads-lengths.txt ${MIN_ALIGNMENT_LENGTH_READ_READ} ${MAX_ERROR_RATE_READ_READ} > ${CONVERSION_FILE}
java -cp ${REVANT_BIN} de.mpi_cbg.revant.factorize.PAF2LAshowReadRepeat ${INPUT_DIR}/reads-repeats.paf $( wc -l < ${INPUT_DIR}/reads-repeats.paf ) ${CONVERSION_FILE} ${INPUT_DIR}/LAshow-reads-repeats.txt ${MIN_ALIGNMENT_LENGTH_READ_REPEAT} ${MAX_ERROR_RATE_READ_REPEAT}