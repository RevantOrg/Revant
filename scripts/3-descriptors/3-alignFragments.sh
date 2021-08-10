#!/bin/bash
#
# Computes all fragment-fragment and fragment-reference alignments.
#
# Remark: the script assumes environment variable $PATH$ to contain the path of all
# DAZZLER DB and DALIGNER executables.
#
# --------------------------------- CONFIGURATION ----------------------------------------
# This is the only section of the script that needs to be customized.
#
PROJECT_DIR=$1
# DALIGNER
IDENTITY="0.7"  # Must be >=0.7
MIN_ALIGNMENT_LENGTH="500"
N_THREADS="4"
MEMORY="16"  # GB
# ----------------------------------------------------------------------------------------




ROOT_DIR="${PROJECT_DIR}/step1/finalOutput/step4/step5"
FRAGMENTS_DIR="${ROOT_DIR}/fragments-strings"
REFERENCES_DIR="${ROOT_DIR}/references-strings"
FRAGMENTS_LIST="${ROOT_DIR}/list-fragments.txt"
REFERENCES_LIST="${ROOT_DIR}/list-references.txt"
ALIGNMENTS_DIR="${ROOT_DIR}/fragments-strings-alignments"
OUTPUT_PREFIX="test-basin"

rm -rf ${ALIGNMENTS_DIR}
mkdir ${ALIGNMENTS_DIR}

# Computing all fragment-fragment alignments
ls ${FRAGMENTS_DIR}/fragments-*.txt > ${FRAGMENTS_LIST}
for INPUT_FILE in $(cat ${FRAGMENTS_LIST}); do
    BASE_NAME=$(basename ${INPUT_FILE} .txt)
    OUTPUT_DIR="${ALIGNMENTS_DIR}/${OUTPUT_PREFIX}-${BASE_NAME}"
    mkdir ${OUTPUT_DIR}
    cp ${FRAGMENTS_DIR}/${BASE_NAME}.txt ${OUTPUT_DIR}/${BASE_NAME}.fasta
    cd ${OUTPUT_DIR}
    fasta2DB ${BASE_NAME}.db ${BASE_NAME}.fasta
    DBsplit ${BASE_NAME}.db
    DBdump -rh ${BASE_NAME}.db > output-DBdump.txt
    daligner -T${N_THREADS} -M${MEMORY} -e${IDENTITY} -l${MIN_ALIGNMENT_LENGTH} ${BASE_NAME}.db ${BASE_NAME}.db
    LAshow ${BASE_NAME}.db ${BASE_NAME}.${BASE_NAME}.las > LAshow.txt
    cd ..
done

# Computing all fragment-reference alignments
ls ${REFERENCES_DIR}/reference-*.txt > ${REFERENCES_LIST}
for INPUT_FILE in $(cat ${REFERENCES_LIST}); do
    BASE_NAME=$(basename ${INPUT_FILE} .txt)
    ID=${BASE_NAME#"reference-"}
    OUTPUT_DIR="${ALIGNMENTS_DIR}/${OUTPUT_PREFIX}-fragments-${ID}"
    if [ -d ${OUTPUT_DIR} ]; then
        cp ${REFERENCES_DIR}/${BASE_NAME}.txt ${OUTPUT_DIR}/${BASE_NAME}.fasta
        cd ${OUTPUT_DIR}
        fasta2DB ${BASE_NAME}.db ${BASE_NAME}.fasta
        daligner -T${N_THREADS} -M${MEMORY} -e${IDENTITY} -l${MIN_ALIGNMENT_LENGTH} fragments-${ID}.db ${BASE_NAME}.db
        LAshow fragments-${ID}.db ${BASE_NAME}.db fragments-${ID}.${BASE_NAME}.las > LAshow-fragments-reference.txt
        cd ..
	else
		# No fragment exists. This is not necessarily an error.
		:
    fi
done
