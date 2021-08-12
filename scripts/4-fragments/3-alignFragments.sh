#!/bin/bash
#
# Computes all alignments between the fragments and the new references built by the
# previous script.
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




ROOT_DIR="${PROJECT_DIR}/step1/finalOutput/step4/step5/fragments-strings-alignments"
FRAGMENTS_DIR="${ROOT_DIR}/fragments-strings-new"
REFERENCES_DIR="${ROOT_DIR}/references-strings-new"
NEW_ALIGNMENTS_DIR="${ROOT_DIR}/fragments-strings-alignments-new"
FRAGMENTS_LIST="${ROOT_DIR}/list-fragments.txt"
REFERENCES_LIST="${ROOT_DIR}/list-references.txt"
OUTPUT_PREFIX="test-basin"
LENGTHS_SUBDIR="lengths"

# Moving length files
rm -rf ${FRAGMENTS_DIR}/${LENGTHS_SUBDIR}/
mkdir ${FRAGMENTS_DIR}/${LENGTHS_SUBDIR}/
mv ${FRAGMENTS_DIR}/*-lengths.txt ${FRAGMENTS_DIR}/${LENGTHS_SUBDIR}/
rm -rf ${REFERENCES_DIR}/${LENGTHS_SUBDIR}/
mkdir ${REFERENCES_DIR}/${LENGTHS_SUBDIR}
mv ${REFERENCES_DIR}/*-lengths.txt ${REFERENCES_DIR}/${LENGTHS_SUBDIR}/

# Building fragment DBs
rm -rf ${NEW_ALIGNMENTS_DIR}
mkdir ${NEW_ALIGNMENTS_DIR}
ls ${FRAGMENTS_DIR}/fragments-*.txt > ${FRAGMENTS_LIST}
rm -rf ${ROOT_DIR}/${OUTPUT_PREFIX}*
for INPUT_FILE in $(cat ${FRAGMENTS_LIST}); do
    BASE_NAME=$(basename ${INPUT_FILE} .txt)
    OUTPUT_DIR="${NEW_ALIGNMENTS_DIR}/${OUTPUT_PREFIX}-${BASE_NAME}"
    mkdir ${OUTPUT_DIR}
    cp ${FRAGMENTS_DIR}/${BASE_NAME}.txt ${OUTPUT_DIR}/${BASE_NAME}.fasta
    cd ${OUTPUT_DIR}
    fasta2DB ${BASE_NAME}.db ${BASE_NAME}.fasta
    DBsplit ${BASE_NAME}.db
    cd ..
done

# Computing all fragment-reference alignments
ls ${REFERENCES_DIR}/reference-*.txt > ${REFERENCES_LIST}
for INPUT_FILE in $(cat ${REFERENCES_LIST}); do
    BASE_NAME=$(basename ${INPUT_FILE} .txt)
    ID=${BASE_NAME#"reference-"}
    OUTPUT_DIR="${NEW_ALIGNMENTS_DIR}/${OUTPUT_PREFIX}-fragments-${ID}"
    if [ -d ${OUTPUT_DIR} ]; then
        cp ${REFERENCES_DIR}/${BASE_NAME}.txt ${OUTPUT_DIR}/${BASE_NAME}.fasta
        cd ${OUTPUT_DIR}
        fasta2DB ${BASE_NAME}.db ${BASE_NAME}.fasta
        daligner -T${N_THREADS} -M${MEMORY} -e${IDENTITY} -l${MIN_ALIGNMENT_LENGTH} fragments-${ID}.db ${BASE_NAME}.db
        LAshow fragments-${ID}.db ${BASE_NAME}.db fragments-${ID}.${BASE_NAME}.las > LAshow-fragments-reference.txt
        cd ..
	else
		# Repeats with no fragment files are not processed.
		:
    fi
done

# Moving back length files
mv ${FRAGMENTS_DIR}/${LENGTHS_SUBDIR}/*-lengths.txt ${FRAGMENTS_DIR}
rm -rf ${FRAGMENTS_DIR}/${LENGTHS_SUBDIR}/
mv ${REFERENCES_DIR}/${LENGTHS_SUBDIR}/*-lengths.txt ${REFERENCES_DIR}
rm -rf ${REFERENCES_DIR}/${LENGTHS_SUBDIR}/
