#!/bin/bash
#
# Checks the consistency of consensus strings by aligning them to: (1) the fragments used
# to compute each consensus; (2) every other consensus; (3) the blocks of the database
# used to infer the consensi.
#
# Remark: the script assumes that consensus strings use "_" instead of "-" in their
# filenames.
#
# Remark: the script assumes environment variable $PATH$ to contain the path of all
# DAZZLER DB and DALIGNER executables.
#
# --------------------------------- CONFIGURATION ----------------------------------------
# This is the only section of the script that needs to be customized.
#
FRAGMENTS_DIR=$1  # Contains DAZZLER DBs
CONSENSI_DIR=$2  # Contains strings
OUTPUT_ROOT=$3
CONSENSUS_SUFFIX=".filt.daccord.consensus"  # Excluding ".fasta"
CONSENSI_LIST="list-consensi.txt"
OUTPUT_PREFIX="test-basin"
# Full DB with blocks
BLOCK_DB=$4  # Path of the .db file
FIRST_BLOCK=$5
LAST_BLOCK=$6
# DALIGNER
IDENTITY_FF="0.7"  # Fragment-fragment. Must be >=0.7.
IDENTITY_CC="0.85"  # Consensus-consensus
MIN_ALIGNMENT_LENGTH="500"
N_THREADS="8"
MEMORY="16"  # GB
# ----------------------------------------------------------------------------------------




# Fragment-consensus alignments
ls ${CONSENSI_DIR}/*${CONSENSUS_SUFFIX}.fasta > ${OUTPUT_ROOT}/${CONSENSI_LIST}
for INPUT_FILE in $(cat ${OUTPUT_ROOT}/${CONSENSI_LIST}); do
    BASE_NAME=$(basename ${INPUT_FILE} ${CONSENSUS_SUFFIX}.fasta)
    ID=${BASE_NAME//[_]/-}
    OUTPUT_DIR="${OUTPUT_PREFIX}-fragments-${ID}"
    if [ -d ${OUTPUT_ROOT}/${OUTPUT_DIR} ]; then
		rm -rf ${OUTPUT_ROOT}/${OUTPUT_DIR}
    fi
    mkdir ${OUTPUT_ROOT}/${OUTPUT_DIR}
    cp ${CONSENSI_DIR}/${INPUT_FILE} ${OUTPUT_ROOT}/${OUTPUT_DIR}
    cd ${OUTPUT_ROOT}/${OUTPUT_DIR}
    fasta2DB ${BASE_NAME}.db ${BASE_NAME}${CONSENSUS_SUFFIX}.fasta
    daligner  -T${N_THREADS} -M${MEMORY} -e${IDENTITY_FF} -l${MIN_ALIGNMENT_LENGTH} ${FRAGMENTS_DIR}/${OUTPUT_DIR}/fragments-${ID}.db ${BASE_NAME}.db
    LAshow ${FRAGMENTS_DIR}/${OUTPUT_DIR}/fragments-${ID}.db ./${BASE_NAME}.db fragments-${ID}.${BASE_NAME}.las > LAshow-fragments-consensus.txt
    cd ..
done

# Consensus-consensus alignments
ALL_CONSENSI="all-consensi"
rm -rf ${OUTPUT_ROOT}/${ALL_CONSENSI}.fasta
for INPUT_FILE in $(cat ${OUTPUT_ROOT}/${CONSENSI_LIST}); do
	cat ${INPUT_FILE} >> ${OUTPUT_ROOT}/${ALL_CONSENSI}.fasta
done
cd ${OUTPUT_ROOT}
fasta2DB ${ALL_CONSENSI}.db ${ALL_CONSENSI}.fasta
daligner -T${MEMORY} -M${N_THREADS} -e${IDENTITY_CC} -l${MIN_ALIGNMENT_LENGTH} ${ALL_CONSENSI}.db ${ALL_CONSENSI}.db
LAshow ${ALL_CONSENSI}.db ${ALL_CONSENSI}.db ${ALL_CONSENSI}.${ALL_CONSENSI}.las > LAshow-consensi-consensi.txt
cd ..

# Alignments between each consensus and blocks of the database.
# Remark: this is sequential just to exemplify; in practice it should be run in parallel
# since there are many repeats and many blocks.
for INPUT_FILE in $(cat ${OUTPUT_ROOT}/${CONSENSI_LIST}); do
	BASE_NAME=$(basename ${INPUT_FILE} ${CONSENSUS_SUFFIX}.fasta)
	ID=${BASE_NAME//[_]/-}
	OUTPUT_DIR="${OUTPUT_ROOT}/${OUTPUT_PREFIX}-fragments-${ID}"
	BLOCK_DB_PREFIX=${BLOCK_DB%.db}
	cd ${OUTPUT_DIR}
	for BLOCK in $(seq ${FIRST_BLOCK} ${LAST_BLOCK}); do
	    daligner -T${N_THREADS} -M${MEMORY} -e${IDENTITY_FF} -l${MIN_ALIGNMENT_LENGTH} ${BLOCK_DB_PREFIX}.${BLOCK}.db ./${BASE_NAME}.db
	done
	LAmerge LAshow-block-consensus.las ${BLOCK_DB_PREFIX}.@${FIRST_BLOCK}-${LAST_BLOCK}.${BASE_NAME}.las
	LAshow ${BLOCK_DB} ./${BASE_NAME}.db LAshow-block-consensus.las > LAshow-block-consensus.txt
	cd ..
done
