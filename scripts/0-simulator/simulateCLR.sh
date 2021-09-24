#!/bin/bash
#
GENOME_DIR=$1
GENOME_FILE=$2  # In Fasta format
OUTPUT_DIR=$3
COVERAGE=$4
MIN_READ_LENGTH=$5
AVG_READ_LENGTH=$6
MAX_READ_LENGTH=$7


PASS_DIR="/Users/ramseysnow/git/PaSS"
PBSIM_DIR="/Users/ramseysnow/git/pbsim-1.0.3/src"
#SimLoRD is assumed to have been installed via python
FASTQSIM_DIR="/Users/ramseysnow/git/FASTQsim_v2.0"
PBSIM2_DIR="/Users/ramseysnow/git/pbsim2/src"
LONGISLND_DIR="/Users/ramseysnow/git/longislnd"
SILICO_DIR="/Users/ramseysnow/git/SiLiCO"
NPBSS_DIR="/Users/ramseysnow/git/NPBSS_Octave-master"


N_THREADS="2"
GENOME_LENGTH=$(wc -c < ${GENOME_DIR}/${GENOME_FILE})
N_READS=$(( (${COVERAGE} * ${GENOME_LENGTH}) / ${AVG_READ_LENGTH} ))


# ----------------------------------------------------------------------------------------
# 2020 - PBSIM2
#
# ${PBSIM2_DIR}/pbsim --depth ${COVERAGE} --length-min ${MIN_READ_LENGTH} --length-max ${MAX_READ_LENGTH} --hmm_model ${PBSIM2_DIR}/../data/P6C4.model --length-mean ${AVG_READ_LENGTH} ${GENOME_DIR}/${GENOME_FILE}
# PBSIM2_OUTPUT_DIR="${OUTPUT_DIR}/PBSIM2"
# rm -rf ${PBSIM2_OUTPUT_DIR}
# mkdir ${PBSIM2_OUTPUT_DIR}
# mv ./sd_* ${PBSIM2_OUTPUT_DIR}


# ----------------------------------------------------------------------------------------
# 2019 - PaSS
#
# PASS_OUTPUT_DIR="${OUTPUT_DIR}/PaSS"
# rm -rf ${PASS_OUTPUT_DIR}
# mkdir ${PASS_OUTPUT_DIR}
# # Building index
# perl ${PASS_DIR}/pacbio_mkindex.pl ${GENOME_DIR}/${GENOME_FILE} ${PASS_OUTPUT_DIR}
# PASS_MACHINE="pacbio_RS"  # CLR
# #PASS_MACHINE="pacbio_sequel"  # CCS
# PASS_ERROR_MODEL="${PASS_DIR}/E.coli/ecoli.config"  # RS
# #PASS_ERROR_MODEL="${PASS_DIR}/C.elegan/elegan.config"  # RS
# #PASS_ERROR_MODEL="${PASS_DIR}/Arabidopsis/arab.config"  # Sequel
# #PASS_ERROR_MODEL="${PASS_DIR}/sim.config"  # Example dataset
# N_THREADS="2"
# ${PASS_DIR}/PaSS -list ${PASS_OUTPUT_DIR}/percentage.txt -index ${PASS_OUTPUT_DIR}/index -m ${PASS_MACHINE} -c ${PASS_ERROR_MODEL} -r ${N_READS} -t ${N_THREADS} -o ${PASS_OUTPUT_DIR}/reads-PaSS -d


# ----------------------------------------------------------------------------------------
# 2016 - SimLoRD
#
# SIMLORD_OUTPUT_DIR="${OUTPUT_DIR}/SimLoRD/"
# rm -rf ${SIMLORD_OUTPUT_DIR}
# mkdir ${SIMLORD_OUTPUT_DIR}
# simlord --read-reference ${GENOME_DIR}/${GENOME_FILE} --coverage ${COVERAGE} --max-passes 1 --min-readlength ${MIN_READ_LENGTH} ${SIMLORD_OUTPUT_DIR}
# mv ${SIMLORD_OUTPUT_DIR}/../SimLoRD* ${SIMLORD_OUTPUT_DIR}



# ----------------------------------------------------------------------------------------
# 2018 - NPBSS 
#
cd ${NPBSS_DIR}
octave --eval "pkg load statistics; NPBSS('${GENOME_DIR}/${GENOME_FILE}','-dep ${COVERAGE}','-min ${MIN_READ_LENGTH}','-max ${MAX_READ_LENGTH}','-len ${AVG_READ_LENGTH}')"
cd -
NPBSS_OUTPUT_DIR="${OUTPUT_DIR}/NPBSS/"
rm -rf ${NPBSS_OUTPUT_DIR}
mkdir ${NPBSS_OUTPUT_DIR}
mv ${GENOME_DIR}/*_CLR.* ${NPBSS_OUTPUT_DIR}









# ------------------------------ SIMULATORS CEMETERY -------------------------------------




# ----------------------------------------------------------------------------------------
# 2016 - LongISLND (Does not ship with a model, cumbersome to create one)
# ${LONGISLND_DIR}/simulate.py --fasta ${GENOME_DIR}/${GENOME_FILE} --model_dir ........


# ----------------------------------------------------------------------------------------
# 2016 - SiLiCO (Says it's successful but produces an empty file)
#
# SILICO_OUTPUT_DIR="${OUTPUT_DIR}/SiLiCO/"
# rm -rf ${SILICO_OUTPUT_DIR}
# mkdir ${SILICO_OUTPUT_DIR}
# python3 ${SILICO_DIR}/SiLiCO.py -i ${GENOME_DIR}/${GENOME_FILE} -o ${SILICO_OUTPUT_DIR} -m ${AVG_READ_LENGTH} -c ${COVERAGE} --pacbio


# ----------------------------------------------------------------------------------------
# 2014 - FASTQSim (Extremely slow, does not work with a reference of length 10k)
#
# FASTQSIM_OUTPUT_DIR="${OUTPUT_DIR}/FASTQSim/"
# rm -rf ${FASTQSIM_OUTPUT_DIR}
# mkdir ${FASTQSIM_OUTPUT_DIR}
# cd ${FASTQSIM_DIR}
# python ./src/pathoSpike.py -nobackground -platform pacbio -source ${COVERAGE} ${GENOME_DIR}/${GENOME_FILE} false -o ${FASTQSIM_OUTPUT_DIR} -threads ${N_THREADS}
# cd -


# ----------------------------------------------------------------------------------------
# 2012 - PBSIM (Superseded by PBSIM2)
#
# PBSIM_DATA_TYPE="CLR"
# #PBSIM_DATA_TYPE="CCS"
# PBSIM_ERROR_MODEL="${PBSIM_DIR}/../data/model_qc_clr"
# ${PBSIM_DIR}/pbsim --data-type ${PBSIM_DATA_TYPE} --depth ${COVERAGE} --length-min ${MIN_READ_LENGTH} --length-max ${MAX_READ_LENGTH} --model_qc ${PBSIM_ERROR_MODEL} --length-mean ${AVG_READ_LENGTH} ${GENOME_DIR}/${GENOME_FILE}
# PBSIM_OUTPUT_DIR="${OUTPUT_DIR}/PBSIM"
# rm -rf ${PBSIM_OUTPUT_DIR}
# mkdir ${PBSIM_OUTPUT_DIR}
# mv ./sd_* ${PBSIM_OUTPUT_DIR}