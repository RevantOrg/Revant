#!/bin/bash
#
#
#
#
# Remark: the script assumes that environment variable $REVANT_BINARIES$ is already set
# to the directory that contains REVANT's binaries, i.e. to the directory that contains
# the $de/mpi_cbg/revant/...$ subpath.
#
# --------------------------------- CONFIGURATION ----------------------------------------
# This is the only section of the script that needs to be customized.
#
INPUT_DIR=$1
REPEAT_MODEL_DIR=$2
READ_SIMULATOR=$3
MIN_ALIGNMENT_LENGTH="1000"  # Arbitrary
NONREPETITIVE_BPS=$(( ${MIN_ALIGNMENT_LENGTH} / 2 ))
REPETITIVE_NONREPETITIVE_RATIO="2000"  # In the genome
COVERAGE="30"
MIN_READ_LENGTH="10000"
AVG_READ_LENGTH="15000"
MAX_READ_LENGTH="30000"
READ_READ_IDENTITY=".75"
N_THREADS="4"
MAX_MEMORY="10"  # GB
DALIGNER_DIR="/Users/ramseysnow/git/DALIGNER"
DAZZDB_DIR="/Users/ramseysnow/git/DAZZ_DB"
# Read simulators
PBSIM_DIR="/Users/ramseysnow/git/pbsim2/src"  #"/Users/ramseysnow/git/pbsim-1.0.3/src"
PASS_DIR="/Users/ramseysnow/git/PaSS"
NPBSS_DIR="/Users/ramseysnow/git/NPBSS_Octave-master"
# SimLoRD is assumed to have been installed via Python
# REVANT
JAVA_RUNTIME_FLAGS="-Xms2G -Xmx10G"
# ----------------------------------------------------------------------------------------

set -o pipefail; set -e; set -u

# Building the genome
GENOME_FILE="${GENOME_DIR}/genome.fasta"
GENOME_IMAGE="${GENOME_DIR}/genome.png"
READS_FILE="${GENOME_DIR}/reads.fasta"
REPEATS_FILE="${GENOME_DIR}/repeats-used-by-simulator.fasta"
rm -f ${GENOME_FILE} ${GENOME_IMAGE} ${REPEATS_FILE}
java ${JAVA_RUNTIME_FLAGS} -classpath "${REVANT_BINARIES}" de.mpi_cbg.revant.biology.GenomeSimulator -1 ${NONREPETITIVE_BPS} ${REPETITIVE_NONREPETITIVE_RATIO} ${REPEAT_MODEL_DIR} ${GENOME_FILE} null 0 ${REPEATS_FILE} null ${GENOME_IMAGE} 10 null 0 0 0 0
GENOME_LENGTH=$(wc -c < ${GENOME_FILE})
N_READS=$(( (${COVERAGE} * ${GENOME_LENGTH}) / ${AVG_READ_LENGTH} ))

# Simulating reads
if [ ${READ_SIMULATOR} = "pbsim" ]; then
	# 2020 - PBSIM2
	${PBSIM_DIR}/pbsim --depth ${COVERAGE} --length-min ${MIN_READ_LENGTH} --length-max ${MAX_READ_LENGTH} --hmm_model ${PBSIM_DIR}/../data/P6C4.model --length-mean ${AVG_READ_LENGTH} ${GENOME_FILE}
	mv ./sd_* ${INPUT_DIR}
	cat ${INPUT_DIR}/sd_0001.fastq | awk '{if(NR%4==1) {printf(">%s\n",substr($0,2));} else if(NR%4==2) print;}' > ${READS_FILE}
	rm -f ${INPUT_DIR}/sd_0001.fastq ${INPUT_DIR}/sd_0001.maf
elif [ ${READ_SIMULATOR} = "pass" ]; then
	# 2019 - PaSS
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
	:
elif [ ${READ_SIMULATOR} = "npbss" ]; then
	# 2018 - NPBSS
	# cd ${NPBSS_DIR}
	# octave --eval "pkg load statistics; NPBSS('${GENOME_DIR}/${GENOME_FILE}','-dep ${COVERAGE}','-min ${MIN_READ_LENGTH}','-max ${MAX_READ_LENGTH}','-len ${AVG_READ_LENGTH}')"
	# cd -
	# NPBSS_OUTPUT_DIR="${OUTPUT_DIR}/NPBSS/"
	# rm -rf ${NPBSS_OUTPUT_DIR}
	# mkdir ${NPBSS_OUTPUT_DIR}
	# mv ${GENOME_DIR}/*_CLR.* ${NPBSS_OUTPUT_DIR}
	:
elif [ ${READ_SIMULATOR} = "simlord" ]; then
	# 2016 - SimLoRD
	# SIMLORD_OUTPUT_DIR="${OUTPUT_DIR}/SimLoRD/"
	# rm -rf ${SIMLORD_OUTPUT_DIR}
	# mkdir ${SIMLORD_OUTPUT_DIR}
	# simlord --read-reference ${GENOME_DIR}/${GENOME_FILE} --coverage ${COVERAGE} --max-passes 1 --min-readlength ${MIN_READ_LENGTH} ${SIMLORD_OUTPUT_DIR}
	# mv ${SIMLORD_OUTPUT_DIR}/../SimLoRD* ${SIMLORD_OUTPUT_DIR}
	:
else
	echo "Simulator not supported"
	exit 1
fi

# Computing alignments
READS_DB="${INPUT_DIR}/reads.db"
${DAZZDB_DIR}/fasta2DB ${READS_DB} ${READS_FILE}
daligner -T${N_THREADS} -M${MAX_MEMORY} -e${READ_READ_IDENTITY} -l${MIN_ALIGNMENT_LENGTH} ${READS_DB} ${READS_DB}
for FILE in $(find ${INPUT_DIR} -name "*.las"); do
	LAsort ${FILE}
	rm -f ${FILE}
done
LAmerge ${INPUT_DIR}/reads.reads.las ${INPUT_DIR}/*.S.las







# ------------------------------ SIMULATORS CEMETERY -------------------------------------
#FASTQSIM_DIR="/Users/ramseysnow/git/FASTQsim_v2.0"
#LONGISLND_DIR="/Users/ramseysnow/git/longislnd"
#SILICO_DIR="/Users/ramseysnow/git/SiLiCO"

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
