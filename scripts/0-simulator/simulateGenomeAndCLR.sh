#!/bin/bash
#
# Simulates a repetitive haploid genome from a Repbase file, samples reads from it, and
# computes all read-read alignments and all read-repeat alignments.
#
# Remark: the script assumes that environment variable $REVANT_BINARIES$ is already set
# to the directory that contains REVANT's binaries, i.e. to the directory that contains
# the $de/mpi_cbg/revant/...$ subpath.
#
# --------------------------------- CONFIGURATION ----------------------------------------
# This is the only section of the script that needs to be customized.
#
OUTPUT_DIR=$1
REPEAT_MODEL_DIR=$2
READ_SIMULATOR=$3  # One of: pbsim, pass, npbss, simlord
# Genome
NONREPETITIVE_BPS="500"
REPETITIVE_NONREPETITIVE_RATIO="2000"  # In the genome
# Reads
COVERAGE="30"
MIN_READ_LENGTH="18000" #"10000"
AVG_READ_LENGTH="30000" #"15000"
MAX_READ_LENGTH="35000" #"30000"
# Alignments
MIN_ALIGNMENT_LENGTH_READ_READ="500"
READ_READ_IDENTITY=".75"
MIN_ALIGNMENT_LENGTH_READ_REPEAT="500"
READ_REPEAT_IDENTITY=".85"
MAX_MEMORY="12"  # In GB
N_THREADS="4"
# Programs
DALIGNER_DIR="/Users/ramseysnow/git/DALIGNER"
DAZZDB_DIR="/Users/ramseysnow/git/DAZZ_DB"
PBSIM_DIR="/Users/ramseysnow/git/pbsim2/src"  #"/Users/ramseysnow/git/pbsim-1.0.3/src"
PASS_DIR="/Users/ramseysnow/git/PaSS"
NPBSS_DIR="/Users/ramseysnow/git/NPBSS_Octave-master"
# SimLoRD is assumed to have been installed via Python
# REVANT
JAVA_RUNTIME_FLAGS="-Xms2G -Xmx10G"
# ----------------------------------------------------------------------------------------

set -o pipefail; set -e; set -u

# Building the genome
GENOME_FILE="${OUTPUT_DIR}/genome.fasta"
GENOME_IMAGE="${OUTPUT_DIR}/genome.png"
OCCURRENCE_LENGTH_HISTOGRAM="${OUTPUT_DIR}/genome-occurrences-lengths.txt"
READS_FILE="${OUTPUT_DIR}/reads.fasta"
REPEATS_FILE="${OUTPUT_DIR}/repeats-used-by-simulator.fasta"
REPEATS_FILE_LENGTHS="${OUTPUT_DIR}/repeats-lengths.txt"
REPEATS_FILE_ISPERIODIC="${OUTPUT_DIR}/repeats-isPeriodic.txt"
rm -f ${GENOME_FILE} ${GENOME_IMAGE} ${REPEATS_FILE} ${REPEATS_FILE_LENGTHS} ${REPEATS_FILE_ISPERIODIC}
java ${JAVA_RUNTIME_FLAGS} -classpath "${REVANT_BINARIES}" de.mpi_cbg.revant.biology.GenomeSimulator -1 ${NONREPETITIVE_BPS} ${REPETITIVE_NONREPETITIVE_RATIO} ${REPEAT_MODEL_DIR} ${GENOME_FILE} null 0 ${REPEATS_FILE} ${REPEATS_FILE_LENGTHS} ${REPEATS_FILE_ISPERIODIC} null ${GENOME_IMAGE} ${AVG_READ_LENGTH} null 0 0 0 0 ${OCCURRENCE_LENGTH_HISTOGRAM}
GENOME_LENGTH=$(wc -c < ${GENOME_FILE})
N_READS=$(( (${COVERAGE} * ${GENOME_LENGTH}) / ${AVG_READ_LENGTH} ))

# Simulating reads
TMP_READS="${OUTPUT_DIR}/reads.txt"
rm -f ${TMP_READS}
if [ ${READ_SIMULATOR} = "pbsim" ]; then
	# 2020 - PBSIM2
	${PBSIM_DIR}/pbsim --depth ${COVERAGE} --length-min ${MIN_READ_LENGTH} --length-max ${MAX_READ_LENGTH} --hmm_model ${PBSIM_DIR}/../data/P6C4.model --length-mean ${AVG_READ_LENGTH} ${GENOME_FILE}
	mv ./sd_* ${OUTPUT_DIR}
	cat ${OUTPUT_DIR}/sd_0001.fastq | awk '{if(NR%4==1) {printf(">%s\n",substr($0,2));} else if(NR%4==2) print;}' > ${TMP_READS}
	rm -f ${OUTPUT_DIR}/sd_0001.fastq ${OUTPUT_DIR}/sd_0001.maf ${OUTPUT_DIR}/sd_0001.ref
elif [ ${READ_SIMULATOR} = "pass" ]; then
	# 2019 - PaSS
	PASS_OUTPUT_DIR="${OUTPUT_DIR}/PaSS"
	rm -rf ${PASS_OUTPUT_DIR}
	mkdir ${PASS_OUTPUT_DIR}
	# Building index
	perl ${PASS_DIR}/pacbio_mkindex.pl ${GENOME_FILE} ${PASS_OUTPUT_DIR}
	PASS_MACHINE="pacbio_RS"  # CLR
	#PASS_MACHINE="pacbio_sequel"  # CCS
	PASS_ERROR_MODEL="${PASS_DIR}/E.coli/ecoli.config"  # RS
	#PASS_ERROR_MODEL="${PASS_DIR}/C.elegan/elegan.config"  # RS
	#PASS_ERROR_MODEL="${PASS_DIR}/Arabidopsis/arab.config"  # Sequel
	#PASS_ERROR_MODEL="${PASS_DIR}/sim.config"  # Example dataset
	${PASS_DIR}/PaSS -list ${PASS_OUTPUT_DIR}/percentage.txt -index ${PASS_OUTPUT_DIR}/index -m ${PASS_MACHINE} -c ${PASS_ERROR_MODEL} -r ${N_READS} -t ${N_THREADS} -o ${PASS_OUTPUT_DIR}/reads-PaSS -d
	cat ${PASS_OUTPUT_DIR}/reads-PaSS.fq | awk '{if(NR%4==1) {printf(">%s\n",substr($0,2));} else if(NR%4==2) print;}' > ${TMP_READS}
	rm -rf ${PASS_OUTPUT_DIR}/../reftemp/
	rm -rf ${PASS_OUTPUT_DIR}
elif [ ${READ_SIMULATOR} = "npbss" ]; then
	# 2018 - NPBSS
	# cd ${NPBSS_DIR}
	# octave --eval "pkg load statistics; NPBSS('${OUTPUT_DIR}/${GENOME_FILE}','-dep ${COVERAGE}','-min ${MIN_READ_LENGTH}','-max ${MAX_READ_LENGTH}','-len ${AVG_READ_LENGTH}')"
	# cd -
	# NPBSS_OUTPUT_DIR="${OUTPUT_DIR}/NPBSS/"
	# rm -rf ${NPBSS_OUTPUT_DIR}
	# mkdir ${NPBSS_OUTPUT_DIR}
	# mv ${OUTPUT_DIR}/*_CLR.* ${NPBSS_OUTPUT_DIR}
	:
elif [ ${READ_SIMULATOR} = "simlord" ]; then
	# 2016 - SimLoRD
	# SIMLORD_OUTPUT_DIR="${OUTPUT_DIR}/SimLoRD/"
	# rm -rf ${SIMLORD_OUTPUT_DIR}
	# mkdir ${SIMLORD_OUTPUT_DIR}
	# simlord --read-reference ${OUTPUT_DIR}/${GENOME_FILE} --coverage ${COVERAGE} --max-passes 1 --min-readlength ${MIN_READ_LENGTH} ${SIMLORD_OUTPUT_DIR}
	# mv ${SIMLORD_OUTPUT_DIR}/../SimLoRD* ${SIMLORD_OUTPUT_DIR}
	:
else
	echo "Simulator not supported"
	exit 1
fi
rm -f ${READS_FILE}
java ${JAVA_RUNTIME_FLAGS} -classpath "${REVANT_BINARIES}" de.mpi_cbg.revant.biology.Repbase2Pacbio ${TMP_READS} ${READS_FILE}
rm -f ${TMP_READS}
seq 0 $(( ($(wc -l < ${READS_FILE}) / 2) - 1 )) > ${OUTPUT_DIR}/reads-ids.txt
java ${JAVA_RUNTIME_FLAGS} -classpath "${REVANT_BINARIES}" de.mpi_cbg.revant.util.Fasta2Lengths ${READS_FILE} > ${OUTPUT_DIR}/reads-lengths.txt

# Computing read-read alignments
READS_DB="${OUTPUT_DIR}/reads.db"
DBrm ${READS_DB} && echo 0 || echo 1
${DAZZDB_DIR}/fasta2DB ${READS_DB} ${READS_FILE}
DBsplit -x16 ${READS_DB}
daligner -T${N_THREADS} -M${MAX_MEMORY} -e${READ_READ_IDENTITY} -l${MIN_ALIGNMENT_LENGTH_READ_READ} ${READS_DB} ${READS_DB}
mv ./*.las ${OUTPUT_DIR}
LAsort ${OUTPUT_DIR}/reads.reads.las
rm -f ${OUTPUT_DIR}/reads.reads.las
LAshow ${READS_DB} ${OUTPUT_DIR}/reads.reads.S.las > ${OUTPUT_DIR}/LAshow-reads-reads.txt
rm -f ${OUTPUT_DIR}/reads.reads.S.las

# Computing read-repeat alignments
REPEATS_DB="${OUTPUT_DIR}/repeats.db"
DBrm ${REPEATS_DB} && echo 0 || echo 1
${DAZZDB_DIR}/fasta2DB ${REPEATS_DB} ${REPEATS_FILE}
daligner -T${N_THREADS} -M0 -t100000000 -e${READ_REPEAT_IDENTITY} -l${MIN_ALIGNMENT_LENGTH_READ_REPEAT} ${REPEATS_DB} ${READS_DB}
mv ./*.las ${OUTPUT_DIR}
REPEATS_FILE_PREFIX=$(basename ${REPEATS_DB} .db)
LAsort ${OUTPUT_DIR}/reads.${REPEATS_FILE_PREFIX}.las
rm -f ${OUTPUT_DIR}/reads.${REPEATS_FILE_PREFIX}.las
LAshow ${READS_DB} ${REPEATS_DB} ${OUTPUT_DIR}/reads.${REPEATS_FILE_PREFIX}.S.las > ${OUTPUT_DIR}/LAshow-reads-repeats.txt
rm -f ${OUTPUT_DIR}/*.las









# ------------------------------ SIMULATORS CEMETERY -------------------------------------
#FASTQSIM_DIR="/Users/ramseysnow/git/FASTQsim_v2.0"
#LONGISLND_DIR="/Users/ramseysnow/git/longislnd"
#SILICO_DIR="/Users/ramseysnow/git/SiLiCO"

# ----------------------------------------------------------------------------------------
# 2016 - LongISLND (Does not ship with a model, cumbersome to create one)
# ${LONGISLND_DIR}/simulate.py --fasta ${OUTPUT_DIR}/${GENOME_FILE} --model_dir ........

# ----------------------------------------------------------------------------------------
# 2016 - SiLiCO (Says it's successful but produces an empty file)
#
# SILICO_OUTPUT_DIR="${OUTPUT_DIR}/SiLiCO/"
# rm -rf ${SILICO_OUTPUT_DIR}
# mkdir ${SILICO_OUTPUT_DIR}
# python3 ${SILICO_DIR}/SiLiCO.py -i ${OUTPUT_DIR}/${GENOME_FILE} -o ${SILICO_OUTPUT_DIR} -m ${AVG_READ_LENGTH} -c ${COVERAGE} --pacbio

# ----------------------------------------------------------------------------------------
# 2014 - FASTQSim (Extremely slow, does not work with a reference of length 10k)
#
# FASTQSIM_OUTPUT_DIR="${OUTPUT_DIR}/FASTQSim/"
# rm -rf ${FASTQSIM_OUTPUT_DIR}
# mkdir ${FASTQSIM_OUTPUT_DIR}
# cd ${FASTQSIM_DIR}
# python ./src/pathoSpike.py -nobackground -platform pacbio -source ${COVERAGE} ${OUTPUT_DIR}/${GENOME_FILE} false -o ${FASTQSIM_OUTPUT_DIR} -threads ${N_THREADS}
# cd -

# ----------------------------------------------------------------------------------------
# 2012 - PBSIM (Superseded by PBSIM2)
#
# PBSIM_DATA_TYPE="CLR"
# #PBSIM_DATA_TYPE="CCS"
# PBSIM_ERROR_MODEL="${PBSIM_DIR}/../data/model_qc_clr"
# ${PBSIM_DIR}/pbsim --data-type ${PBSIM_DATA_TYPE} --depth ${COVERAGE} --length-min ${MIN_READ_LENGTH} --length-max ${MAX_READ_LENGTH} --model_qc ${PBSIM_ERROR_MODEL} --length-mean ${AVG_READ_LENGTH} ${OUTPUT_DIR}/${GENOME_FILE}
# PBSIM_OUTPUT_DIR="${OUTPUT_DIR}/PBSIM"
# rm -rf ${PBSIM_OUTPUT_DIR}
# mkdir ${PBSIM_OUTPUT_DIR}
# mv ./sd_* ${PBSIM_OUTPUT_DIR}
