#!/bin/bash
#
# Given the directory of a project, the script processes every connected component of the
# interval graph and creates repeat descriptor files for it (these are stored in files
# $basin-*-*-*.txt$).
#
# Remark: the script assumes that environment variable $REVANT_BINARIES$ is already set
# to the directory that contains REVANT's binaries, i.e. to the directory that contains
# the $de/mpi_cbg/revant/...$ subpath.
#
# --------------------------------- CONFIGURATION ----------------------------------------
# This is the only section of the script that needs to be customized.
#
PROJECT_DIR=$1
MIN_ALIGNMENT_LENGTH="1500"
MEDIAN_READ_LENGTH="28005"
GENOME_COVERAGE="1"
N_THREADS="1"  # Connected components of the interval graph can be processed in parallel
# REVANT
JAVA_RUNTIME_FLAGS="-Xms2G -Xmx10G"
# ----------------------------------------------------------------------------------------




GRAPHS_DIR="${PROJECT_DIR}/step1"
TAGS_DIR="${GRAPHS_DIR}/finalOutput"
STEP4_DIR="${TAGS_DIR}/step4"
STEP5_DIR="${STEP4_DIR}/step5"
QUALITY_THRESHOLDS_FILE="${PROJECT_DIR}/qualityThresholds.txt"
KERNEL_LENGTHS_LOG="${GRAPHS_DIR}/kernelLengths.txt"
QUALITIES_SUFFIX="reads-phred.damar"
if [ ! -e "${PROJECT_DIR}/reads-phred.damar" ]; then
	QUALITIES_SUFFIX="reads-phred.dbdump"
fi
QUALITIES_FILE="${PROJECT_DIR}/${QUALITIES_SUFFIX}"
THREAD_PREFIX="thread-"




# ------------------------------------ A THREAD ------------------------------------------
function run() {
	local FILES_LIST=$1
	for FILE1 in $(cat ${FILES_LIST}); do
		local ID1=$(basename ${FILE1})
		local ID1=${ID1%.*}
		local N_READS=$(wc -l < ${GRAPHS_DIR}/${ID1}-reads-ids.txt | tr -d ' ')
		mkdir ${GRAPHS_DIR}/${ID1}-clusters/
		java ${JAVA_RUNTIME_FLAGS} -classpath "${REVANT_BINARIES}" de.mpi_cbg.revant.intervalgraph.IntervalGraphStep2 ${GRAPHS_DIR}/${ID1}.graph ${GRAPHS_DIR}/${ID1}-clusters 1 3 ${MIN_ALIGNMENT_LENGTH} 1 ${TAGS_DIR} ${ID1} ${N_READS}
		local EXIT_STATUS=$?
		if [ ${EXIT_STATUS} -ne 0 ]; then
		    echo "The following command returned error ${EXIT_STATUS}:"
			echo "java ${JAVA_RUNTIME_FLAGS} -classpath "${REVANT_BINARIES}" de.mpi_cbg.revant.intervalgraph.IntervalGraphStep2 ${GRAPHS_DIR}/${ID1}.graph ${GRAPHS_DIR}/${ID1}-clusters 1 3 ${MIN_ALIGNMENT_LENGTH} 1 ${TAGS_DIR} ${ID1} ${N_READS}"
			exit ${EXIT_STATUS}
		fi
		if [ -z "$(ls ${GRAPHS_DIR}/${ID1}-clusters/*.graph 2> /dev/null)" ]; then
			cp ${GRAPHS_DIR}/${ID1}.graph ${GRAPHS_DIR}/${ID1}-clusters/${ID1}.graph
			cp "${GRAPHS_DIR}/${ID1}-reads-ids.txt" "${GRAPHS_DIR}/${ID1}-clusters/${ID1}-reads-ids.txt"
			cp "${GRAPHS_DIR}/${ID1}-reads-lengths.txt" "${GRAPHS_DIR}/${ID1}-clusters/${ID1}-reads-lengths.txt"
			cp "${GRAPHS_DIR}/${ID1}-${QUALITIES_SUFFIX}" "${GRAPHS_DIR}/${ID1}-clusters/${ID1}-${QUALITIES_SUFFIX}"
			cp "${GRAPHS_DIR}/${ID1}-reads-shortPeriod.txt" "${GRAPHS_DIR}/${ID1}-clusters/${ID1}-reads-shortPeriod.txt"
			cp "${GRAPHS_DIR}/${ID1}-alignments.txt" "${GRAPHS_DIR}/${ID1}-clusters/${ID1}-alignments.txt"
			cp "${GRAPHS_DIR}/${ID1}-alignments-shortPeriod.txt" "${GRAPHS_DIR}/${ID1}-clusters/${ID1}-alignments-shortPeriod.txt"
		fi
		for FILE2 in $(find ${GRAPHS_DIR}/${ID1}-clusters -maxdepth 1 -type f -name "*.graph"); do
			local ID2=$(basename ${FILE2})
			local ID2=${ID2%.*}
			local N_READS=$(wc -l < ${GRAPHS_DIR}/${ID1}-clusters/${ID2}-reads-ids.txt | tr -d ' ')
			local N_ALIGNMENTS=$(wc -l < ${GRAPHS_DIR}/${ID1}-clusters/${ID2}-alignments.txt | tr -d ' ')
			java ${JAVA_RUNTIME_FLAGS} -classpath "${REVANT_BINARIES}" de.mpi_cbg.revant.intervalgraph.IntervalGraphStep3 ${GRAPHS_DIR}/${ID1}-clusters ${ID2} 1 3 ${MIN_ALIGNMENT_LENGTH} 50000 ${N_READS} ${PROJECT_DIR}/qualityThresholds.txt ${N_ALIGNMENTS} ${TAGS_DIR} ${ID1} ${KERNEL_LENGTHS_LOG}
			local EXIT_STATUS=$?
			if [ ${EXIT_STATUS} -ne 0 ]; then
			    echo "The following command returned error ${EXIT_STATUS}:"
				echo "java ${JAVA_RUNTIME_FLAGS} -classpath "${REVANT_BINARIES}" de.mpi_cbg.revant.intervalgraph.IntervalGraphStep3 ${GRAPHS_DIR}/${ID1}-clusters ${ID2} 1 3 ${MIN_ALIGNMENT_LENGTH} 50000 ${N_READS} ${PROJECT_DIR}/qualityThresholds.txt ${N_ALIGNMENTS} ${TAGS_DIR} ${ID1} ${KERNEL_LENGTHS_LOG}"
				exit ${EXIT_STATUS}
			fi
		done
		local N_ALIGNMENTS=$(wc -l < ${GRAPHS_DIR}/${ID1}-alignments.txt | tr -d ' ')
		local N_READS=$(wc -l < ${GRAPHS_DIR}/${ID1}-reads-ids.txt | tr -d ' ')
		java ${JAVA_RUNTIME_FLAGS} -classpath "${REVANT_BINARIES}" de.mpi_cbg.revant.intervalgraph.IntervalGraphStep4 ${GRAPHS_DIR} ${ID1} ${N_ALIGNMENTS} ${N_READS} ${MIN_ALIGNMENT_LENGTH} ${PROJECT_DIR}/qualityThresholds.txt 1 3
		local EXIT_STATUS=$?
		if [ ${EXIT_STATUS} -eq 255 ]; then
			if [ -n "$(ls ${TAGS_DIR}/basin-${ID1}-* 2> /dev/null)" ]; then
				cp ${TAGS_DIR}/basin-${ID1}-* ${STEP4_DIR}/
			fi
			if [ -n "$(ls ${TAGS_DIR}/tags-${ID1}-* 2> /dev/null)" ]; then
				cp ${TAGS_DIR}/tags-${ID1}-* ${STEP4_DIR}/
			fi
			if [ -e ${TAGS_DIR}/tags-${ID1}.txt ]; then
				cp ${TAGS_DIR}/tags-${ID1}.txt ${STEP4_DIR}/
			fi
		elif [ ${EXIT_STATUS} -eq 0 ]; then
			# Files ${TAGS_DIR}/tags-${ID1}.txt and ${TAGS_DIR}/tags-${ID1}-X.txt for some
			# integer X should not be copied, since their intervals already belong to the
			# output of $IntervalGraphStep4$.
			:
		else
		    echo "The following command returned error ${EXIT_STATUS}:"
			echo "java ${JAVA_RUNTIME_FLAGS} -classpath "${REVANT_BINARIES}" de.mpi_cbg.revant.intervalgraph.IntervalGraphStep4 ${GRAPHS_DIR} ${ID1} ${N_ALIGNMENTS} ${N_READS} ${MIN_ALIGNMENT_LENGTH} ${PROJECT_DIR}/qualityThresholds.txt 1 3"
			exit ${EXIT_STATUS}
		fi
	done
}




# ---------------------------------- MAIN PROGRAM ----------------------------------------

# Cleaning up
find ${GRAPHS_DIR} -maxdepth 1 -type d -name "*-clusters" -exec rm -r {} +
find ${TAGS_DIR} -maxdepth 1 -type f -name "basin-*.txt" -exec rm -r {} +
find ${TAGS_DIR} -maxdepth 1 -type f -name "tags-[0-9]{1,}*.txt" -exec rm -r {} +
rm -f ${TAGS_DIR}/allTags.txt  # Does not delete $tags-root.txt$ if it exists.
rm -rf ${STEP4_DIR}
mkdir ${STEP4_DIR}
rm -f ${KERNEL_LENGTHS_LOG}

# Running parallel jobs
TMP_FILE_1="${GRAPHS_DIR}/tmp-1.txt"
find ${GRAPHS_DIR} -maxdepth 1 -type f -name "*.graph" > ${TMP_FILE_1}
N_FILES=$(wc -l < ${TMP_FILE_1} | tr -d ' ')
FILES_PER_THREAD=$(( ${N_FILES} / ${N_THREADS} ))
split -l ${FILES_PER_THREAD} ${TMP_FILE_1} ${THREAD_PREFIX}
for THREAD in $(ls ${THREAD_PREFIX}* 2> /dev/null); do
	run ${THREAD} &
done
wait
rm -f ${TMP_FILE_1} ${THREAD_PREFIX}*
cp ${TAGS_DIR}/tags-root.txt ${STEP4_DIR}/

# Merging all tag files
TMP_FILE="${TAGS_DIR}/tmp.txt"
rm -f ${TMP_FILE}
find "${TAGS_DIR}" -type f -maxdepth 1 -name "tags-*.txt" -exec cat {} + >> ${TMP_FILE}
sort -m -t , -k 1,1n -k 2,2n -k 3,3n ${TMP_FILE} > ${TAGS_DIR}/allTags.txt
rm -f ${TMP_FILE}
find "${STEP4_DIR}" -type f -maxdepth 1 -name "tags-*.txt" -exec cat {} + >> ${TMP_FILE}
sort -m -t , -k 1,1n -k 2,2n -k 3,3n ${TMP_FILE} > ${STEP4_DIR}/allTags.txt
rm -f ${TMP_FILE}

# Step 5
rm -rf ${STEP5_DIR}
mkdir ${STEP5_DIR}
N_READS=$(wc -l < ${PROJECT_DIR}/reads-lengths.txt | tr -d ' ')
java ${JAVA_RUNTIME_FLAGS} -classpath "${REVANT_BINARIES}" de.mpi_cbg.revant.intervalgraph.IntervalGraphStep5 ${GRAPHS_DIR} ${GENOME_COVERAGE} ${QUALITY_THRESHOLDS_FILE} ${N_READS} ${MIN_ALIGNMENT_LENGTH} ${QUALITIES_FILE} ${MEDIAN_READ_LENGTH}
java ${JAVA_RUNTIME_FLAGS} -classpath "${REVANT_BINARIES}" de.mpi_cbg.revant.graphics.PrintBasinDescriptorTree ${PROJECT_DIR}/step1 ${MIN_ALIGNMENT_LENGTH} 0 > "${STEP4_DIR}/basins.dot"