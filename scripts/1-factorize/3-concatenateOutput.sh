#!/bin/bash
#
# Given the directory of a project, the script concatenates the factorization output from
# multiple chunks.
#
# --------------------------------- CONFIGURATION ----------------------------------------
# This is the only section of the script that needs to be customized.
#
PROJECT_DIR=$1
LAST_CHUNK=$2  # Zero-based
# ----------------------------------------------------------------------------------------




CONCATENATION=${PROJECT_DIR}/intervals-periodic-all.txt
/bin/rm -f ${CONCATENATION}
for FILE in $(seq 0 ${LAST_CHUNK})
do 
   /bin/cat ${PROJECT_DIR}/intervals-periodic-${FILE}.txt >> ${CONCATENATION}
done
CONCATENATION=${PROJECT_DIR}/intervals-dense-all.txt
/bin/rm -f ${CONCATENATION}
for FILE in $(seq 0 ${LAST_CHUNK})
do 
   /bin/cat ${PROJECT_DIR}/intervals-dense-${FILE}.txt >> ${CONCATENATION}
done
CONCATENATION=${PROJECT_DIR}/intervals-alignments-all.txt
/bin/rm -f ${CONCATENATION}
for FILE in $(seq 0 ${LAST_CHUNK})
do
   /bin/cat ${PROJECT_DIR}/intervals-alignments-${FILE}.txt >> ${CONCATENATION}
done
CONCATENATION=${PROJECT_DIR}/connection-all.txt
/bin/rm -f ${CONCATENATION}
for FILE in $(seq 0 ${LAST_CHUNK})
do
   /bin/cat ${PROJECT_DIR}/connection-${FILE}.txt >> ${CONCATENATION}
done 
