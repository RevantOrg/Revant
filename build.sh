#!/bin/sh
#
# Compiles the project
#
SOURCE_DIR="./src"
SOURCE_DIR_NESTED="${SOURCE_DIR}/de/mpi_cbg/revant"
LIBRARIES="./lib/commons-statistics-distribution-1.0.jar"
BUILD_DIR="./bin"
DEBUGGING_INFO="1"  # 0/1


COMPILER_FLAGS="-deprecation"    # "-deprecation -Xlint"
if [ ${DEBUGGING_INFO} -eq 0 ]; then
	COMPILER_FLAGS="${COMPILER_FLAGS} -g:none"
fi
if [ ! -e ${BUILD_DIR} ]; then
	mkdir ${BUILD_DIR}
fi
#rm -rf ${BUILD_DIR}
#mkdir ${BUILD_DIR}
javac ${COMPILER_FLAGS} -d ${BUILD_DIR} ${SOURCE_DIR_NESTED}/util/*.java
javac ${COMPILER_FLAGS} -cp ${BUILD_DIR} -d ${BUILD_DIR} ${SOURCE_DIR_NESTED}/factorize/*.java
javac ${COMPILER_FLAGS} -cp ${BUILD_DIR} -d ${BUILD_DIR} ${SOURCE_DIR_NESTED}/intervalgraph/*.java
javac ${COMPILER_FLAGS} -cp ${BUILD_DIR} -d ${BUILD_DIR} ${SOURCE_DIR_NESTED}/fragments/*.java
javac ${COMPILER_FLAGS} -cp ${BUILD_DIR} -d ${BUILD_DIR} ${SOURCE_DIR_NESTED}/consensus/*.java
javac ${COMPILER_FLAGS} -cp ${BUILD_DIR} -d ${BUILD_DIR} ${SOURCE_DIR_NESTED}/graphics/*.java
javac ${COMPILER_FLAGS} -cp ${BUILD_DIR} -d ${BUILD_DIR} ${SOURCE_DIR_NESTED}/biology/*.java
javac ${COMPILER_FLAGS} -cp ${BUILD_DIR}:${LIBRARIES} -d ${BUILD_DIR} ${SOURCE_DIR_NESTED}/apps/*.java
