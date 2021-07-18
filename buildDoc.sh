#!/bin/sh
#
# Builds the documentation of the project.
#
SOURCE_DIR="./src"
DOCS_DIR="./docs"


if [ ! -e ${DOCS_DIR} ]; then
	mkdir ${DOCS_DIR}
fi
javadoc -quiet -sourcepath ${SOURCE_DIR} -d ${DOCS_DIR} de.mpi_cbg.revant.util de.mpi_cbg.revant.factorize de.mpi_cbg.revant.intervalgraph \
												 	    de.mpi_cbg.revant.fragments de.mpi_cbg.revant.consensus de.mpi_cbg.revant.graphics \
												 		de.mpi_cbg.revant.biology