<img align="center" src="./logo.png" width="297" height="35"/>
De novo repeat inference from long reads.


**Remark: This is just a first draft. More detailed instructions will follow.**


## Compiling the project

This project is written in Java, so you need a Java development kit (any version) for compiling it. Java was chosen because it takes segfaults and memory leaks out of the way, which was very useful while developing the first prototype of a complex project. The code uses few Java-specific features and can be converted to C with little effort.

After cloning this repo, you can compile it by running the `build.sh` script in the main directory.

## External dependencies

* [DAZZ_DB](https://github.com/thegenemyers/DAZZ_DB) 
* [DALIGNER](https://github.com/thegenemyers/DALIGNER)
* [d'accord](https://gitlab.com/german.tischler/daccord)

Please refer to each external dependency for instructions on compiling it.



# Running the pipeline

**Remark: A detailed, step-by-step example with a toy dataset will follow.**

## Setting up a new project

All input and intermediate files used by REVANT are currently human-readable text files: this helps with debugging, and is easy to drop in the future in the unlikely case IO becomes a bottleneck. 

To set up a new project, create a new emtpy directory and put the following text files in a subdirectory called `input`.

* `LAshow.txt`: all pairwise alignments between the reads in a random subset that covers 1X of the genome. This is the output of the `LAshow` tool from DALIGNER.
* `DBdump-lengths.txt`: length of each read. This is the output of the `DBdump -rh` tool from DAZZ_DB.
* `qualityThresholds.txt`: text file containing quality thresholds for deciding low-quality regions. Example in `/scripts/1-factorize/qualityThresholds.txt`.
* `qtrack.txt`: an estimate of intrinsic quality value for each read. This can be computed by aligning the 1X subset against, say, a distinct 10X random subset, and by running the `DASqv` tool from DALIGNER.
* Optionally, one or several text files that each represent a DAZZ_DB track: the union of all such tracks will be discarded by the pipeline. This is the output of the `DBdump -r -mTrackName` tool from DAZZ_DB. This is useful for inferring repeats incrementally, e.g. a track could store the parts of the reads that align to an existing repeat library.

## Running the scripts

The `scripts` directory contains shell scripts that automate most steps of the pipeline. This is just the first draft of a very primitive scripting system, and must be improved. Current issues:

* For simplicity, the current scripts are designed for a single machine rather than for a cluster. A cluster version of the whole pipeline will follow.
* The scripts require several input variables from the user, and such variables are currently set inside the scripts themselves. So, when running the pipeline on a new dataset, you currently have to duplicate the `scripts` directory and to customize its scripts as needed. A better solution will follow.
* Scripts are numbered according to the order in which they should be executed. There is currently no scheduling system for running all scripts automatically.

The pipeline is organized in the following steps, that are reflected both by the `scripts` directory and by the source code. For more details about each step, please refer to the comments inside each script.

### Read factorization

Given the alignments between a read and all other reads in the 1X subset, this step marks some intervals of the read as belonging to a specific repeat type. This is embarrassingly parallel and should be performed on several chunks of the alignments file at the same time.



# Acknowledgements

The main parts of this project were developed while the author was affiliated with [MPI-CBG](https://www.mpi-cbg.de/home/) and [CSBD](http://www.csbdresden.de).
