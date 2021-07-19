<img align="center" src="./logo.png" width="297" height="35"/>

A pipeline for de novo repeat inference from long reads.

**Remark: This is just a first draft of the readme file. More detailed instructions will follow.**


## Compiling the project

This project is written in Java, so you need a Java development kit (any version) for compiling it. Java was chosen because it takes segfaults and memory leaks out of the way, which was very useful for developing the first prototype of a complex project. The code uses very few Java-specific features and can be converted to C with little effort.

After cloning this repo, you can compile it by running the `build.sh` script in the main directory.

## External dependencies

* [DAZZ_DB](https://github.com/thegenemyers/DAZZ_DB) 
* [DALIGNER](https://github.com/thegenemyers/DALIGNER)
* [d'accord](https://gitlab.com/german.tischler/daccord)

Please refer to each external dependency for instructions on compiling it.



# Running the pipeline

**Remark: A detailed, step-by-step example with a toy dataset will follow.**

To set up a new project, create a new emtpy directory and put the following input files in a subdirectory called `input`. All input and intermediate files used by REVANT are currently human-readable text files: this is for debugging and might be dropped in the future.

* `LAshow.txt`: all pairwise alignments between the reads in a random subset that covers 1X of the genome. This is the output of the `LAshow` tool from DALIGNER.
* `DBdump-lengths.txt`: 


The `scripts` directory contains shell scripts that automate most steps of the pipeline. This is just a first draft of a scripting system and must be improved. Current issues:

* For simplicity, the current scripts are designed for a single machine rather than for a cluster. A cluster version of the whole pipeline will follow.
* The scripts require several input variables from the user, and such variables are currently set inside the scripts themselves. So, when running the pipeline on a new dataset, you currently have to duplicate the `scripts` directory and to customize its scripts as needed. A better solution will follow.
* Scripts are numbered according to the order in which they should be executed. There is currently no scheduling system for running all scripts automatically.

The pipeline is organized in the following steps. Both the `scripts` directory and the source code are organized around these steps. 

## Read factorization

This is the process that 
