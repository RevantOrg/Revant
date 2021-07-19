<img align="center" src="./logo.png" width="297" height="35"/>

A pipeline for de novo repeat inference from long reads.

**Remark: This is just a first draft of the readme file. More detailed instructions will follow.**


## Compiling the project

This project is written in Java, so you need a Java development kit (any version) for compiling it. After cloning this repo, you can compile it by running the `build.sh` script in the main directory.

## External dependencies

* [DAZZ_DB](https://github.com/thegenemyers/DAZZ_DB) 
* [DALIGNER](https://github.com/thegenemyers/DALIGNER)
* [d'accord](https://gitlab.com/german.tischler/daccord)

Please refer to each external dependency for instructions on compiling it.



# Running the pipeline

The `scripts` directory contains Bash scripts that automate most steps of the pipeline. For simplicity, the current scripts are designed for a single machine rather than for a cluster: a cluster version of the whole pipeline will follow. The scripts require several input variables from the user, and such variables are currently set inside the scripts themselves. So, when running the pipeline on a new project, please duplicate the `scripts` directory and customize its scripts as needed. 
