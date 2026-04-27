
# Description
The methods used in this module are described in `docs/methods.pdf`.

Insert brief description of module. This is a Snakemake workflow example.



# Quickstart

Quickstart for deploying this pipeline on a local cluster.

## 1. Set up the environment

Install most of the packages via conda:
```bash
# the repo directory
REPO_MODULE="${HOME}/repo/path/to/this/module"

# install environment using Conda
conda env create --name example_module_2 --file ${REPO_MODULE}/environment.yml

# activate the new Conda environment
source activate example_module_2

# to update environment file:
#conda env export --no-builds | grep -v prefix | grep -v name > environment.yml
```


## 2. Prepare the input files

The Cromwell pipeline takes as input:
1. File1.
2. File2.

Prepare the File1.
```bash
python "${REPO_MODULE}/scripts/prepare_file1_script.py" \
    --input_file input-file1.tsv > input-file1.tsv
```

## 3. Set up and run Snakemake

*
```bash
# make a directory for a demo
mkdir snkmk-example_module_2
cd snkmk-example_module_2

# set up the directory with the required data dir
mkdir data
cp ${REPO_MODULE}/data/demo/* data

# set up the config file inside the working dir
# change json file to fit the analysis
cp ${REPO_MODULE}/example_runtime_setup/params-input.json .

# now run snakemake in dryrun
snakemake --snakefile ${REPO_MODULE}/example_module_2.snk --configfile params-input.json --dryrun

# now run snakemake
snakemake --snakefile ${REPO_MODULE}/example_module_2.snk --configfile params-input.json --printshellcmds
```

The above jobs may take a little while, it would be much faster if we utilized a cluster to run the jobs.

Below is a demo for LSF, you can easily tweak the config file for other systems.

```bash
# now run snakemake
snakemake --snakefile ${REPO_MODULE}/example_module_2.snk --configfile params-input.json --printshellcmds --latency-wait 600 --jobs 999 --cluster-config config_cluster.json --cluster 'bsub -g {cluster.group} -M {cluster.memory} -R {cluster.resources} -J {cluster.name} -oo {cluster.output} -eo {cluster.error}'

# alternatively load settings from the config file
#
# set up the config file inside the working dir
# change json files to fit the system
cp ${REPO_MODULE}/../../lib/snakemake/configs/cluster-lsf.json params-cluster.json

snakemake --snakefile ${REPO_MODULE}/example_module_2.snk --configfile params-input.json --printshellcmds --latency-wait 600 --jobs 999 --cluster-config params-cluster.json --cluster ${REPO_MODULE}/../../lib/snakemake/wrappers/lsf.py
```

See example runtime setup in the [example_runtime_setup](example_runtime_setup) dir.


# Pipeline

## Core pipeline scripts

* `example_module_2.snk`: Performs the core analysis.

# Additional notes

## Analysis 1 of example_module_2

* Additional notes.
