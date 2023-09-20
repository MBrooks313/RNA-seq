#!/bin/bash
#SBATCH --time=10:00:00


#########################################
# Author: Matthew J Brooks
# Last update: Sept 19, 2023 by Matthew J Brooks
# Run with: sbatch submit.sh
#########################################

NOW=$(date +"%Y%m%d_%H%M%S")

module load python/3.9

mkdir -p logs

snakemake \
	--snakefile Indexes.py \
	--jobname '{rulename}.{jobid}' \
	--rerun-incomplete \
	--nolock \
	--verbose \
	-k -p \
	-j 3000 \
	--stats rna_pipeline.stats \
	--cluster "sbatch --mail-type=FAIL -o logs/{params.rulename}.%j.o {params.batch}" \
	>& rna_pipeline_${NOW}.log


## DRY Run with Print out the shell commands that will be executed
# snakemake --directory . --snakefile Indexes.py --dryrun -p -r
