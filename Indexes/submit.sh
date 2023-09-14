#!/bin/bash
#SBATCH --job-name="submit"
#SBATCH --mail-type=FAIL
#SBATCH --output=logs/snakemake.%j.o
#SBATCH --cpus-per-task=1
#SBATCH --mem=1g
#SBATCH --time=04:00:00


#########################################
# Author: Matthew J Brooks
# Last update: Sept 14, 2023 by Matthew J Brooks
# Run with: sbatch submit.sh
#########################################

module load python/3.9

mkdir -p logs

snakemake \
	--snakefile Indexes.py \
	--jobname "{rulename}.{jobid}" \
	--verbose -j \
	-k -p
	-j 3000
	--stats Indexes_snakefile_{jobid}.stats \
	--latency-wait 180 \
	--timestamp \
	--rerun-incomplete \
	--cores 300 \
	--cluster "sbatch --mail-type=FAIL -o logs/{params.rulename}.%j.o {params.batch}" --out logs/job_%j.out" \
	>& logs/Indexes_$SLURM_JOB_ID_snakefile.log


## DRY Run with Print out the shell commands that will be executed
# snakemake --directory . --snakefile Indexes_v2.py --dryrun -p -r