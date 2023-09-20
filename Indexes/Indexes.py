#######################################
# This Indexes.py Snakefile creates indexes necessary for RNA-seq analysis
# Written by Matthew Brooks on Mar 30, 2017
# Modified on Sept 19, 2023
#######################################


#######################################
# Import config file and modules needed
#######################################

import json


#############################################################
# List of directories needed and end point files for analysis
#############################################################

#Import the species specific files from the config.py file
configfile: "config.json"

STARDIR = config["star_outdir"]
KALLDIR = config["kall_outdir"]
NAME = config["name"]


# End point files
STAR = STARDIR + 'Log.out'
KAL = KALLDIR + NAME



##############################
# Snakemake rules for analysis
##############################

localrules: all


rule all:
    input: STAR, KAL
    params: batch = config["job_all"]
    

rule star:
    output: STAR
    version: config["star"]
    log:    "logs/star.log"   
    params: 
            rulename = "star",
            batch = config["job_star"],
            outdir = config["star_outdir"],
            ref = config["ref"],
            gtf = config["gtf"]
    shell: """ \
    module load STAR/{version} || exit 1;
    mkdir -p {params.outdir};
    STAR \
    --runThreadN ${{SLURM_CPUS_ON_NODE}} \
    --runMode genomeGenerate \
    --genomeDir {params.outdir} \
    --genomeFastaFiles {params.ref} \
    --sjdbGTFfile {params.gtf} \
    --sjdbOverhang 124 \
    --genomeSAindexNbases 11 \
    """


rule kallisto:
    output: KAL
    version: config["kallisto"]
    log:    "logs/kallisto.log"   
    params: 
            rulename = "kallisto",
            batch = config["job_kallisto"],
            outdir = config["kall_outdir"],
            name = config["name"],
            trans = config["trans"]
    shell: """ \
    module load kallisto; \
    mkdir -p {params.outdir};
    kallisto index -i {params.outdir}/{params.name} {params.trans} \
    """
