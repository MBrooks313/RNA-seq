#######################################
# This Indexes.py Snakefile creates indexes necessary for RNA-seq analysis
# Written by Matthew Brooks on Mar 30, 2017
# Modified on Sept 14, 2023
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
basedir = config["basedir"]
name = config["name"]

#These are not to be changed
DIRS = ['STAR/', 'STAR/124base/', 'kallisto/', 'logs/']
STAR = basedir + 'STAR/124base/' + 'SAindex'
KAL = basedir + 'kallisto/' + name



##############################
# Snakemake rules for analysis
##############################

localrules: dirs


rule all:
    input: STAR, KAL, DIRS
    params: batch = config["job_all"]
    

rule dirs:
    output: DIRS
    log:    "logs/dirs.log"   
    params: 
            rulename = "dirs",
            batch = config["job_all"]
    shell:  "cd {basedir}; mkdir -p "+' '.join(DIRS)

rule star:
    input: dir = 'STAR/124base/'
    output: STAR
    version: config["star"]
    log:    "logs/star.log"   
    params: 
            rulename = "star",
            batch = config["job_star"],
            ref = config["ref"],
            gtf = config["gtf"]
    shell: """ \
    cd {basedir}; \
    module load STAR/{version} || exit 1;
    STAR \
    --runThreadN ${{SLURM_CPUS_ON_NODE}} \
    --runMode genomeGenerate \
    --genomeDir {input.dir} \
    --genomeFastaFiles {params.ref} \
    --sjdbGTFfile {params.gtf} \
    --sjdbOverhang 124 \
    --genomeSAindexNbases 11 \
    """


rule kallisto:
    input: dir = 'kallisto/'
    output: KAL
    version: config["kallisto"]
    log:    "logs/kallisto.log"   
    params: 
            rulename = "kallisto",
            batch = config["job_kallisto"],
            name = config["name"],
            trans = config["trans"]
    shell: """ \
    cd {basedir}; \
    module load kallisto; \
    kallisto index -i {input.dir}{params.name} {params.trans} \
    """
