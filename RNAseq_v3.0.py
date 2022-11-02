#######################################
# This is an RNAseq analysis snakemake script.
# Written my Matthew J. Brooks on December 17th, 2020
# This runs on an HPC running SLURM
#######################################



#######################################
# Import config file and modules needed
#######################################

# Import modules
import glob
import os
import json
import pandas as pd

# Snakemake Base location
try:
	BASE_DIR=os.environ['BASE_DIR']
except KeyError:
	print("I can not locate your BASE_DIR directory.")
	pass

# Import configs
configfile: BASE_DIR + "/rna_config.json"
FASTQS = BASE_DIR + "/FQ"

meta_file = config["meta_file"]
fq_col = config["fastq_column"]
meta = pd.read_csv(BASE_DIR + "/" + meta_file)
pe = config["pe"]
fq_R1_end = config["fastq_R1_end"]
fq_R2_end = config["fastq_R2_end"]
base_R1_end = fq_R1_end.split(".")[0]
base_R2_end = fq_R2_end.split(".")[0]

########################################
# Import sample names from the FQ folder
########################################

# Define sample base names
#SAMPLES = [os.path.basename(fname).split('.')[0] for fname in glob.glob(FASTQS + '/*.R1.fastq.gz')]
SAMPLES = meta[fq_col].tolist()

# Define read IDs
if pe:
    READS = [base_R1_end, base_R2_end]
else:
    READS = [base_R1_end]


#############################################################
# List of directories needed and end point files for analysis
#############################################################

# DIRS = ['trim/','BAMS/','star/','kallisto/','fastqc/','logs/']
FQC = expand("fastqc/{sample}{read}_fastqc.html", sample=SAMPLES, read=READS)
KAL = expand("kallisto/{sample}/abundance.h5", sample=SAMPLES)
IDX = expand("BAMS/{sample}.bam.bai", sample=SAMPLES)
QCGB = ["QC/GeneBody/rseqc.Coverage.heatMap.pdf"]
QCF = ["QC/Multi_FQ/Fastqc_multiqc_report.html"]
QCS = ["QC/Multi_STAR/STAR_multiqc_report.html"]
QCR = expand("count_rRNA/{sample}/rRNA.txt", sample=SAMPLES)

##############################
# Snakemake rules for analysis
##############################

localrules: all

rule all:
        input:  QCR + KAL + FQC + IDX + QCGB + QCF + QCS
        params:
                batch = config["job_all"]


rule trimmo:
    """
    This rule trims the Illumina adapters from the fastq files
    """
    input:
            R1 = FASTQS + "/{sample}" + fq_R1_end,
            R2 = FASTQS + "/{sample}" + fq_R2_end,
            adapter = config["trimmoadapt"]
    output:
            forw = ("trim/{sample}" + fq_R1_end),
            rev = ("trim/{sample}" + fq_R2_end),
            for_un = temp("trim/unpaired_{sample}" + fq_R1_end),
            rev_un = temp("trim/unpaired_{sample}" + fq_R2_end)
    log:    "logs/trim.{sample}.log"
    version: config["trimmo"]
    params:
            rulename = "trimmo",
            batch = config["job_trimmo"]
    shell: """
    module load trimmomatic/{version} || exit 1;
    mkdir -p trim;
    java -jar $TRIMMOJAR PE -threads ${{SLURM_CPUS_ON_NODE}} \
    {input.R1} {input.R2} \
    {output.forw} {output.for_un} \
    {output.rev} {output.rev_un} \
    ILLUMINACLIP:{input.adapter}:2:30:10:1:TRUE \
    """


rule kallisto:
   input:
            R1 = "trim/{sample}" + fq_R1_end,
            R2 = "trim/{sample}" + fq_R2_end
   output:
            "kallisto/{sample}/abundance.h5",
            "kallisto/{sample}/abundance.tsv",
            "kallisto/{sample}/run_info.json",
            dir = directory("kallisto/{sample}")
   log:    "logs/kallisto.{sample}.log"
   version: config["kallisto"]
   params:
           rulename = "kallisto",
           batch = config["job_kallisto"],
           ref = config["kallidx"]
   shell: """
   module load kallisto/{version} || exit 1;
   mkdir -p kallisto;
   kallisto quant \
   -i {params.ref} \
   -b 30 \
   -o {output.dir} \
   {input.R1} {input.R2}
   """


rule star:
   input:
            R1 = "trim/{sample}" + fq_R1_end,
            R2 = "trim/{sample}" + fq_R2_end
   output:
            "star/{sample}.starAligned.out.bam",
            dir = directory("star/{sample}.star/")
   log:    "logs/star.{sample}.log"
   version: config["star"]
   params:
           rulename = "star",
           batch = config["job_star"],
           ref = config["staridx"]
   shell: """
   module load STAR/{version} || exit 1;
   mkdir -p {output.dir};
   STAR --runThreadN ${{SLURM_CPUS_ON_NODE}} \
   --genomeDir {params.ref} \
   --readFilesIn {input.R1} {input.R2} \
   --outFileNamePrefix {output.dir} \
   --readFilesCommand zcat \
   --outSAMtype BAM Unsorted \
   --outSAMprimaryFlag AllBestScore \
   --outReadsUnmapped Fastx \
   --twopassMode Basic \
   --genomeSAindexNbases 11 \
   --outFilterType BySJout \
   --outFilterMultimapNmax 20 \
   --alignSJoverhangMin 8 \
   --alignSJDBoverhangMin 1 \
   --outFilterMismatchNmax 999 \
   --alignIntronMin 20 \
   --alignIntronMax 1000000 \
   --alignMatesGapMax 1000000 \
   """


rule sort_index:
    input: "star/{sample}.starAligned.out.bam"
    output:
            bam = "BAMS/{sample}.bam",
            bai = "BAMS/{sample}.bam.bai"
    log:    "logs/sort_index.{sample}.log"
    version: config["samtools"]
    params:
            rulename = "sort_index",
            batch = config["job_sort_index"]
    shell: """
    module load samtools/{version} || exit 1;
    mkdir -p BAMS;
    samtools sort \
    -o {output.bam} \
    -T {output.bam} \
    -@ ${{SLURM_CPUS_ON_NODE}} \
    {input};
    samtools index {output.bam};
    """


rule genebody_cov:
    input:  expand("BAMS/{sample}.bam", sample = SAMPLES)
    output: "QC/GeneBody/rseqc.Coverage.heatMap.pdf",
            dir = directory("QC/GeneBody/")
    log:    "logs/genebody_cov.log"
    version: config["rseqc"]
    params:
            rulename = "genebody_cov",
            batch = config["job_genebody_cov"],
            ref = config["hkgidx"]
    shell:  """
    module load rseqc/{version} || exit 1;
    mkdir -p {output.dir};
    geneBody_coverage.py \
    -r {params.ref} \
    -f "pdf" \
    -i BAMS \
    -o {output.dir}/rseqc
    """


rule fastqc:
    input:
            R1 = "trim/{sample}" + fq_R1_end,
            R2 = "trim/{sample}" + fq_R2_end
    output:
            "fastqc/{sample}" + base_R1_end + "_fastqc.html",
            "fastqc/{sample}" + base_R2_end + "_fastqc.html"
    log:    "logs/fastqc.{sample}.log"
    version: config["fastqc"]
    params:
            rulename = "fastqc",
            batch = config["job_fastqc"]
    shell:  """
    module load fastqc/{version} || exit 1;
    mkdir -p fastqc;
    fastqc -o fastqc {input.R1} {input.R2}
    """


rule multiqc_fastq:
    input:  expand("fastqc/{sample}"+ base_R1_end + "_fastqc.html", sample = SAMPLES),
            expand("fastqc/{sample}" + base_R1_end + "_fastqc.html", sample = SAMPLES)
    output: "QC/Multi_FQ/Fastqc_multiqc_report.html",
            dir = directory("QC/Multi_FQ/")
    log:    "logs/multiqc_fastq.log"
    version: config["multiqc"]
    params:
            rulename = "multiqc_fastq",
            batch = config["job_multiqc"]
    shell:  """
    module load multiqc/{version} || exit 1;
    mkdir -p QC/Multi_FQ/;
    multiqc --title "Fastqc" -o "QC/Multi_FQ/" fastqc/
    """


rule multiqc_star:
    input:  expand("star/{sample}.starAligned.out.bam", sample = SAMPLES)
    output: "QC/Multi_STAR/STAR_multiqc_report.html",
            dir = directory("QC/Multi_STAR/")
    log:    "logs/multiqc_star.log"
    version: config["multiqc"]
    params:
            rulename = "multiqc_star",
            batch = config["job_multiqc"]
    shell:  """
    module load multiqc/{version} || exit 1;
    mkdir -p QC/Multi_STAR/;
    multiqc --title "STAR" -o QC/Multi_STAR/ star/
    """


rule count_rRNA:
    input:  bam = "BAMS/{sample}.bam"
    output: ct = "count_rRNA/{sample}/rRNA.txt",
            dir = directory("count_rRNA/{sample}")
    log:    "logs/count_rRNA.{sample}.log"
    version:    config["samtools"]
    params:
            rulename = "count_rRNA",
            batch = config["job_count_rRNA"],
            ref = config["rRNAidx"]
    shell:  """
    module load samtools/{version} || exit 1;
    mkdir -p {output.dir};
    samtools view -c -L {params.ref} {input} > {output.ct}
    """


# rule :
#         input:
#         output:
#         log:
#         version:
#         params:
#                 rulename =
#                 batch =
#         shell:  """
#     """
