# RNA-seq
This is the basic bulk RNA sequencing pipeline.

The pipeline was based off the NNRL lab pipeline and is adapted to run on HPCs running SLURM.

This requires the snakefile RNAseq_v3.0.py, rna_config.json, meta.csv and starts using the rna_submit_snakemake.sh script.

### Usage: 
### sbatch --time=24:00:00 rna_submit_snakemake.sh


## Edit the following files ##
meta.csv : add the base FQ names (minus the [R1|R2]_ending.fastq.gz)

rna_submit_snakemake.sh : change "<path/to/base_dir>" line 19

rna_config.json : change path to indexes, adapters, and software versions...main file to edit

## Index Information ##

STAR and Kallisto indexes must be created prior to running the main RNAseq_v3.0.py snakefile analysis

Ensembl assembly/annotation can be downloaded here: [Ensembl](https://ftp.ensembl.org/pub/current_fasta/)

NCBI assembly/amnnotation can be downloaded here: [NCBI](https://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Mus_musculus/)

* NEEDED files prior to creating indexes:
    + Primary assembly FASTA : Mus_musculus.GRCm39.dna.primary_assembly.fa.gz
    + cDNA FASTA : concatenated cDNA and ncRNA FASTA (Mus_musculus.GRCm39.cdna.all.fa.gz, Mus_musculus.GRCm39.ncrna.fa.gz)
    + GTF : Mus_musculus.GRCm39.109.gtf.gz


## Edit when making indexes ##
config.json : edit the path information for the destination indexes

*Make sure the version for STAR and kallisto are the same as needed in rna_config.json for the RNAseq_v3.0.py snakefile*


