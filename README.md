# RNA-seq
This is the basic bulk RNA sequencing pipeline.

The pipeline was based off the NNRL lab pipeline.
This pipeline adapted to run on HPCs running SLURM
This requires the snakefile RNAseq_v3.0.py, rna_config.json, meta.csv and starts using the rna_submit_snakemake.sh script

# Run with: sbatch --time=24:00:00 rna_submit_snakemake.sh


## Edit the following files ##
meta.csv : add the base FQ names (minus the [R1|R2]_ending.fastq.gz)

rna_submit_snakemake.sh : change "<path/to/base_dir>" line 19

rna_config.json : change path to indexes, adapters, and software versions...main file to edit


