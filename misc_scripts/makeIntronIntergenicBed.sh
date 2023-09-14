#!/bin/bash

#:'
#adapted from :
#http://crazyhottommy.blogspot.com/2013/05/find-exons-introns-and-intergenic.html
#'

#Load modules
module load bedtools;

#Generate exon bed file, sort, and merge to get exon.bed
awk 'BEGIN{OFS="\t";} $3=="exon" {print $1,$4-1,$5}' < Mus_musculus.GRCm38.94.gtf > MmENSv94_exon.bed
sortBed -i MmENSv94_exon.bed > MmENSv94_exon_tmp.bed
mv -f MmENSv94_exon_tmp.bed MmENSv94_exon.bed
mergeBed -i MmENSv94_exon.bed > MmENSv94_exon_merge.bed

#Define gene coordinates and subtract exons to get intron.bed
awk 'BEGIN{OFS="\t";} $3=="gene" {print $1,$4-1,$5}' < Mus_musculus.GRCm38.94.gtf > MmENSv94_gene.bed
sortBed -i MmENSv94_gene.bed > MmENSv94_gene_tmp.bed
mv -f MmENSv94_gene_tmp.bed MmENSv94_gene.bed
subtractBed -a MmENSv94_gene.bed -b MmENSv94_exon_merge.bed > MmENSv94_intron.bed

#Make integenic.bed
sortBed -i MmENSv94_gene.bed -g MmENSv94_chrLength.txt > MmENSv94_gene_genSort.bed
complementBed -i MmENSv94_gene_genSort.bed -g MmENSv94_chrLength.txt > MmENSv94_intergenic.bed
