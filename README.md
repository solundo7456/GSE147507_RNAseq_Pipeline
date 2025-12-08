# GSE147507_RNAseq_Pipeline

📌 RNA-seq Analysis of NHBE Cells – GSE147507

This repository contains a full reproducible workflow for processing and analyzing RNA-seq data from Primary Human Lung Epithelium (NHBE) infected with SARS-CoV-2, IAV, and Mock controls.
The dataset originates from the public study GSE147507.

🧬 1. Project Overview

Goal:
To perform quality control, trimming, quantification, and differential gene expression analysis for the NHBE samples only.

Conditions included:

Mock

SARS-CoV-2 infection

Influenza A Virus (IAV) infection

Files produced:

gene_count_matrix.csv

transcript_count_matrix.csv

Processed metadata file

QC reports

R scripts for DESeq2 analysis

🔽 3. Downloading the Data
Using SRA Toolkit
# Download all FASTQs from the run list
cat SRR_Acc_List.txt | while read srr; do
    fasterq-dump $srr -O fastq/ --split-files
done

Filtering NHBE samples
grep -F -f nhbe_sra_list.txt SraRunTable.csv > metadata/SraRunTable_filtered.csv

🔧 4. Quality Control
Run fastp
fastp -i sample_1.fastq.gz -I sample_2.fastq.gz \
      -o trimmed/sample_1.fq.gz -O trimmed/sample_2.fq.gz \
      -h qc/sample_fastp.html

🎯 5. Quantification (Salmon)
Index building
salmon index -t transcripts.fa -i salmon_index

Run quantification
salmon quant -i salmon_index -l A \
             -1 trimmed/sample_1.fq.gz \
             -2 trimmed/sample_2.fq.gz \
             -p 8 -o counts/sample/

📊 6. Create Count Matrices

After quantification:

gene_count_matrix.csv

transcript_count_matrix.csv

Generated using:

tximport(...)

📈 7. Differential Expression (DESeq2)
Example R code:
dds <- DESeqDataSetFromTximport(txi, metadata, ~ treatment)
dds <- DESeq(dds)
res <- results(dds, contrast = c("treatment", "SARSCoV2", "Mock"))


Outputs include:

Normalized counts

PCA plots

Volcano plots

DEG tables

📦 8. Requirements
Software

SRA Toolkit

fastp

Salmon / Kallisto

R (≥4.0)

R packages:

DESeq2

tximport

tidyverse

pheatmap

📝 9. Citation

Please cite the original study when using this data:

Blanco-Melo et al., Imbalanced host response to SARS-CoV-2 drives development of COVID-19. Cell (2020)
