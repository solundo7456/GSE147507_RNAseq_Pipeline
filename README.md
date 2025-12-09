# RNA-Seq Analysis of SARS-CoV-2 Infection in Primary Human Lung Epithelium (NHBE)

This project analyzes RNA-sequencing data from GSE147507, focusing specifically on primary human bronchial epithelial (NHBE) cells under two conditions:

Mock-treated

SARS-CoV-2-infected (24 h)

The full workflow goes from raw FASTQ files → QC → trimming → alignment → transcript quantification → differential expression → pathway analysis (GO/KEGG).

🔬 1. Data Acquisition
1.1 Download metadata

Downloaded GSE147507 metadata from GEO.

Filtered SraRunTable to include only NHBE samples.

Extracted Run accessions (SRR numbers) for downloading.

1.2 Download FASTQ files

Used SRA-tools:

prefetch SRRXXXXXX
fasterq-dump SRRXXXXXX --split-files --outdir fastq/

🧪 2. Quality Control
2.1 FastQC (raw reads)

Initial QC to inspect read quality:

fastqc fastq/*.fastq -o qc/raw/

2.2 fastp (trimming + QC)

Used to:

remove adapters

trim low-quality bases

generate HTML reports

fastp \
  -i sample.fastq \
  -o trimmed/sample_trimmed.fastq \
  -q 20 \
  -t 1 \
  -u 30 \
  -h qc/sample_fastp.html \
  -j qc/sample_fastp.json

2.3 FastQC (post-trim)
fastqc trimmed/*.fastq -o qc/trimmed/

🧬 3. Reference Genome Setup
3.1 Download reference (GRCh38 + GTF)

Obtained from Ensembl or GENCODE.

3.2 Build HISAT2 Index
hisat2-build genome.fa hisat2_index/genome

🧲 4. Alignment to the Genome

Aligned cleaned FASTQ reads using HISAT2:
hisat2 -p 8 \
  -x hisat2_index/genome \
  -U trimmed/${sample}_trimmed.fastq \
  -S bam/${sample}.sam

📦 5. BAM Processing (SAMtools)

Converted, sorted, and indexed alignments:

samtools view -bS bam/${sample}.sam > bam/${sample}.bam

samtools sort bam/${sample}.bam -o bam/${sample}_sorted.bam

samtools index bam/${sample}_sorted.bam

🧫 6. Transcript Assembly & Quantification (StringTie)

Used StringTie to compute transcript-level abundances.

6.1 Quantification per sample
stringtie sample_sorted.bam \
  -G annotation.gtf \
  -o stringtie/sample.gtf \
  -A stringtie/gene_abund.tab \
  -e -B


gene_abund.tab → gene TPMs

sample.gtf → transcript structures

-B → Ballgown-compatible output

6.2 Optional: Merge transcripts
stringtie --merge -G annotation.gtf \
  -o merged.gtf mergelist.txt

6.3 Re-quantify using merged reference
stringtie sample_sorted.bam \
  -G merged.gtf -e -o sample_merged.gtf

📊 7. Gene-Level Counting

Although StringTie gives TPM/FPKM, differential expression requires raw counts.

Used featureCounts:

featureCounts -T 8 \
  -a annotation.gtf \
  -o gene_count_matrix.txt \
  bam/*.bam


Generated:

gene_count_matrix.csv

transcript_count_matrix.csv (from StringTie)

📈 8. Differential Expression Analysis (DESeq2)
Load count matrix + metadata
Run DESeq2
Visualizations

PCA plot

Heatmaps

Volcano plot

Sample distance plot

8.4 DEG Export
write.csv(res, "DEG_results.csv")

🧭 9. Functional & Pathway Analysis

Used R packages:

clusterProfiler

org.Hs.eg.db

enrichplot

9.1 GO Enrichment
9.2 KEGG Pathways

Identified:

Antiviral response pathways

Interferon signaling

Cytokine signaling

Pattern recognition receptor activation



Reference:
Genome Reference Consortium (2013). GRCh38 – Human genome assembly.
Available at: https://www.ncbi.nlm.nih.gov/grc/human

Leinonen, R., Sugawara, H., Shumway, M. (2011). The Sequence Read Archive. Nucleic Acids Research, 39, D19–D21.
https://doi.org/10.1093/nar/gkq1019
