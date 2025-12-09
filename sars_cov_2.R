
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("DESeq2")
library(GEOquery)
library(limma)
library(statmod)
library(dplyr)
library(tibble)
library(DESeq2)
library(ggplot2)


#creating table
metadata = read.csv('metadata.csv')
head(metadata)
str(metadata)

#subset
metadata2 <- metadata %>% select(Run,treatment,time_point)%>%
                       mutate(treatment= ifelse(grepl('Mock',treatment),"Mock",
                             ifelse(grepl('SARS',treatment),"Sars-infected",treatment)))
print(metadata2)
class(metadata2)
# making metadata dseq ready
sampleData <- metadata2 %>% rename(sampleName = Run,
                                   condition = treatment,
                                   time_point = time_point) %>%
              mutate(condition = factor(condition)) %>%
              column_to_rownames(var='sampleName')
print(sampleData)

#import count data
getwd()
path="~/Downloads/bioinformatics/sars_cov_2/project_GSE147507/count/actual_counts/gene_count_matrix.csv"
count_data = read.csv(path, header=TRUE)
head(count_data)
# Move 'gene_id' to row names and remove the column
rownames(count_data) <- count_data$gene_id
count_data <- count_data[ , -1]  # remove the gene_id column
head(count_data)


#confirming we have all samples
all(colnames(count_data) == rownames(sampleData))
colnames(count_data)
rownames(sampleData)
# Keep only samples present in count_data
sampleData <- sampleData[colnames(count_data), ]
all(colnames(count_data) == rownames(sampleData))

#quality check for count data
dim(count_data)
summary(count_data)
#filter remove genes with low expression
keep <- rowSums(count_data>10) >=5
counts_filtered <- count_data[keep,]
dim(counts_filtered)

# Total counts per sample
library_sizes <- colSums(counts_filtered)
summary(library_sizes)
hist(log10(library_sizes), main='library sizes (log10)', xlab='log10(Total counts)')
# Quick barplot
barplot(library_sizes, las=2, main="Library Sizes", col="skyblue")

library_sizes
# Convert to a data frame and join with metadata
library_sizes_df <- data.frame(
  sample = names(library_sizes),
  library_size = library_sizes
) %>%
  left_join(
    sampleData %>% mutate(sample = rownames(sampleData)),
    by = "sample"
  )
library_sizes_df
#plot
ggplot(library_sizes_df, aes(x = condition, y = library_size, fill = condition)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label=sample), vjust=-0.5, size=3) + # optional: show sample IDs
  theme_bw() +
  labs(title = "Library Sizes by Condition",
       x = "Condition",
       y = "Total Counts") +
  theme(legend.position = "none")


#creating dseq object
dds <- DESeqDataSetFromMatrix(
      countData = count_data,
      colData = sampleData,
      design = ~condition)

#dseq model
dds <- DESeq(dds)
plotDispEsts(dds)
res <- results(dds, alpha = 0.05)

#number of significate genes
sum(res$padj < 0.05, na.rm = TRUE)  # number of DE genes at FDR 5%

#MA-plot
plotMA(res, main="DESeq2 MA-plot")

#volcano plot



BiocManager::install("EnhancedVolcano")

library(EnhancedVolcano)
EnhancedVolcano(res,
                lab = rownames(res),
                x = 'log2FoldChange',
                y = 'padj',
                pCutoff = 0.05,
                FCcutoff = 1,
                title = "Positive vs Negative")
#Heatmap 
library(pheatmap)
top_genes <- head(order(res$padj), 50)
# top 50 DE genes
mat <- assay(vst(dds))[top_genes, ]
colnames(meta2)

pheatmap(mat, scale="row", annotation_col = meta2[, c("infection_status","gender")])
write.csv(as.data.frame(res), "DE_results_Positive_vs_Negative.csv")

# 2. Filter significant genes
# -----------------------------
sig_genes <- rownames(res)[which(res$padj < 0.05 & abs(res$log2FoldChange) > 1)]
length(sig_genes)
head(sig_genes)
sig_genes
# 3. Map gene symbols to Entrez IDs
# -----------------------------
entrez_ids <- bitr(sig_genes,
                   fromType = "SYMBOL",
                   toType = "ENTREZID",
                   OrgDb = org.Hs.eg.db)
head(entrez_ids)

# 4. GO enrichment (Biological Process)
# -----------------------------
ego <- enrichGO(gene         = entrez_ids$ENTREZID,
                OrgDb        = org.Hs.eg.db,
                keyType      = "ENTREZID",
                ont          = "BP",
                pAdjustMethod= "BH",
                pvalueCutoff = 0.05,
                qvalueCutoff = 0.05)

ego_1 <- enrichGO(gene         = entrez_ids$ENTREZID,
                  OrgDb        = org.Hs.eg.db,
                  keyType      = "ENTREZID",
                  ont          = "MF",
                  pAdjustMethod= "BH",
                  pvalueCutoff = 0.05,
                  qvalueCutoff = 0.05)
# Top GO terms
head(ego)
head(ego_1)
# Dotplot
dotplot(ego, showCategory=20) + ggtitle("GO Biological Process Enrichment")

# 5. KEGG pathway enrichment
# -----------------------------
options(timeout = 600) 
ekegg <- enrichKEGG(gene         = entrez_ids$ENTREZID,
                    organism     = "hsa",
                    pvalueCutoff = 0.05)

# Top KEGG pathways
head(ekegg)
dotplot(ekegg, showCategory=20) + ggtitle("KEGG Pathway Enrichment")

# -----------------------------
# 6. Reactome pathway enrichment
# -----------------------------
ereact <- enrichPathway(gene=entrez_ids$ENTREZID,
                        organism="human",
                        pvalueCutoff=0.05,
                        readable=TRUE)

# Top Reactome pathways
head(ereact)
dotplot(ereact, showCategory=20) + ggtitle("Reactome Pathway Enrichment")
