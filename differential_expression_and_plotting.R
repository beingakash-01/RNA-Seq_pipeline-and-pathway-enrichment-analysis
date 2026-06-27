library(GEOquery)
library(tidyverse)
library(dplyr)
library(DESeq2)
library(ggplot2)

gse <- getGEO("GSE52778", GSEMatrix = TRUE)

metadata <- pData(phenoData(gse[[1]]))

metadata <- metadata[ , c(2, 49)]
metadata <- metadata%>%
  rename("treatment:ch1" = "treatment")

rownames(metadata) <- NULL 

metadata <- metadata%>%
  column_to_rownames(var = "geo_accession")

metadata <- metadata[-c(3, 4, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16), ]


final_countmatrix <- read.table("D:/Projects/RNA_seqpipeline/feature_counts/final_count_matrix.txt", header = TRUE, sep = "\t")

final_countmatrix <- final_countmatrix[, c(1, 7, 8, 9, 10)]
final_countmatrix <- final_countmatrix%>%
  rename("aligned_SRR1039508.bam" = "GSM1275862",
         "aligned_SRR1039509.bam" = "GSM1275863",
         "aligned_SRR1039512.bam" = "GSM1275866",
         "aligned_SRR1039513.bam" = "GSM1275867")

final_countmatrix <- final_countmatrix%>%
  column_to_rownames(var = "Geneid")
 
write.table(final_countmatrix, file = "counts_data.csv", sep = ',', col.names = T, row.names = T, quote = F)
write.table(metadata, file = "meta_data.csv", sep = ',', col.names = T, row.names = T, quote = F)


#final_countmatrix <- final_countmatrix[, c(2, 5)]
  

all(colnames(final_countmatrix) %in% rownames(metadata))
all(colnames(final_countmatrix) == rownames(metadata))

#metadata <- metadata[-c(1, 3, 4, 6)]

DESeqmatrix <- DESeqDataSetFromMatrix(countData = final_countmatrix,
                                      colData = metadata,
                                      design = ~ treatment) 

keep <- rowSums(counts(DESeqmatrix)) >= 10
matrix <- DESeqmatrix[keep,]                
matrix


matrix$variation <- relevel(matrix$variation, ref = 'Dexamethasone')

#running DESeq2
matrix <- DESeq(matrix)
res <- results(matrix, contrast = c('treatment', 'Untreated', 'Dexamethasone')) 
res

summary(res)
resultsNames(matrix)

results(matrix, contrast = c('treatment', 'Untreated', 'Dexamethasone'))

plotMA(res)

res_df <- as.data.frame(res)
res_df <- res_df[!is.na(res_df$padj), ]
res_df$diffexpressed <- "Not Significant"
res_df$diffexpressed[res_df$log2FoldChange > 1 & res_df$padj < 0.05] <- "Up-regulated"
res_df$diffexpressed[res_df$log2FoldChange < 1 & res_df$padj > 0.05] <- "Down-regulated"

ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj), color = diffexpressed )) +
  geom_point(size = 1.5, alpha = 0.8)+
  scale_color_manual(values = c("Up-regulated" = "blue",
                                "Not Significant" = "grey",
                                "Down-regulated" = "red"))+
  geom_vline(xintercept = c(-1, 1), col = "black", linetype = "dashed")+
  geom_hline(yintercept = -log10(0.05), col = "black", linetype = "dashed")+
  theme_minimal()+
  labs(title = "Volcano plot : Dexamethasone vs Untreated",
       x = "Log2fold change",
       y = "-Log10 adjusted P-value",
       color = "Gene Status")

#--------------------------------------------------------------------------------------------------------------------------------
#PHASE 2: Filtering significant genes (upregulated and downregulted) and generating heatmap based on assay data of the sifnificant genes

count_matrix_ids <- final_countmatrix%>%
  rownames_to_column(var = "Ensembl_ids")

count_matrix_ids <- count_matrix_ids[1]

write.table(count_matrix_ids, file = "Ensemble_gene_ids.csv", sep = ',', col.names = T, row.names = T, quote = F)

up_regulated <- rownames_to_column(up_regulated, var = "Ensembl_ids")
sig_genes <- rownames_to_column(sig_genes, var = "Ensembl_ids")

library(biomaRt)

listEnsembl()
ensemble <- useEnsembl(biomart = "genes")
datasets <- listDatasets(ensemble)

ensembl.con <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")


attributes_ensembl <- listAttributes(ensembl.con)
filters <- listFilters(ensembl.con)

gene_ids <- getBM(attributes = c("ensembl_gene_id", "external_gene_name"),
                  filters = "ensembl_gene_id",
                  values = up_regulated$Ensembl_ids,
                  mart = ensembl.con)

sig_genes <- res_df %>%
  filter(padj < 0.05 & abs(log2FoldChange) > 1) %>%
  arrange(padj)

up_regulated <- sig_genes %>%
  filter(log2FoldChange > 1)

down_regulated <- sig_genes %>%
  filter(log2FoldChange < -1)

print(paste("Total Significant DEGs:", nrow(sig_genes)))
print(paste("Upregulated:", nrow(up_regulated)))
print(paste("Downregulated:", nrow(down_regulated)))


write.csv(sig_genes, file = "significant_deg_results.csv", row.names = FALSE)

# Save just the gene IDs/Symbols for enrichment analysis (e.g., DAVID, g:Profiler)
write_table(up_regulated, file = "up_regulated_genes.csv", sep = ',', col.names = T, row.names = T, quote = F)

up_reg_genes <- up_regulated$Ensembl_ids
write.table(up_regulated$Ensembl_ids, file = "up_regulated_genes.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)

down_reg_genes <- down_regulated[1]
down_regulated <- down_regulated%>%
  rownames_to_column(var = "Ensembl_ids")

write.table(down_regulated$Ensembl_ids, file = "down_regulated_genes.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(down_regulated, file = "down_regulated_genes.csv", sep = ',', col.names = T, row.names = T, quote = F)

library(pheatmap)

vst <- vst(matrix, blind = FALSE)
assay_data <- assay(vst)

sig_gene_ids <- sig_genes$Ensembl_ids
heatmap_matrix <- assay_data[sig_gene_ids, ]

annotation_col <- as.data.frame(colData(vst)[, c("treatment"), drop=FALSE])

#generating the heatmap.
pheatmap(heatmap_matrix, 
         scale = "row",                      # Scale genes by row (Z-score)
         clustering_distance_rows = "euclidean", 
         clustering_distance_cols = "euclidean",
         clustering_method = "complete",     # Standard hierarchical clustering
         annotation_col = annotation_col,    # Adds the top color bar for conditions
         show_rownames = FALSE,              # Turn off if list is too long
         show_colnames = TRUE,               # Show sample names
         main = "Clustered Heatmap of Top Significant DEGs",
         color = colorRampPalette(c("navy", "white", "firebrick3"))(50))

#-------------------------------------------------------------------------------------------------------------------------------
#GO(gene ontology) enrichment analysis

library(clusterProfiler)
library(org.Hs.eg.db)

# 1. Run GO Enrichment for Upregulated Genes
# (Assuming your gene IDs are ENSEMBL IDs. Change keyType if using Gene Symbols)
go_enrich_up <- enrichGO(gene          = up_regulated$Ensembl_id,
                         OrgDb         = org.Hs.eg.db,
                         keyType       = "ENSEMBL",
                         ont           = "BP",           # BP = Biological Process (can also use CC or MF)
                         pAdjustMethod = "BH",           # Benjamini-Hochberg adjustment
                         pvalueCutoff  = 0.05,
                         qvalueCutoff  = 0.05)

# 2. Visualize the top enriched pathways using a dotplot
dotplot(go_enrich_up, showCategory = 15, title = "GO Biological Processes - Upregulated Genes")

