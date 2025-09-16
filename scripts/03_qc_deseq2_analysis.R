############################################################
# Script: 02_qc_and_DE.R
# Goal: Quick QC + first DESeq2 run for TCGA-BRCA subset
############################################################

# -------------------------------
# Step 1: Load Libraries
# -------------------------------
# Core libraries
library(DESeq2)
library(tidyverse)

# QC and plotting
library(ggplot2)
library(pheatmap)
library(RColorBrewer)

# -------------------------------
# Step 2: Load Data
# -------------------------------
# Load counts and metadata
counts <- read.csv("data/counts_sub_head20.csv", row.names = 1)
coldata <- read.csv("data/coldata_br_ca.csv", row.names = 1)

# Check alignment of samples
all(rownames(coldata) == colnames(counts))

# If FALSE, inspect mismatches
colnames(counts)
rownames(coldata)

# Fix sample naming mismatch (replace '-' with '.')
rownames(coldata) <- gsub("-", ".", rownames(coldata))

# Re-check alignment
all(rownames(coldata) == colnames(counts))

# Make condition a factor with explicit order
coldata$condition <- factor(coldata$condition, levels = c("Normal", "Tumor"))

# -------------------------------
# Step 3: Create DESeq2 Dataset
# -------------------------------
dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = coldata,
                              design = ~ condition)

# -------------------------------
# Step 4: Quick QC
# -------------------------------
# 1. Filter out low counts
dds <- dds[rowSums(counts(dds)) > 10, ]

# 2. Variance stabilizing transformation
vsd <- varianceStabilizingTransformation(dds, blind = TRUE)

# 3. PCA Plot
p <- plotPCA(vsd, intgroup = "condition")
ggsave("figures/pca_plot.pdf", plot = p, width = 6, height = 5)

# Save PCA coordinates
pca_df <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
write.csv(pca_df, "data/processed/pca_coords_vsd.csv", row.names = TRUE)

# 4. Sample distances heatmap
sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- vsd$condition

pdf("figures/sample_dist_heatmap.pdf", width = 8, height = 6)
pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         col = colorRampPalette(rev(brewer.pal(9, "Blues")))(255))
dev.off()

# -------------------------------
# Step 5: Run Differential Expression
# -------------------------------
dds <- DESeq(dds)

# Extract results
res <- results(dds, contrast = c("condition", "Tumor", "Normal"))

# Order by adjusted p-value
res <- res[order(res$padj), ]

# Save DE results table
write.csv(as.data.frame(res),
          "data/processed/DESeq2_results.csv",
          row.names = TRUE)

# Quick peek at top 10
head(res, 10)

# -------------------------------
# Step 6: Volcano Plot
# -------------------------------
res_df <- as.data.frame(res)

vplot <- ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(alpha = 0.4) +
  geom_hline(yintercept = -log10(0.05), color = "red", linetype = "dashed") +
  geom_vline(xintercept = c(-1,1), color = "blue", linetype = "dashed") +
  theme_minimal()

# Save volcano plot
ggsave("figures/volcano_plot.pdf", plot = vplot, width = 7, height = 6)

# -------------------------------
# Step 7: Heatmap of Top 50 DEGs
# -------------------------------
# Select top 50 DEGs (by adjusted p-value)
top50_genes <- rownames(head(res, 50))

# Save list of top 50 genes
write.csv(top50_genes,
          "data/processed/top50_DEGs.csv",
          row.names = FALSE)

# Extract variance-stabilized expression for those genes
mat <- assay(vsd)[top50_genes, ]

# Z-score normalize across samples (for better contrast in heatmap)
mat <- t(scale(t(mat)))

# Create heatmap and save
pdf("figures/heatmap_top50.pdf", width = 8, height = 10)
pheatmap(mat,
         show_rownames = TRUE,
         show_colnames = FALSE,
         annotation_col = as.data.frame(coldata["condition"]),
         fontsize_row = 6,
         clustering_distance_rows = "correlation",
         clustering_distance_cols = "correlation",
         col = colorRampPalette(rev(brewer.pal(9, "RdBu")))(255))
dev.off()

# -------------------------------
# Step 8: Save Normalized Counts
# -------------------------------
norm_counts <- counts(dds, normalized = TRUE)
write.csv(as.data.frame(norm_counts),
          "data/processed/normalized_counts.csv",
          row.names = TRUE)

############################################################
# End of Stage 2: You now have:
# - QC plots: PCA + sample distance heatmap
# - DE results: DESeq2_results.csv
# - Volcano plot: volcano_plot.pdf
# - Heatmap of top 50 DEGs: heatmap_top50.pdf
# - List of top 50 DEGs: top50_DEGs.csv
# - Normalized counts table: normalized_counts.csv
############################################################
