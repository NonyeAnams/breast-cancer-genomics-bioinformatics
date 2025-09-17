############################################################
# Stage 1: Data Retrieval and Preparation (TCGA-BRCA)
# Goal: Download RNA-seq STAR counts (10 Tumor, 10 Normal)
#       Clean gene IDs and build counts + metadata matrices
############################################################

# ---- Load libraries ----
library(TCGAbiolinks)
library(SummarizedExperiment)
library(tidyverse)

# ---- Step 1: Query all BRCA RNA-seq STAR counts ----
query_all <- GDCquery(
  project       = "TCGA-BRCA",
  data.category = "Transcriptome Profiling",
  data.type     = "Gene Expression Quantification",
  workflow.type = "STAR - Counts",
  sample.type   = c("Primary Tumor", "Solid Tissue Normal")
)

# Preview available samples
results_all <- getResults(query_all)
head(results_all)

# ---- Step 2: Separate Tumor vs Normal barcodes ----
tumor_samples  <- TCGAquery_SampleTypes(barcode = results_all$cases, typesample = "TP")
normal_samples <- TCGAquery_SampleTypes(barcode = results_all$cases, typesample = "NT")

cat("Available tumor samples:", length(tumor_samples), "\n")
cat("Available normal samples:", length(normal_samples), "\n")

# ---- Step 3: Randomly select a subset (10 tumor + 10 normal) ----
set.seed(123)  # reproducibility
tumor_subset  <- sample(tumor_samples, 10)
normal_subset <- sample(normal_samples, 10)
selected_samples <- c(tumor_subset, normal_subset)

# Save sample info
sample_info <- data.frame(
  barcode = selected_samples,
  type    = c(rep("Tumor", length(tumor_subset)),
              rep("Normal", length(normal_subset)))
)
write.csv(sample_info, "data/selected_barcodes.csv", row.names = FALSE)

# ---- Step 4: Query only the selected samples ----
query_subset <- GDCquery(
  project       = "TCGA-BRCA",
  data.category = "Transcriptome Profiling",
  data.type     = "Gene Expression Quantification",
  workflow.type = "STAR - Counts",
  barcode       = selected_samples
)

# Download data for selected samples
GDCdownload(query_subset, files.per.chunk = 20)

# ---- Step 5: Prepare data into SummarizedExperiment ----
data_se <- GDCprepare(query_subset)

# Extract raw counts and metadata
counts_raw <- assay(data_se)               # gene counts
meta <- as.data.frame(colData(data_se))    # sample metadata

# ---- Step 6: Clean Ensembl IDs (remove version suffix) ----
rownames(counts_raw) <- gsub("\\..*$","", rownames(counts_raw))
counts <- counts_raw   # if duplicates exist, aggregate later

# ---- Step 7: Identify Tumor vs Normal barcodes ----
barcodes        <- colnames(counts)
tumor_barcodes  <- TCGAquery_SampleTypes(barcodes, typesample = "TP")
normal_barcodes <- TCGAquery_SampleTypes(barcodes, typesample = "NT")

cat("Total samples in subset:", ncol(counts), "\n")
cat("Tumor samples (TP):", length(tumor_barcodes), "\n")
cat("Normal samples (NT):", length(normal_barcodes), "\n")

# ---- Step 8: Build filtered counts matrix + coldata ----
keep_barcodes <- c(tumor_barcodes, normal_barcodes)
counts_sub    <- counts[, keep_barcodes]

sample_type <- ifelse(colnames(counts_sub) %in% tumor_barcodes, "Tumor", "Normal")
coldata <- data.frame(
  sample    = colnames(counts_sub),
  condition = sample_type,
  stringsAsFactors = FALSE
)
rownames(coldata) <- coldata$sample

# ---- Step 9: Save outputs ----
# counts 
write.csv(head(counts_sub,20), "data/counts_sub_head20.csv", row.names = TRUE)

saveRDS(counts_sub, "data/counts_sub.rds")
write.csv(counts_sub, "data/counts_sub.csv", row.names = TRUE)

# Metadata for selected subset
meta_subset <- getResults(query_subset)
write.csv(meta_subset, "data/tcga_brca_expression_metadata.csv", row.names = FALSE)

# Condition labels
write.csv(coldata, "data/coldata_br_ca.csv", row.names = FALSE)

############################################################
# End of Stage 1: You now have:
# - counts_sub: gene counts (20 samples, cleaned)
# - coldata: tumor vs normal labels
# - metadata CSVs for reference
############################################################
