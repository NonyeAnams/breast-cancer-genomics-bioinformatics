# Breast Cancer Genomics Bioinformatics Project

This is an independent bioinformatics project analyzing **TCGA breast cancer (BRCA) RNA-seq data**.  
The goal is to identify differentially expressed genes between tumor and normal samples, perform pathway enrichment analysis, and explore network-level interactions.

---

## Objectives
- Retrieve and preprocess BRCA RNA-seq data from TCGA (via TCGAbiolinks).
- Perform quality control and normalization.
- Conduct differential expression analysis (DESeq2).
- Identify significantly enriched pathways and biological processes.
- Explore protein–protein interaction networks (STRING, Cytoscape).
- Summarize findings in a mini-report.

---

## Repository Structure

breast-cancer-genomics-bioinformatics/
│
├── data/ # metadata and small processed CSVs (no raw TCGA data)
├── scripts/ # R scripts (organized by stage)
├── results/ # plots, tables, top DEGs, enrichment results
├── report/ # mini-report (RMarkdown, Word, or PDF)
├── README.md # project overview (this file)
└── .gitignore # ignore rules for large/raw files

---

## Workflow
1. **Setup** → install R packages (01_setup_install_packages.R).  
2. **Data download** → TCGA-BRCA HTSeq counts (02_download_data.R).  
3. **QC + DESeq2 analysis** (03_qc_deseq2_analysis.R).  
4. **Pathway enrichment** (04_pathway_enrichment.R).  
5. **Network analysis** (05_network_analysis.R).  
6. **Figures & reporting** (06_report_plots.R).  

---

## Tools and Packages
- **R/Bioconductor:** TCGAbiolinks, DESeq2, clusterProfiler, fgsea, EnhancedVolcano, pheatmap.  
- **Databases:** TCGA, KEGG, GO, STRING.  
- **Visualization:** ggplot2, Cytoscape.  

---

## Notes
- Raw TCGA data is not uploaded here due to size restrictions.  
- This repo contains **scripts, processed subsets, results, and report materials**.  