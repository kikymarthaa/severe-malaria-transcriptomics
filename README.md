# Transcriptomics Analysis of GSE117613

## Overview

This project presents a transcriptomic analysis of the public microarray dataset **GSE117613** to identify Differentially Expressed Genes (DEGs) between **Cerebral Malaria (CM)** and **Severe Malarial Anemia (SMA)**.

The analysis was performed using **R** and **Bioconductor**, including differential gene expression analysis, functional enrichment analysis, and visualization.

---

## Dataset

- **Accession:** GSE117613
- **Platform:** Illumina HumanHT-12 V4.0 Expression BeadChip
- **Database:** NCBI Gene Expression Omnibus (GEO)

Comparison:

- Cerebral Malaria (CM)
- Severe Malarial Anemia (SMA)

---

## Analysis Workflow

1. Download dataset from GEO
2. Data preprocessing and log2 transformation
3. Differential expression analysis using **limma**
4. Gene annotation
5. Quality control visualization
   - Boxplot
   - Density Plot
   - UMAP
6. Differential expression visualization
   - Volcano Plot
   - Heatmap (Top 50 DEGs)
7. Functional enrichment analysis
   - Gene Ontology (BP, CC, MF)
   - KEGG Pathway

---

## Repository Structure

```text
.
├── data/
├── figures/
├── results/
├── script/
│   └── script.R
├── report.md
└── README.md
```

---

## Software

- R
- GEOquery
- limma
- clusterProfiler
- enrichplot
- ggplot2
- pheatmap
- org.Hs.eg.db
- illuminaHumanv4.db
- umap

---

## Outputs

### Figures

- Boxplot
- Density Plot
- UMAP
- Volcano Plot
- Heatmap of Top 50 DEGs
- GO Biological Process
- GO Cellular Component
- GO Molecular Function
- KEGG Pathway

### Results

- Differentially Expressed Genes (DEGs)
- GO enrichment results
- KEGG enrichment results

---

## Author

**Kiky Martha Ariesaka**

June 2026