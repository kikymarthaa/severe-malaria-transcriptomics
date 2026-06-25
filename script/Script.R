############################################################
# Transcriptomics Analysis using GSE117613
# Comparison : Cerebral Malaria (CM) vs Severe Malarial Anemia (SMA)
# Platform   : Illumina HumanHT-12 V4.0 Expression BeadChip
# Dataset    : GSE117613
# Method     : limma + GO + KEGG enrichment
#
# Author     : Kiky Martha Ariesaka
# Date       : June 2026
############################################################
rm(list = ls())

options(stringsAsFactors = FALSE)

set.seed(123)

############################################################
# Create output folders
############################################################

if (!dir.exists("figures"))
  dir.create("figures")

if (!dir.exists("results"))
  dir.create("results")

############################################################
# PART A. INSTALL PACKAGES
############################################################

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

packages_bioc <- c(
  "GEOquery",
  "limma",
  "illuminaHumanv4.db",
  "AnnotationDbi",
  "clusterProfiler",
  "org.Hs.eg.db",
  "enrichplot"
)

for (pkg in packages_bioc) {
  if (!requireNamespace(pkg, quietly = TRUE))
    BiocManager::install(pkg, ask = FALSE, update = FALSE)
}

packages_cran <- c(
  "ggplot2",
  "dplyr",
  "pheatmap",
  "umap"
)

for (pkg in packages_cran) {
  if (!requireNamespace(pkg, quietly = TRUE))
    install.packages(pkg)
}

############################################################
# PART B. LOAD LIBRARIES
############################################################

library(GEOquery)
library(limma)
library(ggplot2)
library(dplyr)
library(pheatmap)
library(AnnotationDbi)
library(illuminaHumanv4.db)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(umap)

cat("Packages loaded successfully.\n")

############################################################
# PART C. DOWNLOAD DATASET
############################################################

gset <- getGEO(
  "GSE117613",
  GSEMatrix = TRUE,
  AnnotGPL = TRUE
)[[1]]

############################################################
# PART D. EXPRESSION MATRIX
############################################################

ex <- exprs(gset)

############################################################
# PART E. LOG2 TRANSFORMATION (IF NEEDED)
############################################################

qx <- as.numeric(
  quantile(
    ex,
    c(0, 0.25, 0.5, 0.75, 0.99, 1),
    na.rm = TRUE
  )
)

LogTransform <- (
  qx[5] > 100 ||
    (qx[6] - qx[1] > 50 && qx[2] > 0)
)

if (LogTransform) {
  ex[ex <= 0] <- NA
  ex <- log2(ex)
}

############################################################
# PART F. CHECK DATASET
############################################################

cat("Number of genes   :", nrow(ex), "\n")
cat("Number of samples :", ncol(ex), "\n")

############################################################
# PART G. SAMPLE INFORMATION
############################################################

pheno <- pData(gset)

head(pheno)

############################################################
# PART H. SAMPLE GROUPING
############################################################

group_info <- as.character(
  pData(gset)$source_name_ch1
)

table(group_info)

groups <- make.names(group_info)

factor_groups <- factor(groups)

levels(factor_groups)

############################################################
# PART I. KEEP ONLY CM AND SMA SAMPLES
############################################################

keep <- factor_groups %in% c(
  "cerebral.milaria",
  "severe.malaral.anemia"
)

ex2 <- ex[, keep]

group2 <- factor(
  factor_groups[keep]
)

levels(group2) <- c(
  "CM",
  "SMA"
)

table(group2)

############################################################
# PART J. DESIGN MATRIX
############################################################

design <- model.matrix(~0 + group2)

colnames(design) <- c(
  "CM",
  "SMA"
)

design

############################################################
# PART K. CONTRAST MATRIX
############################################################

contrast_matrix <- makeContrasts(
  CM_vs_SMA = CM - SMA,
  levels = design
)

contrast_matrix

############################################################
# PART L. DIFFERENTIAL EXPRESSION ANALYSIS
############################################################

fit <- lmFit(
  ex2,
  design
)

fit2 <- contrasts.fit(
  fit,
  contrast_matrix
)

fit2 <- eBayes(fit2)

############################################################
# PART M. EXTRACT DIFFERENTIALLY EXPRESSED GENES
############################################################

deg_CM_vs_SMA_all <- topTable(
  fit2,
  coef = "CM_vs_SMA",
  number = Inf,
  adjust.method = "BH",
  sort.by = "P"
)

deg_CM_vs_SMA_all$PROBEID <- rownames(
  deg_CM_vs_SMA_all
)

deg_CM_vs_SMA <- subset(
  deg_CM_vs_SMA_all,
  adj.P.Val < 0.05 &
    abs(logFC) > 0.5
)

cat(
  "Number of significant DEGs :",
  nrow(deg_CM_vs_SMA),
  "\n"
)

############################################################
# PART N. GENE ANNOTATION
############################################################

probe_ids <- rownames(ex2)

gene_annotation <- AnnotationDbi::select(
  illuminaHumanv4.db,
  keys = probe_ids,
  columns = c(
    "SYMBOL",
    "GENENAME"
  ),
  keytype = "PROBEID"
)

deg_CM_vs_SMA_all$PROBEID <- rownames(
  deg_CM_vs_SMA_all
)

deg_CM_vs_SMA_all <- merge(
  deg_CM_vs_SMA_all,
  gene_annotation,
  by = "PROBEID",
  all.x = TRUE
)

deg_CM_vs_SMA <- merge(
  deg_CM_vs_SMA,
  gene_annotation,
  by = "PROBEID",
  all.x = TRUE
)

head(
  deg_CM_vs_SMA[
    ,
    c(
      "PROBEID",
      "SYMBOL",
      "GENENAME"
    )
  ]
)

############################################################
# PART O.1 BOXPLOT DISTRIBUSI NILAI EKSPRESI
############################################################

group_colors <- as.numeric(group2)

boxplot(
  ex2,
  col = group_colors,
  las = 2,
  outline = FALSE,
  main = "Boxplot Distribusi Nilai Ekspresi per Sampel",
  ylab = "Expression Value (log2)"
)

legend(
  "topright",
  legend = levels(group2),
  fill = unique(group_colors),
  cex = 0.8
)

############################################################
# PART O.2 DISTRIBUSI NILAI EKSPRESI (DENSITY PLOT)
############################################################

expr_long <- data.frame(
  Expression = as.vector(ex2),
  Group = rep(group2, each = nrow(ex2))
)

density_plot <- ggplot(
  expr_long,
  aes(
    x = Expression,
    color = Group
  )
) +
  geom_density(
    linewidth = 1
  ) +
  theme_minimal() +
  labs(
    title = "Distribusi Nilai Ekspresi Gen (CM vs SMA)",
    x = "Expression Value (log2)",
    y = "Density"
  )

density_plot

ggsave(
  filename = "figures/Density_CM_vs_SMA.png",
  plot = density_plot,
  width = 8,
  height = 6,
  dpi = 300
)

############################################################
# PART O.3 UMAP (VISUALISASI DIMENSI RENDAH)
############################################################

umap_input <- t(ex2)

umap_result <- umap(
  umap_input
)

umap_df <- data.frame(
  UMAP1 = umap_result$layout[, 1],
  UMAP2 = umap_result$layout[, 2],
  Group = group2
)

umap_plot <- ggplot(
  umap_df,
  aes(
    x = UMAP1,
    y = UMAP2,
    color = Group
  )
) +
  geom_point(
    size = 3,
    alpha = 0.8
  ) +
  theme_minimal() +
  labs(
    title = "UMAP Plot (CM vs SMA)",
    x = "UMAP1",
    y = "UMAP2"
  )

umap_plot

ggsave(
  filename = "figures/UMAP_CM_vs_SMA.png",
  plot = umap_plot,
  width = 8,
  height = 6,
  dpi = 300
)

############################################################
# PART P.1 VOLCANO PLOT (CM vs SMA)
############################################################

volcano_data <- deg_CM_vs_SMA_all

volcano_data$status <- "Not Significant"

volcano_data$status[
  volcano_data$logFC > 0.5 &
    volcano_data$adj.P.Val < 0.05
] <- "Upregulated"

volcano_data$status[
  volcano_data$logFC < -0.5 &
    volcano_data$adj.P.Val < 0.05
] <- "Downregulated"

volcano_plot <- ggplot(
  volcano_data,
  aes(
    x = logFC,
    y = -log10(adj.P.Val),
    color = status
  )
) +
  geom_point(
    size = 1.8,
    alpha = 0.7
  ) +
  scale_color_manual(
    values = c(
      "Upregulated" = "red",
      "Downregulated" = "blue",
      "Not Significant" = "grey70"
    )
  ) +
  geom_vline(
    xintercept = c(-0.5, 0.5),
    linetype = "dashed"
  ) +
  geom_hline(
    yintercept = -log10(0.05),
    linetype = "dashed"
  ) +
  labs(
    title = "Volcano Plot: CM vs SMA",
    x = "log2 Fold Change",
    y = "-log10 Adjusted P-value",
    color = ""
  ) +
  theme_classic(
    base_size = 14
  )

volcano_plot

ggsave(
  filename = "figures/Volcano_CM_vs_SMA.png",
  plot = volcano_plot,
  width = 8,
  height = 6,
  dpi = 300
)

############################################################
# PART P.2 TOP 50 DIFFERENTIALLY EXPRESSED GENES
############################################################

top50 <- deg_CM_vs_SMA_all %>%
  arrange(adj.P.Val) %>%
  slice(1:50)

mat_heatmap <- ex2[
  top50$PROBEID,
]

gene_label <- ifelse(
  is.na(top50$SYMBOL) |
    top50$SYMBOL == "",
  top50$PROBEID,
  top50$SYMBOL
)

rownames(mat_heatmap) <- gene_label

mat_heatmap <- mat_heatmap[
  rowSums(is.na(mat_heatmap)) == 0,
]

gene_variance <- apply(
  mat_heatmap,
  1,
  var
)

mat_heatmap <- mat_heatmap[
  gene_variance > 0,
]

annotation_col <- data.frame(
  Group = group2
)

rownames(annotation_col) <- colnames(mat_heatmap)

pheatmap(
  mat_heatmap,
  scale = "row",
  annotation_col = annotation_col,
  show_colnames = FALSE,
  show_rownames = TRUE,
  fontsize_row = 8,
  border_color = NA,
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method = "complete",
  main = "Top 50 DEG CM vs SMA"
)

############################################################
# PART Q.1 KONVERSI GENE SYMBOL MENJADI ENTREZ ID
############################################################

gene_list <- unique(
  na.omit(deg_CM_vs_SMA$SYMBOL)
)

gene.df <- bitr(
  gene_list,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)

head(gene.df)

############################################################
# PART Q.2 GENE ONTOLOGY BIOLOGICAL PROCESS (BP)
############################################################

ego_bp <- enrichGO(
  gene = gene.df$ENTREZID,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.20,
  readable = TRUE
)

bp_plot <- dotplot(
  ego_bp,
  showCategory = 15
) +
  ggtitle("GO Biological Process")

bp_plot

ggsave(
  filename = "figures/GO_BP_CM_vs_SMA.png",
  plot = bp_plot,
  width = 8,
  height = 6,
  dpi = 300
)

write.csv(
  as.data.frame(ego_bp),
  "results/GO_BP_CM_vs_SMA.csv",
  row.names = FALSE
)

############################################################
# PART Q.3 GENE ONTOLOGY CELLULAR COMPONENT (CC)
############################################################

ego_cc <- enrichGO(
  gene = gene.df$ENTREZID,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "CC",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.20,
  readable = TRUE
)

cc_plot <- dotplot(
  ego_cc,
  showCategory = 15
) +
  ggtitle("GO Cellular Component")

cc_plot

ggsave(
  filename = "figures/GO_CC_CM_vs_SMA.png",
  plot = cc_plot,
  width = 8,
  height = 6,
  dpi = 300
)

write.csv(
  as.data.frame(ego_cc),
  "results/GO_CC_CM_vs_SMA.csv",
  row.names = FALSE
)

############################################################
# PART Q.4 GENE ONTOLOGY MOLECULAR FUNCTION (MF)
############################################################

ego_mf <- enrichGO(
  gene = gene.df$ENTREZID,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "MF",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.20,
  readable = TRUE
)

mf_plot <- dotplot(
  ego_mf,
  showCategory = 15
) +
  ggtitle("GO Molecular Function")

mf_plot

ggsave(
  filename = "figures/GO_MF_CM_vs_SMA.png",
  plot = mf_plot,
  width = 8,
  height = 6,
  dpi = 300
)

write.csv(
  as.data.frame(ego_mf),
  "results/GO_MF_CM_vs_SMA.csv",
  row.names = FALSE
)

############################################################
# PART R.1 KEGG PATHWAY ENRICHMENT
############################################################

kegg <- enrichKEGG(
  gene = gene.df$ENTREZID,
  organism = "hsa",
  keyType = "ncbi-geneid",
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH"
)

head(as.data.frame(kegg))

kegg_plot <- dotplot(
  kegg,
  showCategory = 15
) +
  ggtitle("KEGG Pathway Enrichment")

kegg_plot

ggsave(
  filename = "figures/KEGG_CM_vs_SMA.png",
  plot = kegg_plot,
  width = 8,
  height = 6,
  dpi = 300
)

write.csv(
  as.data.frame(kegg),
  "results/KEGG_CM_vs_SMA.csv",
  row.names = FALSE
)

############################################################
# PART R.2 MENYIMPAN HASIL ANALISIS
############################################################

write.csv(
  deg_CM_vs_SMA,
  "results/DEG_CM_vs_SMA.csv",
  row.names = FALSE
)

write.csv(
  deg_CM_vs_SMA_all,
  "results/DEG_CM_vs_SMA_ALL.csv",
  row.names = FALSE
)

message(
  "Seluruh analisis transcriptomics GSE117613 (CM vs SMA) telah selesai."
)