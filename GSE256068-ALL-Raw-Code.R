# Load required packages
library(GEOquery)
library(limma)
library(dplyr)
library(ggplot2)
library(VennDiagram)
library(tibble)
library(org.Hs.eg.db)
library(clusterProfiler)
library(enrichplot)

# ==================== Data Download and Preprocessing ====================

# Set working directory and download dataset
setwd("C:\\Users\\xieru\\Desktop\\小凡凡\\GEO")
cat("Downloading GEO dataset...\n")
gset <- getGEO("GSE256068", destdir = ".", getGPL = TRUE)

# Extract expression set and sample information
eset <- if (is.list(gset)) gset[[1]] else gset
sample_info <- pData(eset)

# Process sample grouping
sample_info$group <- ifelse(grepl("Control", sample_info$characteristics_ch1.3), 
                            "Control", "Seizures")

# ==================== Expression Data Processing ====================

# Read TPM data and perform gene annotation
expr_data <- read.table("Norm_counts_TPM.tsv", header = TRUE, row.names = 1, sep = "\t")
gene_ids <- rownames(expr_data)

# Gene annotation using org.Hs.eg.db
gene_anno <- AnnotationDbi::select(org.Hs.eg.db, 
                                   keys = gene_ids, 
                                   columns = c("SYMBOL"), 
                                   keytype = "ENTREZID")

# Remove duplicates and merge with expression data
gene_anno_unique <- gene_anno[!duplicated(gene_anno$ENTREZID), ]
expr_data_anno <- merge(gene_anno_unique, 
                        data.frame(ENTREZID = rownames(expr_data), expr_data), 
                        by = "ENTREZID", all.y = TRUE)

# Filter samples for Control and Seizures groups
sample_mapping <- sample_info[, c("geo_accession", "group")]
target_samples <- sample_mapping[sample_mapping$group %in% c("Control", "Seizures"), ]
available_samples <- intersect(target_samples$geo_accession, colnames(expr_data_anno))

expr_filtered <- expr_data_anno[, c("SYMBOL", available_samples)]
group_vector <- setNames(target_samples$group[match(available_samples, target_samples$geo_accession)], 
                         available_samples)

# ==================== Data Quality Control ====================

# Data preprocessing
expr_matrix <- as.matrix(expr_filtered[, -1])
gene_symbols <- expr_filtered$SYMBOL

# Filter low-expression genes
keep <- rowSums(expr_matrix > 1) >= 3
expr_matrix_filtered <- expr_matrix[keep, ]
gene_symbols_filtered <- gene_symbols[keep]

# Log2 transformation if needed
max_value <- max(expr_matrix_filtered)
if(max_value > 50) {
  expr_log <- log2(expr_matrix_filtered + 1)
  transformation_applied <- TRUE
} else {
  expr_log <- expr_matrix_filtered
  transformation_applied <- FALSE
}

# ==================== Data Visualization ====================

# PCA analysis
expr_for_pca <- t(expr_log)
pca_result <- prcomp(expr_for_pca, scale. = TRUE)
variance_explained <- summary(pca_result)$importance[2, 1:2] * 100

# Prepare PCA data
pca_data <- data.frame(
  Sample = rownames(expr_for_pca),
  PC1 = pca_result$x[, 1],
  PC2 = pca_result$x[, 2],
  Group = group_vector[rownames(expr_for_pca)]
)

# Create PCA plot
pca_plot <- ggplot(pca_data, aes(x = PC1, y = PC2, color = Group)) +
  geom_point(size = 4, alpha = 0.8) +
  scale_color_manual(values = c("Control" = "#2E86AB", "Seizures" = "#A23B72")) +
  labs(
    title = "PCA Analysis of Gene Expression Data",
    x = paste0("PC1 (", round(variance_explained[1], 1), "%)"),
    y = paste0("PC2 (", round(variance_explained[2], 1), "%)")
  ) +
  theme_bw()

print(pca_plot)

# ==================== Differential Expression Analysis ====================

# Create design matrix
group_factor <- factor(group_vector, levels = c("Control", "Seizures"))
design <- model.matrix(~0 + group_factor)
colnames(design) <- c("Control", "Seizures")

# Create contrast matrix
contrast_matrix <- makeContrasts(
  Seizures_vs_Control = Seizures - Control,
  levels = design
)

# Perform limma analysis
fit <- lmFit(expr_log, design)
fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2)

# Get differential expression results
deg_results <- topTable(fit2, coef = "Seizures_vs_Control", number = Inf)
deg_results$ENTREZID <- rownames(deg_results)
deg_results <- merge(deg_results, 
                     data.frame(ENTREZID = rownames(expr_log), 
                                SYMBOL = gene_symbols_filtered),
                     by = "ENTREZID")

# Classify differential genes
logFC_threshold <- 0.58
pvalue_threshold <- 0.05

deg_results$Change <- "Not Significant"
deg_results$Change[deg_results$logFC > logFC_threshold & deg_results$P.Value < pvalue_threshold] <- "Up-regulated"
deg_results$Change[deg_results$logFC < -logFC_threshold & deg_results$P.Value < pvalue_threshold] <- "Down-regulated"

# ==================== Enrichment Analysis ====================

# Get significant genes for enrichment
sig_genes <- deg_results[deg_results$Change != "Not Significant", ]
sig_entrez <- sig_genes$ENTREZID[!is.na(sig_genes$ENTREZID)]

# GO enrichment analysis
go_bp <- enrichGO(
  gene = sig_entrez,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  readable = TRUE
)

# KEGG pathway analysis
kegg <- enrichKEGG(
  gene = sig_entrez,
  organism = 'hsa',
  pvalueCutoff = 0.05
)

if(nrow(kegg) > 0) {
  kegg <- setReadable(kegg, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
}

# Create enrichment plots
if(nrow(go_bp) > 0) {
  p_go <- dotplot(go_bp, showCategory = 15, title = "GO Biological Process")
  print(p_go)
}

if(nrow(kegg) > 0) {
  p_kegg <- dotplot(kegg, showCategory = 15, title = "KEGG Pathway")
  print(p_kegg)
}

# ==================== Gene Expression Boxplots ====================

# Function to create boxplots for specific genes
create_gene_boxplot <- function(gene_symbol, expr_data, sample_groups) {
  
  # Find gene index
  gene_idx <- which(gene_symbols_filtered == gene_symbol)
  if(length(gene_idx) == 0) {
    cat("Warning: Gene", gene_symbol, "not found\n")
    return(NULL)
  }
  
  # Extract gene expression
  gene_expr <- expr_data[gene_idx[1], ]
  
  # Create plot data
  plot_data <- data.frame(
    Sample = names(gene_expr),
    Expression = as.numeric(gene_expr),
    Group = sample_groups[names(gene_expr)]
  )
  
  # Perform t-test
  control_expr <- plot_data$Expression[plot_data$Group == "Control"]
  seizures_expr <- plot_data$Expression[plot_data$Group == "Seizures"]
  
  if(length(control_expr) > 0 & length(seizures_expr) > 0) {
    t_test <- t.test(seizures_expr, control_expr)
    p_value <- t_test$p.value
    p_text <- ifelse(p_value < 0.001, "p < 0.001", 
                     paste0("p = ", sprintf("%.3f", p_value)))
  } else {
    p_text <- "p = NA"
  }
  
  # Create boxplot
  p <- ggplot(plot_data, aes(x = Group, y = Expression, fill = Group)) +
    geom_boxplot(alpha = 0.7, outlier.shape = NA) +
    geom_jitter(width = 0.2, size = 2.5, alpha = 0.8) +
    scale_fill_manual(values = c("Control" = "#2E86AB", "Seizures" = "#A23B72")) +
    annotate("text", x = 1.5, y = max(plot_data$Expression) * 1.05, 
             label = p_text, size = 4, fontface = "bold") +
    labs(
      title = paste0("Expression of ", gene_symbol),
      x = "Group",
      y = if(transformation_applied) "Log2(TPM + 1)" else "TPM"
    ) +
    theme_bw() +
    theme(legend.position = "none")
  
  return(p)
}

# Plot specific genes
target_genes <- c("GPX4", "TNF", "IL6", "IL1B", "SLC7A11")
for(gene in target_genes) {
  p <- create_gene_boxplot(gene, expr_log, group_vector)
  if(!is.null(p)) print(p)
}

# ==================== Save Results ====================

# Save expression matrix
write.csv(deg_results, "differential_expression_results.csv", row.names = FALSE)

# Save plots
ggsave("PCA_analysis.png", pca_plot, width = 10, height = 8, dpi = 300)
if(exists("p_go") && !is.null(p_go)) {
  ggsave("GO_enrichment.png", p_go, width = 12, height = 8, dpi = 300)
}
if(exists("p_kegg") && !is.null(p_kegg)) {
  ggsave("KEGG_enrichment.png", p_kegg, width = 12, height = 8, dpi = 300)
}

cat("Analysis completed successfully!\n")
