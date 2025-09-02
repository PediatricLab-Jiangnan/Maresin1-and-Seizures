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

cat("Dataset basic information:\n")
cat("- Sample count:", nrow(sample_info), "\n")
cat("- Available columns:", paste(colnames(sample_info)[1:5], collapse = ", "), "...\n")

# Process sample grouping
sample_info$group <- ifelse(grepl("Control", sample_info$characteristics_ch1.3), 
                            "Control", "Seizures")

group_counts <- table(sample_info$group)
cat("Group statistics:\n")
print(group_counts)

# ==================== Expression Data Processing ====================

# Read TPM data and perform gene annotation
cat("\nReading TPM expression data...\n")
expr_data <- read.table("Norm_counts_TPM.tsv", header = TRUE, row.names = 1, sep = "\t")
cat("Original data dimensions:", nrow(expr_data), "genes ×", ncol(expr_data), "samples\n")

# Gene annotation
cat("Performing gene annotation...\n")
gene_ids <- rownames(expr_data)
cat("Gene ID examples:", paste(gene_ids[1:5], collapse = ", "), "\n")

# Annotate genes using org.Hs.eg.db
gene_anno <- AnnotationDbi::select(org.Hs.eg.db, 
                                   keys = gene_ids, 
                                   columns = c("SYMBOL"), 
                                   keytype = "ENTREZID")

# Remove duplicates and merge with expression data
gene_anno_unique <- gene_anno[!duplicated(gene_anno$ENTREZID), ]
expr_data_anno <- expr_data
expr_data_anno$ENTREZID <- rownames(expr_data_anno)

# Add gene symbols
expr_data_anno <- merge(gene_anno_unique, expr_data_anno, 
                        by.x = "ENTREZID", by.y = "ENTREZID", 
                        all.y = TRUE, sort = FALSE)
rownames(expr_data_anno) <- expr_data_anno$ENTREZID

# Get sample columns (excluding ENTREZID and SYMBOL)
sample_columns <- setdiff(colnames(expr_data_anno), c("ENTREZID", "SYMBOL"))

# Filter samples for Control and Seizures groups
cat("\nFiltering target samples...\n")
sample_mapping <- sample_info[, c("geo_accession", "group")]
control_seizures_samples <- sample_mapping[sample_mapping$group %in% c("Control", "Seizures"), ]

# Check sample matching
available_samples <- intersect(control_seizures_samples$geo_accession, sample_columns)
cat("Successfully matched samples:", length(available_samples), "\n")

# Filter expression data
expr_filtered <- expr_data_anno[, c("SYMBOL", available_samples)]
final_groups <- control_seizures_samples[control_seizures_samples$geo_accession %in% available_samples, ]
group_vector <- setNames(final_groups$group, final_groups$geo_accession)

cat("Final sample grouping:\n")
print(table(group_vector))

# ==================== Data Quality Control ====================

# Data preprocessing
cat("\nData preprocessing...\n")
expr_matrix <- as.matrix(expr_filtered[, -1])
gene_symbols <- expr_filtered$SYMBOL

# Filter low-expression genes
keep <- rowSums(expr_matrix > 1) >= 3
expr_matrix_filtered <- expr_matrix[keep, ]
gene_symbols_filtered <- gene_symbols[keep]
cat("Genes retained after filtering:", nrow(expr_matrix_filtered), "\n")

# Check data distribution and log transformation
cat("\nChecking data distribution...\n")
max_value <- max(expr_matrix_filtered)
cat("Data maximum value:", round(max_value, 2), "\n")

# Log2 transformation if needed
if(max_value > 50) {
  cat("Performing Log2 transformation...\n")
  expr_log <- log2(expr_matrix_filtered + 1)
  transformation_applied <- TRUE
} else {
  cat("Using original data...\n")
  expr_log <- expr_matrix_filtered
  transformation_applied <- FALSE
}

# Ensure sample order matching
sample_groups <- group_vector[colnames(expr_log)]

# ==================== Data Visualization ====================

# Data distribution visualization
cat("\nGenerating data distribution plots...\n")
set.seed(123)
sample_genes <- sample(nrow(expr_matrix_filtered), min(1000, nrow(expr_matrix_filtered)))

# Create distribution plots
par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))

# Original data distribution
sample_data <- as.vector(expr_matrix_filtered[sample_genes, ])
hist(sample_data, breaks = 50, 
     main = "Original TPM Data Distribution", 
     xlab = "TPM Values", ylab = "Frequency",
     col = "lightblue", border = "white")

# Log-transformed distribution
hist(log2(sample_data + 1), breaks = 50, 
     main = "Log2 Transformed Distribution", 
     xlab = "Log2(TPM+1)", ylab = "Frequency",
     col = "lightcoral", border = "white")

# PCA analysis
cat("\nPCA analysis...\n")
expr_for_pca <- t(expr_log)
pca_result <- prcomp(expr_for_pca, scale. = TRUE)
variance_explained <- summary(pca_result)$importance[2, 1:2] * 100

cat("PC1 explained variance:", round(variance_explained[1], 2), "%\n")
cat("PC2 explained variance:", round(variance_explained[2], 2), "%\n")

# Prepare PCA data
pca_data <- data.frame(
  Sample = rownames(expr_for_pca),
  PC1 = pca_result$x[, 1],
  PC2 = pca_result$x[, 2],
  Group = sample_groups,
  stringsAsFactors = FALSE
)

# Create PCA plot
pca_plot <- ggplot(pca_data, aes(x = PC1, y = PC2, color = Group, fill = Group)) +
  geom_point(size = 4, alpha = 0.8, stroke = 0.5) +
  labs(
    title = "PCA Analysis of Gene Expression Data",
    subtitle = paste0("Total samples: ", nrow(pca_data), 
                      " | Genes analyzed: ", nrow(expr_log)),
    x = paste0("PC1 (", round(variance_explained[1], 1), "%)"),
    y = paste0("PC2 (", round(variance_explained[2], 1), "%)"),
    color = "Group",
    fill = "Group"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 12, color = "gray60"),
    axis.title = element_text(size = 12, face = "bold"),
    axis.text = element_text(size = 10),
    legend.position = "bottom",
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "gray90", size = 0.3)
  ) +
  scale_color_manual(values = c("Control" = "#2E86AB", "Seizures" = "#A23B72")) +
  scale_fill_manual(values = c("Control" = "#2E86AB", "Seizures" = "#A23B72"))

print(pca_plot)

# ==================== Differential Expression Analysis ====================

# Create design matrix
cat("\nStarting differential expression analysis...\n")
group_factor <- factor(sample_groups, levels = c("Control", "Seizures"))
design <- model.matrix(~0 + group_factor)
colnames(design) <- c("Control", "Seizures")

# Create contrast matrix
contrast_matrix <- makeContrasts(
  Seizures_vs_Control = Seizures - Control,
  levels = design
)

cat("Contrast matrix:\n")
print(contrast_matrix)

# Perform limma analysis
fit <- lmFit(expr_log, design)
fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2)

# Get differential expression results
deg_results <- topTable(fit2, coef = "Seizures_vs_Control", number = Inf, sort.by = "P")

# Add gene symbols
deg_results$ENTREZID <- rownames(deg_results)
deg_results <- merge(deg_results, 
                     data.frame(ENTREZID = rownames(expr_log), 
                                SYMBOL = gene_symbols_filtered),
                     by = "ENTREZID", all.x = TRUE)

# Remove genes without symbols
deg_results <- deg_results[!is.na(deg_results$SYMBOL), ]
deg_results <- deg_results[order(deg_results$P.Value), ]

cat("Differential expression analysis completed\n")
cat("Total genes:", nrow(deg_results), "\n")

# ==================== Differential Gene Classification ====================

# Set thresholds
logFC_threshold <- 0.58      # |log2FC| >= 0.58
pvalue_threshold <- 0.05     # P-value < 0.05

# Add differential expression classification
deg_results$Change <- "Not Significant"
deg_results$Change[deg_results$logFC > logFC_threshold & deg_results$P.Value < pvalue_threshold] <- "Up-regulated"
deg_results$Change[deg_results$logFC < -logFC_threshold & deg_results$P.Value < pvalue_threshold] <- "Down-regulated"

# Statistics
cat("\n=== Differential Gene Statistics (P < 0.05, |logFC| >= 0.58) ===\n")
change_summary <- table(deg_results$Change)
print(change_summary)

# Show top differential genes
cat("\n=== Top 10 Up-regulated Genes ===\n")
top_up <- deg_results[deg_results$Change == "Up-regulated", ][1:10, c("SYMBOL", "logFC", "P.Value", "adj.P.Val")]
print(top_up)

cat("\n=== Top 10 Down-regulated Genes ===\n")
top_down <- deg_results[deg_results$Change == "Down-regulated", ][1:10, c("SYMBOL", "logFC", "P.Value", "adj.P.Val")]
print(top_down)

# ==================== Enrichment Analysis ====================

# Get significant genes for enrichment
sig_genes <- deg_results[deg_results$Change != "Not Significant", ]
sig_entrez <- sig_genes$ENTREZID[!is.na(sig_genes$ENTREZID)]

cat("=== GO Enrichment Analysis ===\n")

# GO enrichment analysis - Biological Process
go_bp <- enrichGO(
  gene = sig_entrez,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2,
  readable = TRUE
)

# GO enrichment analysis - Molecular Function
go_mf <- enrichGO(
  gene = sig_entrez,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "MF",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2,
  readable = TRUE
)

# GO enrichment analysis - Cellular Component
go_cc <- enrichGO(
  gene = sig_entrez,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "CC",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2,
  readable = TRUE
)

# KEGG pathway analysis
cat("=== KEGG Enrichment Analysis ===\n")
kegg <- enrichKEGG(
  gene = sig_entrez,
  organism = 'hsa',
  keyType = "kegg",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2
)

# Convert to gene symbols
if(nrow(kegg) > 0) {
  kegg <- setReadable(kegg, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
}

# Create enrichment plots
cat("=== Creating Enrichment Plots ===\n")

# GO BP dotplot
if(nrow(go_bp) > 0) {
  p_go_bp <- dotplot(go_bp, showCategory = 10, title = "GO Biological Process") +
    theme(axis.text.y = element_text(size = 10))
  print(p_go_bp)
}

# GO MF dotplot
if(nrow(go_mf) > 0) {
  p_go_mf <- dotplot(go_mf, showCategory = 10, title = "GO Molecular Function") +
    theme(axis.text.y = element_text(size = 10))
  print(p_go_mf)
}

# GO CC dotplot
if(nrow(go_cc) > 0) {
  p_go_cc <- dotplot(go_cc, showCategory = 10, title = "GO Cellular Component") +
    theme(axis.text.y = element_text(size = 10))
  print(p_go_cc)
}

# KEGG dotplot
if(nrow(kegg) > 0) {
  p_kegg <- dotplot(kegg, showCategory = 15, title = "KEGG Pathway Enrichment") +
    theme(axis.text.y = element_text(size = 10))
  print(p_kegg)
} else {
  cat("No significant KEGG pathway enrichment results\n")
}

# ==================== Gene Expression Boxplots ====================

# Function to create boxplots for specific genes
create_gene_boxplot <- function(gene_symbol, expr_data, sample_groups, title_suffix = "") {
  
  # Check if gene exists
  gene_idx <- which(gene_symbols_filtered == gene_symbol)
  if(length(gene_idx) == 0) {
    cat("Warning: Gene", gene_symbol, "not found in data\n")
    return(NULL)
  }
  
  # Extract gene expression data
  gene_expr <- expr_data[gene_idx[1], ]  # Take first if duplicated
  
  # Create plot data frame
  plot_data <- data.frame(
    Sample = names(gene_expr),
    Expression = as.numeric(gene_expr),
    Group = sample_groups[names(gene_expr)],
    stringsAsFactors = FALSE
  )
  
  # Remove NA values
  plot_data <- plot_data[complete.cases(plot_data), ]
  
  # Perform t-test
  control_expr <- plot_data$Expression[plot_data$Group == "Control"]
  seizures_expr <- plot_data$Expression[plot_data$Group == "Seizures"]
  
  if(length(control_expr) > 0 & length(seizures_expr) > 0) {
    t_test_result <- t.test(seizures_expr, control_expr)
    p_value <- t_test_result$p.value
    
    # Format p-value
    if(p_value < 0.001) {
      p_text <- "p < 0.001"
    } else if(p_value < 0.01) {
      p_text <- paste0("p = ", sprintf("%.3f", p_value))
    } else {
      p_text <- paste0("p = ", sprintf("%.3f", p_value))
    }
    
    # Calculate fold change
    mean_control <- mean(control_expr, na.rm = TRUE)
    mean_seizures <- mean(seizures_expr, na.rm = TRUE)
    fold_change <- mean_seizures - mean_control  # Log-transformed data
    
  } else {
    p_text <- "p = NA"
    fold_change <- NA
  }
  
  # Calculate statistics for display
  stats_summary <- plot_data %>%
    group_by(Group) %>%
    summarise(
      n = n(),
      mean = mean(Expression, na.rm = TRUE),
      sd = sd(Expression, na.rm = TRUE),
      .groups = 'drop'
    )
  
  cat("Gene", gene_symbol, "statistics:\n")
  print(stats_summary)
  cat("Fold Change (log2):", round(fold_change, 3), "\n")
  cat("P-value:", p_text, "\n\n")
  
  # Create boxplot
  p <- ggplot(plot_data, aes(x = Group, y = Expression, fill = Group)) +
    # Add boxplot
    geom_boxplot(alpha = 0.7, outlier.shape = NA, width = 0.6) +
    # Add scatter points
    geom_jitter(width = 0.2, size = 2.5, alpha = 0.8, 
                aes(color = Group), stroke = 0.3) +
    
    # Set colors
    scale_fill_manual(values = c("Control" = "#2E86AB", "Seizures" = "#A23B72")) +
    scale_color_manual(values = c("Control" = "#1B5E7A", "Seizures" = "#7A1E4F")) +
    
    # Add statistical significance annotation
    annotate("text", x = 1.5, y = max(plot_data$Expression) * 1.05, 
             label = p_text, size = 4, fontface = "bold") +
    
    # Set titles and labels
    labs(
      title = paste0("Expression of ", gene_symbol, title_suffix),
      subtitle = paste0("n(Control) = ", sum(plot_data$Group == "Control"), 
                        ", n(Seizures) = ", sum(plot_data$Group == "Seizures")),
      x = "Group",
      y = if(transformation_applied) "Log2(TPM + 1)" else "TPM",
      fill = "Group",
      color = "Group"
    ) +
    
    # Beautify theme
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 11, color = "gray60"),
      axis.title = element_text(size = 12, face = "bold"),
      axis.text = element_text(size = 11),
      legend.position = "none",  # Hide legend since x-axis shows grouping
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.major.y = element_line(color = "gray90", size = 0.3),
      plot.margin = margin(15, 15, 15, 15)
    ) +
    
    # Set y-axis range
    coord_cartesian(ylim = c(min(plot_data$Expression) * 0.95, 
                             max(plot_data$Expression) * 1.15))
  
  return(p)
}

# Plot specific genes
cat("=== Creating Gene Expression Boxplots ===\n")

# Define target genes
target_genes <- c("GPX4", "TNF", "IL6", "IL1B", "SLC7A11", "FTH1", "NFE2L2", "NFKB1", "NFKB2")
genes_to_plot <- target_genes

cat("Genes to plot:", paste(genes_to_plot, collapse = ", "), "\n")

# Create individual gene plots
gene_plots <- list()
valid_genes <- c()

for(gene in genes_to_plot) {
  plot_obj <- create_gene_boxplot(gene, expr_log, sample_groups)
  if(!is.null(plot_obj)) {
    gene_plots[[gene]] <- plot_obj
    valid_genes <- c(valid_genes, gene)
    print(plot_obj)
    
    # Save individual gene plot
    ggsave(paste0("boxplot_", gene, ".png"), plot_obj, 
           width = 8, height = 6, dpi = 300)
    ggsave(paste0("boxplot_", gene, ".pdf"), plot_obj, 
           width = 8, height = 6)
  }
}

# ==================== Save Results ====================

cat("=== Saving Results ===\n")

# Save expression matrix
expr_output <- data.frame(
  ENTREZID = rownames(expr_log),
  SYMBOL = gene_symbols_filtered,
  expr_log,
  stringsAsFactors = FALSE
)

output_filename <- paste0("expression_matrix_", 
                          if(transformation_applied) "log2_transformed" else "original",
                          ".csv")
write.csv(expr_output, file = output_filename, row.names = FALSE)

# Save differential expression results
write.csv(deg_results, "differential_expression_results.csv", row.names = FALSE)

# Save PCA plot
ggsave("PCA_analysis.png", plot = pca_plot, width = 10, height = 8, dpi = 300)
ggsave("PCA_analysis.pdf", plot = pca_plot, width = 10, height = 8)

# Save sample information
sample_summary <- data.frame(
  Sample_ID = names(sample_groups),
  Group = sample_groups,
  stringsAsFactors = FALSE
)
write.csv(sample_summary, file = "sample_information.csv", row.names = FALSE)

# Save enrichment analysis results
if(nrow(go_bp) > 0) {
  write.csv(as.data.frame(go_bp), "GO_BP_enrichment.csv", row.names = FALSE)
  ggsave("GO_BP_dotplot.png", p_go_bp, width = 12, height = 8, dpi = 300)
  ggsave("GO_BP_dotplot.pdf", p_go_bp, width = 12, height = 8)
}

if(nrow(go_mf) > 0) {
  write.csv(as.data.frame(go_mf), "GO_MF_enrichment.csv", row.names = FALSE)
  ggsave("GO_MF_dotplot.png", p_go_mf, width = 12, height = 8, dpi = 300)
  ggsave("GO_MF_dotplot.pdf", p_go_mf, width = 12, height = 8)
}

if(nrow(go_cc) > 0) {
  write.csv(as.data.frame(go_cc), "GO_CC_enrichment.csv", row.names = FALSE)
}

if(nrow(kegg) > 0) {
  write.csv(as.data.frame(kegg), "KEGG_enrichment.csv", row.names = FALSE)
  ggsave("KEGG_dotplot.png", p_kegg, width = 12, height = 8, dpi = 300)
  ggsave("KEGG_dotplot.pdf", p_kegg, width = 12, height = 8)
}

# Enrichment analysis summary
cat("=== Enrichment Analysis Summary ===\n")
cat("GO BP significant terms:", nrow(go_bp), "\n")
cat("GO MF significant terms:", nrow(go_mf), "\n")
cat("GO CC significant terms:", nrow(go_cc), "\n")
cat("KEGG significant pathways:", nrow(kegg), "\n")

# Show top results
if(nrow(go_bp) > 0) {
  cat("\n=== Top 5 GO BP Terms ===\n")
  print(as.data.frame(go_bp)[1:min(5, nrow(go_bp)), c("Description", "Count", "p.adjust")])
}

if(nrow(kegg) > 0) {
  cat("\n=== Top 5 KEGG Pathways ===\n")
  print(as.data.frame(kegg)[1:min(5, nrow(kegg)), c("Description", "Count", "p.adjust")])
}

cat("\nAnalysis completed! All results and plots have been saved.\n")

