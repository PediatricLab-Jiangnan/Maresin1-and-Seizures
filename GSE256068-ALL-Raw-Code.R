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
library(pheatmap)
library(RColorBrewer)
library(gridExtra)
library(grid)
library(ggrepel)
library(MASS)

# ==================== Data Acquisition and Preprocessing ====================

# Step 1: Set working directory and download dataset
setwd("C:\\Users\\xieru\\Desktop\\姜小凡新一轮\\GEO")
cat("Downloading GEO dataset...\n")
gset <- getGEO("GSE256068", destdir = ".", getGPL = TRUE)

# Step 2: Extract expression set object and sample information
cat("Extracting dataset information...\n")
eset <- if (is.list(gset)) gset[[1]] else gset
sample_info <- pData(eset)

cat("Basic dataset information:\n")
cat("- Number of samples:", nrow(sample_info), "\n")
cat("- Available columns:", paste(colnames(sample_info)[1:min(5, ncol(sample_info))], collapse = ", "), "...\n")

# Step 3: Process sample grouping information
cat("\nProcessing sample groups...\n")
sample_info$group <- ifelse(grepl("Control", sample_info$characteristics_ch1.3), 
                            "Control", 
                            ifelse(grepl("TLE-HS", sample_info$characteristics_ch1.3), 
                                   "TLE", "Unknown"))

group_counts <- table(sample_info$group)
cat("Group statistics:\n")
print(group_counts)

# ==================== Expression Data Processing ====================

# Step 4: Read TPM data and perform gene annotation
cat("\nReading TPM expression data...\n")
expr_data <- read.table("Norm_counts_TPM.tsv", header = TRUE, row.names = 1, sep = "\t")
cat("Original data dimensions:", nrow(expr_data), "genes ×", ncol(expr_data), "samples\n")

# Gene annotation
cat("Performing gene annotation...\n")
gene_ids <- rownames(expr_data)
gene_anno <- AnnotationDbi::select(org.Hs.eg.db, 
                                   keys = gene_ids, 
                                   columns = c("SYMBOL"), 
                                   keytype = "ENTREZID")
gene_anno <- gene_anno[!duplicated(gene_anno$ENTREZID), ]

# Merge annotation information
expr_data_anno <- merge(gene_anno, expr_data, by.x = "ENTREZID", by.y = "row.names", all.y = TRUE)
rownames(expr_data_anno) <- expr_data_anno$ENTREZID
expr_data_anno <- expr_data_anno[, -1]

# Step 5: Filter Control and TLE group samples
cat("\nFiltering target samples...\n")
sample_mapping <- sample_info[, c("geo_accession", "group")]
control_tle_samples <- sample_mapping[sample_mapping$group %in% c("Control", "TLE"), ]
available_samples <- intersect(control_tle_samples$geo_accession, colnames(expr_data_anno))

# Filter expression data
expr_filtered <- expr_data_anno[, c("SYMBOL", available_samples)]
cat("Filtered data dimensions:", nrow(expr_filtered), "genes ×", length(available_samples), "samples\n")

# Create group vector
final_groups <- control_tle_samples[control_tle_samples$geo_accession %in% available_samples, ]
group_vector <- setNames(final_groups$group, final_groups$geo_accession)

cat("Final sample grouping:\n")
print(table(group_vector))

# ==================== Data Quality Control and Transformation ====================

# Step 6: Data preprocessing
cat("\nData preprocessing...\n")
expr_matrix <- as.matrix(expr_filtered[, -1])
gene_symbols <- expr_filtered$SYMBOL

# Filter low expression genes (keep genes expressed in at least 3 samples)
keep <- rowSums(expr_matrix > 1) >= 3
expr_matrix_filtered <- expr_matrix[keep, ]
gene_symbols_filtered <- gene_symbols[keep]
cat("After filtering, retained", nrow(expr_matrix_filtered), "genes\n")

# Step 7: Data distribution check and Log transformation
cat("\nChecking data distribution...\n")
max_value <- max(expr_matrix_filtered, na.rm = TRUE)
cat("Maximum data value:", round(max_value, 2), "\n")

# Determine if log transformation is needed
if(max_value > 50) {
  cat("Performing Log2 transformation...\n")
  expr_log <- log2(expr_matrix_filtered + 1)
  transformation_applied <- TRUE
} else {
  cat("Using original data...\n")
  expr_log <- expr_matrix_filtered
  transformation_applied <- FALSE
}

# Ensure sample order matches
sample_groups <- group_vector[colnames(expr_log)]

# ==================== Data Visualization ====================

# Step 8: Improved data distribution visualization
cat("\nGenerating data distribution plots...\n")
set.seed(123)
sample_genes <- sample(nrow(expr_matrix_filtered), min(1000, nrow(expr_matrix_filtered)))

# Create more aesthetic distribution plots
par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))

# Original data distribution
sample_data <- as.vector(expr_matrix_filtered[sample_genes, ])
hist(sample_data, breaks = 50, 
     main = "Original TPM Data Distribution", 
     xlab = "TPM Values", ylab = "Frequency",
     col = "lightblue", border = "white")

# Log-transformed distribution
hist(log2(sample_data + 1), breaks = 50, 
     main = "Log2-transformed Distribution", 
     xlab = "Log2(TPM+1)", ylab = "Frequency",
     col = "lightcoral", border = "white")

# Step 9: Optimized PCA analysis
cat("\nPCA analysis...\n")

# Perform PCA
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

# Step 10: Draw optimized PCA plot with confidence ellipses
cat("Drawing PCA plot...\n")

# Function to calculate confidence ellipse data
stat_ellipse_data <- function(data, level = 0.95) {
  ellipse_data <- data %>%
    group_by(Group) %>%
    do({
      if(nrow(.) >= 3) {  # Need at least 3 points to draw ellipse
        # Calculate covariance matrix
        cov_matrix <- cov(cbind(.$PC1, .$PC2))
        center <- c(mean(.$PC1), mean(.$PC2))
        
        # Calculate ellipse parameters
        eigenvals <- eigen(cov_matrix)$values
        eigenvecs <- eigen(cov_matrix)$vectors
        
        # Calculate confidence interval
        chisq_val <- qchisq(level, df = 2)
        
        # Generate ellipse points
        angles <- seq(0, 2*pi, length.out = 100)
        ellipse_points <- matrix(0, nrow = 100, ncol = 2)
        
        for(i in 1:100) {
          point <- c(cos(angles[i]), sin(angles[i]))
          ellipse_points[i, ] <- center + sqrt(chisq_val) * 
            (eigenvecs %*% diag(sqrt(eigenvals)) %*% point)
        }
        
        data.frame(
          PC1 = ellipse_points[, 1],
          PC2 = ellipse_points[, 2]
        )
      } else {
        data.frame(PC1 = numeric(0), PC2 = numeric(0))
      }
    })
  
  return(ellipse_data)
}

# Calculate ellipse data
ellipse_data <- stat_ellipse_data(pca_data, level = 0.95)

# Draw improved PCA plot
pca_plot <- ggplot(pca_data, aes(x = PC1, y = PC2, color = Group, fill = Group)) +
  # Add confidence ellipses
  geom_polygon(data = ellipse_data, 
               aes(x = PC1, y = PC2, fill = Group), 
               alpha = 0.2, color = NA) +
  # Add ellipse borders
  geom_path(data = ellipse_data, 
            aes(x = PC1, y = PC2, color = Group), 
            size = 1, alpha = 0.8) +
  # Add data points
  geom_point(size = 4, alpha = 0.8, stroke = 0.5) +
  # Set titles and labels
  labs(
    title = "PCA Analysis of Gene Expression Data",
    subtitle = paste0("Total samples: ", nrow(pca_data), 
                      " | Genes analyzed: ", nrow(expr_log)),
    x = paste0("PC1 (", round(variance_explained[1], 1), "%)"),
    y = paste0("PC2 (", round(variance_explained[2], 1), "%)"),
    color = "Group",
    fill = "Group"
  ) +
  # Beautify theme
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
    panel.grid.major = element_line(color = "gray90", size = 0.3),
    plot.margin = margin(20, 20, 20, 20)
  ) +
  # Set colors
  scale_color_manual(values = c("Control" = "#2E86AB", "TLE" = "#A23B72")) +
  scale_fill_manual(values = c("Control" = "#2E86AB", "TLE" = "#A23B72")) +
  # Add guides
  guides(
    color = guide_legend(override.aes = list(size = 5)),
    fill = guide_legend(override.aes = list(alpha = 0.3))
  )

# Display PCA plot
print(pca_plot)

# ==================== Data Output ====================

# Step 11: Output processed data
cat("\nOutputting processed data...\n")

# Create output matrix
expr_output <- data.frame(
  ENTREZID = rownames(expr_log),
  SYMBOL = gene_symbols_filtered,
  expr_log,
  stringsAsFactors = FALSE
)

# Save expression matrix
output_filename <- paste0("expression_matrix_", 
                          if(transformation_applied) "log2_transformed" else "original",
                          ".csv")
write.csv(expr_output, file = output_filename, row.names = FALSE)

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

# ==================== Differential Expression Analysis ====================

# Step 12: Differential expression analysis
cat("\nStarting differential expression analysis...\n")

# 1. Create group factor
group_factor <- factor(sample_groups, levels = c("Control", "TLE"))

# 2. Build design matrix
design <- model.matrix(~0 + group_factor)
colnames(design) <- c("Control", "TLE")

# 3. Create contrast matrix
contrast_matrix <- makeContrasts(
  TLE_vs_Control = TLE - Control,
  levels = design
)

cat("Contrast matrix:\n")
print(contrast_matrix)

# 4. Use limma for differential analysis
fit <- lmFit(expr_log, design)
fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2)

# 5. Get differential expression results
deg_results <- topTable(fit2, coef = "TLE_vs_Control", number = Inf, sort.by = "P")

# 6. Add gene symbols
deg_results$ENTREZID <- rownames(deg_results)
deg_results <- merge(deg_results, 
                     data.frame(ENTREZID = rownames(expr_log), 
                                SYMBOL = gene_symbols_filtered),
                     by = "ENTREZID", all.x = TRUE)

# Remove genes without symbols
deg_results <- deg_results[!is.na(deg_results$SYMBOL), ]
deg_results <- deg_results[order(deg_results$P.Value), ]

cat("Differential expression analysis completed\n")
cat("Total number of genes:", nrow(deg_results), "\n")

# ==================== Differential Gene Statistics ====================

# Step 13: Differential gene classification and statistics
cat("\nDifferential gene statistics...\n")

# Set thresholds
logFC_threshold <- 0.58      # |log2FC| >= 0.58 (approximately 1.5-fold change)
pvalue_threshold <- 0.05     # P-value < 0.05
padj_threshold <- 0.05       # Adjusted P-value < 0.05

# Add differential expression classification
deg_results$Change <- "Not Significant"
deg_results$Change[deg_results$logFC > logFC_threshold & deg_results$P.Value < pvalue_threshold] <- "Up-regulated"
deg_results$Change[deg_results$logFC < -logFC_threshold & deg_results$P.Value < pvalue_threshold] <- "Down-regulated"

# Add strict classification based on adjusted p-value
deg_results$Change_strict <- "Not Significant"
deg_results$Change_strict[deg_results$logFC > logFC_threshold & deg_results$adj.P.Val < padj_threshold] <- "Up-regulated"
deg_results$Change_strict[deg_results$logFC < -logFC_threshold & deg_results$adj.P.Val < padj_threshold] <- "Down-regulated"

# Statistics results
cat("=== Differential Gene Statistics (P-value < 0.05, |logFC| >= 0.58) ===\n")
change_summary <- table(deg_results$Change)
print(change_summary)

cat("\n=== Differential Gene Statistics (Adjusted P-value < 0.05, |logFC| >= 0.58) ===\n")
change_strict_summary <- table(deg_results$Change_strict)
print(change_strict_summary)

# Display top differential genes
cat("\n=== Top 10 Up-regulated Genes ===\n")
top_up <- deg_results[deg_results$Change == "Up-regulated", ][1:min(10, sum(deg_results$Change == "Up-regulated")), 
                                                              c("SYMBOL", "logFC", "P.Value", "adj.P.Val")]
print(top_up)

cat("\n=== Top 10 Down-regulated Genes ===\n")
top_down <- deg_results[deg_results$Change == "Down-regulated", ][1:min(10, sum(deg_results$Change == "Down-regulated")), 
                                                                  c("SYMBOL", "logFC", "P.Value", "adj.P.Val")]
print(top_down)

# Save differential expression results
write.csv(deg_results, "differential_expression_results.csv", row.names = FALSE)

# ==================== Enrichment Analysis ====================

# Extract significant genes for enrichment analysis
sig_genes <- deg_results[deg_results$Change != "Not Significant", ]
sig_entrez <- sig_genes$ENTREZID[!is.na(sig_genes$ENTREZID)]

cat("=== Enrichment Analysis ===\n")
cat("Number of significant genes for enrichment:", length(sig_entrez), "\n")

# Check if there are enough genes for enrichment analysis
if(length(sig_entrez) < 5) {
  cat("Warning: Too few significant genes for enrichment analysis. Using relaxed criteria.\n")
  # Use more relaxed criteria
  relaxed_genes <- deg_results[abs(deg_results$logFC) > 0.3 & deg_results$P.Value < 0.1, ]
  sig_entrez <- relaxed_genes$ENTREZID[!is.na(relaxed_genes$ENTREZID)]
  cat("Using", length(sig_entrez), "genes with relaxed criteria\n")
}

# Proceed only if we have enough genes
if(length(sig_entrez) >= 5) {
  
  # Step 1: GO enrichment analysis
  cat("=== GO Enrichment Analysis ===\n")
  
  # Biological Process
  go_bp <- enrichGO(
    gene          = sig_entrez,
    OrgDb         = org.Hs.eg.db,
    keyType       = "ENTREZID",
    ont           = "BP",           # Biological Process
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    qvalueCutoff  = 0.2,
    readable      = TRUE
  )
  
  # Molecular Function
  go_mf <- enrichGO(
    gene          = sig_entrez,
    OrgDb         = org.Hs.eg.db,
    keyType       = "ENTREZID",
    ont           = "MF",           # Molecular Function
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    qvalueCutoff  = 0.2,
    readable      = TRUE
  )
  
  # Cellular Component
  go_cc <- enrichGO(
    gene          = sig_entrez,
    OrgDb         = org.Hs.eg.db,
    keyType       = "ENTREZID",
    ont           = "CC",           # Cellular Component
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    qvalueCutoff  = 0.2,
    readable      = TRUE
  )
  
  # All GO categories combined
  go_all <- enrichGO(
    gene          = sig_entrez,
    OrgDb         = org.Hs.eg.db,
    keyType       = "ENTREZID",
    ont           = "ALL",          # All categories
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    qvalueCutoff  = 0.2,
    readable      = TRUE
  )
  
  # Step 2: KEGG pathway enrichment analysis
  cat("=== KEGG Pathway Enrichment Analysis ===\n")
  
  kegg <- enrichKEGG(
    gene         = sig_entrez,
    organism     = 'hsa',
    keyType      = "kegg",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.2
  )
  
  # Convert to gene symbols
  if(nrow(kegg) > 0) {
    kegg <- setReadable(kegg, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
  }
  
  # Step 3: Custom dotplot function for better visualization
  create_custom_dotplot <- function(enrich_result, title, top_n = 15) {
    if(nrow(enrich_result) == 0) {
      cat("No significant terms found for", title, "\n")
      return(NULL)
    }
    
    # Prepare data
    plot_data <- as.data.frame(enrich_result)
    plot_data <- plot_data[1:min(top_n, nrow(plot_data)), ]
    
    # Calculate gene ratio
    plot_data$GeneRatio_num <- sapply(plot_data$GeneRatio, function(x) {
      nums <- as.numeric(unlist(strsplit(x, "/")))
      return(nums[1]/nums[2])
    })
    
    # Truncate long descriptions
    plot_data$Description_short <- ifelse(
      nchar(plot_data$Description) > 50,
      paste0(substr(plot_data$Description, 1, 47), "..."),
      plot_data$Description
    )
    
    # Sort by p.adjust
    plot_data <- plot_data[order(plot_data$p.adjust), ]
    plot_data$Description_short <- factor(plot_data$Description_short, 
                                          levels = rev(plot_data$Description_short))
    
    # Create dotplot
    p <- ggplot(plot_data, aes(x = GeneRatio_num, y = Description_short)) +
      geom_point(aes(size = Count, color = p.adjust), alpha = 0.8) +
      scale_color_gradient(low = "#E31A1C", high = "#1F78B4", 
                           name = "Adjusted\nP-value",
                           trans = "log10",
                           labels = scales::trans_format("log10", scales::math_format(10^.x))) +
      scale_size_continuous(name = "Gene\nCount", range = c(3, 10)) +
      labs(
        title = title,
        x = "Gene Ratio",
        y = ""
      ) +
      theme_bw() +
      theme(
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10),
        axis.title.x = element_text(size = 12, face = "bold"),
        legend.title = element_text(size = 10, face = "bold"),
        legend.text = element_text(size = 9),
        panel.grid.minor = element_blank(),
        panel.grid.major.y = element_line(color = "gray90", size = 0.5),
        panel.grid.major.x = element_line(color = "gray90", size = 0.5)
      ) +
      guides(
        color = guide_colorbar(order = 1),
        size = guide_legend(order = 2)
      )
    
    return(p)
  }
  
  # Step 4: Generate and save enrichment plots
  cat("=== Generating Enrichment Plots ===\n")
  
  # GO BP dotplot
  if(nrow(go_bp) > 0) {
    p_custom_bp <- create_custom_dotplot(go_bp, "GO Biological Process Enrichment", 15)
    if(!is.null(p_custom_bp)) {
      print(p_custom_bp)
      ggsave("GO_BP_dotplot.png", p_custom_bp, width = 12, height = 8, dpi = 300)
      ggsave("GO_BP_dotplot.pdf", p_custom_bp, width = 12, height = 8)
    }
    write.csv(as.data.frame(go_bp), "GO_BP_enrichment.csv", row.names = FALSE)
  }
  
  # GO MF dotplot
  if(nrow(go_mf) > 0) {
    p_custom_mf <- create_custom_dotplot(go_mf, "GO Molecular Function Enrichment", 15)
    if(!is.null(p_custom_mf)) {
      print(p_custom_mf)
      ggsave("GO_MF_dotplot.png", p_custom_mf, width = 12, height = 8, dpi = 300)
      ggsave("GO_MF_dotplot.pdf", p_custom_mf, width = 12, height = 8)
    }
    write.csv(as.data.frame(go_mf), "GO_MF_enrichment.csv", row.names = FALSE)
  }
  
  # GO CC results (save data only)
  if(nrow(go_cc) > 0) {
    write.csv(as.data.frame(go_cc), "GO_CC_enrichment.csv", row.names = FALSE)
  }
  
  # KEGG dotplot
  if(nrow(kegg) > 0) {
    p_custom_kegg <- create_custom_dotplot(kegg, "KEGG Pathway Enrichment", 15)
    if(!is.null(p_custom_kegg)) {
      print(p_custom_kegg)
      ggsave("KEGG_dotplot.png", p_custom_kegg, width = 12, height = 8, dpi = 300)
      ggsave("KEGG_dotplot.pdf", p_custom_kegg, width = 12, height = 8)
    }
    write.csv(as.data.frame(kegg), "KEGG_enrichment.csv", row.names = FALSE)
  }
  
  # Step 5: Enrichment analysis summary
  cat("=== Enrichment Analysis Summary ===\n")
  cat("GO BP significant terms:", nrow(go_bp), "\n")
  cat("GO MF significant terms:", nrow(go_mf), "\n")
  cat("GO CC significant terms:", nrow(go_cc), "\n")
  cat("KEGG significant pathways:", nrow(kegg), "\n")
  
  # Display top results
  if(nrow(go_bp) > 0) {
    cat("\n=== Top 5 GO BP Terms ===\n")
    print(as.data.frame(go_bp)[1:min(5, nrow(go_bp)), c("Description", "Count", "p.adjust")])
  }
  
  if(nrow(kegg) > 0) {
    cat("\n=== Top 5 KEGG Pathways ===\n")
    print(as.data.frame(kegg)[1:min(5, nrow(kegg)), c("Description", "Count", "p.adjust")])
  }
  
} else {
  cat("Insufficient genes for enrichment analysis. Skipping enrichment steps.\n")
}

# ==================== Gene-specific Boxplots ====================

# Step 1: Define target gene list
cat("=== Drawing Gene-specific Boxplots ===\n")

# Method 1: Manually specify genes of interest (gene symbols)
target_genes <- c("GPX4", "TNF", "IL6", "IL1B", "SLC7A11", "FTH1", "NFE2L2", "NFKB1", "NFKB2")

# Method 2: Or select top differential genes
if(nrow(deg_results[deg_results$Change == "Up-regulated", ]) > 0) {
  top_up_genes <- deg_results[deg_results$Change == "Up-regulated", "SYMBOL"][1:min(5, sum(deg_results$Change == "Up-regulated"))]
  top_up_genes <- top_up_genes[!is.na(top_up_genes)]
} else {
  top_up_genes <- c()
}

if(nrow(deg_results[deg_results$Change == "Down-regulated", ]) > 0) {
  top_down_genes <- deg_results[deg_results$Change == "Down-regulated", "SYMBOL"][1:min(5, sum(deg_results$Change == "Down-regulated"))]
  top_down_genes <- top_down_genes[!is.na(top_down_genes)]
} else {
  top_down_genes <- c()
}

top_deg_genes <- c(top_up_genes, top_down_genes)

# Choose which gene list to use
genes_to_plot <- target_genes  # or use top_deg_genes

cat("Genes to be plotted:", paste(genes_to_plot, collapse = ", "), "\n")

# Step 2: Create plotting function
create_gene_boxplot <- function(gene_symbol, expr_data, sample_groups, title_suffix = "") {
  
  # Check if gene exists
  gene_idx <- which(gene_symbols_filtered == gene_symbol)
  if(length(gene_idx) == 0) {
    cat("Warning: Gene", gene_symbol, "not found in data\n")
    return(NULL)
  }
  
  # Extract gene expression data
  gene_expr <- expr_data[gene_idx[1], ]  # Take first if duplicated
  
  # Create plotting data frame
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
  tle_expr <- plot_data$Expression[plot_data$Group == "TLE"]
  
  if(length(control_expr) > 0 & length(tle_expr) > 0) {
    t_test_result <- t.test(tle_expr, control_expr)
    p_value <- t_test_result$p.value
    
    # Format p-value
    if(p_value < 0.001) {
      p_text <- "p < 0.001"
    } else if(p_value < 0.01) {
      p_text <- paste0("p = ", sprintf("%.3f", p_value))
    } else {
      p_text <- paste0("p = ", sprintf("%.3f", p_value))
    }
     # Calculate fold change (log2 scale for transformed data)
    mean_control <- mean(control_expr, na.rm = TRUE)
    mean_tle <- mean(tle_expr, na.rm = TRUE)
    fold_change <- mean_tle - mean_control  # For log-transformed data
    
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
    scale_fill_manual(values = c("Control" = "#2E86AB", "TLE" = "#A23B72")) +
    scale_color_manual(values = c("Control" = "#1B5E7A", "TLE" = "#7A1E4F")) +
    
    # Add statistical significance annotation
    annotate("text", x = 1.5, y = max(plot_data$Expression) * 1.05, 
             label = p_text, size = 4, fontface = "bold") +
    
    # Set titles and labels
    labs(
      title = paste0("Expression of ", gene_symbol, title_suffix),
      subtitle = paste0("n(Control) = ", sum(plot_data$Group == "Control"), 
                        ", n(TLE) = ", sum(plot_data$Group == "TLE")),
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
      legend.position = "none",  # Hide legend since x-axis shows groups
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

# Step 3: Batch generate individual gene plots
cat("Drawing individual gene boxplots...\n")

# Store all plots
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

# Step 4: Create combined plot
if(length(gene_plots) > 1) {
  cat("Creating combined plot...\n")
  
  # Calculate grid layout
  n_plots <- length(gene_plots)
  if(n_plots <= 4) {
    ncol <- 2
    nrow <- ceiling(n_plots / 2)
  } else if(n_plots <= 9) {
    ncol <- 3
    nrow <- ceiling(n_plots / 3)
  } else {
    ncol <- 4
    nrow <- ceiling(n_plots / 4)
  }
  
  # Create combined plot
  combined_plot <- do.call(grid.arrange, c(gene_plots, ncol = ncol, nrow = nrow))
  
  # Save combined plot
  ggsave("combined_gene_boxplots.png", combined_plot, 
         width = ncol * 6, height = nrow * 5, dpi = 300)
  ggsave("combined_gene_boxplots.pdf", combined_plot, 
         width = ncol * 6, height = nrow * 5)
}

# Step 5: Create gene expression heatmap
cat("Creating gene expression heatmap...\n")

if(length(valid_genes) > 1) {
  # Extract target gene expression data
  target_indices <- which(gene_symbols_filtered %in% valid_genes)
  heatmap_data <- expr_log[target_indices, ]
  rownames(heatmap_data) <- gene_symbols_filtered[target_indices]
  
  # Prepare annotation information
  col_annotation <- data.frame(
    Group = sample_groups[colnames(heatmap_data)],
    row.names = colnames(heatmap_data)
  )
  
  # Set annotation colors
  ann_colors <- list(
    Group = c("Control" = "#2E86AB", "TLE" = "#A23B72")
  )
  
  # Draw heatmap
  pheatmap(
    heatmap_data,
    scale = "row",  # Row-wise standardization
    clustering_distance_rows = "euclidean",
    clustering_distance_cols = "euclidean",
    clustering_method = "complete",
    annotation_col = col_annotation,
    annotation_colors = ann_colors,
    color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
    border_color = NA,
    fontsize = 10,
    fontsize_row = 10,
    fontsize_col = 8,
    main = paste0("Expression Heatmap of Selected Genes\n(n=", length(valid_genes), " genes)"),
    filename = "gene_expression_heatmap.png",
    width = 10,
    height = 8
  )
  
  pheatmap(
    heatmap_data,
    scale = "row",
    clustering_distance_rows = "euclidean",
    clustering_distance_cols = "euclidean",
    clustering_method = "complete",
    annotation_col = col_annotation,
    annotation_colors = ann_colors,
    color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
    border_color = NA,
    fontsize = 10,
    fontsize_row = 10,
    fontsize_col = 8,
    main = paste0("Expression Heatmap of Selected Genes\n(n=", length(valid_genes), " genes)"),
    filename = "gene_expression_heatmap.pdf",
    width = 10,
    height = 8
  )
}

# Step 6: Create gene expression statistics table
cat("Creating gene expression statistics table...\n")

# Calculate statistics for each gene in both groups
gene_stats_list <- list()

for(gene in valid_genes) {
  gene_idx <- which(gene_symbols_filtered == gene)[1]
  gene_expr <- expr_log[gene_idx, ]
  
  plot_data <- data.frame(
    Expression = as.numeric(gene_expr),
    Group = sample_groups[names(gene_expr)]
  )
  plot_data <- plot_data[complete.cases(plot_data), ]
  
  # Calculate statistics
  control_expr <- plot_data$Expression[plot_data$Group == "Control"]
  tle_expr <- plot_data$Expression[plot_data$Group == "TLE"]
  
  if(length(control_expr) > 0 & length(tle_expr) > 0) {
    t_test <- t.test(tle_expr, control_expr)
    
    gene_stats_list[[gene]] <- data.frame(
      Gene = gene,
      Control_n = length(control_expr),
      Control_mean = mean(control_expr),
      Control_sd = sd(control_expr),
      TLE_n = length(tle_expr),
      TLE_mean = mean(tle_expr),
      TLE_sd = sd(tle_expr),
      Log2FC = mean(tle_expr) - mean(control_expr),
      P_value = t_test$p.value,
      CI_lower = t_test$conf.int[1],
      CI_upper = t_test$conf.int[2],
      stringsAsFactors = FALSE
    )
  }
}

# Merge statistical results
if(length(gene_stats_list) > 0) {
  gene_stats_df <- do.call(rbind, gene_stats_list)
  gene_stats_df$P_adjusted <- p.adjust(gene_stats_df$P_value, method = "BH")
  gene_stats_df$Significance <- ifelse(gene_stats_df$P_value < 0.001, "***",
                                      ifelse(gene_stats_df$P_value < 0.01, "**",
                                            ifelse(gene_stats_df$P_value < 0.05, "*", "ns")))
  
  # Format numeric values
  numeric_cols <- c("Control_mean", "Control_sd", "TLE_mean", "TLE_sd", 
                   "Log2FC", "P_value", "P_adjusted", "CI_lower", "CI_upper")
  gene_stats_df[numeric_cols] <- lapply(gene_stats_df[numeric_cols], function(x) round(x, 4))
  
  # Save statistics table
  write.csv(gene_stats_df, "selected_genes_statistics.csv", row.names = FALSE)
  
  cat("=== Selected Gene Statistics Summary ===\n")
  print(gene_stats_df[, c("Gene", "Log2FC", "P_value", "P_adjusted", "Significance")])
}

# Step 7: Create volcano plot highlighting selected genes
cat("Creating volcano plot with selected genes highlighted...\n")

# Prepare volcano plot data
volcano_data <- deg_results
volcano_data$log10_pval <- -log10(volcano_data$P.Value)
volcano_data$Highlight <- ifelse(volcano_data$SYMBOL %in% valid_genes, "Selected", "Other")

# Draw volcano plot
volcano_plot <- ggplot(volcano_data, aes(x = logFC, y = log10_pval)) +
  # Add all points
  geom_point(data = volcano_data[volcano_data$Highlight == "Other", ], 
             aes(color = Change), alpha = 0.6, size = 1) +
  # Highlight selected genes
  geom_point(data = volcano_data[volcano_data$Highlight == "Selected", ], 
             color = "red", size = 3, alpha = 0.8) +
  # Add gene labels
  geom_text_repel(
    data = volcano_data[volcano_data$Highlight == "Selected", ],
    aes(label = SYMBOL),
    color = "red",
    fontface = "bold",
    size = 4,
    box.padding = 0.5,
    point.padding = 0.3,
    segment.color = "red",
    segment.alpha = 0.6
  ) +
  # Add threshold lines
  geom_vline(xintercept = c(-logFC_threshold, logFC_threshold), 
             linetype = "dashed", color = "gray60") +
  geom_hline(yintercept = -log10(pvalue_threshold), 
             linetype = "dashed", color = "gray60") +
  # Set colors
  scale_color_manual(values = c("Up-regulated" = "#E31A1C", 
                               "Down-regulated" = "#1F78B4", 
                               "Not Significant" = "gray70")) +
  # Set titles and labels
  labs(
    title = "Volcano Plot with Selected Genes Highlighted",
    subtitle = paste0("Selected genes: ", paste(valid_genes, collapse = ", ")),
    x = "Log2 Fold Change (TLE vs Control)",
    y = "-Log10(P-value)",
    color = "Regulation"
  ) +
  # Beautify theme
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 10, color = "gray60"),
    axis.title = element_text(size = 12, face = "bold"),
    axis.text = element_text(size = 10),
    legend.position = "bottom",
    legend.title = element_text(size = 11, face = "bold"),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "gray90", size = 0.3)
  )

print(volcano_plot)

# Save volcano plot
ggsave("volcano_plot_with_selected_genes.png", volcano_plot, 
       width = 12, height = 10, dpi = 300)
ggsave("volcano_plot_with_selected_genes.pdf", volcano_plot, 
       width = 12, height = 10)

# ==================== Final Summary and Cleanup ====================

# Step 8: Final analysis summary
cat("=== Gene-specific Analysis Summary ===\n")
cat("Successfully analyzed genes:", length(valid_genes), "\n")
cat("Analyzed genes:", paste(valid_genes, collapse = ", "), "\n")
cat("Generated files:\n")
cat("- Individual gene boxplots: boxplot_[gene_name].png/pdf\n")
if(length(gene_plots) > 1) {
  cat("- Combined boxplots: combined_gene_boxplots.png/pdf\n")
}
if(length(valid_genes) > 1) {
  cat("- Expression heatmap: gene_expression_heatmap.png/pdf\n")
}
cat("- Statistics summary table: selected_genes_statistics.csv\n")
cat("- Highlighted volcano plot: volcano_plot_with_selected_genes.png/pdf\n")

# Create comprehensive analysis summary
cat("\n=== Complete Analysis Summary ===\n")
cat("Dataset: GSE256068\n")
cat("Total samples processed:", length(sample_groups), "\n")
cat("Sample groups:", paste(names(table(sample_groups)), "=", table(sample_groups), collapse = ", "), "\n")
cat("Total genes after filtering:", nrow(expr_log), "\n")
cat("Differential genes (p < 0.05, |logFC| >= 0.58):", sum(deg_results$Change != "Not Significant"), "\n")
cat("  - Up-regulated:", sum(deg_results$Change == "Up-regulated"), "\n")
cat("  - Down-regulated:", sum(deg_results$Change == "Down-regulated"), "\n")

if(exists("go_bp") && nrow(go_bp) > 0) {
  cat("GO BP enriched terms:", nrow(go_bp), "\n")
}
if(exists("go_mf") && nrow(go_mf) > 0) {
  cat("GO MF enriched terms:", nrow(go_mf), "\n")
}
if(exists("go_cc") && nrow(go_cc) > 0) {
  cat("GO CC enriched terms:", nrow(go_cc), "\n")
}
if(exists("kegg") && nrow(kegg) > 0) {
  cat("KEGG enriched pathways:", nrow(kegg), "\n")
}

cat("Specifically analyzed genes:", length(valid_genes), "\n")

# Save session info for reproducibility
session_info <- sessionInfo()
capture.output(print(session_info), file = "session_info.txt")

cat("\n=== Analysis Complete ===\n")
cat("All results have been saved to the working directory.\n")
cat("Session information saved to: session_info.txt\n")
cat("Analysis completed successfully!\n")

# Reset graphics parameters
par(mfrow = c(1, 1))

# Optional: Clean up large objects to free memory
# rm(expr_data, expr_data_anno, expr_matrix, expr_matrix_filtered)
# gc()  # Garbage collection

cat("\nGene-specific boxplot analysis completed!\n")
