#!/usr/bin/env Rscript

# DESeq2 differential expression: HGSOC versus fallopian-tube epithelium (FTE)
# Project: GSE190688
# Input count matrix must contain gene IDs in the first column and one raw-count
# column for each sample. Counts must be non-negative integers, not TPM/FPKM.

suppressPackageStartupMessages(library(DESeq2))

project_dir <- path.expand("~/GSE190688")
count_file <- file.path(project_dir, "counts", "GSE190688_gene_count_matrix.tsv")
metadata_file <- file.path(project_dir, "metadata", "sample_metadata.csv")
out_dir <- file.path(project_dir, "deseq2_results")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(count_file)) stop("Count matrix not found: ", count_file)
if (!file.exists(metadata_file)) stop("Sample metadata not found: ", metadata_file)

# 1. Read input files
count_df <- read.delim(count_file, row.names = 1, check.names = FALSE,
                       stringsAsFactors = FALSE)
metadata <- read.csv(metadata_file, stringsAsFactors = FALSE, check.names = FALSE)

if (!all(c("sample", "condition") %in% colnames(metadata))) {
  stop("Metadata must contain columns named 'sample' and 'condition'.")
}
if (anyDuplicated(rownames(count_df))) stop("Duplicated gene IDs found in count matrix.")
if (anyDuplicated(metadata$sample)) stop("Duplicated sample IDs found in metadata.")
if (!setequal(colnames(count_df), metadata$sample)) {
  stop("Sample columns in count matrix do not match metadata$sample. Check sample IDs.")
}

# Put metadata in exactly the same order as count-matrix columns.
metadata <- metadata[match(colnames(count_df), metadata$sample), , drop = FALSE]
rownames(metadata) <- metadata$sample

count_matrix <- as.matrix(count_df)
suppressWarnings(storage.mode(count_matrix) <- "numeric")
if (anyNA(count_matrix)) stop("Missing or non-numeric values found in count matrix.")
if (any(!is.finite(count_matrix)) || any(count_matrix < 0)) {
  stop("Counts must be finite and non-negative.")
}
if (any(count_matrix != round(count_matrix))) {
  stop("DESeq2 requires raw integer counts; non-integer values were found.")
}
count_matrix <- round(count_matrix)

# Set FTE as the reference level so positive log2 fold changes mean higher in HGSOC.
if (!all(c("FTE", "HGSOC") %in% unique(metadata$condition))) {
  stop("Condition column must include the exact labels 'FTE' and 'HGSOC'.")
}
metadata$condition <- factor(metadata$condition, levels = c("FTE", "HGSOC"))

# Save the actual sample/condition mapping used by this analysis.
write.csv(metadata, file.path(out_dir, "sample_metadata_used.csv"), row.names = FALSE)

# 2. Filter very low-count genes: at least 10 reads in at least 2 samples.
keep <- rowSums(count_matrix >= 10) >= 2
count_filtered <- count_matrix[keep, , drop = FALSE]
if (nrow(count_filtered) == 0) stop("No genes remain after low-count filtering.")

# 3. Build and fit the DESeq2 model.
dds <- DESeqDataSetFromMatrix(
  countData = count_filtered,
  colData = metadata,
  design = ~ condition
)
dds <- DESeq(dds)

# Explicit contrast: HGSOC relative to FTE.
res <- results(dds, contrast = c("condition", "HGSOC", "FTE"), alpha = 0.05)
res_df <- as.data.frame(res)
res_df$gene_id <- rownames(res_df)
res_df <- res_df[, c("gene_id", setdiff(colnames(res_df), "gene_id"))]
res_df <- res_df[order(res_df$padj, na.last = TRUE), , drop = FALSE]

# 4. Export all DESeq2 results and direction-specific significant gene lists.
write.csv(res_df,
          file.path(out_dir, "HGSOC_vs_FTE_DESeq2_results.csv"),
          row.names = FALSE, na = "")

sig_df <- subset(res_df, !is.na(padj) & padj < 0.05)
up_df <- subset(sig_df, log2FoldChange > 0)
down_df <- subset(sig_df, log2FoldChange < 0)

write.csv(sig_df, file.path(out_dir, "Significant_DEGs_padj_lt_0.05.csv"),
          row.names = FALSE, na = "")
write.csv(up_df, file.path(out_dir, "Upregulated_in_HGSOC.csv"),
          row.names = FALSE, na = "")
write.csv(down_df, file.path(out_dir, "Downregulated_in_HGSOC_up_in_FTE.csv"),
          row.names = FALSE, na = "")

# Export normalized counts for expression plots or candidate-gene review.
norm_counts <- counts(dds, normalized = TRUE)
write.csv(data.frame(gene_id = rownames(norm_counts), norm_counts, check.names = FALSE),
          file.path(out_dir, "normalized_counts.csv"), row.names = FALSE)

# 5. Save a compact run summary.
summary_lines <- c(
  "GSE190688 DESeq2 analysis summary",
  paste0("Samples: ", ncol(dds)),
  paste0("FTE samples: ", sum(metadata$condition == "FTE")),
  paste0("HGSOC samples: ", sum(metadata$condition == "HGSOC")),
  paste0("Genes in input matrix: ", nrow(count_matrix)),
  paste0("Genes retained after filter: ", nrow(dds)),
  paste0("Genes with non-NA adjusted p-values: ", sum(!is.na(res$padj))),
  paste0("Significant DEGs (padj < 0.05): ", nrow(sig_df)),
  paste0("Upregulated in HGSOC (log2FC > 0, padj < 0.05): ", nrow(up_df)),
  paste0("Downregulated in HGSOC / higher in FTE (log2FC < 0, padj < 0.05): ", nrow(down_df)),
  paste0("Significant DEGs with |log2FC| >= 1: ",
         sum(!is.na(res$padj) & res$padj < 0.05 & abs(res$log2FoldChange) >= 1)),
  "Direction: positive log2FoldChange = higher in HGSOC; negative = higher in FTE.",
  "Significance threshold: adjusted p-value (Benjamini-Hochberg) < 0.05."
)
writeLines(summary_lines, file.path(out_dir, "DESeq2_summary.txt"))

# 6. Standard DESeq2 diagnostic plots.
pdf(file.path(out_dir, "MA_plot.pdf"), width = 7, height = 6)
plotMA(res, alpha = 0.05, main = "HGSOC vs FTE", ylim = c(-6, 6))
dev.off()

pdf(file.path(out_dir, "dispersion_plot.pdf"), width = 7, height = 6)
plotDispEsts(dds, main = "DESeq2 dispersion estimates")
dev.off()

# Variance-stabilized values are for visualization, not for DESeq2 model fitting.
vsd <- vst(dds, blind = FALSE)
pdf(file.path(out_dir, "PCA_by_condition.pdf"), width = 7, height = 6)
print(plotPCA(vsd, intgroup = "condition") + ggplot2::ggtitle("PCA by condition: HGSOC vs FTE"))
dev.off()

# Volcano plot; significant genes are highlighted by direction.
volcano_df <- res_df[!is.na(res_df$padj) & !is.na(res_df$log2FoldChange), , drop = FALSE]
volcano_df$minus_log10_padj <- -log10(pmax(volcano_df$padj, .Machine$double.xmin))
volcano_df$direction <- "Not significant"
volcano_df$direction[volcano_df$padj < 0.05 & volcano_df$log2FoldChange > 0] <- "Higher in HGSOC"
volcano_df$direction[volcano_df$padj < 0.05 & volcano_df$log2FoldChange < 0] <- "Higher in FTE"
pdf(file.path(out_dir, "Volcano_plot.pdf"), width = 7, height = 6)
plot(volcano_df$log2FoldChange, volcano_df$minus_log10_padj,
     pch = 20, cex = 0.45,
     col = ifelse(volcano_df$direction == "Higher in HGSOC", "firebrick",
           ifelse(volcano_df$direction == "Higher in FTE", "royalblue", "grey70")),
     xlab = "log2 fold change (HGSOC / FTE)", ylab = "-log10 adjusted p-value",
     main = "Differential expression: HGSOC vs FTE")
abline(v = c(-1, 1), lty = 2, col = "grey40")
abline(h = -log10(0.05), lty = 2, col = "grey40")
legend("topright", legend = c("Higher in HGSOC", "Higher in FTE", "Not significant"),
       col = c("firebrick", "royalblue", "grey70"), pch = 20, cex = 0.75, bty = "n")
dev.off()

# Heatmap of the top 30 significant genes using row-wise z-scores.
if (nrow(sig_df) >= 2) {
  top_ids <- head(sig_df$gene_id, min(30, nrow(sig_df)))
  mat <- assay(vsd)[top_ids, , drop = FALSE]
  mat_z <- t(scale(t(mat)))
  mat_z[!is.finite(mat_z)] <- 0
  pdf(file.path(out_dir, "Top30_DEG_heatmap.pdf"), width = 10, height = 9)
  heatmap(mat_z, scale = "none", margins = c(8, 8),
          main = "Top significant DEGs (row-wise z-scores)",
          xlab = "Samples", ylab = "Genes")
  dev.off()
}

cat(paste(summary_lines, collapse = "\n"), "\n\nOutputs saved in:\n", out_dir, "\n")
