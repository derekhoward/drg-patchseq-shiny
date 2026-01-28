#!/usr/bin/env Rscript

library(Seurat)
library(Matrix)

input_seurat <- "data/pig_patchseq.RDS"
input_meta <- "data/patchseq-meta.csv"
output_dir <- "data"

# Set to FALSE if gene-level aggregation is too heavy for your machine.
compute_gene_medians <- TRUE

exclude_meta_cols <- c(
  "cell_id", "labels", "orig.ident", "Animal", "age", "weight",
  "batch", "species", "PC1"
)

starts_with_exclude <- "ADS_"

seu <- readRDS(input_seurat)
patch_meta <- read.csv(input_meta, stringsAsFactors = FALSE)

meta <- seu@meta.data
if ("cell_id" %in% colnames(meta)) {
  meta$cell_id <- as.character(meta$cell_id)
} else {
  meta$cell_id <- rownames(meta)
}

patch_meta$cell_id <- as.character(patch_meta$cell_id)

labels <- patch_meta$labels[match(meta$cell_id, patch_meta$cell_id)]
meta$labels <- labels

missing_labels <- sum(is.na(meta$labels))
if (missing_labels > 0) {
  message("Warning: missing labels for ", missing_labels, " cells.")
}

is_numeric <- vapply(patch_meta, is.numeric, logical(1))
numeric_cols <- names(patch_meta)[is_numeric]
ephys_cols <- setdiff(numeric_cols, exclude_meta_cols)
ephys_cols <- ephys_cols[!startsWith(ephys_cols, starts_with_exclude)]
ephys_cols <- ephys_cols[ephys_cols != "Sinusscore_new"]

ephys_df <- patch_meta[match(meta$cell_id, patch_meta$cell_id), ephys_cols, drop = FALSE]
if ("Sinusscore_old" %in% colnames(ephys_df)) {
  colnames(ephys_df)[colnames(ephys_df) == "Sinusscore_old"] <- "Sinus score"
}
ephys_cols <- colnames(ephys_df)
cell_meta <- data.frame(
  cell_id = meta$cell_id,
  labels = meta$labels,
  ephys_df,
  check.names = FALSE
)

assay_name <- DefaultAssay(seu)
counts <- GetAssayData(seu, assay = assay_name, slot = "counts")
lib_sizes <- Matrix::colSums(counts)
lib_sizes[lib_sizes == 0] <- NA_real_

expr_mat_cpm <- Matrix::t(Matrix::t(counts) / lib_sizes) * 1e6

gene_list <- rownames(expr_mat_cpm)
efeat_list <- ephys_cols

agg_ephys <- aggregate(
  cell_meta[, ephys_cols, drop = FALSE],
  by = list(labels = cell_meta$labels),
  FUN = median,
  na.rm = TRUE
)

agg_gene_cpm <- NULL
if (compute_gene_medians) {
  if (!requireNamespace("matrixStats", quietly = TRUE)) {
    stop("matrixStats is required for gene-level median aggregation.")
  }

  label_levels <- unique(cell_meta$labels)
  agg_gene_cpm <- lapply(label_levels, function(lbl) {
    cols <- which(cell_meta$labels == lbl)
    if (length(cols) == 0) {
      rep(NA_real_, nrow(expr_mat_cpm))
    } else {
      matrixStats::rowMedians(as.matrix(expr_mat_cpm[, cols, drop = FALSE]), na.rm = TRUE)
    }
  })

  agg_gene_cpm <- do.call(cbind, agg_gene_cpm)
  colnames(agg_gene_cpm) <- label_levels
  rownames(agg_gene_cpm) <- gene_list
}

saveRDS(expr_mat_cpm, file.path(output_dir, "expr_mat_cpm.rds"))
saveRDS(cell_meta, file.path(output_dir, "cell_meta.rds"))
saveRDS(agg_ephys, file.path(output_dir, "agg_ephys.rds"))
saveRDS(agg_gene_cpm, file.path(output_dir, "agg_gene_cpm_median.rds"))
saveRDS(gene_list, file.path(output_dir, "gene_list.rds"))
saveRDS(efeat_list, file.path(output_dir, "efeat_list.rds"))

message("Done. Wrote app data to ", output_dir)
