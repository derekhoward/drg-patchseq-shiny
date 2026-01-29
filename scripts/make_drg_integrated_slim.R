#!/usr/bin/env Rscript

library(Seurat)

args <- commandArgs(trailingOnly = TRUE)
input_path <- if (length(args) >= 1) args[1] else "data/drg_integrated.RDS"
output_path <- if (length(args) >= 2) args[2] else "data/drg_integrated_slim.RDS"

if (!file.exists(input_path)) {
  stop("Input file not found: ", input_path)
}
if (!file.exists("data/gene_list.rds")) {
  stop("Required file not found: data/gene_list.rds")
}

obj <- readRDS(input_path)
gene_list <- readRDS("data/gene_list.rds")

if (!"integrated" %in% names(obj@assays)) {
  stop("Assay 'integrated' not found in the Seurat object.")
}
if (!"umap" %in% names(obj@reductions)) {
  stop("UMAP reduction not found in the Seurat object.")
}

DefaultAssay(obj) <- "integrated"
keep_genes <- intersect(gene_list, rownames(obj@assays[["integrated"]]))

diet_args <- list(
  object = obj,
  assays = "integrated",
  features = keep_genes,
  counts = FALSE,
  data = TRUE,
  scale.data = FALSE
)

diet_formals <- names(formals(Seurat::DietSeurat))
if ("reductions" %in% diet_formals) {
  diet_args$reductions <- "umap"
} else {
  diet_args$dimreducs <- "umap"
}

obj <- do.call(Seurat::DietSeurat, diet_args)
obj@meta.data <- obj@meta.data[, c("labels", "labels.p", "dataset"), drop = FALSE]

saveRDS(obj, output_path)
message("Wrote slim Seurat object to ", output_path)
