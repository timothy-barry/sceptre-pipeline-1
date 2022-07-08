#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
precomp_submats <- lapply(args, function(fp) readRDS(fp))
combined_precomp_mat <- data.table::rbindlist(precomp_submats)
rm(precomp_submats)
ids <- combined_precomp_mat$id
combined_precomp_mat$id <- NULL
combined_precomp_mat_t <- data.table::as.data.table(t(combined_precomp_mat))
rm(combined_precomp_mat)
colnames(combined_precomp_mat_t) <- ids
fst::write_fst(x = combined_precomp_mat_t, path = "precomputation_matrix.fst")
