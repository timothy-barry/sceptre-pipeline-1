#!/usr/bin/env Rscript

# Get CL args
args <- commandArgs(trailingOnly = TRUE)
n_args <- length(args)
result_file_name <- args[1]
pairs_fp <- args[2]
raw_results <- args[seq(3, n_args)]

# combine raw results and save
results <- lapply(X = raw_results, FUN = readRDS) |> data.table::rbindlist()
pairs <- readRDS(pairs_fp)
data.table::setDF(results)
out <- dplyr::left_join(x = pairs, y = results, by = c("gene_id", "grna_group")) |>
  dplyr::relocate(gene_id, grna_group, p_value, z_value)
saveRDS(object = out, file = result_file_name)
