#!/usr/bin/env Rscript

# Get CL args
args <- commandArgs(trailingOnly = TRUE)
n_args <- length(args)
result_file_name <- args[1]
raw_results <- args[seq(2, n_args)]

# combine raw results and save
out <- lapply(X = raw_results, FUN = readRDS) |> dplyr::bind_rows()
saveRDS(object = out, file = result_file_name)
