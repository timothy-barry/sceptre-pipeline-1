#!/usr/bin/env Rscript

#######
# SETUP
#######
# get the command line args
args <- commandArgs(trailingOnly = TRUE)
multimodal_metadata_fp <- args[1]
gene_odm_fp <- args[2]
gRNA_odm_fp <- args[3]
other_args <- args[seq(4, length(args))]

# process "other_args"
l_other_args <- length(other_args)
idx_list <- vector(mode = "list", length = 4)
for (i in seq(1, 4)) idx_list[[i]] <- seq(from =  1 + l_other_args/4 * (i - 1), length.out = l_other_args/4)
gene_ids <- other_args[idx_list[[1]]]
gRNA_ids <- other_args[idx_list[[2]]]
gene_precomp_fps <- other_args[idx_list[[3]]]
gRNA_precomp_fps <- other_args[idx_list[[4]]]
n_pairs <- length(gene_ids)

# library ondisc
library(ondisc)

##################################
# PREPARE DATA AND SET HYPERPARAMS
##################################
mm_odm <- read_multimodal_odm(c(gene_odm_fp, gRNA_odm_fp), multimodal_metadata_fp)
gene_odm <- mm_odm |> get_modality("gene")
gRNA_odm <- mm_odm |> get_modality("gRNA")
threshold <- get_modality(mm_odm, "gRNA")@misc$threshold
side <- get_modality(mm_odm, "gRNA")@misc$side
rm(mm_odm)

##########################
# ITERATE AND APPLY METHOD
##########################
out_l <- vector(mode = "list", length = n_pairs)
for (i in seq(1, n_pairs)) {
  # only load gene data if necessary
  if (i == 1 || gene_ids[i] != gene_ids[i - 1]) {
    gene_precomp <- readRDS(gene_precomp_fps[i])
    gene_expressions <- as.numeric(gene_odm[[gene_precomp$id,]])
  }
  # only load gRNA data if necessary
  if (i == 1 || gRNA_ids[i] != gRNA_ids[i - 1]) {
    gRNA_precomp <- readRDS(gRNA_precomp_fps[i])
    gRNA_indicators <- load_thresholded_and_grouped_gRNA(covariate_odm = gRNA_odm,
                                                         gRNA_group = gRNA_precomp$id,
                                                         threshold = threshold) |> as.integer()
  }
  # run the dCRT; append gene id and gRNA group to the output
  res <- sceptre:::run_sceptre_using_precomp_fast(expressions = gene_expressions,
                                                  gRNA_indicators = gRNA_indicators,
                                                  gRNA_precomp = gRNA_precomp$resamples,
                                                  side = side,
                                                  gene_precomp_size = gene_precomp$gene_precomp_size,
                                                  gene_precomp_offsets = gene_precomp$gene_precomp_offsets,
                                                  full_output = FALSE)
  res <- res |> dplyr::mutate(gene_id = gene_precomp$id, gRNA_group = gRNA_precomp$id)
  out_l[[i]] <- res
}

out <- dplyr::mutate_at(do.call(rbind, out_l), c("gene_id", "gRNA_group"), factor)
saveRDS(object = out, file = "raw_result.rds")
