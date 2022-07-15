#!/usr/bin/env Rscript

#######
# SETUP
#######
library(ondisc)

# get the command line args
args <- commandArgs(trailingOnly = TRUE)
multimodal_metadata_fp <- args[1] # multimodal metadata fp
gene_odm_fp <- args[2] # gene odm fp
grna_odm_fp <- args[3] # grna odm fp
gene_precomp_matrix_fp <- args[4] # gene precomp matrix fp
grna_precomp_matrix_fp <- args[5] # grna precomp matrix fp
pair_to_pod_id_map_fp <- args[6] # pair to pod id map fp
pod_id <- as.integer(args[7]) # the pod id
threshold <- as.integer(args[8]) # threshold
B <- as.integer(args[9]) # B
side <- args[10] # sidedness of test
full_output <- as.logical(args[11]) # full_output

# obtain the pairs to analyze and prepare data
pairs_to_analyze <- readRDS(pair_to_pod_id_map_fp) |>
  dplyr::filter(pod_id == !!pod_id) |>
  dplyr::select(gene_id, grna_group)
gene_ids <- as.character(pairs_to_analyze$gene_id)
grna_group_ids <- as.character(pairs_to_analyze$grna_group)
n_pairs <- nrow(pairs_to_analyze)
rm(pairs_to_analyze)
mm_odm <- read_multimodal_odm(c(gene_odm_fp, grna_odm_fp), multimodal_metadata_fp)
gene_odm <- mm_odm |> get_modality("gene")
grna_odm <- mm_odm |> get_modality("grna")
global_cell_covariates <- mm_odm |> get_cell_covariates() |> as.matrix()
rm(mm_odm)


###########################################################
# Iterate over pairs and apply pairwise test of association
###########################################################
binomial_obj <- binomial()
out_l <- vector(mode = "list", length = n_pairs)
for (i in seq(1, n_pairs)) {
  gene_id <- gene_ids[i]
  grna_group_id <- grna_group_ids[i]
  print(paste0("Testing gene ", gene_id, " and gRNA group ", grna_group_id, "."))
  # Load gene expressions and fitted means (if necessary)
  if (i == 1 || gene_ids[i] != gene_ids[i - 1]) {
    # get fitted means and theta (if applicable)
    gene_precomp <- fst::read_fst(gene_precomp_matrix_fp, gene_id) |> dplyr::pull()
    gene_fitted_coefs <- gene_precomp[-length(gene_precomp)]
    gene_theta <- gene_precomp[length(gene_precomp)]
    gene_fitted_linear_components <- (global_cell_covariates %*% gene_fitted_coefs)[,1]
    # get expressions
    gene_expressions <- as.numeric(gene_odm[[gene_id,]])
  } 
  # Load grna indicators and fitted means (if necessary)
  if (i == 1 || grna_group_ids[i] != grna_group_ids[i - 1]) {
    # get fitted means
    grna_group_fitted_coefs <- fst::read_fst(grna_precomp_matrix_fp, grna_group_id) |> dplyr::pull()
    grna_group_fitted_means <- binomial_obj$linkinv((global_cell_covariates %*% grna_group_fitted_coefs)[,1]) 
    # get indicators
    grna_group_indicators <- ondisc::load_thresholded_and_grouped_grna(covariate_odm = grna_odm,
                                                                       grna_group = grna_group_id,
                                                                       threshold = threshold) |> as.integer()
    # sample the sparse matrix of grna presences/absences
    synthetic_indicator_matrix <- sceptre:::generate_synthetic_grna_data(fitted_probs = grna_group_fitted_means, B = B)
  }
  # carry out the dCRT test
  out <- sceptre:::run_sceptre_using_precomp_fast(expressions = gene_expressions,
                                                  gRNA_indicators = grna_group_indicators,
                                                  gRNA_precomp = synthetic_indicator_matrix, 
                                                  side = side,
                                                  gene_precomp_size = gene_theta,
                                                  gene_precomp_offsets = gene_fitted_linear_components,
                                                  full_output = full_output) |>
    dplyr::mutate(gene_id = gene_id, grna_group = grna_group_id)
  # add to output list
  out <- data.table::setDT(out) |>
    dplyr::mutate(gene_id = factor(gene_id), grna_group = factor(grna_group))
  out_l[[i]] <- out
}

raw_result <- data.table::rbindlist(out_l)
saveRDS(object = raw_result, file = "raw_result.rds")