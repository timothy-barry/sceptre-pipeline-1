#!/usr/bin/env Rscript

#######
# SETUP
#######

# get the command line args
args <- commandArgs(trailingOnly = TRUE)
modality_name <- args[1]
multimodal_metadata_fp <- args[2]
gene_odm_fp <- args[3]
gRNA_odm_fp <- args[4]
ids <- args[seq(5, length(args))]

# load ondisc
library(ondisc)

##############################
# PREPARE DATA AND HYPERPARAMS
##############################
mm_odm <- read_multimodal_odm(c(gene_odm_fp, gRNA_odm_fp), multimodal_metadata_fp)
threshold <- get_modality(mm_odm, "gRNA")@misc$threshold
B <- get_modality(mm_odm, "gRNA")@misc$B
modality_odm <- get_modality(mm_odm, modality_name)
global_cell_covariates <- get_cell_covariates(mm_odm)
rm(mm_odm)

####################################
# ITERATE AND PERFORM PRECOMPUTATION
####################################
for (id in ids) {
  if (modality_name == "gene") {
    expressions <- as.numeric(modality_odm[[id,]])
    precomp <- sceptre:::run_gene_precomputation(expressions = expressions,
                                                 covariate_matrix = global_cell_covariates,
                                                 gene_precomp_size = NULL) 
  } else { # gRNA
    indicators <- load_thresholded_and_grouped_gRNA(covariate_odm = modality_odm,
                                                    gRNA_group = id,
                                                    threshold = threshold) |> as.integer()
    resamples <- sceptre:::run_gRNA_precomputation(gRNA_indicators = indicators,
                                                   covariate_matrix = global_cell_covariates,
                                                   B = B,
                                                   seed = 1)
    precomp <- list(resamples = resamples)
  }
  precomp$id <- id
  to_save_fp <- paste0(id, ".rds")
  saveRDS(precomp, to_save_fp)
}
