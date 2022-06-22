#!/usr/bin/env Rscript

#######
# SETUP
#######

# get the command line args
args <- commandArgs(trailingOnly = TRUE)
multimodal_metadata_fp <- args[1]
gene_odm_fp <- args[2]
gRNA_odm_fp <- args[3]
pair_fp <- args[4]

# load ondisc
library(ondisc)

# create the multimodal odm
mm_odm <- read_multimodal_odm(c(gene_odm_fp, gRNA_odm_fp), multimodal_metadata_fp)

# read the pairs
pairs <- readRDS(pair_fp)

########
# CHECKS
########

modality_names <- names(mm_odm@modalities)
# 1. "gene" must be in modality_names
if (!("gene" %in% modality_names)) {
  stop("`gene` must be a modality name in the multimodal ondisc matrix.")
}

# 2. "gRNA" must be in modality_names
if (!("gRNA" %in% modality_names)) {
  stop("`gRNA` must be a modality name in the multimodal ondisc matrix.")
}

# 3. "gRNA_group" must be a column of the gRNA feature covariate matrix
gRNA_feature_covariates <- mm_odm |>
  get_modality("gRNA") |>
  get_feature_covariates()
if (!("gRNA_group" %in% colnames(gRNA_feature_covariates))) {
  stop("The feature covariates matrix of the `gRNA` ondisc matrix must have a column called `gRNA_group`.")
}

# 4. "gene_id" must be a column of pairs
if (!("gene_id" %in% colnames(pairs))) {
  stop("The pairs data frame must contain a column called `gene_id`.")
}

# 5. "gRNA_group" must be a column of pairs
if (!("gRNA_group" %in% colnames(pairs))) {
  stop("The pairs data frame must contain a column called `gRNA_group`.")
}

# 6. check that the gene IDs in the pairs data frame are a subset of the gene IDs of the gene ODM
odm_gene_ids <- mm_odm |>
  get_modality("gene") |>
  get_feature_ids()
pairs_gene_ids <- as.character(pairs$gene_id)
if (!all(pairs_gene_ids %in% odm_gene_ids)) {
  stop("The gene IDs in the pairs data frame are not a subset of the gene IDs in the gene ondisc matrix.")
}

# 7. check that the gRNA groups in the pairs data frame are a subset of the gRNA groups in the gRNA ODM
odm_gRNA_groups <- gRNA_feature_covariates$gRNA_group
pairs_gRNA_groups <- as.character(pairs$gRNA_group)
if (!all(pairs_gRNA_groups %in% odm_gRNA_groups)) {
  stop("The gRNA groups in the pairs data frame are not a subset of the gRNA groups in the gRNA ondisc matrix.")
}

# 8. check for weird numbers in the global cell covariate matrix
global_cell_covariates <- mm_odm |> get_cell_covariates()
if (ncol(global_cell_covariates) == 0) {
  stop("The global cell covariate matrix of the multimodal ondisc matrix should contain at least one column.")
}
for (col_name in colnames(global_cell_covariates)) {
  vect <- global_cell_covariates[[col_name]]
  if (any(vect == -Inf) || any(vect == Inf) || any(is.na(vect))) {
    stop(paste0("The column `", col_name, "` of the global cell covariate matrix contains entries that are either NA, -Inf, or Inf. Remove these entries (by, for example, removing the corresponding cells from the multimodal ondisc matrix)."))
  }
}

################################
# PRINT GENE IDS and GRNA GROUPS
################################
if (TRUE) {
  set.seed(4)
  pairs <- pairs |> dplyr::filter(gene_id %in% c("ENSG00000135018", "ENSG00000155592", "ENSG00000197937") &
                                  gRNA_group %in% c("random_9", "scrambled_21", "random_22"))
}

write_vector <- function(file_name, vector) {
  file_con <- file(file_name)
  writeLines(as.character(vector), file_con)
  close(file_con)
}
gene_id_v <- as.character(pairs$gene_id)
gRNA_group_v <- as.character(pairs$gRNA_group)

write_vector("gene_ids.txt", unique(gene_id_v))
write_vector("grna_groups.txt", unique(gRNA_group_v))
write_vector("pairs.txt", paste0(gene_id_v, " ", gRNA_group_v))
