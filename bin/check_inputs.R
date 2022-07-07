#!/usr/bin/env Rscript

#######
# SETUP
#######

# get the command line args
args <- commandArgs(trailingOnly = TRUE)
multimodal_metadata_fp <- args[1] # multimodal metadata fp
gene_odm_fp <- args[2] # gene ODM backing file
grna_odm_fp <- args[3] # grna ODM backing file
pair_fp <- args[4] # pairs df
form <- args[5] # formula string
threshold <- as.integer(args[6]) # threshold
B <- as.integer(args[7]) # B
side <- args[8] # sidedness
n_pairs_to_sample <- as.integer(args[9]) # n pairs
gene_modality_name <- args[10] # gene modality name
grna_modality_name <- args[11] # grna modality name
gene_pod_size <- max(as.integer(args[12]), 2) # gene pod size
grna_group_pod_size <- max(as.integer(args[13]), 2) # grna group pod size
pair_pod_size <- max(as.integer(args[14]), 2) # pair pod size

# load ondisc
library(ondisc)

# create the multimodal odm
mm_odm <- read_multimodal_odm(c(gene_odm_fp, grna_odm_fp), multimodal_metadata_fp)

# read the pairs
pairs <- readRDS(pair_fp)


########
# CHECKS
########
modality_names <- names(mm_odm@modalities)
# 1. gene_modality_name must be in modality_names; set gene modality name to `gene`
if (!(gene_modality_name %in% modality_names)) {
  stop(paste0("`", gene_modality_name, "` must be a modality name in the multimodal ondisc matrix."))
} else {
  names(mm_odm@modalities)[names(mm_odm@modalities) == gene_modality_name] <- "gene"
}

# 2. "grna" must be in modality_names; set grna modality name to `grna`
if (!(grna_modality_name %in% modality_names)) {
  stop(paste0("`", grna_modality_name, "` must be a modality name in the multimodal ondisc matrix."))
} else {
  names(mm_odm@modalities)[names(mm_odm@modalities) == grna_modality_name] <- "grna"
}

# 3. "grna_group" must be a column of the grna feature covariate matrix
grna_feature_covariates <- mm_odm |> get_modality("grna") |> get_feature_covariates()
if (!("grna_group" %in% colnames(grna_feature_covariates))) {
  stop("The feature covariates matrix of the `grna` ondisc matrix must have a column called `grna_group`.")
}

# 4. "gene_id" must be a column of pairs
if (!("gene_id" %in% colnames(pairs))) {
  stop("The pairs data frame must contain a column called `gene_id`.")
}

# 5. "grna_group" must be a column of pairs
if (!("grna_group" %in% colnames(pairs))) {
  stop("The pairs data frame must contain a column called `grna_group`.")
}

# 6. check that the gene IDs in the pairs data frame are a subset of the gene IDs of the gene ODM
odm_gene_ids <- mm_odm |> get_modality("gene") |> get_feature_ids()
pairs_gene_ids <- as.character(pairs$gene_id)
if (!all(pairs_gene_ids %in% odm_gene_ids)) {
  stop("The gene IDs in the pairs data frame are not a subset of the gene IDs in the gene ondisc matrix.")
}

# 7. check that the grna groups in the pairs data frame are a subset of the grna groups in the grna ODM
odm_grna_groups <- grna_feature_covariates$grna_group
pairs_grna_groups <- as.character(pairs$grna_group)
if (!all(pairs_grna_groups %in% odm_grna_groups)) {
  stop("The grna groups in the pairs data frame are not a subset of the grna groups in the grna ondisc matrix.")
}


###############
# UPDATE MM ODM
###############
# 1. add threshold, B, and side to misc
mm_odm@modalities$grna@misc$threshold <- threshold
mm_odm@modalities$grna@misc$B <- B
mm_odm@modalities$grna@misc$side <- side
# 2. apply a formula object to the global cell covariates
if (form != "NA") {
  if (grepl("offset", form)) stop("Offsets are not currently supported in formulas.")
  form <- paste0(form, "+0")
  global_cell_covariates <- mm_odm |> get_cell_covariates()
  global_cell_covariates_new <- model.matrix(object = as.formula(form), data = global_cell_covariates) |> as.data.frame()
  mm_odm@global_cell_covariates <- global_cell_covariates_new
}

# Check for weird numbers in the global cell covariate matrix
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


#################################################################
# OUTPUT GENE IDS, GRNA GROUPS, PAIRS, AND UPDATED MULTIMODAL ODM
#################################################################
set.seed(4)
if (!(n_pairs_to_sample == 0)) pairs <- pairs |> dplyr::sample_n(n_pairs_to_sample)
pairs <- dplyr::arrange(pairs, gene_id, grna_group)

# obtain unique genes and grna groups
unique_genes <- data.frame(gene_id = unique(as.character(pairs$gene_id)))
unique_grna_groups <- data.frame(grna_group = unique(as.character(pairs$grna_group)))

# map the unique genes, unique grna groups, and gene-grna group pairs to their pod id
unique_genes$pod_id <- as.integer(cut(seq(1, nrow(unique_genes)), gene_pod_size))
unique_grna_groups$pod_id <- as.integer(cut(seq(1, nrow(unique_grna_groups)), grna_group_pod_size))
pairs$pod_id <- as.integer(cut(seq(1, nrow(pairs)), pair_pod_size))

# write vector
write_vector <- function(file_name, vector) {
  file_con <- file(file_name)
  writeLines(as.character(vector), file_con)
  close(file_con)
}

# save several outputs:
# (i) the gene to pod id map
# (ii) the grna group to pod id map
# (iii) the pair to pod id map
# (iv) the gene pod ids
# (v) the grna pod ids
# (vi) the pair pod ids
# (vii) the new multimodal matrix
saveRDS(unique_genes, "gene_to_pod_id_map.rds")
saveRDS(unique_grna_groups, "grna_group_to_pod_id_map.rds")
saveRDS(pairs, "pairs_to_pod_id_map.rds")
write_vector("gene_pods.txt", unique(unique_genes$pod_id))
write_vector("grna_group_pods.txt", unique(unique_grna_groups$pod_id))
write_vector("pair_pods.txt", unique(pairs$pod_id))
save_multimodal_odm(multimodal_odm = mm_odm, "mm_odm_new.rds")
