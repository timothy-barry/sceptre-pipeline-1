#!/bin/bash

# Limit NF driver to 2 GB memory
export NXF_OPTS="-Xms500M -Xmx2G"

#############
# INPUT FILES
#############
source ~/.research_config
gasp_data_dir=$LOCAL_SCEPTRE2_DATA_DIR"data/gasperini/at_scale/"
# i) multimodal metadata file
multimodal_metadata_fp=$gasp_data_dir"multimodal_metadata.rds"
# ii) gene ODM
gene_odm_fp=$gasp_data_dir"gene/matrix.odm"
# iii) grna ODM
grna_odm_fp=$gasp_data_dir"grna_expression/matrix.odm"
# iv) gene-grna group pairs
pair_fp=$gasp_data_dir"neg_control_pairs.rds"

##############
# OUTPUT FILE:
##############
result_fp=$PWD"/sceptre_result.rds"

###############
# OPTIONAL ARGS
###############
# formula, threshold, B, side, n_pairs_to_sample, gene_pod_size, grna_group_pod_size, pair_pod_size are optional args
formula="~p_mito+batch+log(gene_n_nonzero)+log(gene_n_umis)+log(grna_expression_n_nonzero)+log(grna_expression_n_umis)"
gene_pod_size=5
grna_group_pod_size=5
pair_pod_size=10
n_pairs_to_sample=20

########################
# invoke the NF pipeline
########################
# nextflow pull timothy-barry/sceptre-pipeline
# nextflow run timothy-barry/sceptre-pipeline -r main \
nextflow run main.nf \
 --multimodal_metadata_fp $multimodal_metadata_fp \
 --gene_odm_fp $gene_odm_fp \
 --grna_odm_fp $grna_odm_fp \
 --pair_fp $pair_fp \
 --result_fp $result_fp \
 --formula $formula \
 --gene_pod_size $gene_pod_size \
 --grna_group_pod_size $grna_group_pod_size \
 --pair_pod_size $pair_pod_size \
 --grna_modality_name "grna_expression" \
 --n_pairs_to_sample $n_pairs_to_sample \
 --full_output "true" \
 --result_fp "$PWD/sceptre_result_crt.rds" \
 --inference_method "crt" \
 -resume
