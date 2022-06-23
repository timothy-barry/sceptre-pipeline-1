#!/bin/bash

# Limit NF driver to 2 GB memory
export NXF_OPTS="-Xms500M -Xmx2G"

#############
# INPUT FILES
#############
source ~/.research_config
# i) multimodal metadata file
multimodal_metadata_fp=$LOCAL_GASPERINI_2019_DATA_DIR"at-scale/processed/multimodal/multimodal_metadata.rds"
# ii) gene ODM
gene_odm_fp=$LOCAL_GASPERINI_2019_DATA_DIR"at-scale/processed/gene/gasp_scale_gene_expressions.odm"
# iii) gRNA ODM
gRNA_odm_fp=$LOCAL_GASPERINI_2019_DATA_DIR"at-scale/processed/gRNA_ungrouped/gasp_scale_gRNA_counts_ungrouped.odm"
# iv) gene-gRNA group pairs
pair_fp=$LOCAL_GASPERINI_2019_DATA_DIR"at-scale/processed/multimodal/pairs.rds"

##############
# OUTPUT FILE:
##############
result_fp=$PWD"/sceptre_result.rds"

###############
# OPTIONAL ARGS
###############
# formula, threshold, B, side, n_pairs_to_sample, gene_pod_size, gRNA_group_pod_size, pair_pod_size are optional args
formula="~gene_p_mito+gene_batch+log(gene_n_nonzero)+log(gene_n_umis)+log(gRNA_n_nonzero)+log(gRNA_n_umis)"
gene_pod_size=5
gRNA_group_pod_size=5
pair_pod_size=10
n_pairs_to_sample=25

########################
# invoke the NF pipeline
########################
nextflow pull timothy-barry/sceptre-pipeline
nextflow run timothy-barry/sceptre-pipeline -r main \
 --multimodal_metadata_fp $multimodal_metadata_fp \
 --gene_odm_fp $gene_odm_fp \
 --gRNA_odm_fp $gRNA_odm_fp \
 --pair_fp $pair_fp \
 --result_fp $result_fp \
 --formula $formula \
 --gene_pod_size $gene_pod_size \
 --gRNA_group_pod_size $gRNA_group_pod_size \
 --pair_pod_size $pair_pod_size \
 --n_pairs_to_sample $n_pairs_to_sample \
