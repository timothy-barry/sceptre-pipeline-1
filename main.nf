process check_inputs {
  debug true

  input:
  path multimodal_metadata_fp from params.multimodal_metadata_fp
  path gene_odm_fp from params.gene_odm_fp
  path gRNA_odm_fp from params.gRNA_odm_fp
  path pair_fp from params.pair_fp

  //output:
  //path "gRNA_groups.txt" into gRNA_groups_raw_ch
  //path "gene_ids.txt" into gene_ids raw_ch

  """
  check_inputs.R $multimodal_metadata_fp $gene_odm_fp $gRNA_odm_fp $pair_fp
  """
}
