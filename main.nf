// Use DSL 2
nextflow.enable.dsl = 2

// Define optional parameters
params.formula = "NA"
params.threshold = 3
params.B = 1000
params.side = "both"
params.gene_pod_size = 50
params.grna_group_pod_size = 50
params.pair_pod_size = 100
params.n_pairs_to_sample = 0
params.result_fp = "$PWD/sceptre_result.rds"
params.gene_modality_name = "gene"
params.grna_modality_name = "grna"
params.inference_method = "crt" // one of `crt` and `gcm`

// Mild command line argument processing
File out_f = new File(params.result_fp)
result_file_name = out_f.getName()
result_dir = out_f.getAbsoluteFile().getParent()
formula = params.formula.replaceAll("\\(", "\\\\(").replaceAll("\\)", "\\\\)")

// PROCESS 1: Check inputs; output the list of gene IDs and grna groups
process check_inputs {
  time "5m"
  memory "4 GB"
  debug true

  input:
  path "multimodal_metadata_fp"
  path "gene_odm_fp"
  path "grna_odm_fp"
  path "pair_fp"

  output:
  path "gene_to_pod_id_map.rds", emit: gene_to_pod_id_map_ch
  path "grna_group_to_pod_id_map.rds", emit: grna_group_to_pod_id_map_ch
  path "pair_to_pod_id_map.rds", emit: pair_to_pod_id_map_ch
  path "gene_pods.txt", emit: gene_pods_ch
  path "grna_group_pods.txt", emit: grna_group_pods_ch
  path "pair_pods.txt", emit: pair_pods_ch
  path "mm_odm_new.rds", emit: multimodal_metadata_ch

  """
  check_inputs.R $multimodal_metadata_fp $gene_odm_fp $grna_odm_fp $pair_fp $formula $params.threshold $params.B $params.side $params.n_pairs_to_sample $params.gene_modality_name $params.grna_modality_name $params.gene_pod_size $params.grna_group_pod_size $params.pair_pod_size
  """
}

// PROCESS 2: Perform gene precomputation
process perform_gene_precomputation {
  time { 30.s * params.gene_pod_size }
  memory "2 GB"

  input:
  path "multimodal_metadata_fp"
  path "gene_odm_fp"
  path "gene_to_pod_id_map"
  val gene_pod_id

  output:
  path "precomp_sub_matrix.rds"

  """
  perform_precomputation.R "gene" $multimodal_metadata_fp $gene_odm_fp $gene_to_pod_id_map $gene_pod_id $params.threshold
  """
}

// PROCESS 3: Perform grna precomputation
process perform_grna_precomputation {
  time { 30.s * params.gene_pod_size }
  memory "2 GB"

  input:
  path "multimodal_metadata_fp"
  path "grna_odm_fp"
  path "grna_group_to_pod_id_map"
  val grna_group_pod_id

  output:
  path "precomp_sub_matrix.rds"

  """
  perform_precomputation.R "grna" $multimodal_metadata_fp $grna_odm_fp $grna_group_to_pod_id_map $grna_group_pod_id $params.threshold
  """
}


// PROCESS 6: Perform pairwise association test
process perform_pairwise_association_test {
  time { 30.s * params.pair_pod_size }
  memory "5 GB"

  input:
  path "multimodal_metadata_fp"
  path "gene_odm_fp"
  path "grna_odm_fp"
  path "gene_precomp_matrix_fp"
  path "grna_precomp_matrix_fp"
  path "pair_to_pod_id_map"
  val pair_pod_id

  output:
  path "raw_result.rds"

  """
  run_association_test.R $multimodal_metadata_fp $gene_odm_fp $grna_odm_fp $gene_precomp_matrix_fp $grna_precomp_matrix_fp $pair_to_pod_id_map $pair_pod_id $params.threshold $params.B $params.side
  """
}

// PROCESS 7: Combine results
process combine_results {
  debug true

  time "5m"
  memory "2 GB"

  publishDir result_dir, mode: "copy"

  output:
  path "$result_file_name"

  input:
  path "raw_result"
  path pair_fp

  """
  combine_results.R $result_file_name $pair_fp raw_result*
  """
}

// IMPORT SHARED PROCESSES FROM shared_processes.nf
include { join_precomputations as join_gene_precomputations } from "./shared_processes.nf"
include { join_precomputations as join_grna_precomputations } from "./shared_processes.nf"

// DEFINE WORKFLOW
workflow {
  // Step 1: Check inputs for correctness; output channels for gene IDs, grna groups, and pairs
  check_inputs(params.multimodal_metadata_fp,
               params.gene_odm_fp,
               params.grna_odm_fp,
               params.pair_fp)

  // Step 2: Clean up the output channels of the first process
  gene_pods = check_inputs.out.gene_pods_ch.splitText().map{it.trim()}
  grna_groups_pods = check_inputs.out.grna_group_pods_ch.splitText().map{it.trim()}
  pair_pods = check_inputs.out.pair_pods_ch.splitText().map{it.trim()}
  gene_to_pod_id_map = check_inputs.out.gene_to_pod_id_map_ch
  grna_group_to_pod_id_map = check_inputs.out.grna_group_to_pod_id_map_ch
  pair_to_pod_id_map_ch = check_inputs.out.pair_to_pod_id_map_ch
  multimodal_metadata_ch = check_inputs.out.multimodal_metadata_ch

  // Step 3: Perform the precomputation on genes
  perform_gene_precomputation(multimodal_metadata_ch,
                              params.gene_odm_fp,
                              gene_to_pod_id_map,
                              gene_pods)
  gene_precomp_ch_raw = perform_gene_precomputation.out

  // Step 4: Perform the precomputation on grna groups
  perform_grna_precomputation(
    multimodal_metadata_ch,
    params.grna_odm_fp,
    grna_group_to_pod_id_map,
    grna_groups_pods)
  grna_precomp_ch_raw = perform_grna_precomputation.out

  // Step 5: Combine the gene and grna group precomputations into a single .fst file
  join_gene_precomputations(gene_precomp_ch_raw.collect())
  gene_precomp_matrix = join_gene_precomputations.out
  join_grna_precomputations(grna_precomp_ch_raw.collect())
  grna_precomp_matrix = join_grna_precomputations.out

  // Step 6: Iterate over the pair pods, performing the pairwise association tests
  perform_pairwise_association_test(
    multimodal_metadata_ch,
    params.gene_odm_fp,
    params.grna_odm_fp,
    gene_precomp_matrix,
    grna_precomp_matrix,
    pair_to_pod_id_map_ch,
    pair_pods)

  // Step 7: Gather the results
  combine_results(perform_pairwise_association_test.out.collect(),
                  params.pair_fp)
}
