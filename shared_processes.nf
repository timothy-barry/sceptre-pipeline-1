process join_precomputations {
  input:
  path "precomp_submatrix"

  output:
  path "precomputation_matrix.fst"

  """
  join_precomputations.R precomp_submatrix*
  """
}
