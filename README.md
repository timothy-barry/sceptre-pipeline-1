
Tim Barry

<!-- README.md is generated from README.Rmd. Please edit that file -->

# `sceptre` Nextflow pipeline

This repository contains the `sceptre` Nextflow pipeline. The `sceptre`
Nextflow pipeline is a command line utility that facilitates running
`sceptre` (i) out-of-core on a laptop or desktop or (ii) in a
distributed fashion on a cluster or cloud. The `sceptre` Nextflow
pipeline is highly scalable and memory-efficient; therefore, we
recommend using the `sceptre` Nextflow pipeline instead of the R
interface to `sceptre` in most cases.

# Requirements

-   Install [Nextflow](https://www.nextflow.io/)
-   Install the `sceptre` and `ondisc` packages:

<!-- -->

    # in R
    install.packages("devtools")
    setRepositories(ind = 1:4)
    devtools::install_github("katsevich-lab/sceptre")
    devtools::install_github("timothy-barry/ondisc")

-   Install the the `sceptre` Nextflow pipeline.

<!-- -->

    # on command line 
    nextflow pull timothy-barry/sceptre-pipeline

# Pipeline arguments

We describe the `sceptre` pipeline arguments here. As an example we use
data from the
[paper](https://www.sciencedirect.com/science/article/pii/S009286741831554X)
“A genome-wide framework for mapping gene regulation via cellular
genetic screens” by Gasperini et al., 2019.

## Input files

Four input files are required: (i) the multimodal ondisc matrix metadata
file (`multimodal_metadata_fp`), (ii) the backing .odm file of the gene
ondisc matrix (`gene_odm_fp`), (iii) the backing .odm file of the gRNA
ondisc matrix (`grna_odm_fp`), and (iv) a data frame containing the set
of gene-gRNA group pairs to analyze (`pair_fp`). On my (Tim’s) machine
these files are located in the following places:

``` r
# in R
gasp_offsite_dir <- paste0(.get_config_path("LOCAL_SCEPTRE2_DATA_DIR"), "data/gasperini/at_scale/")
multimodal_metadata_fp <- paste0(gasp_offsite_dir, "multimodal_metadata.rds")
gene_odm_fp <- paste0(gasp_offsite_dir, "gene/matrix.odm")
grna_odm_fp <- paste0(gasp_offsite_dir, "grna_expression/matrix.odm")
pair_fp <- paste0(gasp_offsite_dir, "neg_control_pairs.rds")
```

The multimodal ondisc matrix should satisfy the following conditions.

-   The multimodal ondisc matrix should have modalities for the gene and
    gRNA data.

``` r
# in R
library(ondisc)
crispr_experiment <- read_multimodal_odm(odm_fps = c(gene_odm_fp, grna_odm_fp),
                                         multimodal_metadata_fp = multimodal_metadata_fp)
gene_modality <- get_modality(crispr_experiment, "gene")
grna_modality <- get_modality(crispr_experiment, "grna_expression")
```

-   The gRNA modality should be an integer-valued matrix of gRNA
    expressions or (less commonly) a logical matrix of gRNA-to-cell
    assignments. The feature covariate matrix of the gRNA modality
    should contain a column called `grna_group` indicating the “group”
    to which each gRNA belongs. Typically, targeting gRNAs are grouped
    according to the site that they target, and non-targeting gRNAs are
    grouped randomly into sets of size two or three.

``` r
# in R
grna_modality
#> A covariate_ondisc_matrix with the following components:
#>  An integer-valued ondisc_matrix with 13182 features and 206271 cells.
#>  A cell covariate matrix with columns n_nonzero, n_umis.
#>  A feature covariate matrix with columns mean_expression, n_nonzero, grna_group, target_type, target, target_gene.
grna_modality |>
  get_feature_covariates() |>
  head()
#>                      mean_expression n_nonzero   grna_group target_type
#> AAACCGCTCCCGAGCACGGG      0.08524339      1450 SH3BGRL3_TSS    gene_tss
#> AAATAGTGGGAAGATTCGTG      0.02863151       556 MTRNR2L8_TSS    gene_tss
#> AACACACCACGGAGGAGTGG      0.06421350       987   FAM83A_TSS    gene_tss
#> AACAGCCCGGCCGGCCAAGG      0.07196465      1192   ZNF593_TSS    gene_tss
#> AACGAGAGACTGCTTGCTGG      0.03283749       693   ATPIF1_TSS    gene_tss
#> AACGGCTCGGAAGCCTAGGG      0.07587158      1251    TIPRL_TSS    gene_tss
#>                                        target     target_gene
#> AAACCGCTCCCGAGCACGGG   chr1:26605667-26605668 ENSG00000142669
#> AAATAGTGGGAAGATTCGTG  chr11:10530735-10530736 ENSG00000255823
#> AACACACCACGGAGGAGTGG chr8:124191287-124191288 ENSG00000147689
#> AACAGCCCGGCCGGCCAAGG   chr1:26496362-26496363 ENSG00000142684
#> AACGAGAGACTGCTTGCTGG   chr1:28562620-28562621 ENSG00000130770
#> AACGGCTCGGAAGCCTAGGG chr1:168148171-168148172 ENSG00000143155
```

-   The gene modality should be an integer-valued matrix of gene
    expressions.

``` r
gene_modality
#> A covariate_ondisc_matrix with the following components:
#>  An integer-valued ondisc_matrix with 12129 features and 206271 cells.
#>  A cell covariate matrix with columns n_nonzero, n_umis, p_mito, batch.
#>  A feature covariate matrix with columns mean_expression, coef_of_variation, n_nonzero.
```

Next, the data frame containing the pairs to analyze should contain
columns `gene_id` and `grna_group`. Additional columns are permitted but
ignored.

``` r
# in R
pairs_df <- readRDS(pair_fp)
head(pairs_df)
#>           gene_id grna_group
#> 1 ENSG00000238009  random_24
#> 2 ENSG00000237683  random_24
#> 3 ENSG00000228463  random_24
#> 4 ENSG00000237094  random_24
#> 5 ENSG00000235373  random_24
#> 6 ENSG00000228327  random_24
```

The gene IDs within the `gene_id` column should be a subset of the
feature IDs of the gene modality; meanwhile, the grna groups within the
`grna_group` column should be a subset of the entries of the
`grna_group` column of the feature covariate matrix of the gRNA
modality.

``` r
# in R
all(pairs_df$gene_id %in% get_feature_ids(gene_modality)) # gene ID check
#> [1] TRUE
all(pairs_df$grna_group  %in% get_feature_covariates(grna_modality)$grna_group) # gRNA group check
#> [1] TRUE
```

## Output file

The output file path (`result_fp`) specifies the location to write the
results.

## Additional arguments

The `sceptre` pipeline accepts several additional arguments, ordered
here from most important to least important.

-   `formula`: an R formula (stored as a string) specifying the
    covariates to adjust for in the analysis. The variables in the
    string are assumed to be columns of the global cell covariate
    matrix. An example formula is as follows:

<!-- -->

    # in bash
    formula="~gene_p_mito+gene_batch+log(gene_n_nonzero)+log(gene_n_umis)+log(grna_n_nonzero)+log(grna_n_umis)"

The variables in this formula (`gene_p_mito`, `gene_batch`,
`gene_n_nonzero`, `gene_n_umis`, `grna_n_nonzero`, and `grna_n_umis`)
are columns of the global cell covariate matrix:

``` r
# in R
crispr_experiment |> get_cell_covariates() |> head() 
#>                                  gene_n_nonzero gene_n_umis
#> AAACCTGAGAGGTACC-1_1A_1_SI-GA-E2           3549       17566
#> AAACCTGAGTCAATAG-1_1A_1_SI-GA-E2           2543        8917
#> AAACCTGCAAACAACA-1_1A_1_SI-GA-E2           3191       14626
#> AAACCTGCACTTCTGC-1_1A_1_SI-GA-E2           4539       22783
#> AAACCTGCATGTAGTC-1_1A_1_SI-GA-E2           2605       10124
#> AAACCTGGTAGCGCAA-1_1A_1_SI-GA-E2           2187        9743
#>                                  grna_assignment_n_nonzero
#> AAACCTGAGAGGTACC-1_1A_1_SI-GA-E2                        68
#> AAACCTGAGTCAATAG-1_1A_1_SI-GA-E2                        31
#> AAACCTGCAAACAACA-1_1A_1_SI-GA-E2                        63
#> AAACCTGCACTTCTGC-1_1A_1_SI-GA-E2                        41
#> AAACCTGCATGTAGTC-1_1A_1_SI-GA-E2                        38
#> AAACCTGGTAGCGCAA-1_1A_1_SI-GA-E2                        58
#>                                  grna_expression_n_nonzero
#> AAACCTGAGAGGTACC-1_1A_1_SI-GA-E2                        84
#> AAACCTGAGTCAATAG-1_1A_1_SI-GA-E2                        44
#> AAACCTGCAAACAACA-1_1A_1_SI-GA-E2                        82
#> AAACCTGCACTTCTGC-1_1A_1_SI-GA-E2                        56
#> AAACCTGCATGTAGTC-1_1A_1_SI-GA-E2                        71
#> AAACCTGGTAGCGCAA-1_1A_1_SI-GA-E2                        79
#>                                  grna_expression_n_umis        batch
#> AAACCTGAGAGGTACC-1_1A_1_SI-GA-E2                    994 prep_batch_1
#> AAACCTGAGTCAATAG-1_1A_1_SI-GA-E2                    347 prep_batch_1
#> AAACCTGCAAACAACA-1_1A_1_SI-GA-E2                    930 prep_batch_1
#> AAACCTGCACTTCTGC-1_1A_1_SI-GA-E2                    579 prep_batch_1
#> AAACCTGCATGTAGTC-1_1A_1_SI-GA-E2                   1098 prep_batch_1
#> AAACCTGGTAGCGCAA-1_1A_1_SI-GA-E2                   1276 prep_batch_1
#>                                       p_mito
#> AAACCTGAGAGGTACC-1_1A_1_SI-GA-E2 0.058786706
#> AAACCTGAGTCAATAG-1_1A_1_SI-GA-E2 0.036086518
#> AAACCTGCAAACAACA-1_1A_1_SI-GA-E2 0.069823051
#> AAACCTGCACTTCTGC-1_1A_1_SI-GA-E2 0.026186508
#> AAACCTGCATGTAGTC-1_1A_1_SI-GA-E2 0.007991318
#> AAACCTGGTAGCGCAA-1_1A_1_SI-GA-E2 0.022356681
```

The default behavior is to adjust for all (untransformed) variables
stored within the global cell covariate matrix.

**Note**: The formula string should contain **no** white spaces (e.g.,
spaces, tabs, etc.).

-   `gene_modality_name`: the name of the gene modality within the
    multimodal ondisc matrix (“gene”, in the example above)

-   `grna_modality_name`: the name of the gRNA modality within the
    multimodal ondisc matrix (“grna_expression”, in the example above)

-   `threshold`: the threshold to use to assign gRNAs to cells. For a
    given cell and gRNA, if the UMI count of the gRNA within the cell is
    equal to or greater than the threshold, then the gRNA is taken to be
    present in the cell. The default is 3.

-   `B`: the number of resamples to draw for the conditional
    randomization test. The default is 1000.

-   `side`: the sidedness of the test, one of “left”, “right”, or
    “both”. The default is “both”.

-   `gene_pod_size`, `grna_group_pod_size`, and `pair_pod_size`:
    parameters that control the amount of parallelization. At a high
    level the pipeline works as follows: first, all genes are regressed
    onto the technical factors; next, all gRNA groups are regressed onto
    the technical factors; finally, all pairs of genes and grna groups
    (as specified in the pairs data frame) are tested for association.
    `gene_pod_size` (resp., `grna_group_pod_size`) is the number of
    genes (resp., gRNA groups) to regress onto the technical factors in
    a given Nextflow process. Meanwhile, `pair_pod_size` is the number
    of gene-gRNA group pairs to test for association in a given Nextflow
    process. The default values are 200, 200, and 5000, respectively.

-   `n_pairs_to_sample`: the number of randomly-selected gene-gRNA group
    pairs on which to run `sceptre`. This parameter is for debugging and
    testing purposes; often, it is useful to run the pipeline on (say)
    25 randomly-selected pairs to ensure that the pipeline is set up
    correctly. The default is to run the pipeline on the entire set of
    pairs, i.e., to not subsample at all.

-   `full_output`: “true” or “false” (default “false”); output the
    resampled test statistics alongside the p-value, z-score, and log
    fold change?

-   `inference_method`: (**EXPERIMENTAL**) “crt” or “gcm” (default
    “crt”); perform inference using the vanilla conditional
    randomization test (“crt”) or the faster (but possibly less
    accurate) generalized covariance measure (“gcm”)?

## Invoking the pipeline

An example invocation script is available in the `sceptre` Nextflow
repository. Git clone the repository.

    # on command line
    git clone git@github.com:timothy-barry/sceptre-pipeline.git
    cd sceptre-pipeline

The bash script `example_launch.sh` contains an example invocation of
the pipeline. We can launch this script via a call to `bash` on a
laptop/desktop, `qsub` on a cluster running a Sun Gride Engine
scheduler, `sbatch` on a cluster running a SLURM scheduler, etc.

    # on command line
    bash example_launch.sh # laptop
    qsub example_launch.sh # sun grid engine
    sbatch example_launch.sh # slurm

# Interpreting the results

The results data frame contains columns `gene_id`, `grna_group`,
`p_value`, `z_value`, and `log_fold_change`. The results can be analyzed
using the tips outlined in Step 7 (“Analyze the Results”) of the
[in-memory sceptre
tutorial](https://katsevich-lab.github.io/sceptre/articles/using_sceptre_v2.html).
