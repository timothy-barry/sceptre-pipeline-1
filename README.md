
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
    devtools::install_github("katsevich-lab/sceptre")
    devtools::install_github("timothy-barry/ondisc")

-   Install the the `sceptre` nextflow pipeline.

<!-- -->

    # on command line 
    nextflow pull timothy-barry/sceptre-pipeline

# Pipeline arguments

We describe the `sceptre` pipeline arguments here.

## Input files

Four input files are required: (i) the multimodal ondisc matrix metadata
file (`multimodal_metadata_fp`), (ii) the backing .odm file of the gene
ondisc matrix (`gene_odm_fp`), (iii) the backing .odm file of the gRNA
ondisc matrix (`gRNA_odm_fp`), and (iv) a data frame containing the set
of gene-gRNA group pairs to analyze (`pair_fp`). On my (Tim’s) machine
these files are located in the following places:

``` r
# in R
gasp_offsite_dir <- .get_config_path("LOCAL_GASPERINI_2019_DATA_DIR")
multimodal_metadata_fp <- paste0(gasp_offsite_dir, "at-scale/processed/multimodal/multimodal_metadata.rds")
gene_odm_fp <- paste0(gasp_offsite_dir, "at-scale/processed/gene/gasp_scale_gene_expressions.odm")
gRNA_odm_fp <- paste0(gasp_offsite_dir, "at-scale/processed/gRNA_ungrouped/gasp_scale_gRNA_counts_ungrouped.odm")
pair_fp <- paste0(gasp_offsite_dir, "at-scale/processed/multimodal/pairs.rds")
```

The multimodal ondisc matrix should satisfy the following conditions.

-   The multimodal ondisc matrix should have modalities named “gene” and
    “gRNA”.

``` r
# in R
library(ondisc)
crispr_experiment <- read_multimodal_odm(odm_fps = c(gene_odm_fp, gRNA_odm_fp),
                                         multimodal_metadata_fp = multimodal_metadata_fp)
gene_modality <- get_modality(crispr_experiment, "gene")
gRNA_modality <- get_modality(crispr_experiment, "gRNA")
```

-   The gRNA modality should be an integer-valued matrix of gRNA
    expressions or (less commonly) a logical matrix of gRNA-to-cell
    assignments. The feature covariate matrix of the gRNA modality
    should contain a column called `gRNA_group` indicating the “group”
    to which each gRNA belongs. Typically, targeting gRNAs are grouped
    according to the site that they target, and non-targeting gRNAs are
    grouped randomly into sets of size two or three.

``` r
# in R
gRNA_modality
#> A covariate_ondisc_matrix with the following components:
#>  An integer-valued ondisc_matrix with 13189 features and 207320 cells.
#>  A cell covariate matrix with columns n_nonzero, n_umis.
#>  A feature covariate matrix with columns mean_expression, n_nonzero, gRNA_group.
gRNA_modality |>
  get_feature_covariates() |>
  head()
#>                      mean_expression n_nonzero   gRNA_group
#> AAACCGCTCCCGAGCACGGG      0.08524339      1450 SH3BGRL3_TSS
#> AAATAGTGGGAAGATTCGTG      0.02863151       556 MTRNR2L8_TSS
#> AACACACCACGGAGGAGTGG      0.06421350       987   FAM83A_TSS
#> AACAGCCCGGCCGGCCAAGG      0.07196465      1192   ZNF593_TSS
#> AACGAGAGACTGCTTGCTGG      0.03283749       693   ATPIF1_TSS
#> AACGGCTCGGAAGCCTAGGG      0.07587158      1251    TIPRL_TSS
```

-   The gene modality should be an integer-valued matrix of gene
    expressions.

``` r
gene_modality
#> A covariate_ondisc_matrix with the following components:
#>  An integer-valued ondisc_matrix with 13135 features and 207320 cells.
#>  A cell covariate matrix with columns n_nonzero, n_umis, p_mito, batch.
#>  A feature covariate matrix with columns mean_expression, coef_of_variation, n_nonzero.
```

Next, the data frame containing the pairs to analyze should contain
columns `gene_id` and `gRNA_group`. Additional columns are permitted but
ignored.

``` r
# in R
pairs_df <- readRDS(pair_fp)
head(pairs_df)
#>           gene_id gRNA_group site_type
#> 1 ENSG00000008256   ACTB_TSS       TSS
#> 2 ENSG00000011275   ACTB_TSS       TSS
#> 3 ENSG00000075618   ACTB_TSS       TSS
#> 4 ENSG00000075624   ACTB_TSS   selfTSS
#> 5 ENSG00000086232   ACTB_TSS       TSS
#> 6 ENSG00000106305   ACTB_TSS       TSS
```

The gene IDs within the `gene_id` column should be a subset of the
feature IDs of the gene modality; meanwhile, the gRNA groups within the
`gRNA_group` column should be a subset of the entries of the
`gRNA_group` column of the feature covariate matrix of the gRNA
modality.

``` r
# in R
all(pairs_df$gene_id %in% get_feature_ids(gene_modality)) # gene ID check
#> [1] TRUE
all(pairs_df$gRNA_group  %in% get_feature_covariates(gRNA_modality)$gRNA_group) # gRNA group check
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
    formula="~gene_p_mito+gene_batch+log(gene_n_nonzero)+log(gene_n_umis)+log(gRNA_n_nonzero)+log(gRNA_n_umis)"

The variables in this formula (`gene_p_mito`, `gene_batch`,
`gene_n_nonzero`, `gene_n_umis`, `gRNA_n_nonzero`, and `gRNA_n_umis`)
are columns of the global cell covariate matrix:

``` r
# in R
crispr_experiment |> get_cell_covariates() |> head() 
#>                                  gene_n_nonzero gene_n_umis gene_p_mito
#> AAACCTGAGAGGTACC-1_1A_1_SI-GA-E2           3549       17566 0.058786706
#> AAACCTGAGTCAATAG-1_1A_1_SI-GA-E2           2543        8917 0.036086518
#> AAACCTGCAAACAACA-1_1A_1_SI-GA-E2           3191       14626 0.069823051
#> AAACCTGCACTTCTGC-1_1A_1_SI-GA-E2           4539       22783 0.026186508
#> AAACCTGCATGTAGTC-1_1A_1_SI-GA-E2           2605       10124 0.007991318
#> AAACCTGGTAGCGCAA-1_1A_1_SI-GA-E2           2187        9743 0.022356681
#>                                    gene_batch gRNA_n_nonzero gRNA_n_umis
#> AAACCTGAGAGGTACC-1_1A_1_SI-GA-E2 prep_batch_1             84         994
#> AAACCTGAGTCAATAG-1_1A_1_SI-GA-E2 prep_batch_1             44         347
#> AAACCTGCAAACAACA-1_1A_1_SI-GA-E2 prep_batch_1             82         930
#> AAACCTGCACTTCTGC-1_1A_1_SI-GA-E2 prep_batch_1             56         579
#> AAACCTGCATGTAGTC-1_1A_1_SI-GA-E2 prep_batch_1             71        1098
#> AAACCTGGTAGCGCAA-1_1A_1_SI-GA-E2 prep_batch_1             79        1276
```

The default behavior is to adjust for all (untransformed) variables
stored within the global cell covariate matrix.

**Note**: The formula string should contain **no** white spaces (e.g.,
spaces, tabs, etc.).

-   `threshold`: the threshold to use to assign gRNAs to cells. For a
    given cell and gRNA, if the UMI count of the gRNA within the cell is
    equal to or greater than the threshold, then the gRNA is taken to be
    present in the cell. The default is 3.

-   `B`: the number of resamples to draw for the conditional
    randomization test. The default is 1000.

-   `side`: the sidedness of the test, one of “left”, “right”, or
    “both”. The default is “both”.

-   `n_pairs_to_sample`: the number of randomly-selected gene-gRNA group
    pairs on which to run `sceptre`. This parameter is for debugging
    purposes; often, it is useful to run the pipeline on (say) 25
    randomly-selected pairs to ensure that the pipeline is set up
    correctly. The default is to run the pipeline on the entireset of
    pairs, i.e., to not subsample at all.

-   `gene_pod_size`, `gRNA_group_pod_size`, and `pair_pod_size`:
    parameters that control the amount of parallelization. At a high
    level the pipeline works as follows: first, all genes are regressed
    onto the technical factors; next, all gRNA groups are regressed onto
    the technical factors; finally, all pairs of genes and gRNA groups
    (as specified in the pairs data frame) are tested for association.
    `gene_pod_size` (resp., `gRNA_group_pod_size`) is the number of
    genes (resp., gRNA groups) to regress onto the technical factors in
    a given Nextflow process. Meanwhile, `pair_pod_size` is the number
    of gene-gRNA group pairs to test for association in a given Nextflow
    process. The default values are 50, 50, and 100, respectively.

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

The results data frame contains the columns
