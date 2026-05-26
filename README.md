# MSIClust

MSIClust is a spatially adaptive fuzzy clustering workflow for mass
spectrometry imaging (MSI) pixel data. It clusters spectra with `vsclust` while
allowing each pixel to receive a local fuzzifier based on its neighborhood
correlation.

## Start Here

- [MSIClust_vignette.md](MSIClust_vignette.md): GitHub-rendered walkthrough with
  example output and figures.
- [MSIClust_vignette.Rmd](MSIClust_vignette.Rmd): Source vignette.
- [MSIClust_function_map.md](MSIClust_function_map.md): Short reference for the
  helper functions.
- [MSIClust_helpers.R](MSIClust_helpers.R): Standalone implementation used by
  the vignette.
- [example_data/Synthetic_MSIClust](example_data/Synthetic_MSIClust): Synthetic
  imzML example data used by the vignette.

## Minimal Workflow

```r
source("MSIClust_helpers.R")

prep <- prepare_msiclust_input(
  msi_df = msi_df,
  normalize_method = "tic",
  feature_standardize = "sd",
  neighbor_radius = 1,
  cor_cores = 8,
  cor_scale = 25
)

msiclust_res <- run_msiclust(
  prep = prep,
  nclust = 3,
  iter_max = 100,
  min_membership = 0.5,
  seed = 1
)

clustered_df <- append_msiclust_labels(prep, msiclust_res)
plot_msiclust_map(clustered_df, cluster_col = "MSIClust_cluster")
```

Input data should contain one row per pixel with `x`, `y`, and numeric `mz_...`
feature columns.
