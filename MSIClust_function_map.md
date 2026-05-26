# MSIClust Core Function Map

This is a compact lookup for the standalone MSIClust helper script:
`MSIClust_helpers.R`.

## Main Entry Points

| Function | Role |
|---|---|
| `prepare_msiclust_input()` | Prepare an MSI pixel data frame for MSIClust. |
| `run_msiclust()` | Run MSIClust on a prepared object. |
| `append_msiclust_labels()` | Add MSIClust labels back to the prepared pixel table. |
| `run_msiclust_workflow()` | One-call wrapper that prepares, clusters, and returns `clustered_df`. |
| `plot_msiclust_map()` | Plot spatial MSIClust labels. |

## Input Helpers

| Function | Role |
|---|---|
| `make_msi_dataframe()` | Convert a binned Cardinal object into `runNames`, `x`, `y`, and `mz_...` columns. |
| `msiclust_signal_columns()` | Locate the `mz_` feature columns. |
| `msiclust_feature_mz_values()` | Parse numeric m/z values from feature names. |
| `msiclust_make_unique_signal_columns()` | Make duplicated `mz_` names unique before clustering. |

## Preparation Helpers

| Function | Role |
|---|---|
| `msiclust_compute_neighbor_cor()` | Compute local spectral correlation between each pixel and its neighbors. |
| `msiclust_normalize_pixels()` | Apply pixel-wise TIC, median, or RMS normalization. |
| `msiclust_standardize_features()` | Apply optional feature-wise SD scaling or z-scoring. |

## MSIClust Internals

| Object or step | Meaning |
|---|---|
| `avg_corr_neighbors` | Average local spectral correlation for each pixel. |
| `inv_cor` | `1 - avg_corr_neighbors`. |
| `inv_cor_scaled` | `inv_cor * cor_scale`; input to `determine_fuzz()`. |
| `vsclust::determine_fuzz()` | Converts scaled inverse correlations into per-pixel fuzzifiers. |
| `vsclust::vsclust_algorithm()` | Runs fuzzy clustering using the per-pixel fuzzifier vector. |
| `msiclust_labels_from_membership()` | Applies `min_membership` and creates `No_cluster` labels. |

## Plotting Helpers

| Function | Role |
|---|---|
| `msiclust_cluster_palette()` | Build cluster colors while reserving grey for `No_cluster`. |
| `plot_msiclust_map()` | Draw the spatial cluster raster. |

## Compatibility Aliases

`MSIClust_helpers.R` also defines aliases matching names used in older local
scripts, such as `prepare_msi_clustering_input()`, `run_msiclust_method()`,
`plot_cluster_map()`, `compute_neighbor_cor()`, and `normalize_pixels()`.

