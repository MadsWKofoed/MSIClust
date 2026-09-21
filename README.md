# MSIClust

MSIClust is a spatially adaptive fuzzy clustering workflow for mass
spectrometry imaging (MSI) pixel data. It clusters spectra with
[`vsclust`](https://bioconductor.org/packages/vsclust) while allowing each pixel
to receive a local fuzzifier based on its neighborhood correlation. Pixels that
disagree with their neighbors are clustered more softly, and pixels whose best
cluster membership is too weak are labelled `No_cluster`.

![MSIClust cluster map of the synthetic example data](MSIClust_vignette_files/figure-gfm/synthetic-example-plot-1.png)

*MSIClust on the bundled synthetic example: four spatial regions recovered from
the spectra alone, with ambiguous border pixels left as `No_cluster` (grey).*

## Installation

MSIClust is a set of R functions in a single script, so there is nothing to
install beyond its dependencies. Clone or download this repository and start R
**from the repository root**, because the examples use relative paths.

```r
install.packages(c("BiocManager", "ggplot2", "matrixStats", "RColorBrewer"))
BiocManager::install(c("vsclust", "Cardinal"))
```

| Package | Needed for |
|---|---|
| `vsclust` | Required. Runs the fuzzy clustering. |
| `ggplot2` | Required for `plot_msiclust_map()`. |
| `Cardinal` | Reading imzML files and `make_msi_dataframe()`. Attach it with `library(Cardinal)` before calling `make_msi_dataframe()`. Not needed if you already have a pixel data frame. |
| `matrixStats`, `RColorBrewer` | Optional. Faster row maxima and nicer cluster colors. |
| `CardinalIO` | Only to regenerate the synthetic example data. |

Tested with R 4.4.2, `vsclust` 1.8.0, `Cardinal` 3.8.3 and `ggplot2` 3.5.1.

## Quick Start: Run the Bundled Example

This runs MSIClust on the synthetic imzML dataset in
[`example_data/Synthetic_MSIClust`](example_data/Synthetic_MSIClust) (45 x 35
pixels, 120 m/z features, four simulated regions). It takes a few seconds.

```r
source("MSIClust_helpers.R")
library(Cardinal)   # must be attached (not just Cardinal::) for make_msi_dataframe()

msi <- readImzML(
  "example_data/Synthetic_MSIClust/synthetic_msiclust.imzML",
  memory = FALSE,
  check = FALSE
)
msi_df <- make_msi_dataframe(msi)   # one row per pixel: x, y, mz_... columns

prep <- prepare_msiclust_input(
  msi_df = msi_df,
  normalize_method = "tic",
  feature_standardize = "sd",
  neighbor_radius = 1,
  cor_cores = 1,
  cor_scale = 25
)

msiclust_res <- run_msiclust(
  prep = prep,
  nclust = 4,
  iter_max = 100,
  min_membership = 0.5,
  seed = 1
)

clustered_df <- append_msiclust_labels(prep, msiclust_res)
table(clustered_df$MSIClust_cluster)
plot_msiclust_map(clustered_df, cluster_col = "MSIClust_cluster")
```

With `seed = 1` you should get 460 / 431 / 462 / 142 pixels in clusters 1-4 and
80 pixels labelled `No_cluster`.

## Using Your Own Data

MSIClust needs a data frame with **one row per pixel**:

| Column | Required | Description |
|---|:---:|---|
| `x`, `y` | yes | Integer pixel coordinates on a regular grid. |
| `mz_...` | yes | Numeric intensity columns, one per m/z feature (for example `mz_369.3`). |
| `runNames` | no | Optional run or slide identifier. |

If you have a binned Cardinal object, `make_msi_dataframe()` builds this table
for you. Then use the same four calls as above (`prepare_msiclust_input()`,
`run_msiclust()`, `append_msiclust_labels()`, `plot_msiclust_map()`), or the
one-call wrapper:

```r
res <- run_msiclust_workflow(msi_df, nclust = 3, seed = 1)
res$clustered_df   # pixel table with an MSIClust_cluster column
```

Neighbor correlations are computed on the pixel grid, so all pixels of one
image should share a coordinate system. For multi-run data, cluster one run at a
time or offset the coordinates so runs do not overlap.

## Key Parameters

| Parameter | Function | Meaning | Starting value |
|---|---|---|---|
| `nclust` | `run_msiclust()` | Number of clusters. | Dataset dependent |
| `normalize_method` | `prepare_msiclust_input()` | Pixel normalization: `"tic"`, `"median"` or `"rms"`. | `"tic"` |
| `feature_standardize` | `prepare_msiclust_input()` | Feature scaling: `"none"`, `"sd"` or `"zscore"`. | `"sd"` |
| `neighbor_radius` | `prepare_msiclust_input()` | Neighborhood radius (in pixels) for local correlation. | `1` |
| `cor_scale` | `prepare_msiclust_input()` | Scales `1 - correlation` before the fuzzifier is computed. | `25` |
| `min_membership` | `run_msiclust()` | A pixel needs a maximum membership above this value, otherwise it becomes `No_cluster`. | `0.5` |
| `iter_max` | `run_msiclust()` | Maximum clustering iterations. | `100` |
| `cor_cores` | `prepare_msiclust_input()` | Parallel workers for the neighbor correlations, the slowest step on large images. | `1` |

## What MSIClust Returns

`run_msiclust()` returns a list with, among others, `labels` (final labels
including `No_cluster`), `membership` (pixel-by-cluster fuzzy membership
matrix), `max_membership`, `centers`, and `model$fuzzifier` (the per-pixel
fuzzifier). `append_msiclust_labels()` adds the labels back to the pixel table
as a new column. See the vignette for the full list.

## Documentation

- [MSIClust_vignette.md](MSIClust_vignette.md): GitHub-rendered walkthrough with
  example output and figures. Start here for a detailed explanation.
- [MSIClust_vignette.Rmd](MSIClust_vignette.Rmd): Source of the vignette.
- [MSIClust_function_map.md](MSIClust_function_map.md): Short reference for the
  helper functions.
- [MSIClust_helpers.R](MSIClust_helpers.R): The complete implementation.
- [example_data/Synthetic_MSIClust](example_data/Synthetic_MSIClust): Synthetic
  imzML example data and the script that generates it.

## Citation and License

MSIClust builds on the variance-sensitive clustering method in `vsclust`; please
cite it if you use this workflow:

> Schwaemmle V (2024). *vsclust: Feature-based variance-sensitive quantitative
> clustering.* R package. <https://doi.org/10.18129/B9.bioc.vsclust>

This repository does not yet include a license file, so all rights are reserved
by default. Please contact the repository owner before reusing the code.
