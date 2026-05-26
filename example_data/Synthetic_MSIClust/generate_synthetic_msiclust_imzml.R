#!/usr/bin/env Rscript

# Regenerate the synthetic MSIClust example imzML/ibd pair.
# The synthetic region labels are saved only to document how the data were made;
# MSIClust does not use them.

if (!requireNamespace("CardinalIO", quietly = TRUE)) {
  stop("Package 'CardinalIO' is required to write the synthetic imzML file.", call. = FALSE)
}

args <- commandArgs(trailingOnly = TRUE)
cmd_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", cmd_args, value = TRUE)
script_dir <- if (length(file_arg) == 1L) {
  dirname(normalizePath(sub("^--file=", "", file_arg), mustWork = FALSE))
} else {
  getwd()
}

output_dir <- if (length(args) >= 1L) {
  normalizePath(args[[1]], mustWork = FALSE)
} else {
  script_dir
}

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

set.seed(20260524)

stem <- "synthetic_msiclust"
nx <- 45L
ny <- 35L
nmz <- 120L

positions <- expand.grid(
  x = seq_len(nx),
  y = seq_len(ny)
)

mz <- seq(400, 900, length.out = nmz)

peak <- function(center, width, height) {
  height * exp(-0.5 * ((mz - center) / width)^2)
}

shared_profile <- peak(520, 11, 18) + peak(705, 15, 14)
matrix_profile <- peak(430, 12, 24) + peak(790, 15, 20)
region_a_profile <- peak(470, 10, 85) + peak(615, 12, 70) + peak(760, 14, 55)
region_b_profile <- peak(505, 9, 80) + peak(675, 11, 78) + peak(835, 16, 62)
region_c_profile <- peak(585, 12, 92) + peak(725, 13, 58) + peak(875, 14, 74)
baseline <- 3 + 0.004 * (mz - min(mz))

x <- positions$x
y <- positions$y

w_a <- exp(-0.5 * (((x - 14) / 7)^2 + ((y - 15) / 9)^2))
w_b <- exp(-0.5 * (((x - 33) / 7)^2 + ((y - 13) / 8)^2))
w_c <- exp(-0.5 * ((y - 28) / 4.5)^2) *
  stats::plogis((x - 8) / 2) *
  stats::plogis((39 - x) / 2)
w_matrix <- rep(0.35, length(x))

weights <- cbind(
  matrix = w_matrix,
  region_a = w_a,
  region_b = w_b,
  region_c = w_c
)
weights <- weights / rowSums(weights)

synthetic_region <- colnames(weights)[max.col(weights, ties.method = "first")]
synthetic_region <- factor(
  synthetic_region,
  levels = c("matrix", "region_a", "region_b", "region_c")
)

intensity <- matrix(NA_real_, nrow = nmz, ncol = nrow(positions))
for (i in seq_len(nrow(positions))) {
  profile <- baseline + shared_profile +
    weights[i, "matrix"] * matrix_profile +
    weights[i, "region_a"] * region_a_profile +
    weights[i, "region_b"] * region_b_profile +
    weights[i, "region_c"] * region_c_profile

  spatial_gradient <- 1 + 0.10 * (x[i] / nx) - 0.06 * (y[i] / ny)
  pixel_scale <- stats::rlnorm(1L, meanlog = 0, sdlog = 0.05)
  noise <- stats::rnorm(nmz, mean = 0, sd = 1.2)
  intensity[, i] <- pmax(0, profile * spatial_gradient * pixel_scale + noise)
}

meta <- CardinalIO::ImzMeta(
  spectrumType = "MS1 spectrum",
  spectrumRepresentation = "profile",
  instrumentModel = "LTQ FT Ultra",
  ionSource = "matrix-assisted laser desorption ionization",
  analyzer = "time-of-flight",
  detectorType = "microchannel plate detector"
)

imzml_path <- file.path(output_dir, paste0(stem, ".imzML"))
ibd_path <- file.path(output_dir, paste0(stem, ".ibd"))
design_path <- file.path(output_dir, paste0(stem, "_design.csv"))
summary_path <- file.path(output_dir, paste0(stem, "_summary.txt"))

unlink(c(imzml_path, ibd_path, design_path, summary_path), force = TRUE)

CardinalIO::writeImzML(
  meta,
  file = imzml_path,
  positions = positions,
  mz = mz,
  intensity = intensity
)

design <- data.frame(
  x = positions$x,
  y = positions$y,
  synthetic_region = synthetic_region,
  stringsAsFactors = FALSE
)

utils::write.csv(design, design_path, row.names = FALSE)

summary_lines <- c(
  "Synthetic MSIClust imzML example",
  paste("seed:", 20260524),
  paste("grid:", paste0(nx, "x", ny)),
  paste("pixels:", nrow(positions)),
  paste("m/z features:", length(mz)),
  paste("m/z range:", paste(round(range(mz), 4), collapse = " - ")),
  "regions: matrix, region_a, region_b, region_c",
  "note: synthetic_region labels document the simulated design and are not used by MSIClust"
)

writeLines(summary_lines, summary_path)

message("Wrote: ", imzml_path)
message("Wrote: ", ibd_path)
message("Wrote: ", design_path)
message("Wrote: ", summary_path)
