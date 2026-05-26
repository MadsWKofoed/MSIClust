#!/usr/bin/env Rscript

# Standalone helpers for running MSIClust on MSI pixel data frames.
# Expected input: one row per pixel, x/y coordinate columns, and numeric mz_*
# feature columns.

.msiclust_require <- function(pkg, purpose = NULL) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    msg <- paste0("Package '", pkg, "' is required")
    if (!is.null(purpose) && nzchar(purpose)) {
      msg <- paste0(msg, " ", purpose)
    }
    stop(msg, ".", call. = FALSE)
  }
  invisible(TRUE)
}

.msiclust_row_maxs <- function(x) {
  if (requireNamespace("matrixStats", quietly = TRUE)) {
    return(matrixStats::rowMaxs(x))
  }
  apply(x, 1L, max, na.rm = TRUE)
}

msiclust_signal_columns <- function(msi_df, pattern = "^mz_") {
  signal_cols <- grep(pattern, names(msi_df), value = TRUE)
  if (length(signal_cols) == 0L) {
    stop("No feature columns found. Expected columns named like 'mz_...'.", call. = FALSE)
  }
  signal_cols
}

msiclust_feature_mz_values <- function(feature_names) {
  clean_names <- sub("__dup[0-9]+$", "", feature_names)
  suppressWarnings(as.numeric(sub("^mz_", "", clean_names)))
}

msiclust_make_unique_signal_columns <- function(msi_df,
                                                signal_cols = msiclust_signal_columns(msi_df),
                                                pattern = "^mz_") {
  signal_idx <- which(names(msi_df) %in% signal_cols & grepl(pattern, names(msi_df)))
  if (length(signal_idx) == 0L) {
    stop("No feature columns found. Expected columns named like 'mz_...'.", call. = FALSE)
  }

  unique_signal_cols <- make.unique(names(msi_df)[signal_idx], sep = "__dup")
  names(msi_df)[signal_idx] <- unique_signal_cols

  list(
    data = msi_df,
    signal_cols = unique_signal_cols,
    signal_idx = signal_idx
  )
}

make_msi_dataframe <- function(msi_data_binned) {
  .msiclust_require("Cardinal", "to convert Cardinal MSI objects")

  msi_matrix <- t(as.matrix(Cardinal::spectra(msi_data_binned)))
  mz_names <- paste0("mz_", Cardinal::mz(msi_data_binned))
  coords <- as.data.frame(Cardinal::coord(msi_data_binned))

  run_name <- tryCatch(Cardinal::runNames(msi_data_binned), error = function(e) "run1")
  pixel_names <- if (length(run_name) == 1L) {
    rep(run_name, nrow(msi_matrix))
  } else if (length(run_name) == nrow(msi_matrix)) {
    run_name
  } else {
    rep(paste(run_name, collapse = ","), nrow(msi_matrix))
  }

  full_df <- data.frame(
    runNames = pixel_names,
    x = coords$x,
    y = coords$y,
    msi_matrix,
    check.names = FALSE
  )

  colnames(full_df) <- c("runNames", "x", "y", mz_names)
  full_df
}

msiclust_normalize_pixels <- function(data,
                                      signal_cols,
                                      spatial_cols = c("x", "y"),
                                      method = c("tic", "median", "rms"),
                                      na.rm = TRUE) {
  method <- match.arg(method)
  missing_cols <- setdiff(signal_cols, names(data))
  if (length(missing_cols) > 0L) {
    stop("Missing signal columns: ", paste(missing_cols, collapse = ", "), call. = FALSE)
  }

  X_signal <- as.matrix(data[, signal_cols, drop = FALSE])
  storage.mode(X_signal) <- "double"
  X_spatial <- data[, intersect(spatial_cols, names(data)), drop = FALSE]

  denom <- switch(
    method,
    tic = rowSums(X_signal, na.rm = na.rm),
    median = apply(X_signal, 1L, stats::median, na.rm = na.rm),
    rms = sqrt(rowMeans(X_signal^2, na.rm = na.rm))
  )

  denom[!is.finite(denom) | denom == 0] <- NA_real_
  X_signal_norm <- sweep(X_signal, 1L, denom, "/")
  cbind(X_spatial, as.data.frame(X_signal_norm, check.names = FALSE))
}

msiclust_standardize_features <- function(data,
                                          signal_cols,
                                          method = c("none", "sd", "zscore")) {
  method <- match.arg(method)
  missing_cols <- setdiff(signal_cols, names(data))
  if (length(missing_cols) > 0L) {
    stop("Missing signal columns: ", paste(missing_cols, collapse = ", "), call. = FALSE)
  }

  if (method == "none") {
    return(list(
      data = data,
      center = setNames(rep(0, length(signal_cols)), signal_cols),
      scale = setNames(rep(1, length(signal_cols)), signal_cols),
      zero_scale_features = character(0)
    ))
  }

  X_signal <- as.matrix(data[, signal_cols, drop = FALSE])
  storage.mode(X_signal) <- "double"

  center <- colMeans(X_signal, na.rm = TRUE)
  scale <- apply(X_signal, 2L, stats::sd, na.rm = TRUE)
  zero_scale_features <- names(scale)[is.na(scale) | scale == 0]
  scale[is.na(scale) | scale == 0] <- 1

  X_scaled <- if (method == "zscore") {
    sweep(X_signal, 2L, center, "-")
  } else {
    X_signal
  }
  X_scaled <- sweep(X_scaled, 2L, scale, "/")

  data[, signal_cols] <- as.data.frame(X_scaled, check.names = FALSE)

  list(
    data = data,
    center = if (method == "zscore") center else setNames(rep(0, length(signal_cols)), signal_cols),
    scale = scale,
    zero_scale_features = zero_scale_features
  )
}

msiclust_compute_neighbor_cor <- function(dat,
                                          x_col = "x",
                                          y_col = "y",
                                          mz_cols = NULL,
                                          r = 1L,
                                          cores = 1L,
                                          weighted = TRUE) {
  if (is.null(mz_cols)) {
    mz_cols <- msiclust_signal_columns(dat)
  }

  required_cols <- c(x_col, y_col, mz_cols)
  missing_cols <- setdiff(required_cols, names(dat))
  if (length(missing_cols) > 0L) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "), call. = FALSE)
  }

  r <- as.integer(r)
  cores <- as.integer(cores)
  if (length(r) != 1L || is.na(r) || r < 1L) {
    stop("r must be a positive integer.", call. = FALSE)
  }
  if (length(cores) != 1L || is.na(cores) || cores < 1L) {
    stop("cores must be a positive integer.", call. = FALSE)
  }

  xy <- as.matrix(dat[, c(x_col, y_col), drop = FALSE])
  intens <- as.matrix(dat[, mz_cols, drop = FALSE])
  storage.mode(intens) <- "double"

  key_vec <- paste(xy[, 1L], xy[, 2L], sep = "_")
  index_lookup <- setNames(seq_len(nrow(dat)), key_vec)

  neighbor_idx_fun <- function(i, xy, r, index_lookup) {
    x <- xy[i, 1L]
    y <- xy[i, 2L]

    grid <- expand.grid(
      xx = (x - r):(x + r),
      yy = (y - r):(y + r)
    )
    coords <- as.matrix(grid)
    coords <- coords[!(coords[, 1L] == x & coords[, 2L] == y), , drop = FALSE]

    idx <- index_lookup[paste(coords[, 1L], coords[, 2L], sep = "_")]
    ok <- !is.na(idx)
    if (!any(ok)) {
      return(list(idx = integer(0), step = integer(0)))
    }

    coords_ok <- coords[ok, , drop = FALSE]
    dx <- abs(coords_ok[, 1L] - x)
    dy <- abs(coords_ok[, 2L] - y)

    list(
      idx = as.integer(idx[ok]),
      step = pmax(dx, dy)
    )
  }

  pixel_cor_fun <- function(i, xy, intens, r, index_lookup, weighted) {
    nei <- neighbor_idx_fun(i, xy, r, index_lookup)
    if (length(nei$idx) == 0L) {
      return(NA_real_)
    }

    v <- intens[i, ]
    mats <- intens[nei$idx, , drop = FALSE]
    cors <- apply(mats, 1L, function(z) stats::cor(v, z, use = "pairwise.complete.obs"))
    ok <- !is.na(cors) & is.finite(cors)

    if (!any(ok)) {
      return(NA_real_)
    }

    if (isTRUE(weighted)) {
      weights <- 1 / nei$step
      ok <- ok & !is.na(weights) & is.finite(weights)
      if (!any(ok)) {
        return(NA_real_)
      }
      return(sum(weights[ok] * cors[ok]) / sum(weights[ok]))
    }

    mean(cors[ok])
  }

  if (cores > 1L) {
    cl <- parallel::makeCluster(cores)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    parallel::clusterExport(
      cl,
      varlist = c("xy", "intens", "index_lookup", "neighbor_idx_fun", "pixel_cor_fun", "r", "weighted"),
      envir = environment()
    )
    return(parallel::parSapply(
      cl,
      X = seq_len(nrow(dat)),
      FUN = pixel_cor_fun,
      xy = xy,
      intens = intens,
      r = r,
      index_lookup = index_lookup,
      weighted = weighted
    ))
  }

  vapply(
    seq_len(nrow(dat)),
    pixel_cor_fun,
    numeric(1L),
    xy = xy,
    intens = intens,
    r = r,
    index_lookup = index_lookup,
    weighted = weighted
  )
}

prepare_msiclust_input <- function(msi_df,
                                   signal_cols = msiclust_signal_columns(msi_df),
                                   x_col = "x",
                                   y_col = "y",
                                   normalize_method = c("tic", "median", "rms"),
                                   feature_standardize = c("none", "sd", "zscore"),
                                   neighbor_radius = 1L,
                                   cor_cores = 1L,
                                   cor_scale = 25,
                                   clamp_negative_cor = FALSE,
                                   weighted_neighbors = TRUE) {
  normalize_method <- match.arg(normalize_method)
  feature_standardize <- match.arg(feature_standardize)

  unique_signal <- msiclust_make_unique_signal_columns(msi_df, signal_cols)
  msi_df <- unique_signal$data
  signal_cols <- unique_signal$signal_cols

  required_cols <- c(x_col, y_col, signal_cols)
  missing_cols <- setdiff(required_cols, names(msi_df))
  if (length(missing_cols) > 0L) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "), call. = FALSE)
  }

  cor_data <- msi_df[, required_cols, drop = FALSE]
  cor_data$avg_corr_neighbors <- msiclust_compute_neighbor_cor(
    dat = cor_data,
    x_col = x_col,
    y_col = y_col,
    mz_cols = signal_cols,
    r = neighbor_radius,
    cores = cor_cores,
    weighted = weighted_neighbors
  )
  cor_data$avg_corr_neighbors[is.nan(cor_data$avg_corr_neighbors)] <- NA_real_

  if (isTRUE(clamp_negative_cor)) {
    cor_data$avg_corr_neighbors[cor_data$avg_corr_neighbors < 0] <- 0
  }

  cor_data$inv_cor <- 1 - cor_data$avg_corr_neighbors
  keep <- !is.na(cor_data$avg_corr_neighbors)

  if (!any(keep)) {
    stop("No pixels with valid neighboring-pixel correlations were found.", call. = FALSE)
  }

  cor_data <- cor_data[keep, , drop = FALSE]
  df_base <- msi_df[keep, , drop = FALSE]
  row.names(cor_data) <- NULL
  row.names(df_base) <- NULL

  normalized_df <- msiclust_normalize_pixels(
    data = cor_data,
    signal_cols = signal_cols,
    spatial_cols = c(x_col, y_col),
    method = normalize_method
  )

  standardized <- msiclust_standardize_features(
    data = normalized_df,
    signal_cols = signal_cols,
    method = feature_standardize
  )
  clustering_df <- standardized$data

  X <- as.matrix(clustering_df[, signal_cols, drop = FALSE])
  storage.mode(X) <- "double"

  list(
    df_base = df_base,
    cor_data = cor_data,
    normalized_df = normalized_df,
    clustering_df = clustering_df,
    X = X,
    signal_cols = signal_cols,
    feature_standardize = feature_standardize,
    feature_center = standardized$center,
    feature_scale = standardized$scale,
    zero_scale_features = standardized$zero_scale_features,
    inv_cor_scaled = cor_data$inv_cor * cor_scale,
    pixel_key = paste(df_base[[x_col]], df_base[[y_col]], sep = "_"),
    x_col = x_col,
    y_col = y_col,
    parameters = list(
      normalize_method = normalize_method,
      feature_standardize = feature_standardize,
      neighbor_radius = neighbor_radius,
      cor_cores = cor_cores,
      cor_scale = cor_scale,
      clamp_negative_cor = clamp_negative_cor,
      weighted_neighbors = weighted_neighbors
    )
  )
}

msiclust_labels_from_membership <- function(cluster,
                                            membership,
                                            min_membership = 0.5,
                                            no_cluster = "No_cluster") {
  max_membership <- .msiclust_row_maxs(membership)
  labels <- ifelse(max_membership > min_membership, as.character(cluster), no_cluster)

  list(
    labels = labels,
    max_membership = max_membership
  )
}

run_msiclust <- function(prep,
                         nclust,
                         iter_max = 100L,
                         min_membership = 0.5,
                         no_cluster = "No_cluster",
                         seed = NULL) {
  .msiclust_require("vsclust", "to run MSIClust")

  if (!is.list(prep) || is.null(prep$X) || is.null(prep$inv_cor_scaled)) {
    stop("prep must be an object returned by prepare_msiclust_input().", call. = FALSE)
  }

  nclust <- as.integer(nclust)
  if (length(nclust) != 1L || is.na(nclust) || nclust < 2L) {
    stop("nclust must be a single integer >= 2.", call. = FALSE)
  }
  if (nrow(prep$X) < nclust) {
    stop("Cannot cluster fewer pixels than nclust.", call. = FALSE)
  }

  if (!is.null(seed)) {
    set.seed(seed)
  }

  start_time <- Sys.time()
  fuzz <- vsclust::determine_fuzz(
    dims = dim(prep$X),
    NClust = nclust,
    Sds = prep$inv_cor_scaled
  )
  alg <- vsclust::vsclust_algorithm(
    prep$X,
    centers = nclust,
    iterMax = iter_max,
    m = fuzz$m
  )

  membership <- alg$membership
  label_info <- msiclust_labels_from_membership(
    cluster = alg$cluster,
    membership = membership,
    min_membership = min_membership,
    no_cluster = no_cluster
  )

  list(
    method = "MSIClust",
    labels = label_info$labels,
    centers = alg$centers,
    center_features = prep$signal_cols,
    raw_cluster = alg$cluster,
    membership = membership,
    max_membership = label_info$max_membership,
    model = list(
      fuzzifier = fuzz$m,
      cluster = alg$cluster,
      membership = membership
    ),
    parameters = c(prep$parameters, list(
      nclust = nclust,
      iter_max = iter_max,
      min_membership = min_membership,
      no_cluster = no_cluster,
      seed = seed
    )),
    runtime_sec = as.numeric(difftime(Sys.time(), start_time, units = "secs"))
  )
}

append_msiclust_labels <- function(prep,
                                   msiclust_res,
                                   cluster_col = "MSIClust_cluster") {
  if (length(msiclust_res$labels) != nrow(prep$df_base)) {
    stop("Number of labels does not match number of prepared pixels.", call. = FALSE)
  }

  df <- prep$df_base
  df[[cluster_col]] <- msiclust_res$labels
  df
}

run_msiclust_workflow <- function(msi_df,
                                  nclust,
                                  signal_cols = msiclust_signal_columns(msi_df),
                                  x_col = "x",
                                  y_col = "y",
                                  normalize_method = c("tic", "median", "rms"),
                                  feature_standardize = c("none", "sd", "zscore"),
                                  neighbor_radius = 1L,
                                  cor_cores = 1L,
                                  cor_scale = 25,
                                  clamp_negative_cor = FALSE,
                                  weighted_neighbors = TRUE,
                                  iter_max = 100L,
                                  min_membership = 0.5,
                                  no_cluster = "No_cluster",
                                  seed = NULL,
                                  cluster_col = "MSIClust_cluster") {
  prep <- prepare_msiclust_input(
    msi_df = msi_df,
    signal_cols = signal_cols,
    x_col = x_col,
    y_col = y_col,
    normalize_method = normalize_method,
    feature_standardize = feature_standardize,
    neighbor_radius = neighbor_radius,
    cor_cores = cor_cores,
    cor_scale = cor_scale,
    clamp_negative_cor = clamp_negative_cor,
    weighted_neighbors = weighted_neighbors
  )

  msiclust_res <- run_msiclust(
    prep = prep,
    nclust = nclust,
    iter_max = iter_max,
    min_membership = min_membership,
    no_cluster = no_cluster,
    seed = seed
  )

  list(
    clustered_df = append_msiclust_labels(prep, msiclust_res, cluster_col = cluster_col),
    msiclust = msiclust_res,
    prep = prep,
    cluster_col = cluster_col
  )
}

msiclust_cluster_palette <- function(clusters,
                                     no_cluster = "No_cluster",
                                     no_cluster_color = "#D9D9D9",
                                     palette = "Set3") {
  clusters <- as.character(clusters)
  valid_clusters <- sort(setdiff(unique(clusters), no_cluster))
  n_valid <- length(valid_clusters)

  if (n_valid == 0L) {
    return(setNames(no_cluster_color, no_cluster))
  }

  cols <- NULL
  if (requireNamespace("RColorBrewer", quietly = TRUE)) {
    palette_names <- c(
      palette,
      "Dark2", "Set1", "Paired", "Accent", "Set2", "Spectral", "RdYlBu"
    )
    palette_names <- unique(palette_names[palette_names %in% rownames(RColorBrewer::brewer.pal.info)])
    cols <- unlist(lapply(palette_names, function(pal) {
      info <- RColorBrewer::brewer.pal.info[pal, ]
      RColorBrewer::brewer.pal(info$maxcolors, pal)
    }), use.names = FALSE)
    cols <- unique(toupper(cols))
  }

  if (is.null(cols) || length(cols) < n_valid) {
    cols <- grDevices::hcl(
      h = seq(15, 375, length.out = n_valid + 1L)[seq_len(n_valid)],
      c = 80,
      l = 55
    )
  }

  cols <- cols[seq_len(n_valid)]
  names(cols) <- valid_clusters
  c(setNames(no_cluster_color, no_cluster), cols)
}

plot_msiclust_map <- function(result_or_df,
                              cluster_col = "MSIClust_cluster",
                              x_col = "x",
                              y_col = "y",
                              no_cluster = "No_cluster",
                              no_cluster_color = "#D9D9D9",
                              title = "MSIClust",
                              palette = "Set3") {
  .msiclust_require("ggplot2", "to plot MSIClust maps")

  df <- if (is.list(result_or_df) && "clustered_df" %in% names(result_or_df)) {
    result_or_df$clustered_df
  } else {
    result_or_df
  }

  required_cols <- c(x_col, y_col, cluster_col)
  missing_cols <- setdiff(required_cols, names(df))
  if (length(missing_cols) > 0L) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "), call. = FALSE)
  }

  df_plot <- data.frame(
    x = df[[x_col]],
    y = df[[y_col]],
    cluster = as.character(df[[cluster_col]]),
    stringsAsFactors = FALSE
  )

  present_clusters <- unique(df_plot$cluster)
  present_clusters <- c(
    if (no_cluster %in% present_clusters) no_cluster,
    sort(setdiff(present_clusters, no_cluster))
  )
  df_plot$cluster <- factor(df_plot$cluster, levels = present_clusters)

  all_colors <- msiclust_cluster_palette(
    clusters = present_clusters,
    no_cluster = no_cluster,
    no_cluster_color = no_cluster_color,
    palette = palette
  )
  all_colors <- all_colors[present_clusters]

  ggplot2::ggplot(df_plot, ggplot2::aes(x = .data$x, y = .data$y, fill = .data$cluster)) +
    ggplot2::geom_raster() +
    ggplot2::scale_fill_manual(values = all_colors, drop = FALSE) +
    ggplot2::coord_equal() +
    ggplot2::scale_y_reverse() +
    ggplot2::labs(title = title, x = x_col, y = y_col, fill = "Cluster") +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.line = ggplot2::element_line(color = "black", linewidth = 0.4),
      axis.ticks = ggplot2::element_line(color = "black"),
      panel.grid = ggplot2::element_blank(),
      legend.position = "bottom"
    )
}

# Compatibility aliases matching names used in older local scripts.
msi_signal_columns <- msiclust_signal_columns
msi_feature_mz_values <- msiclust_feature_mz_values
make_unique_signal_columns <- msiclust_make_unique_signal_columns
normalize_pixels <- msiclust_normalize_pixels
standardize_feature_columns <- msiclust_standardize_features
compute_neighbor_cor <- msiclust_compute_neighbor_cor
prepare_msi_clustering_input <- prepare_msiclust_input
run_msiclust_method <- run_msiclust
plot_cluster_map <- plot_msiclust_map

