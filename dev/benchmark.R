# BMIQ density + calibrated-map plots
#
# From package root:
#   source("dev/benchmark.R")
#
# Requires: ggplot2, plotly (interactive HTML with click-to-toggle legend)

stopifnot(requireNamespace("ggplot2", quietly = TRUE))
stopifnot(requireNamespace("plotly", quietly = TRUE))
library(ggplot2)
library(plotly)

devtools::load_all(".", quiet = TRUE)
source("dev/bmiq-calibration-legacy.R", local = TRUE)

# -----------------------------------------------------------------------------
# Data: one sample, gold aligned by probe name (not positional hope)
# -----------------------------------------------------------------------------
data("GPL21145_sample", envir = environment())
data("horvath_goldstandard", envir = environment())

sample_i <- 1L
stopifnot(sample_i >= 1L, sample_i <= nrow(GPL21145_sample))

probe_names <- colnames(GPL21145_sample)
gold_idx <- match(probe_names, as.character(horvath_goldstandard$Name))
if (anyNA(gold_idx)) {
  stop(
    sum(is.na(gold_idx)),
    " sample probes not found in horvath_goldstandard$Name.",
    call. = FALSE
  )
}
gold <- as.numeric(horvath_goldstandard$goldstandard2)[gold_idx]

dat1 <- GPL21145_sample[sample_i, , drop = FALSE]
# Logical mask is positional on the 1-row matrix / gold vector (same probe order).
ok <- is.finite(as.numeric(dat1[1L, ])) &
  is.finite(gold) &
  as.numeric(dat1[1L, ]) >= 0 &
  as.numeric(dat1[1L, ]) <= 1 &
  gold >= 0 &
  gold <= 1

dat1 <- dat1[, ok, drop = FALSE]
gold <- gold[ok]
storage.mode(dat1) <- "double"

sample_id <- rownames(dat1)
if (is.null(sample_id) || !nzchar(sample_id[[1L]])) {
  sample_id <- paste0("sample_", sample_i)
}

# 1-row matrix: as.numeric() is column-major, which is probe order for nrow == 1.
before_v <- as.numeric(dat1)
gold_v <- as.numeric(gold)
stopifnot(
  length(before_v) == ncol(dat1),
  length(gold_v) == ncol(dat1),
  isTRUE(all.equal(before_v, as.numeric(GPL21145_sample[sample_i, ok])))
)

message(
  "Benchmark data: sample_i=",
  sample_i,
  " (",
  sample_id,
  ") x ",
  ncol(dat1),
  " probes; cor(before, gold)=",
  signif(cor(before_v, gold_v), 4)
)

nfit <- 20000L
niter <- 5L

run_bmiq <- function(nL) {
  res <- bmiq_calibration(
    datM = dat1,
    goldstandard.beta = gold_v,
    nL = nL,
    nfit = nfit,
    niter = niter,
    tol = 0.001,
    verbose = FALSE,
    debug = TRUE
  )
  stopifnot(isTRUE(res$success[[1L]]))
  calibrated <- as.numeric(res$calibrated[1L, ])
  stopifnot(length(calibrated) == length(before_v), all(is.finite(calibrated)))
  list(
    calibrated = calibrated,
    thresholds = as.numeric(res$diagnostics$samples[[1L]]$thresholds)
  )
}

# -----------------------------------------------------------------------------
# Calibrate: production nL = 2 / nL = 3, plus legacy nL = 3
# -----------------------------------------------------------------------------
message("Running bmiq_calibration(nL = 2) ...")
nL2 <- run_bmiq(2L)

message("Running bmiq_calibration(nL = 3) ...")
nL3 <- run_bmiq(3L)

message("Running legacy BMIQcalibration.old(nL = 3) ...")
legacy_mat <- BMIQcalibration.old(
  datM = dat1,
  goldstandard.beta = gold_v,
  nL = 3,
  doH = TRUE,
  nfit = nfit,
  niter = niter,
  tol = 0.001
)
stopifnot(nrow(legacy_mat) == 1L, ncol(legacy_mat) == length(before_v))
legacy_v <- as.numeric(legacy_mat[1L, ])
stopifnot(length(legacy_v) == length(before_v), all(is.finite(legacy_v)))

# -----------------------------------------------------------------------------
# Density: raw gold (empirical KDE of gold betas) + sample before / calibrated
# Gold is NOT the fitted mixture density (no dbeta / eta reconstruction).
# -----------------------------------------------------------------------------
sample_series <- list(
  before = before_v,
  `nL=2` = nL2$calibrated,
  `nL=3` = nL3$calibrated,
  legacy = legacy_v
)
stopifnot(all(lengths(sample_series) == length(before_v)))

sample_dens_df <- data.frame(
  value = unlist(sample_series, use.names = FALSE),
  type = factor(
    rep(names(sample_series), times = lengths(sample_series)),
    levels = names(sample_series)
  )
)

# Empirical density of the raw gold-standard beta vector (same probes as sample).
gold_dens_df <- data.frame(
  value = gold_v,
  type = factor("gold (raw)", levels = "gold (raw)")
)

type_cols <- c(
  "gold (raw)" = "#E45756",
  "before" = "#4C78A8",
  "nL=2" = "#F58518",
  "nL=3" = "#54A24B",
  "legacy" = "#B279A2"
)

p_density <- ggplot() +
  geom_density(
    data = sample_dens_df,
    aes(x = value, colour = type, fill = type),
    alpha = 0.10,
    linewidth = 0.8
  ) +
  # Reference: kernel density of raw gold betas only (not EM mixture fit).
  geom_density(
    data = gold_dens_df,
    aes(x = value, colour = type, fill = type),
    alpha = 0.05,
    linewidth = 1.1,
    linetype = "dashed"
  ) +
  # Do not use scale_x limits here — they drop mass and warp the KDE near 0/1.
  coord_cartesian(xlim = c(0, 1)) +
  scale_colour_manual(values = type_cols, breaks = names(type_cols)) +
  scale_fill_manual(values = type_cols, breaks = names(type_cols)) +
  labs(
    title = paste0("Beta density (", sample_id, ", i=", sample_i, ")"),
    subtitle = "raw gold-standard density (empirical) vs sample before / calibrated",
    x = "beta",
    y = "density",
    colour = NULL,
    fill = NULL
  ) +
  theme_bw(base_size = 12) +
  theme(legend.position = "top")

# -----------------------------------------------------------------------------
# Calibrated map: raw beta -> calibrated (dense-grid nearest-neighbor approx)
# -----------------------------------------------------------------------------
approx_map <- function(raw, cal, grid) {
  stopifnot(length(raw) == length(cal))
  o <- order(raw)
  x <- raw[o]
  y <- cal[o]
  ux <- unique(x)
  uy <- vapply(ux, function(xx) mean(y[x == xx]), numeric(1))
  approx(ux, uy, xout = grid, rule = 2, ties = mean)$y
}

grid <- seq(0.001, 0.999, length.out = 500L)
map_df <- data.frame(
  beta = rep(grid, 3L),
  calibrated = c(
    approx_map(before_v, nL2$calibrated, grid),
    approx_map(before_v, nL3$calibrated, grid),
    approx_map(before_v, legacy_v, grid)
  ),
  method = factor(
    rep(c("nL=2", "nL=3", "legacy"), each = length(grid)),
    levels = c("nL=2", "nL=3", "legacy")
  )
)

map_cols <- c(
  "nL=2" = "#F58518",
  "nL=3" = "#54A24B",
  "legacy" = "#B279A2"
)

p_map <- ggplot(map_df, aes(x = beta, y = calibrated, colour = method)) +
  geom_abline(slope = 1, intercept = 0, colour = "grey70", linetype = 3) +
  geom_line(linewidth = 0.7, alpha = 0.9) +
  {
    if (length(nL3$thresholds)) {
      geom_vline(
        xintercept = nL3$thresholds,
        colour = "grey40",
        linetype = 2,
        linewidth = 0.4
      )
    } else {
      NULL
    }
  } +
  coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
  scale_colour_manual(values = map_cols) +
  labs(
    title = paste0(
      "Calibrated map g(beta) (",
      sample_id,
      ", i=",
      sample_i,
      ")"
    ),
    subtitle = "Grey dashed: production nL = 3 density cuts",
    x = "raw beta",
    y = "calibrated beta",
    colour = NULL
  ) +
  theme_bw(base_size = 12) +
  theme(legend.position = "top", aspect.ratio = 1)

print(p_density)
print(p_map)

# -----------------------------------------------------------------------------
# Plotly: same series, click legend entries to toggle traces
# -----------------------------------------------------------------------------
kde_xy <- function(x, n = 512L) {
  d <- stats::density(x, from = 0, to = 1, n = n)
  data.frame(x = d$x, y = d$y)
}

dens_series <- c(
  list(`gold (raw)` = gold_v),
  sample_series
)
dens_linetypes <- c(
  "gold (raw)" = "dash",
  "before" = "solid",
  "nL=2" = "solid",
  "nL=3" = "solid",
  "legacy" = "solid"
)

p_density_plotly <- plot_ly()
for (nm in names(dens_series)) {
  xy <- kde_xy(dens_series[[nm]])
  p_density_plotly <- add_trace(
    p_density_plotly,
    x = xy$x,
    y = xy$y,
    type = "scatter",
    mode = "lines",
    name = nm,
    line = list(
      color = unname(type_cols[[nm]]),
      width = if (nm == "gold (raw)") 2.5 else 1.8,
      dash = unname(dens_linetypes[[nm]])
    ),
    hovertemplate = paste0(
      nm,
      "<br>beta: %{x:.3f}<br>density: %{y:.3f}<extra></extra>"
    )
  )
}
p_density_plotly <- layout(
  p_density_plotly,
  title = list(
    text = paste0(
      "Beta density (",
      sample_id,
      ", i=",
      sample_i,
      ")<br>",
      "<sup>raw gold (empirical) vs sample before / calibrated — ",
      "click legend to toggle</sup>"
    )
  ),
  xaxis = list(title = "beta", range = c(0, 1)),
  yaxis = list(title = "density"),
  legend = list(orientation = "h", y = 1.08),
  hovermode = "x unified"
)

p_map_plotly <- plot_ly()
p_map_plotly <- add_trace(
  p_map_plotly,
  x = c(0, 1),
  y = c(0, 1),
  type = "scatter",
  mode = "lines",
  name = "y = x",
  line = list(color = "grey70", dash = "dot", width = 1),
  hoverinfo = "skip",
  showlegend = TRUE
)
for (nm in levels(map_df$method)) {
  sub <- map_df[map_df$method == nm, , drop = FALSE]
  p_map_plotly <- add_trace(
    p_map_plotly,
    x = sub$beta,
    y = sub$calibrated,
    type = "scatter",
    mode = "lines",
    name = nm,
    line = list(color = unname(map_cols[[nm]]), width = 1.8),
    hovertemplate = paste0(
      nm,
      "<br>raw: %{x:.3f}<br>calibrated: %{y:.3f}<extra></extra>"
    )
  )
}
if (length(nL3$thresholds)) {
  for (i in seq_along(nL3$thresholds)) {
    th <- nL3$thresholds[[i]]
    p_map_plotly <- add_trace(
      p_map_plotly,
      x = c(th, th),
      y = c(0, 1),
      type = "scatter",
      mode = "lines",
      name = if (i == 1L) "nL=3 cuts" else paste0("nL=3 cut ", i),
      line = list(color = "grey40", dash = "dash", width = 1),
      hovertemplate = paste0("threshold: ", signif(th, 4), "<extra></extra>"),
      showlegend = i == 1L,
      legendgroup = "nL3_cuts"
    )
  }
}
p_map_plotly <- layout(
  p_map_plotly,
  title = list(
    text = paste0(
      "Calibrated map g(beta) (",
      sample_id,
      ", i=",
      sample_i,
      ")<br>",
      "<sup>click legend to toggle</sup>"
    )
  ),
  xaxis = list(title = "raw beta", range = c(0, 1), scaleanchor = "y"),
  yaxis = list(title = "calibrated beta", range = c(0, 1)),
  legend = list(orientation = "h", y = 1.08),
  hovermode = "closest"
)

print(p_density_plotly)
print(p_map_plotly)

# -----------------------------------------------------------------------------
# Save
# -----------------------------------------------------------------------------
out_dir <- file.path("dev", "tmp")
if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

res_nl2 <- bmiq_calibration(
  datM = dat1,
  goldstandard.beta = gold_v,
  nL = 2,
  nfit = nfit,
  niter = niter,
  tol = 0.001,
  verbose = FALSE,
  debug = TRUE
)

res_nl3 <- bmiq_calibration(
  datM = dat1,
  goldstandard.beta = gold_v,
  nL = 3,
  nfit = nfit,
  niter = niter,
  tol = 0.001,
  verbose = FALSE,
  debug = TRUE
)

# Drop huge index vectors so the dump stays readable.
strip_random_indices <- function(diagnostics) {
  if (!is.null(diagnostics$gold)) {
    diagnostics$gold$random_indices <- NULL
  }
  if (
    length(diagnostics$samples) >= 1L &&
      !is.null(diagnostics$samples[[1L]])
  ) {
    diagnostics$samples[[1L]]$random_indices <- NULL
  }
  diagnostics
}

diag_nl2 <- strip_random_indices(res_nl2$diagnostics)
diag_nl3 <- strip_random_indices(res_nl3$diagnostics)

sink_diagnostics <- function(diagnostics, path) {
  sink(path)
  on.exit(sink(), add = TRUE)
  print(diagnostics)
}

sink_diagnostics(diag_nl2, file.path(out_dir, "benchmark-diagnostics-nl2.txt"))
sink_diagnostics(diag_nl3, file.path(out_dir, "benchmark-diagnostics-nl3.txt"))
message("Wrote: ", file.path(out_dir, "benchmark-diagnostics-nl2.txt"))
message("Wrote: ", file.path(out_dir, "benchmark-diagnostics-nl3.txt"))

ggsave(
  filename = file.path(out_dir, "benchmark-density.png"),
  plot = p_density,
  width = 8,
  height = 5,
  dpi = 120
)
ggsave(
  filename = file.path(out_dir, "benchmark-map.png"),
  plot = p_map,
  width = 7,
  height = 6,
  dpi = 120
)

htmlwidgets::saveWidget(
  p_density_plotly,
  file = normalizePath(
    file.path(out_dir, "benchmark-density.html"),
    mustWork = FALSE
  ),
  selfcontained = TRUE,
  title = paste0("Beta density — ", sample_id)
)
htmlwidgets::saveWidget(
  p_map_plotly,
  file = normalizePath(
    file.path(out_dir, "benchmark-map.html"),
    mustWork = FALSE
  ),
  selfcontained = TRUE,
  title = paste0("Calibrated map — ", sample_id)
)

message("Saved: ", file.path(out_dir, "benchmark-density.png"))
message("Saved: ", file.path(out_dir, "benchmark-map.png"))
message("Saved: ", file.path(out_dir, "benchmark-density.html"))
message("Saved: ", file.path(out_dir, "benchmark-map.html"))
message("Done.")
