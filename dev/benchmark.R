# =============================================================================
# BMIQ: timing, agreement, and geometry (production only)
#
# Run from the package root:
#   source("dev/benchmark.R")
#
# Requires: bench, ggplot2, tidyr, dplyr (and the package itself via load_all).
# Entry point: BMIQcalibration()  (nL = 3 and nL = 2)
#
# Production path: adaptive blc clipping + weighted-density intersection
# thresholds + continuous H (nL=3); truncated U/M maps joined at gold cut
# (nL=2).
#
# Plots include gold standard (target), sample before, and sample after
# (nL=3 / nL=2 / nL=3 H-off).
# =============================================================================

stopifnot(requireNamespace("bench", quietly = TRUE))
stopifnot(requireNamespace("ggplot2", quietly = TRUE))
stopifnot(requireNamespace("tidyr", quietly = TRUE))
stopifnot(requireNamespace("dplyr", quietly = TRUE))

library(bench)
library(ggplot2)
library(tidyr)
library(dplyr)

devtools::load_all(".", quiet = TRUE)

# -----------------------------------------------------------------------------
# Data: one sample aligned to Horvath gold standard
# -----------------------------------------------------------------------------
data("GPL21145_sample", envir = environment())
data("horvath_goldstandard", envir = environment())

gold <- as.numeric(horvath_goldstandard$goldstandard2)
if (length(gold) != ncol(GPL21145_sample)) {
  idx <- match(colnames(GPL21145_sample), horvath_goldstandard$Name)
  gold <- as.numeric(horvath_goldstandard$goldstandard2[idx])
}

sample_i <- 1L
dat1 <- GPL21145_sample[sample_i, , drop = FALSE]
ok <- is.finite(dat1[1L, ]) &
  is.finite(gold) &
  dat1[1L, ] >= 0 &
  dat1[1L, ] <= 1 &
  gold >= 0 &
  gold <= 1
dat1 <- dat1[, ok, drop = FALSE]
gold <- gold[ok]
storage.mode(dat1) <- "double"

sample_id <- rownames(dat1)
if (is.null(sample_id) || !nzchar(sample_id[[1L]])) {
  sample_id <- paste0("sample_", sample_i)
}
cpg_ids <- colnames(dat1)
if (is.null(cpg_ids)) {
  cpg_ids <- paste0("cg", seq_len(ncol(dat1)))
}

before_v <- as.numeric(dat1[1L, ])
gold_v <- as.numeric(gold)
n_zero <- sum(before_v == 0)
n_one <- sum(before_v == 1)

message(
  "Benchmark data: 1 sample x ",
  ncol(dat1),
  " probes (sample = ",
  sample_id,
  "; exact 0s = ",
  n_zero,
  ", exact 1s = ",
  n_one,
  ")"
)

run_bmiq <- function(nfit, nL = 3L, doH = NULL, debug = FALSE) {
  BMIQcalibration(
    datM = dat1,
    goldstandard.beta = gold_v,
    nL = as.integer(nL),
    doH = doH,
    nfit = as.integer(nfit),
    niter = 5L,
    tol = 0.001,
    verbose = FALSE,
    debug = isTRUE(debug)
  )
}

# =============================================================================
# 1) Timing: nL=3 vs nL=2  (nfit = 10000)
# =============================================================================
message("\n=== bench::mark (nfit = 10000, check = FALSE) ===")

bm <- bench::mark(
  nL3 = run_bmiq(10000L, nL = 3L)$calibrated,
  nL2 = run_bmiq(10000L, nL = 2L)$calibrated,
  check = FALSE,
  iterations = 3,
  memory = TRUE
)

print(bm)
print(summary(bm))

# =============================================================================
# 2) Agreement: gold / before / nL3 / nL2 / nL3_noH  (nfit = 20000)
# =============================================================================
message("\n=== agreement run (nfit = 20000, 1 sample) ===")

new3_res <- run_bmiq(20000L, nL = 3L, debug = TRUE)
new2_res <- run_bmiq(20000L, nL = 2L, debug = TRUE)
# Isolate U compression from H stitching (nL = 3, H off).
new3_noH_res <- run_bmiq(20000L, nL = 3L, doH = FALSE, debug = TRUE)

new_nL3_v <- as.numeric(new3_res$calibrated[1L, ])
new_nL2_v <- as.numeric(new2_res$calibrated[1L, ])
new_nL3_noH_v <- as.numeric(new3_noH_res$calibrated[1L, ])

stopifnot(
  length(before_v) == length(new_nL3_v),
  length(before_v) == length(new_nL2_v),
  length(before_v) == length(new_nL3_noH_v),
  length(before_v) == length(gold_v),
  length(before_v) == length(cpg_ids)
)

make_long_block <- function(value, type) {
  data.frame(
    id = sample_id,
    cpg = cpg_ids,
    value = value,
    type = type,
    stringsAsFactors = FALSE
  )
}

long <- rbind(
  make_long_block(gold_v, "gold"),
  make_long_block(before_v, "before"),
  make_long_block(new_nL3_v, "nL3"),
  make_long_block(new_nL2_v, "nL2"),
  make_long_block(new_nL3_noH_v, "nL3_noH")
)

wide <- long |>
  tidyr::pivot_wider(names_from = type, values_from = value) |>
  dplyr::mutate(
    abs_diff_nL3_nL2 = abs(nL3 - nL2),
    abs_diff_before_gold = abs(before - gold),
    abs_diff_nL3_gold = abs(nL3 - gold),
    abs_diff_nL2_gold = abs(nL2 - gold),
    abs_diff_nL3_noH_gold = abs(nL3_noH - gold),
    delta_nL3 = nL3 - before,
    delta_nL2 = nL2 - before
  )

message("abs(nL3 - nL2) summary:")
print(summary(wide$abs_diff_nL3_nL2))
message("abs(before - gold) summary:")
print(summary(wide$abs_diff_before_gold))
message("abs(nL3 - gold) summary:")
print(summary(wide$abs_diff_nL3_gold))
message("abs(nL2 - gold) summary:")
print(summary(wide$abs_diff_nL2_gold))
message("abs(nL3_noH - gold) summary:")
print(summary(wide$abs_diff_nL3_noH_gold))

summ_diff <- function(x, label) {
  message(
    label,
    ": max = ",
    signif(max(x, na.rm = TRUE), 8),
    "; mean = ",
    signif(mean(x, na.rm = TRUE), 8),
    "; median = ",
    signif(median(x, na.rm = TRUE), 8)
  )
}
summ_diff(wide$abs_diff_nL3_nL2, "max/mean/median |nL3 - nL2|")
summ_diff(wide$abs_diff_before_gold, "max/mean/median |before - gold|")
summ_diff(wide$abs_diff_nL3_gold, "max/mean/median |nL3 - gold|")
summ_diff(wide$abs_diff_nL2_gold, "max/mean/median |nL2 - gold|")
summ_diff(wide$abs_diff_nL3_noH_gold, "max/mean/median |nL3_noH - gold|")

message(
  "cor(*, gold): before=",
  signif(cor(before_v, gold_v), 6),
  "; nL3=",
  signif(cor(new_nL3_v, gold_v), 6),
  "; nL2=",
  signif(cor(new_nL2_v, gold_v), 6),
  "; nL3_noH=",
  signif(cor(new_nL3_noH_v, gold_v), 6)
)
message(
  "cor(before, nL3)=",
  signif(cor(before_v, new_nL3_v), 6),
  "; cor(before, nL2)=",
  signif(cor(before_v, new_nL2_v), 6),
  "; cor(nL3, nL2)=",
  signif(cor(new_nL3_v, new_nL2_v), 6)
)

# =============================================================================
# 3) Geometry checks (component means, thresholds, H anchors, boundary jumps)
# =============================================================================
message("\n=== geometry checks (debug) ===")

geom_report <- function(res, label, raw = before_v, cal = NULL) {
  if (is.null(cal)) {
    cal <- as.numeric(res$calibrated[1L, ])
  }
  d <- res$diagnostics$samples[[1L]]
  mu <- as.numeric(d$component_means)
  thr <- as.numeric(d$thresholds)
  nL <- length(mu)

  cat("\n-- ", label, " --\n", sep = "")
  cat("  component means: ", paste(signif(mu, 5), collapse = ", "), "\n", sep = "")
  cat("  thresholds:      ", paste(signif(thr, 5), collapse = ", "), "\n", sep = "")
  cat(
    "  a shapes:        ",
    paste(signif(as.numeric(d$component_a), 5), collapse = ", "),
    "\n",
    sep = ""
  )
  cat(
    "  b shapes:        ",
    paste(signif(as.numeric(d$component_b), 5), collapse = ", "),
    "\n",
    sep = ""
  )
  cat(
    "  eta:             ",
    paste(signif(as.numeric(d$eta), 5), collapse = ", "),
    "\n",
    sep = ""
  )

  # Thresholds should be weighted-density crossings (recompute from fit).
  a <- as.numeric(d$component_a)
  b <- as.numeric(d$component_b)
  eta <- as.numeric(d$eta)
  expected_thr <- tryCatch(
    thresholdsFromDensityCrossings(
      a = a,
      b = b,
      eta = eta,
      component.means = mu,
      context = "benchmark-recompute"
    ),
    error = function(e) e
  )
  if (inherits(expected_thr, "error")) {
    cat(
      "  density-crossing recompute failed: ",
      conditionMessage(expected_thr),
      "\n",
      sep = ""
    )
  } else {
    thr_ok <- isTRUE(all.equal(thr, expected_thr, tolerance = 1e-8))
    cat("  thresholds == density crossings? ", thr_ok, "\n", sep = "")
    mean_mid <- (mu[-nL] + mu[-1L]) / 2
    cat(
      "  |thr - mean midpoint|: ",
      paste(signif(abs(thr - mean_mid), 4), collapse = ", "),
      "\n",
      sep = ""
    )
  }

  if (nL >= 3L) {
    inside <- mu[2L] > thr[1L] && mu[2L] < thr[2L]
    ordered <- mu[1L] < thr[1L] &&
      thr[1L] < mu[2L] &&
      mu[2L] < thr[2L] &&
      thr[2L] < mu[3L]
    cat(
      "  mu_H inside H band? ",
      inside,
      "  [not required for density MAP]\n",
      sep = ""
    )
    cat("  full U < t1 < H < t2 < M geometry? ", ordered, "\n", sep = "")
  } else {
    cat(
      "  mu_U < t < mu_M? ",
      mu[1L] < thr[1L] && thr[1L] < mu[2L],
      "\n",
      sep = ""
    )
    if (!is.null(d$nl2_gold_threshold)) {
      cat(
        "  nL=2 gold cut t_g: ",
        signif(d$nl2_gold_threshold, 5),
        "\n",
        sep = ""
      )
    }
  }

  # Class labels from thresholds (lower-class equality).
  cls <- rep.int(1L, length(raw))
  for (b in seq_along(thr)) {
    cls[raw > thr[b]] <- b + 1L
  }
  cat(
    "  class counts:     ",
    paste(tabulate(cls, nL), collapse = " / "),
    "\n",
    sep = ""
  )

  # Boundary jumps on observed probes: g(upper class min) - g(lower class max).
  for (b in seq_along(thr)) {
    lower <- which(cls == b)
    upper <- which(cls == b + 1L)
    if (length(lower) && length(upper)) {
      jump <- min(cal[upper]) - max(cal[lower])
      cat(
        "  jump class ",
        b,
        "|",
        b + 1L,
        " (min upper - max lower): ",
        signif(jump, 5),
        "\n",
        sep = ""
      )
    }
  }

  # Monotonicity / displacement of the observed map.
  o <- order(raw, cal)
  mono_viol <- sum(diff(cal[o]) < -1e-12)
  cat("  rank inversions (order by raw then cal): ", mono_viol, "\n", sep = "")
  cat(
    "  max |cal - raw|: ",
    signif(max(abs(cal - raw)), 5),
    "; mean |cal - raw|: ",
    signif(mean(abs(cal - raw)), 5),
    "\n",
    sep = ""
  )
  cat(
    "  max |cal - gold|: ",
    signif(max(abs(cal - gold_v)), 5),
    "; mean |cal - gold|: ",
    signif(mean(abs(cal - gold_v)), 5),
    "\n",
    sep = ""
  )

  h_applied <- res$h.applied[1L]
  cat("  h.applied: ", h_applied, "\n", sep = "")
  if (isTRUE(h_applied)) {
    cat(
      "  H input range:    ",
      paste(signif(d$H_input_range, 5), collapse = " -> "),
      "\n",
      sep = ""
    )
    cat(
      "  H output anchors: ",
      paste(signif(d$H_output_anchors, 5), collapse = " -> "),
      "\n",
      sep = ""
    )
    cat("  H scale (hf):     ", signif(d$H_scale, 5), "\n", sep = "")
    # Continuous stitch: left H anchor should equal max calibrated U.
    u_idx <- which(cls == 1L)
    if (length(u_idx)) {
      cat(
        "  nminH - max(cal U): ",
        signif(d$H_output_anchors[1L] - max(cal[u_idx]), 5),
        "\n",
        sep = ""
      )
    }
  }

  invisible(list(
    means = mu,
    thresholds = thr,
    class = cls,
    h.applied = h_applied,
    diagnostics = d
  ))
}

geom_nL3 <- geom_report(new3_res, "nL=3 (H on)")
geom_nL3_noH <- geom_report(new3_noH_res, "nL=3 (H off)")
geom_nL2 <- geom_report(new2_res, "nL=2 (truncated maps)")

# nL=2: cut gap should be ~0 under truncated maps (continuous join at t_g).
if (length(geom_nL2$thresholds) == 1L) {
  t2 <- geom_nL2$thresholds[1L]
  left <- which(before_v <= t2)
  right <- which(before_v > t2)
  if (length(left) && length(right)) {
    gap <- min(new_nL2_v[right]) - max(new_nL2_v[left])
    message(
      "nL=2 cut gap (min M - max U): ",
      signif(gap, 5),
      if (gap < -1e-8) {
        "  [DOWNWARD — rank inversion risk]"
      } else if (abs(gap) < 1e-4) {
        "  [~continuous]"
      } else {
        "  [non-negative]"
      }
    )
  }
}

# Gold mixture geometry (from nL=3 run diagnostics if present).
if (!is.null(new3_res$diagnostics$gold)) {
  g <- new3_res$diagnostics$gold
  cat("\n-- gold-standard mixture (from nL=3 fit) --\n")
  cat(
    "  component means: ",
    paste(signif(as.numeric(g$component_means), 5), collapse = ", "),
    "\n",
    sep = ""
  )
  cat(
    "  thresholds:      ",
    paste(signif(as.numeric(g$thresholds), 5), collapse = ", "),
    "\n",
    sep = ""
  )
  cat(
    "  a shapes:        ",
    paste(signif(as.numeric(g$component_a), 5), collapse = ", "),
    "\n",
    sep = ""
  )
  cat(
    "  b shapes:        ",
    paste(signif(as.numeric(g$component_b), 5), collapse = ", "),
    "\n",
    sep = ""
  )
  cat(
    "  eta:             ",
    paste(signif(as.numeric(g$eta), 5), collapse = ", "),
    "\n",
    sep = ""
  )
}

# =============================================================================
# 4) Dense-grid mapping function (observed nearest-neighbor step map)
# =============================================================================
# BMIQ is probe-set based; approximate the map as the isotonic nearest-neighbour
# of observed (raw -> cal) pairs so we can plot g(x) on a dense beta grid.
message("\n=== dense-grid map approximation ===")

approx_map <- function(raw, cal, grid) {
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
    approx_map(before_v, new_nL3_v, grid),
    approx_map(before_v, new_nL2_v, grid),
    approx_map(before_v, new_nL3_noH_v, grid)
  ),
  method = factor(
    rep(c("nL3", "nL2", "nL3_noH"), each = length(grid)),
    levels = c("nL3", "nL2", "nL3_noH")
  )
)

# Finite-difference slopes near class boundaries for nL=3.
if (length(geom_nL3$thresholds) >= 1L) {
  thr3 <- geom_nL3$thresholds
  g3 <- approx_map(before_v, new_nL3_v, grid)
  for (t in thr3) {
    i <- which.min(abs(grid - t))
    i_lo <- max(1L, i - 3L)
    i_hi <- min(length(grid), i + 3L)
    slope_left <- (g3[i] - g3[i_lo]) / (grid[i] - grid[i_lo])
    slope_right <- (g3[i_hi] - g3[i]) / (grid[i_hi] - grid[i])
    message(
      "nL3 slope near t=",
      signif(t, 4),
      ": left=",
      signif(slope_left, 4),
      ", right=",
      signif(slope_right, 4),
      ", ratio=",
      signif(slope_right / slope_left, 4)
    )
  }
}

# nL=2 slopes near the single cut.
if (length(geom_nL2$thresholds) == 1L) {
  t2 <- geom_nL2$thresholds[1L]
  g2 <- approx_map(before_v, new_nL2_v, grid)
  i <- which.min(abs(grid - t2))
  i_lo <- max(1L, i - 3L)
  i_hi <- min(length(grid), i + 3L)
  slope_left <- (g2[i] - g2[i_lo]) / (grid[i] - grid[i_lo])
  slope_right <- (g2[i_hi] - g2[i]) / (grid[i_hi] - grid[i])
  message(
    "nL2 slope near t=",
    signif(t2, 4),
    ": left=",
    signif(slope_left, 4),
    ", right=",
    signif(slope_right, 4),
    ", ratio=",
    signif(slope_right / slope_left, 4)
  )
}

# Scatter helper
scatter_theme <- theme_bw(base_size = 12) +
  theme(aspect.ratio = 1)

scatter_pair <- function(data, x, y, title, subtitle = NULL) {
  ggplot(data, aes(x = .data[[x]], y = .data[[y]])) +
    geom_point(alpha = 0.15, size = 0.4) +
    geom_abline(slope = 1, intercept = 0, colour = "red", linetype = 2) +
    coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
    labs(title = title, x = x, y = y, subtitle = subtitle) +
    scatter_theme
}

# Gold target scatters (before / after calibration)
p_before_gold <- scatter_pair(
  wide,
  "gold",
  "before",
  "Gold vs before (sample)",
  paste0("mean |diff| = ", signif(mean(wide$abs_diff_before_gold), 4))
)
p_nL3_gold <- scatter_pair(
  wide,
  "gold",
  "nL3",
  "Gold vs nL=3 (after)",
  paste0("mean |diff| = ", signif(mean(wide$abs_diff_nL3_gold), 4))
)
p_nL2_gold <- scatter_pair(
  wide,
  "gold",
  "nL2",
  "Gold vs nL=2 (after)",
  paste0("mean |diff| = ", signif(mean(wide$abs_diff_nL2_gold), 4))
)

# Sample path scatters
p_before_nL3 <- scatter_pair(wide, "before", "nL3", "Before vs nL=3")
p_before_nL2 <- scatter_pair(wide, "before", "nL2", "Before vs nL=2")
p_nL3_nL2 <- scatter_pair(
  wide,
  "nL3",
  "nL2",
  "nL=3 vs nL=2",
  paste0("mean |diff| = ", signif(mean(wide$abs_diff_nL3_nL2), 4))
)
p_before_nL3_noH <- scatter_pair(
  wide,
  "before",
  "nL3_noH",
  "Before vs nL=3 (H off)"
)

# Diff histograms — emphasize distance to gold
diff_long <- tidyr::pivot_longer(
  wide,
  cols = c(
    abs_diff_before_gold,
    abs_diff_nL3_gold,
    abs_diff_nL2_gold,
    abs_diff_nL3_noH_gold,
    abs_diff_nL3_nL2
  ),
  names_to = "pair",
  values_to = "abs_diff"
)
diff_long$pair <- factor(
  diff_long$pair,
  levels = c(
    "abs_diff_before_gold",
    "abs_diff_nL3_gold",
    "abs_diff_nL2_gold",
    "abs_diff_nL3_noH_gold",
    "abs_diff_nL3_nL2"
  ),
  labels = c(
    "|before - gold|",
    "|nL3 - gold|",
    "|nL2 - gold|",
    "|nL3_noH - gold|",
    "|nL3 - nL2|"
  )
)

p_abs_hist <- ggplot(diff_long, aes(x = abs_diff)) +
  geom_histogram(bins = 60, fill = "grey40", colour = NA) +
  scale_x_continuous(trans = "sqrt") +
  facet_wrap(~pair, scales = "free_y", ncol = 3) +
  labs(
    title = "Absolute differences (vs gold and method pairs)",
    x = "abs(diff) (sqrt scale)",
    y = "count"
  ) +
  theme_bw(base_size = 12)

# Overlaid densities: gold / before / nL3 / nL2
long_plot <- long[long$type %in% c("gold", "before", "nL3", "nL2"), ]
long_plot$type <- factor(
  long_plot$type,
  levels = c("gold", "before", "nL3", "nL2")
)
type_cols <- c(
  gold = "#E45756",
  before = "#4C78A8",
  nL3 = "#54A24B",
  nL2 = "#B279A2"
)

p_density <- ggplot(long_plot, aes(x = value, colour = type, fill = type)) +
  geom_density(alpha = 0.10, linewidth = 0.8) +
  scale_x_continuous(limits = c(0, 1), expand = c(0.01, 0)) +
  scale_colour_manual(values = type_cols) +
  scale_fill_manual(values = type_cols) +
  labs(
    title = "Beta density: gold / before / nL=3 / nL=2",
    x = "beta",
    y = "density",
    colour = NULL,
    fill = NULL
  ) +
  theme_bw(base_size = 12) +
  theme(legend.position = "top")

map_cols <- c(
  nL3 = "#54A24B",
  nL2 = "#B279A2",
  nL3_noH = "#72B7B2"
)

p_map <- ggplot(map_df, aes(x = beta, y = calibrated, colour = method)) +
  geom_abline(slope = 1, intercept = 0, colour = "grey70", linetype = 3) +
  geom_line(linewidth = 0.7, alpha = 0.9) +
  {
    if (length(geom_nL3$thresholds)) {
      geom_vline(
        xintercept = geom_nL3$thresholds,
        colour = "grey40",
        linetype = 2,
        linewidth = 0.4
      )
    } else {
      NULL
    }
  } +
  {
    if (length(geom_nL2$thresholds) == 1L) {
      geom_vline(
        xintercept = geom_nL2$thresholds,
        colour = "#B279A2",
        linetype = 3,
        linewidth = 0.4
      )
    } else {
      NULL
    }
  } +
  coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
  scale_colour_manual(values = map_cols) +
  labs(
    title = "Calibrated map g(beta) on dense grid (NN / linear approx)",
    subtitle = "Grey dashed: nL=3 cuts; purple dotted: nL=2 cut",
    x = "raw beta",
    y = "calibrated beta",
    colour = NULL
  ) +
  theme_bw(base_size = 12) +
  theme(legend.position = "top", aspect.ratio = 1)

if (requireNamespace("patchwork", quietly = TRUE)) {
  p1 <- (p_before_gold | p_nL3_gold | p_nL2_gold) /
    (p_before_nL3 | p_before_nL2 | p_nL3_nL2) /
    (p_density | p_abs_hist) /
    (p_map | p_before_nL3_noH)
  print(p1)
} else {
  print(p_before_gold)
  print(p_nL3_gold)
  print(p_nL2_gold)
  print(p_before_nL3)
  print(p_before_nL2)
  print(p_nL3_nL2)
  print(p_density)
  print(p_abs_hist)
  print(p_map)
  print(p_before_nL3_noH)
  p1 <- list(
    before_gold = p_before_gold,
    nL3_gold = p_nL3_gold,
    nL2_gold = p_nL2_gold,
    before_nL3 = p_before_nL3,
    before_nL2 = p_before_nL2,
    nL3_nL2 = p_nL3_nL2,
    density = p_density,
    abs_hist = p_abs_hist,
    map = p_map,
    before_nL3_noH = p_before_nL3_noH
  )
}

out_dir <- file.path("dev", "tmp")
if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

saveRDS(
  list(
    bench = bm,
    long = long,
    wide = wide,
    map = map_df,
    geometry = list(
      nL3 = geom_nL3,
      nL3_noH = geom_nL3_noH,
      nL2 = geom_nL2,
      gold = new3_res$diagnostics$gold
    ),
    sample_id = sample_id,
    n_zero = n_zero,
    n_one = n_one
  ),
  file = file.path(out_dir, "benchmark-results.rds")
)

ggsave(
  filename = file.path(out_dir, "benchmark-scatters.png"),
  plot = if (inherits(p1, "ggplot") || inherits(p1, "patchwork")) {
    p1
  } else {
    p_density
  },
  width = 14,
  height = 16,
  dpi = 120
)

ggsave(
  filename = file.path(out_dir, "benchmark-map.png"),
  plot = p_map,
  width = 7,
  height = 6,
  dpi = 120
)

message("\nSaved: ", file.path(out_dir, "benchmark-results.rds"))
message("Saved: ", file.path(out_dir, "benchmark-scatters.png"))
message("Saved: ", file.path(out_dir, "benchmark-map.png"))
message("Done.")
