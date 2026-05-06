# -----------------------------------------------------------------------------
# 04_plot_brfss_2023_marginals.R
#
# Create the marginal-effect figure for the BRFSS 2023 real-data application.
#
# Run scripts/05_download_brfss_2023.R, scripts/06_prepare_brfss_2023.R,
# and scripts/07_fit_brfss_2023.R before running this script.
# -----------------------------------------------------------------------------

# Packages --------------------------------------------------------------------

if (!requireNamespace("bnpmccr", quietly = TRUE)) {
  stop(
    "The 'bnpmccr' package is not installed. ",
    "Install it before running this script, or use devtools::load_all('.') ",
    "from the package root during development."
  )
}

if (!requireNamespace("matrixStats", quietly = TRUE)) {
  stop("The 'matrixStats' package is required to run this script.")
}

# Retrieve package utilities. These are expected to be exported by bnpmccr, but
# the fallback also works during package development if they are internal.
get_bnpmccr_function <- function(name) {
  if (exists(name, where = asNamespace("bnpmccr"), inherits = FALSE)) {
    get(name, envir = asNamespace("bnpmccr"), inherits = FALSE)
  } else {
    stop("Could not find function '", name, "' in the bnpmccr namespace.")
  }
}

NaturalCubicBasis <- get_bnpmccr_function("NaturalCubicBasis")
select_knots <- get_bnpmccr_function("select_knots")

# Paths -----------------------------------------------------------------------

processed_rds_file <- file.path(
  "data", "processed", "brfss2023_preprocessed_for_bnpmccr.rds"
)

fit_rds_file <- file.path(
  "results", "brfss_2023", "bnpmccr_brfss2023_fit.rds"
)

settings_rds_file <- file.path(
  "results", "brfss_2023", "bnpmccr_brfss2023_fit_settings.rds"
)

figure_dir <- file.path("figures", "brfss_2023")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

figure_file <- file.path(
  figure_dir, "plot_marginals.pdf"
)

if (!file.exists(processed_rds_file)) {
  stop(
    "Preprocessed BRFSS file was not found: ", processed_rds_file, "\n",
    "Run scripts/06_prepare_brfss_2023.R first."
  )
}

if (!file.exists(fit_rds_file)) {
  stop(
    "BNPMCCR fit object was not found: ", fit_rds_file, "\n",
    "Run scripts/07_fit_brfss_2023.R first."
  )
}

# Load inputs -----------------------------------------------------------------

brfss_data <- readRDS(processed_rds_file)
fit <- readRDS(fit_rds_file)
fit_settings <- if (file.exists(settings_rds_file)) readRDS(settings_rds_file) else list()

Y <- brfss_data$Y
covT <- as.numeric(brfss_data$X)
predX <- as.matrix(brfss_data$X3)
age_raw <- brfss_data$age_raw

if (is.null(age_raw)) {
  stop(
    "The preprocessed RDS object must contain 'age_raw'. ",
    "Regenerate it using scripts/06_prepare_brfss_2023.R."
  )
}

stopifnot(
  nrow(Y) == length(covT),
  nrow(Y) == nrow(predX)
)

# Plot settings ---------------------------------------------------------------

responses <- c(
  "Diabetes", "HighBp", "HighChol", "Stroke",
  "HeartDisease", "Arthritis", "Asthma", "LogBMI"
)

resp_lbl_plotmath <- expression(
  paste("Diabetes", " (", y[i1], ")"),
  paste("HighBp", " (", y[i2], ")"),
  paste("HighChol", " (", y[i3], ")"),
  paste("Stroke", " (", y[i4], ")"),
  paste("HeartDisease", " (", y[i5], ")"),
  paste("Arthritis", " (", y[i6], ")"),
  paste("Asthma", " (", y[i7], ")"),
  paste("LogBMI", " (", y[i8], ")")
)

predictors <- c(
  "Intercept", "Gender", "Race", "Marriage", "Education",
  "NumChild", "Employment", "Income", "OwnHome", "Smoking",
  "HeavyDrink", "Insurance", "PhyRec", "Urban"
)

n_resp <- length(responses)
n_pred <- length(predictors)

stopifnot(
  ncol(Y) == n_resp,
  ncol(predX) + 1L == n_pred
)

num_basis <- if (!is.null(fit_settings$num_basis)) fit_settings$num_basis else 12L
n_grid <- 1000L
x_grid <- seq(min(covT), max(covT), length.out = n_grid)

# Posterior draws -------------------------------------------------------------

extract_beta_draws <- function(fit) {
  if (!is.null(fit$posterior_samples$beta)) {
    return(as.matrix(fit$posterior_samples$beta))
  }
  if (!is.null(fit$B)) {
    return(as.matrix(fit$B))
  }
  stop("Could not find posterior beta draws in 'fit$posterior_samples$beta' or 'fit$B'.")
}

extract_diabetes_cutpoint_trace <- function(fit) {
  if (!is.null(fit$cutpointsave)) {
    return(as.numeric(sapply(fit$cutpointsave, function(x) x[[1]][3])))
  }
  if (!is.null(fit$posterior_samples$cutpoint)) {
    cutpoint <- fit$posterior_samples$cutpoint
    if (is.matrix(cutpoint) || is.data.frame(cutpoint)) {
      cutpoint <- as.matrix(cutpoint)
      if (ncol(cutpoint) >= 3L) return(as.numeric(cutpoint[, 3L]))
    }
  }
  stop(
    "Could not find the diabetes cutpoint trace. ",
    "Expected 'fit$cutpointsave' or 'fit$posterior_samples$cutpoint'."
  )
}

beta_draws <- extract_beta_draws(fit)
cutpoint_trace <- extract_diabetes_cutpoint_trace(fit)

n_draws <- nrow(beta_draws)
posterior_draws <- 1001:3000
posterior_draws <- posterior_draws[posterior_draws <= n_draws]

if (length(posterior_draws) == 0L) {
  stop("No posterior draws are available for the requested range 1001:3000.")
}

cut_use <- cutpoint_trace[posterior_draws]
if (length(cut_use) != length(posterior_draws) || anyNA(cut_use)) {
  stop("The selected diabetes cutpoint trace does not match the selected posterior draws.")
}

# Basis matrix for marginal functions -----------------------------------------

knots <- select_knots(covT, num_basis - 2L)
basis_grid <- as.matrix(NaturalCubicBasis(x_grid, knots))
basis_grid <- sweep(basis_grid, 2L, colMeans(basis_grid), "-")

nb <- ncol(basis_grid)
P <- n_pred
J <- P * n_resp

# Wmat stores the centered spline basis repeated for each predictor.
Wmat <- do.call(cbind, replicate(P, basis_grid, simplify = FALSE))
block_width <- P + ncol(Wmat)

if (ncol(beta_draws) < n_resp * block_width) {
  stop(
    "The beta draw matrix has fewer columns than expected. ",
    "Check num_basis and the fit object structure."
  )
}

# Posterior summaries for marginal functions ---------------------------------

make_ylim <- function(mean_curve, lower_curve, upper_curve) {
  ymax <- max(upper_curve) + 0.5 * (max(upper_curve) - max(mean_curve))
  ymin <- min(lower_curve) - 0.5 * (-min(lower_curve) + min(mean_curve))
  c(ymin = ymin, ymax = ymax)
}

x_list <- vector("list", J)
mean_list <- vector("list", J)
lower_list <- vector("list", J)
upper_list <- vector("list", J)
ymin_vec <- numeric(J)
ymax_vec <- numeric(J)

# Extra curve for the first panel: beta_0(t) - cutpoint.
mean_list_cut <- vector("list", J)
lower_list_cut <- vector("list", J)
upper_list_cut <- vector("list", J)

for (panel_idx in seq_len(J)) {
  predictor_idx <- ((panel_idx - 1L) %% P) + 1L
  response_idx <- ((panel_idx - 1L) %/% P) + 1L

  block_cols <- ((response_idx - 1L) * block_width + 1L):(response_idx * block_width)
  samp <- beta_draws[posterior_draws, block_cols, drop = FALSE]

  intercept_draws <- samp[, predictor_idx]
  coef_cols <- P + seq_len(nb) + (predictor_idx - 1L) * nb
  basis_coef_draws <- samp[, coef_cols, drop = FALSE]

  basis_cols <- ((predictor_idx - 1L) * nb + 1L):(predictor_idx * nb)
  basis_segment <- Wmat[, basis_cols, drop = FALSE]

  beta_grid <- basis_segment %*% t(basis_coef_draws)
  beta_grid <- sweep(beta_grid, 2L, intercept_draws, "+")

  mean_curve <- matrixStats::rowMeans2(beta_grid)
  cred_curve <- matrixStats::rowQuantiles(beta_grid, probs = c(0.025, 0.975))
  lower_curve <- cred_curve[, 1L]
  upper_curve <- cred_curve[, 2L]

  x_list[[panel_idx]] <- x_grid
  mean_list[[panel_idx]] <- mean_curve
  lower_list[[panel_idx]] <- lower_curve
  upper_list[[panel_idx]] <- upper_curve

  ylim <- make_ylim(mean_curve, lower_curve, upper_curve)
  ymin_vec[panel_idx] <- ylim["ymin"]
  ymax_vec[panel_idx] <- ylim["ymax"]

  # In the (Diabetes, Intercept) panel, add the second cutpoint curve.
  if (response_idx == 1L && predictor_idx == 1L) {
    beta_grid_cut <- sweep(beta_grid, 2L, cut_use, "-")

    mean_curve_cut <- matrixStats::rowMeans2(beta_grid_cut)
    cred_curve_cut <- matrixStats::rowQuantiles(beta_grid_cut, probs = c(0.025, 0.975))

    mean_list_cut[[panel_idx]] <- mean_curve_cut
    lower_list_cut[[panel_idx]] <- cred_curve_cut[, 1L]
    upper_list_cut[[panel_idx]] <- cred_curve_cut[, 2L]
  }
}

# Axis helpers ----------------------------------------------------------------

age_ticks <- c(21, 27, 32, 37, 42, 47, 52, 57, 62, 67, 72, 77, 82)
age_ticks_scaled <- (age_ticks - mean(age_raw)) / stats::sd(age_raw)

make_y_ticks <- function(ymin, ymax,
                         pad_frac = 0.06, edge_frac = 0.10,
                         n_pretty = 2, min_ticks = 3,
                         dec_min = 2, dec_max = 6) {
  yr <- ymax - ymin
  if (!is.finite(yr) || yr <= 0) {
    yr <- ifelse(is.finite(ymin) && ymin != 0, abs(ymin), 1)
    ymin <- ymin - 0.5 * yr
    ymax <- ymax + 0.5 * yr
    yr <- ymax - ymin
  }

  ypad <- yr * pad_frac
  yt <- pretty(c(ymin - ypad, ymax + ypad), n = n_pretty)
  yt <- sort(unique(yt))

  inner_lo <- ymin + yr * edge_frac
  inner_hi <- ymax - yr * edge_frac
  if (!is.finite(inner_lo) || !is.finite(inner_hi) || inner_hi <= inner_lo) {
    inner_lo <- ymin
    inner_hi <- ymax
  }

  yt <- yt[yt > inner_lo & yt < inner_hi]
  if (length(yt) < min_ticks) {
    yt <- seq(inner_lo, inner_hi, length.out = min_ticks)
  }

  fix_neg0 <- function(s, dec) {
    pat <- paste0("^-0\\.", paste(rep("0", dec), collapse = ""))
    sub(pat, paste0("0.", paste(rep("0", dec), collapse = "")), s)
  }

  dec_used <- dec_min
  repeat {
    ylab <- sprintf(paste0("%.", dec_used, "f"), yt)
    ylab <- fix_neg0(ylab, dec_used)

    keep_unique <- !duplicated(ylab)
    yt2 <- yt[keep_unique]
    ylab2 <- ylab[keep_unique]

    if (length(yt2) >= min_ticks) {
      yt <- yt2
      ylab <- ylab2
      break
    }

    if (dec_used < dec_max) {
      dec_used <- dec_used + 1L
    } else {
      yt <- seq(inner_lo, inner_hi, length.out = min_ticks)
      ylab <- sprintf(paste0("%.", dec_used, "f"), yt)
      ylab <- fix_neg0(ylab, dec_used)
      keep_unique <- !duplicated(ylab)
      yt <- yt[keep_unique]
      ylab <- ylab[keep_unique]
      break
    }
  }

  list(at = yt, lab = ylab)
}

# Drawing helpers -------------------------------------------------------------

draw_band <- function(x, lower, upper, col) {
  polygon(
    x = c(x, rev(x)),
    y = c(upper, rev(lower)),
    col = col,
    border = NA
  )
}

draw_overlap_band <- function(x, lower, upper, idx, col) {
  rle_idx <- rle(idx)
  ends <- cumsum(rle_idx$lengths)
  starts <- c(1L, head(ends, -1L) + 1L)

  for (ii in seq_along(rle_idx$values)) {
    if (isTRUE(rle_idx$values[ii])) {
      jj <- starts[ii]:ends[ii]
      draw_band(x[jj], lower[jj], upper[jj], col = col)
    }
  }
}

# Plot ------------------------------------------------------------------------

grDevices::pdf(file = figure_file, width = 14, height = 18, useDingbats = FALSE)

old_par <- par(no.readonly = TRUE)
on.exit({
  par(old_par)
  grDevices::dev.off()
}, add = TRUE)

par(
  mfrow = c(n_pred, n_resp),
  mar = c(0.0, 1, 0.0, 1.5),
  oma = c(3.2, 3.2, 3.0, 0.2),
  mgp = c(1.6, 0.2, 0),
  tcl = -0.3,
  cex.axis = 1.1,
  lwd = 0.8
)

pad_frac <- 0.2
edge_frac <- 0.2
min_yticks <- 2L

for (predictor_idx in seq_len(n_pred)) {
  for (response_idx in seq_len(n_resp)) {
    panel_idx <- (response_idx - 1L) * n_pred + predictor_idx

    x <- x_list[[panel_idx]]
    mean_curve <- mean_list[[panel_idx]]
    lower_curve <- lower_list[[panel_idx]]
    upper_curve <- upper_list[[panel_idx]]
    ymin <- ymin_vec[panel_idx]
    ymax <- ymax_vec[panel_idx]

    is_top <- predictor_idx == 1L
    is_bottom <- predictor_idx == n_pred
    is_left <- response_idx == 1L

    # First panel: show two diabetes cutpoint curves and their CI overlap.
    if (panel_idx == 1L && !is.null(mean_list_cut[[1L]])) {
      mean_cut <- mean_list_cut[[1L]]
      lower_cut <- lower_list_cut[[1L]]
      upper_cut <- upper_list_cut[[1L]]

      ymin <- min(lower_curve, lower_cut, na.rm = TRUE)
      ymax <- max(upper_curve, upper_cut, na.rm = TRUE)
      yr <- ymax - ymin
      if (!is.finite(yr) || yr <= 0) yr <- 1
      ymin <- ymin - 0.08 * yr
      ymax <- ymax + 0.08 * yr

      plot(
        x, mean_curve, type = "n",
        xlab = "", ylab = "", main = "",
        ylim = c(ymin, ymax),
        xaxt = "n", yaxt = "n"
      )

      draw_band(x, lower_curve, upper_curve, grDevices::adjustcolor("grey70", alpha.f = 0.55))
      draw_band(x, lower_cut, upper_cut, grDevices::adjustcolor("grey70", alpha.f = 0.55))

      overlap_lower <- pmax(lower_curve, lower_cut)
      overlap_upper <- pmin(upper_curve, upper_cut)
      draw_overlap_band(
        x, overlap_lower, overlap_upper,
        idx = overlap_lower < overlap_upper,
        col = grDevices::adjustcolor("grey50", alpha.f = 0.75)
      )

      if (ymin <= 0 && ymax >= 0) {
        abline(h = 0, lty = "dashed", col = "red", lwd = 1)
      }

      lines(x, mean_curve, lty = 2, col = "blue", lwd = 1.2)
      lines(x, mean_cut, lty = 2, col = "blue", lwd = 1.2)
      box()
    } else {
      plot(
        x, mean_curve, type = "n",
        xlab = "", ylab = "", main = "",
        ylim = c(ymin, ymax),
        xaxt = "n", yaxt = "n"
      )

      draw_band(x, lower_curve, upper_curve, col = "grey85")

      if (ymin <= 0 && ymax >= 0) {
        abline(h = 0, lty = "dashed", col = "red", lwd = 1)
      }

      lines(x, mean_curve, lty = 2, col = "blue", lwd = 1)
      box()
    }

    dec_min_here <- 2L
    dec_max_here <- 6L
    if (response_idx == n_resp && predictor_idx %in% c(2L, 3L, 9L)) {
      dec_min_here <- 3L
      dec_max_here <- 3L
    }

    ytick <- make_y_ticks(
      ymin, ymax,
      pad_frac = pad_frac,
      edge_frac = edge_frac,
      n_pretty = 2,
      min_ticks = min_yticks,
      dec_min = dec_min_here,
      dec_max = dec_max_here
    )

    old_mgp <- par("mgp")
    par(mgp = c(old_mgp[1L], 0.5, old_mgp[3L]))
    axis(2, at = ytick$at, labels = ytick$lab, cex.axis = 1.1, line = 0)
    par(mgp = old_mgp)

    if (is_bottom) {
      old_mgp <- par("mgp")
      par(mgp = c(old_mgp[1L], 0.6, old_mgp[3L]))
      axis(1, at = age_ticks_scaled, labels = age_ticks, cex.axis = 1.2, line = 0)
      par(mgp = old_mgp)
    }

    if (is_left) {
      at_row <- 1 - (predictor_idx - 0.5) / n_pred
      mtext(
        predictors[predictor_idx], side = 2, outer = TRUE,
        at = at_row, line = 1.0, cex = 0.9
      )
    }

    if (is_top) {
      at_col <- (response_idx - 0.5) / n_resp
      mtext(
        resp_lbl_plotmath[response_idx], side = 3, outer = TRUE,
        at = at_col, line = 0.3, cex = 1.0
      )
    }
  }
}

mtext("Age", side = 1, outer = TRUE, line = 2.0, cex = 1.2)

cat("Saved marginal-effect figure:", figure_file, "\n")

# -----------------------------------------------------------------------------
# Correlation plot for BRFSS 2023
# -----------------------------------------------------------------------------

# Required objects from the previous script:
#   brfss_data, output, fit_settings
#   num_basis2, H, burnin, X_new may already exist from marginal plotting.
#
# This block computes the posterior mean and 95% credible intervals of the
# age-varying copula-latent correlations and draws the final correlation plot.

# Packages --------------------------------------------------------------------

library(matrixStats)

# Basic settings ---------------------------------------------------------------

responseType <- brfss_data$responseType
covT <- as.numeric(brfss_data$X)
age_raw <- brfss_data$age_raw

m_resp <- length(responseType)
stopifnot(m_resp == 8L)

if (!exists("num_basis2")) {
  num_basis2 <- fit_settings$num_basis2
}

if (!exists("H")) {
  H <- fit_settings$initH
}

if (!exists("X_new")) {
  X_new <- seq(min(covT), max(covT), length.out = 1000L)
}

n_grid <- length(X_new)

if (!exists("burnin")) {
  n_saved <- nrow(output$posterior_samples$alpha_h)
  burnin <- if (!is.null(fit_settings$burnin) && fit_settings$burnin > 0L) {
    (fit_settings$burnin + 1L):n_saved
  } else {
    seq_len(n_saved)
  }
}

# Extract posterior samples ----------------------------------------------------

alpha_h <- output$posterior_samples$alpha_h
beta_x  <- output$posterior_samples$beta_x

if (!is.null(output$posterior_samples$Rmat)) {
  R_list <- output$posterior_samples$Rmat
} else if (!is.null(output$Rmat)) {
  R_list <- output$Rmat
} else if (!is.null(output$Rsave)) {
  R_list <- output$Rsave
} else {
  stop("Cannot find posterior correlation matrices in output.")
}

alpha_h <- as.matrix(alpha_h)
beta_x  <- as.matrix(beta_x)

# Build PSBP basis on the plotting grid ---------------------------------------

W_x <- as.matrix(
  NaturalCubicBasis(
    X_new,
    select_knots(covT, num_basis2 - 1L)
  )
)

if (ncol(W_x) != num_basis2) {
  stop("ncol(W_x) does not match num_basis2.")
}

for (j in seq_len(ncol(W_x))) {
  W_x[, j] <- W_x[, j] - mean(W_x[, j])
}

# Helper: compute PSBP weights at X_new for one posterior draw -----------------

compute_psbp_weights <- function(iter_id, alpha_h, beta_x, W_x, H) {
  
  n_grid <- nrow(W_x)
  q <- ncol(W_x)
  
  pi_mat <- matrix(NA_real_, nrow = n_grid, ncol = H)
  remaining <- rep(1, n_grid)
  
  for (h in seq_len(H - 1L)) {
    cols_h <- ((h - 1L) * q + 1L):(h * q)
    eta_h <- as.numeric(alpha_h[iter_id, h] + W_x %*% beta_x[iter_id, cols_h])
    
    stick_h <- pnorm(eta_h)
    pi_mat[, h] <- stick_h * remaining
    remaining <- remaining * (1 - stick_h)
  }
  
  pi_mat[, H] <- remaining
  pi_mat
}

# Helper: extract one R matrix block ------------------------------------------

get_R_block <- function(R_list, iter_id, m_resp, H) {
  
  R_iter <- R_list[[iter_id]]
  
  if (is.list(R_iter)) {
    R_mat <- do.call(
      cbind,
      lapply(R_iter, function(Rh) as.vector(as.matrix(Rh)))
    )
  } else {
    R_mat <- as.matrix(R_iter)
  }
  
  if (nrow(R_mat) != m_resp * m_resp) {
    stop("Each R block must have m_resp^2 rows.")
  }
  
  if (ncol(R_mat) < H) {
    stop("Each R block must contain at least H columns/components.")
  }
  
  R_mat[, seq_len(H), drop = FALSE]
}

# Compute posterior correlation curves ----------------------------------------

responses <- c(
  "Diabetes", "HighBP", "HighChol", "Stroke",
  "HeartDis", "Arthritis", "Asthma", "LogBMI"
)

pair_matrix <- t(combn(seq_len(m_resp), 2L))
n_pairs <- nrow(pair_matrix)

# Column-major vector indices for R[response_a, response_b]
pair_vec_idx <- pair_matrix[, 1L] + (pair_matrix[, 2L] - 1L) * m_resp

set.seed(2025)

n_corr_draws <- min(1500L, length(burnin))
draw_ids <- sample(burnin, size = n_corr_draws, replace = FALSE)

corr_array <- array(
  NA_real_,
  dim = c(n_corr_draws, n_grid, n_pairs)
)

for (s in seq_along(draw_ids)) {
  
  iter_id <- draw_ids[s]
  
  pi_mat <- compute_psbp_weights(
    iter_id = iter_id,
    alpha_h = alpha_h,
    beta_x = beta_x,
    W_x = W_x,
    H = H
  )
  
  R_mat <- get_R_block(
    R_list = R_list,
    iter_id = iter_id,
    m_resp = m_resp,
    H = H
  )
  
  # V[, g] = vec(sum_h pi_h(t_g) R_h)
  V <- R_mat %*% t(pi_mat)
  
  corr_array[s, , ] <- t(V[pair_vec_idx, , drop = FALSE])
  
  if (s %% 50L == 0L) {
    cat(s, "posterior correlation draws processed\n")
  }
}

# Summarize posterior correlation curves --------------------------------------

mean_values  <- vector("list", n_pairs)
lower_values <- vector("list", n_pairs)
upper_values <- vector("list", n_pairs)
min_values   <- numeric(n_pairs)
max_values   <- numeric(n_pairs)

for (k in seq_len(n_pairs)) {
  
  corr_k <- corr_array[, , k, drop = FALSE][, , 1L]
  
  mean_k <- colMeans2(corr_k, na.rm = TRUE)
  cred_k <- rowQuantiles(
    t(corr_k),
    probs = c(0.025, 0.975),
    na.rm = TRUE
  )
  
  lower_k <- cred_k[, 1L]
  upper_k <- cred_k[, 2L]
  
  mean_values[[k]]  <- mean_k
  lower_values[[k]] <- lower_k
  upper_values[[k]] <- upper_k
  
  max_values[k] <- max(upper_k) + 0.5 * (max(upper_k) - max(mean_k))
  min_values[k] <- min(lower_k) - 0.5 * (-min(lower_k) + min(mean_k))
}

# Age-axis ticks ---------------------------------------------------------------

xticks_orig <- c(21, 27, 32, 37, 42, 47, 52, 57, 62, 67, 72, 77, 82)
xticks_std <- (xticks_orig - mean(age_raw)) / sd(age_raw)

# Final plot ------------------------------------------------------------------

fig_dir <- file.path("figures", "brfss_2023")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

pair_lbl_plotmath <- apply(pair_matrix, 1L, function(ii) {
  paste0(
    'paste("', responses[ii[1L]], '", " (", y[i', ii[1L], '], ")", ',
    '" vs ", ',
    '"', responses[ii[2L]], '", " (", y[i', ii[2L], '], ")")'
  )
})

pdf(
  file = file.path(fig_dir, "plot_corr.pdf"),
  width = 18,
  height = 11
)

op <- par(
  mfrow    = c(4, 7),
  mar      = c(2.6, 3, 3, 0.5),
  mgp      = c(0.8, 1, 0),
  tcl      = -0.3,
  cex.lab  = 1.5,
  cex.axis = 1.5,
  cex.main = 1.48
)

zero_col <- "red"

for (k in seq_len(n_pairs)) {
  
  x  <- X_new
  m  <- mean_values[[k]]
  lo <- lower_values[[k]]
  hi <- upper_values[[k]]
  
  plot(
    x, m,
    type = "n",
    xlab = NA,
    ylab = NA,
    main = NA,
    axes = FALSE,
    ylim = c(min_values[k], max_values[k])
  )
  
  title(main = parse(text = pair_lbl_plotmath[k]), font.main = 1)
  
  polygon(
    c(x, rev(x)),
    c(hi, rev(lo)),
    col = "grey85",
    border = NA
  )
  
  if (min_values[k] <= 0 && max_values[k] >= 0) {
    abline(h = 0, lty = "dashed", col = zero_col, lwd = 1)
  }
  
  lines(x, m, lwd = 1, lty = 2, col = "blue")
  
  axis(1, at = xticks_std, labels = xticks_orig, cex.axis = 1.5)
  axis(2)
  box()
}

par(op)
dev.off()

cat("Saved correlation plot: ", file.path(fig_dir, "plot_corr.pdf"), "\n")

