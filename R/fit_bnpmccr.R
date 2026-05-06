
.bnpmccr_as_matrix <- function(x) {
  if (is.data.frame(x)) {
    x <- as.matrix(x)
  }
  if (is.vector(x)) {
    x <- matrix(x, ncol = 1L)
  }
  if (!is.matrix(x)) {
    stop("Input must be coercible to a matrix.", call. = FALSE)
  }
  storage.mode(x) <- "double"
  x
}

.bnpmccr_scale_nonbinary <- function(x) {
  out <- x
  for (j in seq_len(ncol(out))) {
    vals <- unique(stats::na.omit(out[, j]))
    if (length(vals) != 2L) {
      out[, j] <- as.numeric(stats::scale(out[, j]))
    }
  }
  out
}

.bnpmccr_center_columns <- function(x) {
  sweep(x, 2L, colMeans(x), FUN = "-")
}

.bnpmccr_stabilize_truncation_bounds <- function(lower, upper) {
  replace_idx <- union(which(upper == lower), which(lower == Inf))
  if (length(replace_idx) > 0L) {
    lower[replace_idx] <- stats::qnorm(0.9999999999999)
    upper[replace_idx] <- stats::qnorm(0.99999999999999)
  }
  list(lower = lower, upper = upper)
}

.bnpmccr_update_g_scales <- function(g, gamma_block, beta_block, varying_design_matrix, num_basis, p, response_index, n) {
  varying_gamma <- gamma_block[(p + 1L):length(gamma_block)]
  varying_beta <- beta_block[(p + 1L):length(beta_block)]

  for (q in seq_len(p)) {
    block_start <- (num_basis - 1L) * (q - 1L) + 1L
    block_end <- (num_basis - 1L) * q
    gamma_q <- varying_gamma[block_start:block_end]
    beta_q <- varying_beta[block_start:block_end]
    fitted_q <- varying_design_matrix[, block_start:block_end, drop = FALSE] %*% beta_q

    g[(response_index - 1L) * p + q] <- 1 / stats::rgamma(
      1,
      shape = sum(gamma_q) / 2 + 0.5,
      rate = as.numeric(crossprod(fitted_q)) / 2 + n / 2
    )
  }

  g
}

.bnpmccr_rinvgamma <- function(n, shape, scale) {
  1 / stats::rgamma(n = n, shape = shape, rate = scale)
}

.bnpmccr_log_prior <- function(omega, M, num_knots) {
  num_knots * log(1 - omega) + log(omega) - log(choose(M, num_knots))
}

.bnpmccr_log_dmvnorm <- function(x, mean, sigma) {
  if (length(x) == 0L) {
    return(0)
  }
  mvtnorm::dmvnorm(as.vector(x), mean = mean, sigma = sigma, log = TRUE)
}

.bnpmccr_log_scale_prior_ratio <- function(proposed_log_value, current_log_value, sd) {
  stats::dnorm(proposed_log_value, mean = 0, sd = sd, log = TRUE) - proposed_log_value -
    stats::dnorm(current_log_value, mean = 0, sd = sd, log = TRUE) + current_log_value
}

.bnpmccr_normal_prior_ratio <- function(proposed_value, current_value, sd) {
  sum(stats::dnorm(proposed_value, mean = 0, sd = sd, log = TRUE)) -
    sum(stats::dnorm(current_value, mean = 0, sd = sd, log = TRUE))
}

.bnpmccr_adapt_positive_tuning <- function(
    current_value,
    acceptance_rate,
    target_rate,
    lower_tolerance,
    upper_tolerance,
    log_step,
    min_value,
    max_value
) {
  updated_log_value <- log(current_value)

  if (acceptance_rate < target_rate - lower_tolerance) {
    updated_log_value <- updated_log_value - log_step
  } else if (acceptance_rate > target_rate + upper_tolerance) {
    updated_log_value <- updated_log_value + log_step
  }

  updated_value <- exp(updated_log_value)
  max(min(updated_value, max_value), min_value)
}

.bnpmccr_adapt_rw_proposal_sd <- function(
    current_sd,
    acceptance_rate,
    target_rate,
    tolerance,
    log_step,
    min_sd,
    max_sd
) {
  .bnpmccr_adapt_positive_tuning(
    current_value = current_sd,
    acceptance_rate = acceptance_rate,
    target_rate = target_rate,
    lower_tolerance = tolerance,
    upper_tolerance = tolerance,
    log_step = log_step,
    min_value = min_sd,
    max_value = max_sd
  )
}

.bnpmccr_adaptation_batch_size <- function(iteration, burnin, batch_size, start_iteration = 1L) {
  if (iteration < start_iteration || iteration > burnin) {
    return(0L)
  }

  iter_since_start <- iteration - start_iteration + 1L
  if (iter_since_start %% batch_size == 0L) {
    return(as.integer(batch_size))
  }
  if (iteration == burnin) {
    return(as.integer(iter_since_start %% batch_size))
  }
  0L
}

.bnpmccr_softmax_rows <- function(x) {
  row_max <- apply(x, 1L, max)
  z <- exp(x - row_max)
  z / rowSums(z)
}

.bnpmccr_require_backend <- function() {
  needed <- c(
    "NaturalCubicBasis",
    "select_knots",
    "sampleGammaCpp2",
    "updateCoefficients2",
    "updateZCurrGamma",
    "updateDiffLike",
    "updateZCurr",
    "compute_Rh",
    "compute_zeta",
    "computeZetaGamma",
    "computeZetaGamma2",
    "update_Z_star",
    "updateDPParameters"
  )

  ns <- asNamespace("bnpmccr")

  missing <- needed[!vapply(
    needed,
    function(f) {
      exists(f, envir = ns, mode = "function", inherits = FALSE)
    },
    logical(1)
  )]

  if (length(missing) > 0L) {
    stop(
      "Missing required backend functions in the package namespace: ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }

  invisible(TRUE)
}

#' Fit Bayesian nonparametric multivariate conditional copula regression
#'
#' Package-oriented refactoring of the original implementation of the
#' Bayesian nonparametric multivariate conditional copula regression model.
#' The function assumes that the spline basis utilities (`NaturalCubicBasis`,
#' `select_knots`) and all backend update functions are already available in
#' the package namespace.
#'
#' @param Y Response matrix with `n` rows and `m` columns.
#' @param covT Numeric covariate used for varying coefficients / PSBP weights.
#' @param predX Predictor matrix (without intercept).
#' @param responseType Character vector of length `ncol(Y)`.
#' @param niter Total number of MCMC iterations.
#' @param burnin Burn-in iterations.
#' @param tau Global prior scale for spline coefficients.
#' @param nTrial Binomial trial sizes, one entry per binomial margin.
#' @param negBinomParInit Initial size parameter(s) for negative-binomial margins.
#' @param AlphaInit Initial shape parameter(s) for gamma margins.
#' @param mean_alpha Mean used for PSBP intercept initialization.
#' @param std Logical; if `TRUE`, non-binary predictors and Gaussian responses
#'   are standardized as in the original code.
#' @param num_basis Number of marginal basis terms per predictor block,
#'   including the linear term.
#' @param initH Truncation level for the PSBP mixture.
#' @param num_basis2 Number of basis terms for the PSBP weight functions.
#' @param tau2 Prior scale for PSBP spline coefficients.
#' @param sdbeta1 Proposal scale for PSBP spline coefficients.
#' @param initAlpha Prior variance for fixed effects.
#' @param gaussian_variance_prior_eps Small positive hyperparameter `\epsilon`
#'   in the Gaussian marginal prior `IG(\epsilon, \epsilon)`.
#' @param log_additional_param_prior_sd Prior standard deviation for the normal
#'   prior on `log(\xi_k)` in gamma and negative-binomial marginals.
#' @param ordinal_threshold_prior_sd Prior standard deviation for the normal
#'   prior on transformed ordinal threshold spacings.
#' @param adapt_ordinal_cutpoint_proposals Logical; if `TRUE`, adapt the
#'   random-walk proposal standard deviation for ordinal cutpoint updates during
#'   burn-in.
#' @param ordinal_cutpoint_adapt_batch Batch size used to compute observed
#'   acceptance rates for ordinal cutpoint proposal adaptation.
#' @param ordinal_cutpoint_adapt_delay Number of initial iterations kept at the
#'   initial ordinal cutpoint proposal scale before adaptation begins.
#' @param ordinal_cutpoint_proposal_sd_init Initial random-walk proposal
#'   standard deviation for transformed ordinal cutpoint spacings.
#' @param ordinal_cutpoint_target_accept Target acceptance rate used when
#'   adapting ordinal cutpoint proposals during burn-in.
#' @param ordinal_cutpoint_accept_tol Symmetric tolerance around the target
#'   acceptance rate before ordinal cutpoint adaptation is triggered.
#' @param ordinal_cutpoint_adapt_log_step Log-scale multiplicative step size
#'   used when adapting the ordinal cutpoint proposal standard deviation.
#' @param ordinal_cutpoint_min_proposal_sd,ordinal_cutpoint_max_proposal_sd
#'   Lower and upper bounds for the adaptive ordinal cutpoint proposal
#'   standard deviation.
#' @param adapt_additional_params Logical; if `TRUE`, adapt the log-scale
#'   random-walk proposal standard deviations for gamma and negative-binomial
#'   additional parameters during burn-in. Binomial margins do not have a
#'   model-specific additional parameter `\xi_k` in the paper, so no
#'   additional-parameter proposal is used for them.
#' @param additional_param_adapt_batch Batch size used to compute observed
#'   acceptance rates for additional-parameter proposal adaptation.
#' @param additional_param_adapt_delay Number of initial iterations kept at the
#'   initial additional-parameter proposal scales before adaptation begins.
#' @param gamma_additional_proposal_sd_init Initial log-scale random-walk
#'   proposal standard deviation for gamma shape parameters.
#' @param negbin_additional_proposal_sd_init Initial log-scale random-walk
#'   proposal standard deviation for negative-binomial size parameters.
#'   Binomial margins do not have a separate additional parameter in the
#'   current model, so no analogous proposal adaptation is applied to them.
#' @param gamma_additional_target_accept Target acceptance rate used when
#'   adapting gamma additional-parameter proposals during burn-in.
#' @param negbin_additional_target_accept Target acceptance rate used when
#'   adapting negative-binomial additional-parameter proposals during burn-in.
#' @param additional_param_accept_tol Symmetric tolerance around the target
#'   acceptance rate before proposal adaptation is triggered.
#' @param additional_param_adapt_log_step Log-scale multiplicative step size
#'   used when adapting additional-parameter proposal standard deviations.
#' @param additional_param_min_proposal_sd,additional_param_max_proposal_sd
#'   Lower and upper bounds for the adaptive proposal standard deviations.
#' @param adapt_coefficient_proposals Logical; if `TRUE`, adapt the proposal
#'   variance parameter used in the Metropolis-Hastings coefficient updates for
#'   gamma, binomial, and negative-binomial marginals during burn-in.
#' @param coefficient_adapt_batch Batch size used to compute observed
#'   acceptance rates for coefficient-proposal adaptation.
#' @param coefficient_adapt_delay Number of initial MCMC iterations kept at the
#'   initial coefficient proposal variance before adaptation begins. If set to
#'   `300`, adaptation starts after iteration 300 and is still stopped at
#'   `burnin`.
#' @param coefficient_adapt_log_step Log-scale multiplicative step size used
#'   when adapting the coefficient proposal variance parameter.
#' @param coefficient_proposal_var_min,coefficient_proposal_var_max Lower and
#'   upper bounds for the adaptive coefficient proposal variance.
#' @param gamma_coefficient_proposal_var_init,binomial_coefficient_proposal_var_init,
#'   negbin_coefficient_proposal_var_init Initial proposal variance values used
#'   in coefficient updates for gamma, binomial, and negative-binomial
#'   marginals, respectively.
#' @param gamma_coefficient_target_accept,binomial_coefficient_target_accept,
#'   negbin_coefficient_target_accept Target acceptance rates for the
#'   burn-in-only coefficient-proposal adaptation by response family.
#' @param gamma_coefficient_accept_tol,binomial_coefficient_accept_tol,
#'   negbin_coefficient_accept_tol Length-2 numeric vectors giving the lower
#'   and upper tolerances around the target acceptance rate. For example,
#'   `c(0.02, 0.05)` means adapt downward below `target - 0.02` and upward
#'   above `target + 0.05`.
#' @param negbin_coefficient_proposal_var_reset Optional variance value applied
#'   once, at the end of the initial delay period, for negative-binomial
#'   coefficient proposals. The default `2.5` preserves the legacy reset; set
#'   this to `NULL` to disable it.
#' @param omegaM,omegaC Knot-selection hyperparameters for the marginal and
#'   copula components.
#' @param verbose Logical; print progress messages when `TRUE`.
#' @param progress_every Progress frequency.
#' @param keep_burnin Logical; if `FALSE`, only post-burn-in draws are returned.
#'
#' @return An object of class `"bnpmccr_fit"`.
#' @export

fit_bnpmccr <- function(
    Y,
    covT,
    predX,
    responseType,
    niter,
    burnin,
    tau = 1,
    nTrial = NULL,
    negBinomParInit = NULL,
    AlphaInit = NULL,
    mean_alpha = 1,
    std = TRUE,
    num_basis,
    initH = 10,
    num_basis2 = 10,
    tau2 = 10,
    sdbeta1 = 0.1,
    initAlpha = 1e10,
    gaussian_variance_prior_eps = 0.001,
    log_additional_param_prior_sd = 1000,
    ordinal_threshold_prior_sd = 1000,
    adapt_ordinal_cutpoint_proposals = TRUE,
    ordinal_cutpoint_adapt_batch = 50L,
    ordinal_cutpoint_adapt_delay = 300L,
    ordinal_cutpoint_proposal_sd_init = sqrt(0.1),
    ordinal_cutpoint_target_accept = 0.30,
    ordinal_cutpoint_accept_tol = 0.05,
    ordinal_cutpoint_adapt_log_step = 0.05,
    ordinal_cutpoint_min_proposal_sd = 0.01,
    ordinal_cutpoint_max_proposal_sd = 2,
    adapt_additional_params = TRUE,
    additional_param_adapt_batch = 50L,
    additional_param_adapt_delay = 300L,
    gamma_additional_proposal_sd_init = 0.2,
    negbin_additional_proposal_sd_init = 0.05,
    gamma_additional_target_accept = 0.35,
    negbin_additional_target_accept = 0.35,
    additional_param_accept_tol = 0.05,
    additional_param_adapt_log_step = 0.1,
    additional_param_min_proposal_sd = 0.01,
    additional_param_max_proposal_sd = 2,
    adapt_coefficient_proposals = TRUE,
    coefficient_adapt_batch = 50L,
    coefficient_adapt_delay = 300L,
    coefficient_adapt_log_step = 0.05,
    coefficient_proposal_var_min = 0.5,
    coefficient_proposal_var_max = 10,
    gamma_coefficient_proposal_var_init = 1.5,
    binomial_coefficient_proposal_var_init = 1.4,
    negbin_coefficient_proposal_var_init = 4,
    gamma_coefficient_target_accept = 0.20,
    gamma_coefficient_accept_tol = c(0.02, 0.05),
    binomial_coefficient_target_accept = 0.07,
    binomial_coefficient_accept_tol = c(0.03, 0.03),
    negbin_coefficient_target_accept = 0.04,
    negbin_coefficient_accept_tol = c(0.01, 0.03),
    negbin_coefficient_proposal_var_reset = 2.5,
    omegaM = 0.1,
    omegaC = 0.1,
    verbose = interactive(),
    progress_every = 100L,
    keep_burnin = TRUE
){

  .bnpmccr_require_backend()

  response_matrix <- .bnpmccr_as_matrix(Y)
  varying_covariate <- as.numeric(covT)
  predictor_matrix <- .bnpmccr_as_matrix(predX)

  n <- nrow(response_matrix)
  m <- ncol(response_matrix)
  if (length(varying_covariate) != n) {
    stop("`covT` must have length equal to nrow(Y).", call. = FALSE)
  }
  if (nrow(predictor_matrix) != n) {
    stop("`predX` must have nrow equal to nrow(Y).", call. = FALSE)
  }
  if (length(responseType) != m) {
    stop("`responseType` must have length equal to ncol(Y).", call. = FALSE)
  }
  allowed_response <- c("Gaussian", "gamma", "binary", "binomial", "negative binomial", "ordinal")
  bad_response <- setdiff(unique(responseType), allowed_response)
  if (length(bad_response) > 0L) {
    stop("Unsupported response type(s): ", paste(bad_response, collapse = ", "), call. = FALSE)
  }
  if (num_basis < 2L) {
    stop("`num_basis` must be at least 2.", call. = FALSE)
  }
  if (num_basis2 < 2L) {
    stop("`num_basis2` must be at least 2.", call. = FALSE)
  }
  if (burnin >= niter) {
    stop("`burnin` must be smaller than `niter`.", call. = FALSE)
  }
  if (!is.numeric(progress_every) || length(progress_every) != 1L || !is.finite(progress_every) || progress_every <= 0) {
    stop("`progress_every` must be a positive scalar.", call. = FALSE)
  }
  progress_every <- as.integer(progress_every)
  if (!is.numeric(gaussian_variance_prior_eps) || length(gaussian_variance_prior_eps) != 1L ||
      !is.finite(gaussian_variance_prior_eps) || gaussian_variance_prior_eps <= 0) {
    stop("`gaussian_variance_prior_eps` must be a positive finite scalar.", call. = FALSE)
  }
  if (!is.numeric(log_additional_param_prior_sd) || length(log_additional_param_prior_sd) != 1L ||
      !is.finite(log_additional_param_prior_sd) || log_additional_param_prior_sd <= 0) {
    stop("`log_additional_param_prior_sd` must be a positive finite scalar.", call. = FALSE)
  }
  if (!is.numeric(ordinal_threshold_prior_sd) || length(ordinal_threshold_prior_sd) != 1L ||
      !is.finite(ordinal_threshold_prior_sd) || ordinal_threshold_prior_sd <= 0) {
    stop("`ordinal_threshold_prior_sd` must be a positive finite scalar.", call. = FALSE)
  }
  if (!is.logical(adapt_ordinal_cutpoint_proposals) || length(adapt_ordinal_cutpoint_proposals) != 1L || is.na(adapt_ordinal_cutpoint_proposals)) {
    stop("`adapt_ordinal_cutpoint_proposals` must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.numeric(ordinal_cutpoint_adapt_batch) || length(ordinal_cutpoint_adapt_batch) != 1L ||
      !is.finite(ordinal_cutpoint_adapt_batch) || ordinal_cutpoint_adapt_batch <= 0) {
    stop("`ordinal_cutpoint_adapt_batch` must be a positive scalar.", call. = FALSE)
  }
  ordinal_cutpoint_adapt_batch <- as.integer(ordinal_cutpoint_adapt_batch)
  if (!is.numeric(ordinal_cutpoint_adapt_delay) || length(ordinal_cutpoint_adapt_delay) != 1L ||
      !is.finite(ordinal_cutpoint_adapt_delay) || ordinal_cutpoint_adapt_delay < 0) {
    stop("`ordinal_cutpoint_adapt_delay` must be a non-negative scalar.", call. = FALSE)
  }
  ordinal_cutpoint_adapt_delay <- as.integer(ordinal_cutpoint_adapt_delay)
  if (!is.numeric(ordinal_cutpoint_proposal_sd_init) || length(ordinal_cutpoint_proposal_sd_init) != 1L ||
      !is.finite(ordinal_cutpoint_proposal_sd_init) || ordinal_cutpoint_proposal_sd_init <= 0) {
    stop("`ordinal_cutpoint_proposal_sd_init` must be a positive finite scalar.", call. = FALSE)
  }
  if (!is.numeric(ordinal_cutpoint_target_accept) || length(ordinal_cutpoint_target_accept) != 1L ||
      !is.finite(ordinal_cutpoint_target_accept) || ordinal_cutpoint_target_accept <= 0 || ordinal_cutpoint_target_accept >= 1) {
    stop("`ordinal_cutpoint_target_accept` must be a scalar in (0, 1).", call. = FALSE)
  }
  if (!is.numeric(ordinal_cutpoint_accept_tol) || length(ordinal_cutpoint_accept_tol) != 1L ||
      !is.finite(ordinal_cutpoint_accept_tol) || ordinal_cutpoint_accept_tol < 0) {
    stop("`ordinal_cutpoint_accept_tol` must be a non-negative finite scalar.", call. = FALSE)
  }
  if (!is.numeric(ordinal_cutpoint_adapt_log_step) || length(ordinal_cutpoint_adapt_log_step) != 1L ||
      !is.finite(ordinal_cutpoint_adapt_log_step) || ordinal_cutpoint_adapt_log_step <= 0) {
    stop("`ordinal_cutpoint_adapt_log_step` must be a positive finite scalar.", call. = FALSE)
  }
  if (!is.numeric(ordinal_cutpoint_min_proposal_sd) || length(ordinal_cutpoint_min_proposal_sd) != 1L ||
      !is.finite(ordinal_cutpoint_min_proposal_sd) || ordinal_cutpoint_min_proposal_sd <= 0) {
    stop("`ordinal_cutpoint_min_proposal_sd` must be a positive finite scalar.", call. = FALSE)
  }
  if (!is.numeric(ordinal_cutpoint_max_proposal_sd) || length(ordinal_cutpoint_max_proposal_sd) != 1L ||
      !is.finite(ordinal_cutpoint_max_proposal_sd) || ordinal_cutpoint_max_proposal_sd <= 0) {
    stop("`ordinal_cutpoint_max_proposal_sd` must be a positive finite scalar.", call. = FALSE)
  }
  if (ordinal_cutpoint_min_proposal_sd >= ordinal_cutpoint_max_proposal_sd) {
    stop("`ordinal_cutpoint_min_proposal_sd` must be smaller than `ordinal_cutpoint_max_proposal_sd`.", call. = FALSE)
  }
  if (!is.logical(adapt_additional_params) || length(adapt_additional_params) != 1L || is.na(adapt_additional_params)) {
    stop("`adapt_additional_params` must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.numeric(additional_param_adapt_batch) || length(additional_param_adapt_batch) != 1L ||
      !is.finite(additional_param_adapt_batch) || additional_param_adapt_batch <= 0) {
    stop("`additional_param_adapt_batch` must be a positive scalar.", call. = FALSE)
  }
  additional_param_adapt_batch <- as.integer(additional_param_adapt_batch)
  if (!is.numeric(additional_param_adapt_delay) || length(additional_param_adapt_delay) != 1L ||
      !is.finite(additional_param_adapt_delay) || additional_param_adapt_delay < 0) {
    stop("`additional_param_adapt_delay` must be a non-negative scalar.", call. = FALSE)
  }
  additional_param_adapt_delay <- as.integer(additional_param_adapt_delay)
  if (!is.numeric(gamma_additional_proposal_sd_init) || length(gamma_additional_proposal_sd_init) != 1L ||
      !is.finite(gamma_additional_proposal_sd_init) || gamma_additional_proposal_sd_init <= 0) {
    stop("`gamma_additional_proposal_sd_init` must be a positive finite scalar.", call. = FALSE)
  }
  if (!is.numeric(negbin_additional_proposal_sd_init) || length(negbin_additional_proposal_sd_init) != 1L ||
      !is.finite(negbin_additional_proposal_sd_init) || negbin_additional_proposal_sd_init <= 0) {
    stop("`negbin_additional_proposal_sd_init` must be a positive finite scalar.", call. = FALSE)
  }
  if (!is.numeric(gamma_additional_target_accept) || length(gamma_additional_target_accept) != 1L ||
      !is.finite(gamma_additional_target_accept) || gamma_additional_target_accept <= 0 || gamma_additional_target_accept >= 1) {
    stop("`gamma_additional_target_accept` must be a scalar in (0, 1).", call. = FALSE)
  }
  if (!is.numeric(negbin_additional_target_accept) || length(negbin_additional_target_accept) != 1L ||
      !is.finite(negbin_additional_target_accept) || negbin_additional_target_accept <= 0 || negbin_additional_target_accept >= 1) {
    stop("`negbin_additional_target_accept` must be a scalar in (0, 1).", call. = FALSE)
  }
  if (!is.numeric(additional_param_accept_tol) || length(additional_param_accept_tol) != 1L ||
      !is.finite(additional_param_accept_tol) || additional_param_accept_tol < 0) {
    stop("`additional_param_accept_tol` must be a non-negative finite scalar.", call. = FALSE)
  }
  if (!is.numeric(additional_param_adapt_log_step) || length(additional_param_adapt_log_step) != 1L ||
      !is.finite(additional_param_adapt_log_step) || additional_param_adapt_log_step <= 0) {
    stop("`additional_param_adapt_log_step` must be a positive finite scalar.", call. = FALSE)
  }
  if (!is.numeric(additional_param_min_proposal_sd) || length(additional_param_min_proposal_sd) != 1L ||
      !is.finite(additional_param_min_proposal_sd) || additional_param_min_proposal_sd <= 0) {
    stop("`additional_param_min_proposal_sd` must be a positive finite scalar.", call. = FALSE)
  }
  if (!is.numeric(additional_param_max_proposal_sd) || length(additional_param_max_proposal_sd) != 1L ||
      !is.finite(additional_param_max_proposal_sd) || additional_param_max_proposal_sd <= 0) {
    stop("`additional_param_max_proposal_sd` must be a positive finite scalar.", call. = FALSE)
  }
  if (additional_param_min_proposal_sd >= additional_param_max_proposal_sd) {
    stop("`additional_param_min_proposal_sd` must be smaller than `additional_param_max_proposal_sd`.", call. = FALSE)
  }
  if (!is.logical(adapt_coefficient_proposals) || length(adapt_coefficient_proposals) != 1L || is.na(adapt_coefficient_proposals)) {
    stop("`adapt_coefficient_proposals` must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.numeric(coefficient_adapt_batch) || length(coefficient_adapt_batch) != 1L ||
      !is.finite(coefficient_adapt_batch) || coefficient_adapt_batch <= 0) {
    stop("`coefficient_adapt_batch` must be a positive scalar.", call. = FALSE)
  }
  coefficient_adapt_batch <- as.integer(coefficient_adapt_batch)
  if (!is.numeric(coefficient_adapt_delay) || length(coefficient_adapt_delay) != 1L ||
      !is.finite(coefficient_adapt_delay) || coefficient_adapt_delay < 0) {
    stop("`coefficient_adapt_delay` must be a non-negative scalar.", call. = FALSE)
  }
  coefficient_adapt_delay <- as.integer(coefficient_adapt_delay)
  if (!is.numeric(coefficient_adapt_log_step) || length(coefficient_adapt_log_step) != 1L ||
      !is.finite(coefficient_adapt_log_step) || coefficient_adapt_log_step <= 0) {
    stop("`coefficient_adapt_log_step` must be a positive finite scalar.", call. = FALSE)
  }
  if (!is.numeric(coefficient_proposal_var_min) || length(coefficient_proposal_var_min) != 1L ||
      !is.finite(coefficient_proposal_var_min) || coefficient_proposal_var_min <= 0) {
    stop("`coefficient_proposal_var_min` must be a positive finite scalar.", call. = FALSE)
  }
  if (!is.numeric(coefficient_proposal_var_max) || length(coefficient_proposal_var_max) != 1L ||
      !is.finite(coefficient_proposal_var_max) || coefficient_proposal_var_max <= 0) {
    stop("`coefficient_proposal_var_max` must be a positive finite scalar.", call. = FALSE)
  }
  if (coefficient_proposal_var_min >= coefficient_proposal_var_max) {
    stop("`coefficient_proposal_var_min` must be smaller than `coefficient_proposal_var_max`.", call. = FALSE)
  }
  if (!is.numeric(gamma_coefficient_proposal_var_init) || length(gamma_coefficient_proposal_var_init) != 1L ||
      !is.finite(gamma_coefficient_proposal_var_init) || gamma_coefficient_proposal_var_init <= 0) {
    stop("`gamma_coefficient_proposal_var_init` must be a positive finite scalar.", call. = FALSE)
  }
  if (!is.numeric(binomial_coefficient_proposal_var_init) || length(binomial_coefficient_proposal_var_init) != 1L ||
      !is.finite(binomial_coefficient_proposal_var_init) || binomial_coefficient_proposal_var_init <= 0) {
    stop("`binomial_coefficient_proposal_var_init` must be a positive finite scalar.", call. = FALSE)
  }
  if (!is.numeric(negbin_coefficient_proposal_var_init) || length(negbin_coefficient_proposal_var_init) != 1L ||
      !is.finite(negbin_coefficient_proposal_var_init) || negbin_coefficient_proposal_var_init <= 0) {
    stop("`negbin_coefficient_proposal_var_init` must be a positive finite scalar.", call. = FALSE)
  }
  if (!is.numeric(gamma_coefficient_target_accept) || length(gamma_coefficient_target_accept) != 1L ||
      !is.finite(gamma_coefficient_target_accept) || gamma_coefficient_target_accept <= 0 || gamma_coefficient_target_accept >= 1) {
    stop("`gamma_coefficient_target_accept` must be a scalar in (0, 1).", call. = FALSE)
  }
  if (!is.numeric(binomial_coefficient_target_accept) || length(binomial_coefficient_target_accept) != 1L ||
      !is.finite(binomial_coefficient_target_accept) || binomial_coefficient_target_accept <= 0 || binomial_coefficient_target_accept >= 1) {
    stop("`binomial_coefficient_target_accept` must be a scalar in (0, 1).", call. = FALSE)
  }
  if (!is.numeric(negbin_coefficient_target_accept) || length(negbin_coefficient_target_accept) != 1L ||
      !is.finite(negbin_coefficient_target_accept) || negbin_coefficient_target_accept <= 0 || negbin_coefficient_target_accept >= 1) {
    stop("`negbin_coefficient_target_accept` must be a scalar in (0, 1).", call. = FALSE)
  }
  validate_coefficient_tol <- function(x, arg_name) {
    if (!is.numeric(x) || length(x) != 2L || any(!is.finite(x)) || any(x < 0)) {
      stop(sprintf("`%s` must be a length-2 non-negative numeric vector.", arg_name), call. = FALSE)
    }
    unname(as.numeric(x))
  }
  gamma_coefficient_accept_tol <- validate_coefficient_tol(gamma_coefficient_accept_tol, "gamma_coefficient_accept_tol")
  binomial_coefficient_accept_tol <- validate_coefficient_tol(binomial_coefficient_accept_tol, "binomial_coefficient_accept_tol")
  negbin_coefficient_accept_tol <- validate_coefficient_tol(negbin_coefficient_accept_tol, "negbin_coefficient_accept_tol")
  if (gamma_coefficient_target_accept - gamma_coefficient_accept_tol[1] <= 0 ||
      gamma_coefficient_target_accept + gamma_coefficient_accept_tol[2] >= 1) {
    stop("`gamma_coefficient_target_accept` and `gamma_coefficient_accept_tol` must define bounds inside (0, 1).", call. = FALSE)
  }
  if (binomial_coefficient_target_accept - binomial_coefficient_accept_tol[1] <= 0 ||
      binomial_coefficient_target_accept + binomial_coefficient_accept_tol[2] >= 1) {
    stop("`binomial_coefficient_target_accept` and `binomial_coefficient_accept_tol` must define bounds inside (0, 1).", call. = FALSE)
  }
  if (negbin_coefficient_target_accept - negbin_coefficient_accept_tol[1] <= 0 ||
      negbin_coefficient_target_accept + negbin_coefficient_accept_tol[2] >= 1) {
    stop("`negbin_coefficient_target_accept` and `negbin_coefficient_accept_tol` must define bounds inside (0, 1).", call. = FALSE)
  }
  if (!is.null(negbin_coefficient_proposal_var_reset) &&
      (!is.numeric(negbin_coefficient_proposal_var_reset) || length(negbin_coefficient_proposal_var_reset) != 1L ||
       !is.finite(negbin_coefficient_proposal_var_reset) || negbin_coefficient_proposal_var_reset <= 0)) {
    stop("`negbin_coefficient_proposal_var_reset` must be NULL or a positive finite scalar.", call. = FALSE)
  }

  log_prior <- .bnpmccr_log_prior
  rinvgamma <- .bnpmccr_rinvgamma
  log_scale_prior_ratio <- .bnpmccr_log_scale_prior_ratio
  normal_prior_ratio <- .bnpmccr_normal_prior_ratio
  adapt_positive_tuning <- .bnpmccr_adapt_positive_tuning
  adapt_rw_proposal_sd <- .bnpmccr_adapt_rw_proposal_sd
  adaptation_batch_size <- .bnpmccr_adaptation_batch_size
  dmvnorm <- function(x, mean, sigma, log = TRUE) {
    .bnpmccr_log_dmvnorm(x, mean, sigma)
  }

  original_response_matrix <- response_matrix
  original_predictor_matrix <- predictor_matrix
  if (std) {
    predictor_matrix <- .bnpmccr_scale_nonbinary(predictor_matrix)
  }
  if (std) {
    gaussian_idx <- which(responseType == "Gaussian")
    if (length(gaussian_idx) > 0L) {
      for (k in gaussian_idx) {
        response_matrix[, k] <- as.numeric(stats::scale(response_matrix[, k]))
      }
    }
  }

  Y <- response_matrix
  X <- varying_covariate
  X3 <- predictor_matrix
  negBinomPar <- negBinomParInit
  H <- initH
  Xmat <- cbind(`(Intercept)` = 1, X3)
  blockPriorInv <- diag(1 / initAlpha, ncol(Xmat))

  S <- niter
  p <- ncol(Xmat)
  mC <- sum(responseType == "Gaussian")

  psbp_basis_matrix <- .bnpmccr_center_columns(
    as.matrix(NaturalCubicBasis(X, select_knots(X, num_basis2 - 1L)))
  )
  marginal_basis_matrix <- .bnpmccr_center_columns(
    as.matrix(NaturalCubicBasis(X, select_knots(X, num_basis - 2L)))
  )

  WXB_blocks <- lapply(seq_len(ncol(Xmat)), function(j) Xmat[, j] * marginal_basis_matrix)
  W <- do.call(cbind, WXB_blocks)
  p2 <- ncol(psbp_basis_matrix)
  totalW <- cbind(Xmat, W)

  g <- rep(tau, p * m)
  negind <- which(responseType == "negative binomial")
  gamind <- which(responseType == "gamma")
  binomind <- which(responseType == "binomial")
  coefficient_proposal_var <- rep(1, m)
  if (length(negind) > 0L) {
    coefficient_proposal_var[negind] <- negbin_coefficient_proposal_var_init
  }
  if (length(binomind) > 0L) {
    coefficient_proposal_var[binomind] <- binomial_coefficient_proposal_var_init
  }
  if (length(gamind) > 0L) {
    coefficient_proposal_var[gamind] <- gamma_coefficient_proposal_var_init
  }
  if (!is.null(negbin_coefficient_proposal_var_reset) && coefficient_adapt_delay == 0L && length(negind) > 0L) {
    coefficient_proposal_var[negind] <- negbin_coefficient_proposal_var_reset
  }
  var_z <- matrix(rep(coefficient_proposal_var, each = n), nrow = n, ncol = m)
  coefficient_accept_batch <- numeric(m)
  coefficient_accept_total <- numeric(m)
  coefficient_proposal_total <- numeric(m)
  coefficient_target_accept <- rep(NA_real_, m)
  coefficient_accept_tol_lower <- rep(NA_real_, m)
  coefficient_accept_tol_upper <- rep(NA_real_, m)
  if (length(gamind) > 0L) {
    coefficient_target_accept[gamind] <- gamma_coefficient_target_accept
    coefficient_accept_tol_lower[gamind] <- gamma_coefficient_accept_tol[1]
    coefficient_accept_tol_upper[gamind] <- gamma_coefficient_accept_tol[2]
  }
  if (length(binomind) > 0L) {
    coefficient_target_accept[binomind] <- binomial_coefficient_target_accept
    coefficient_accept_tol_lower[binomind] <- binomial_coefficient_accept_tol[1]
    coefficient_accept_tol_upper[binomind] <- binomial_coefficient_accept_tol[2]
  }
  if (length(negind) > 0L) {
    coefficient_target_accept[negind] <- negbin_coefficient_target_accept
    coefficient_accept_tol_lower[negind] <- negbin_coefficient_accept_tol[1]
    coefficient_accept_tol_upper[negind] <- negbin_coefficient_accept_tol[2]
  }
  additional_param_accept_count <- numeric(m)
  additional_param_accept_total <- numeric(m)
  additional_param_proposal_total <- numeric(m)
  additional_param_proposal_sd <- rep(NA_real_, m)
  additional_param_target_accept <- rep(NA_real_, m)
  if (length(gamind) > 0L) {
    additional_param_proposal_sd[gamind] <- gamma_additional_proposal_sd_init
    additional_param_target_accept[gamind] <- gamma_additional_target_accept
  }
  if (length(negind) > 0L) {
    additional_param_proposal_sd[negind] <- negbin_additional_proposal_sd_init
    additional_param_target_accept[negind] <- negbin_additional_target_accept
  }
  ordind <- which(responseType == "ordinal")
  ordinal_cutpoint_accept_batch <- numeric(m)
  ordinal_cutpoint_accept_total <- numeric(m)
  ordinal_cutpoint_proposal_total <- numeric(m)
  ordinal_cutpoint_proposal_sd <- rep(NA_real_, m)
  if (length(ordind) > 0L) {
    ordinal_cutpoint_proposal_sd[ordind] <- ordinal_cutpoint_proposal_sd_init
  }

  if (length(negind) == 0L) {
    negBinomPar <- numeric(0L)
  } else if (is.null(negBinomPar)) {
    negBinomPar <- rep(1, length(negind))
  } else if (length(negBinomPar) == 1L) {
    negBinomPar <- rep(negBinomPar, length(negind))
  } else if (length(negBinomPar) != length(negind)) {
    stop("`negBinomParInit` must have length equal to the number of negative-binomial margins.", call. = FALSE)
  }

  if (length(gamind) == 0L) {
    gammaAlpha <- numeric(0L)
  } else if (is.null(AlphaInit)) {
    gammaAlpha <- rep(10, length(gamind))
  } else if (length(AlphaInit) == 1L) {
    gammaAlpha <- rep(AlphaInit, length(gamind))
  } else if (length(AlphaInit) != length(gamind)) {
    stop("`AlphaInit` must have length equal to the number of gamma margins.", call. = FALSE)
  } else {
    gammaAlpha <- AlphaInit
  }

  if (length(binomind) == 0L) {
    nTrial <- numeric(0L)
  } else if (is.null(nTrial)) {
    stop("`nTrial` must be supplied for binomial margins.", call. = FALSE)
  } else if (length(nTrial) == 1L) {
    nTrial <- rep(nTrial, length(binomind))
  } else if (length(nTrial) != length(binomind)) {
    stop("`nTrial` must have length equal to the number of binomial margins.", call. = FALSE)
  }



  posterior_samples <- list()
  posterior_samples$beta <- matrix(NA, S, (ncol(W) + ncol(Xmat))*m )
  posterior_samples$beta_x <- matrix(NA, S, p2*H)
  posterior_samples$gamma <- matrix(NA, S, (ncol(W) + ncol(Xmat))*m )
  posterior_samples$g <- matrix(0, S, p*m )
  posterior_samples$g_h <- matrix(0, S, H)
  posterior_samples$numPI <- matrix(NA, S, H)
  posterior_samples$d <- matrix(NA, S, m)
  posterior_samples$alpha_h <- matrix(NA, S, H)
  posterior_samples$neg <- matrix(NA, S, length(negBinomPar))
  Rsave <- list()
  cutpointsave <- list()
  negparsave <- list()
  gammaparsave <- list()
  gh <- rep(tau2, H)
  dpropormat <- matrix(1,n,m)



  gamma_current <- c(rep(1,ncol(Xmat)),rbinom( ncol(W), 1, 0.1)  )
  for(j in 1:ncol(Xmat)){
    gamma_current[(ncol(Xmat) + 1 + (num_basis-1)*(j-1))] <- 1
  }

  gamma_current <- rep( gamma_current , m)


  psbp_gamma_current <- rbinom(p2 * H, 1, 0.1)
  psbp_gamma_current[1L + (seq_len(H) - 1L) * p2] <- 1L

  beta_current <- rep(0, (ncol(Xmat) + ncol(W)) * m)
  beta_current[gamma_current != 0] <- rnorm(sum(gamma_current != 0), 0, 0.01)

  linear_predictor_current <- totalW %*% matrix(
    beta_current,
    nrow = ncol(totalW),
    ncol = m
  )

  psbp_beta_current <- rep(0, p2 * H)
  psbp_beta_current[psbp_gamma_current != 0] <- rnorm(sum(psbp_gamma_current != 0), 0, sdbeta1)
  psbp_beta_matrix <- matrix(psbp_beta_current, nrow = H, ncol = p2, byrow = TRUE)
  psbp_linear_predictor <- psbp_basis_matrix %*% t(psbp_beta_matrix)

  acceptratio <- rep(0, m)


  D <- list()
  RCurr <- list()
  invRCurr <- list()
  invSigma <- list()
  Sigma <- list()
  Wmat <- list()

  for( h in 1:H ){
    D[[h]] <- diag(rep(1,m))
    RCurr[[h]] <- diag(rep(1,m))
    invRCurr[[h]] <- solve(RCurr[[h]])
    invSigma[[h]] <- diag(1/ diag(D[[h]])) %*% invRCurr[[h]] %*% diag(1/ diag(D[[h]]))
    Sigma[[h]] <- solve(matrix(invSigma[[h]],m,m))
  }

  d<-rep(1,m)
  d2 <- rep(1,m)
  cutPoint <- list()
  nCat<- NULL

  for (k in 1 : m)
  {
    if ((responseType[k] == "Gaussian") | (responseType[k] == "binomial") | (responseType[k] == "negative binomial")| (responseType[k] == "gamma"))
    {
      cutPoint[[k]] <- NA
    }
    if (responseType[k] == "binary")
    {
      cutPoint[[k]] <- c(-Inf, 0, Inf)
    }
    if (responseType[k] == "ordinal")
    {
      nOlvs <- max(Y[, k], na.rm = TRUE) + 1
      nCat[k] <- nOlvs
      cutPoint[[k]] <- rep(NA, nOlvs + 1)
      cutPoint[[k]][c(1, 2, nOlvs + 1)] <- c(-Inf, 0, Inf)
      countCat <- 0
      for(kk in 3 : nOlvs)
      {
        cutPoint[[k]][kk] <- 0.5 + countCat
        countCat <- countCat + 1.5
      }
    }
  }


  latent_z <- matrix(0, n, m)
  if (mC > 0)
  {
    latent_z[, which(responseType == "Gaussian")] <- Y[, which(responseType == "Gaussian")] - linear_predictor_current[, which(responseType == "Gaussian")]
    latent_z[, which(responseType == "Gaussian")] <- latent_z[, which(responseType == "Gaussian")]/sqrt(d[which(responseType == "Gaussian")])
  }

  TU <- TL <- matrix(NA, n, m)
  TUProp <- TLProp <- rep(NA, n)

  latent_linear_predictor <- latent_z + linear_predictor_current
  z_star <- matrix(NA_real_, n, max(H - 1L, 0L))

  component_membership <- rep(1,n)
  num_pi <- rep(0, H)
  num_pi2 <- rep(0,H)
  alpha_h <- rnorm(H, mean_alpha, 0.1)


  dpropor <- matrix(1,n,m)

  for(s in 1: S){
    countbinom = countnegbinom = countgamma= 0


    for(k in 1:m){
      start <- (k - 1) * (ncol(Xmat) + ncol(W)) + 1
      end <- k * (ncol(Xmat) + ncol(W ))
      gammaind <- which(gamma_current[start : end] > 0)

      if(responseType[k]=="Gaussian"||responseType[k]=="binary"||responseType[k]=="ordinal" ){

        Rh <- compute_Rh(invRCurr, component_membership, k)

        zeta <- compute_zeta(latent_z, invRCurr, component_membership, k)

        gamma_current  <- sampleGammaCpp2(gamma_current, totalW, Xmat, W, blockPriorInv,
                                          d, Y, latent_linear_predictor, responseType, g,
                                          num_basis, p, k, start-1, end-1, omegaM, Rh, zeta)

        gammaind <- which(gamma_current[start : end] > 0)
        Wgamma <- totalW[,gammaind]

      }else{

        newGamma <- rep(NA, end-start + 1)
        newGamma[1:p] <- rep(1,p)
        ind1 <- 1
        ind2 <- p*num_basis
        subGamma <- gamma_current[start:end]
        subGammaprop <- rep(0, end-start + 1 )
        subGammaprop[1:p] <- rep(1,p)

        lengthGamma <-  end-start + 1 - p*2
        subGamma2 <- rep(NA, lengthGamma)

        l<- 1
        for(q in (ind1+p):ind2){
          if( (q-p) %% (num_basis-1) != 1 ){
            subGamma2[l] <- subGamma[q]
            l <- l+1
          }
        }
        propj <- sample(1:lengthGamma, 1 )
        propj_value <- sample(c(0,1),1,0.5)
        subGamma2[propj] <- propj_value

        l <- 1
        gammaprop1 <- subGamma2

        for(q in (ind1+p):ind2){
          if( (q-p) %% (num_basis-1) != 1 ){
            subGammaprop[q] <- gammaprop1[l]
            l <- l+1
          }else{
            subGammaprop[q] <- 1
          }
        }



        logpriorgamma <- log_prior(omegaM,  M = lengthGamma , num_knots= sum(gammaprop1) ) -
          log_prior(omegaM, M = lengthGamma, num_knots = sum( gamma_current[start:end]) - 2*p  )


        newGamma <- subGammaprop


        gammaind <- which(gamma_current[start : end] > 0)
        gammastar <- which(newGamma > 0)

        if(length(gammaind)>0){

          Wgamma <- totalW[,gammaind]
          B_k <- beta_current[start:end]
          B_kgamma <- B_k[gammaind]
          WBgamma <- Wgamma %*% as.matrix(B_kgamma)

        } else {
          Wgamma <- matrix(0, nrow = n, ncol = 0)
          B_k <- beta_current[start:end]
          B_kgamma <- numeric(0L)
          WBgamma <- rep(0, n)
        }

        Wgammastar <- totalW[,gammastar]
      }




      Rh <- compute_Rh(invRCurr, component_membership, k)

      if (responseType[k] %in% c("gamma", "binomial", "negative binomial")) {
        var_z[, k] <- coefficient_proposal_var[k]
      }

      dpropormat[,which(responseType=="gamma")] <- dpropor[,which(responseType == "gamma")]
      if(responseType[k] %in% c("Gaussian", "binary", "ordinal")){

        coef_results <- updateCoefficients2(k-1,as.character(responseType[k]), d, latent_z, linear_predictor_current,
                                            invRCurr, component_membership, gamma_current, start-1, end-1,
                                            totalW, Xmat, p, num_basis, g,
                                            blockPriorInv, Rh, Y, var_z, W, rep(0,n),
                                            matrix(0, n, 2), countgamma,  gammaAlpha,
                                            countnegbinom, negBinomPar,
                                            countbinom, nTrial, RCurr, dpropormat, s)

        Bprop <- coef_results[[1]]
        latent_linear_predictor[,k] <- coef_results[[2]]


      }else if(responseType[k] %in% c("binomial", "negative binomial")){

        coef_results <- updateCoefficients2(k-1, as.character(responseType[k]), d, latent_z, linear_predictor_current,
                                            invRCurr, component_membership, gamma_current, start-1, end-1,
                                            totalW, Xmat, p, num_basis, g,
                                            blockPriorInv, Rh, Y, var_z, W, newGamma,
                                            Wgammastar, countgamma, gammaAlpha,
                                            countnegbinom, negBinomPar,
                                            countbinom, nTrial, RCurr, dpropormat, s)

        Bprop <- coef_results[[1]]
        muprop_OLD <- coef_results[[2]]
        inv_V_OLD <- coef_results[[3]]
        inv_bdiagmat_OLD <- coef_results[[4]]
        muprop <- coef_results[[5]]
        inv_V <- coef_results[[6]]
        inv_bdiagmat <- coef_results[[7]]
        latent_linear_predictor[,k] <- coef_results[[8]]
        ZTilde_new <- coef_results[[11]]
        countnegbinom <- coef_results[[12]]
        negBinomPar <- coef_results[[13]]
        countbinom <- coef_results[[14]]
        nTrial <- coef_results[[15]]
        TLProp <- coef_results[[16]]
        TUProp <- coef_results[[17]]

      }else{

        coef_results <- updateCoefficients2(k-1, as.character(responseType[k]), d, latent_z, linear_predictor_current,
                                            invRCurr, component_membership, gamma_current, start-1, end-1,
                                            totalW, Xmat, p, num_basis, g,
                                            blockPriorInv, Rh, Y, var_z, W, newGamma,
                                            Wgammastar, countgamma, gammaAlpha,
                                            countnegbinom, negBinomPar,
                                            countbinom, nTrial, RCurr, dpropormat, s)
        Bprop <- coef_results[[1]]
        muprop_OLD <- coef_results[[2]]
        inv_V_OLD <- coef_results[[3]]
        inv_bdiagmat_OLD <- coef_results[[4]]
        muprop <- coef_results[[5]]
        inv_V <- coef_results[[6]]
        inv_bdiagmat <- coef_results[[7]]
        latent_linear_predictor[,k] <- coef_results[[8]]
        countgamma <- coef_results[[10]]
        ZTilde_new <- coef_results[[11]]

      }


      if(responseType[k] %in% c("Gaussian", "binary","ordinal")){
        if (length(gammaind) > 0)
        {
          linear_predictor_proposed<-(as.matrix(totalW[, gammaind]) %*% (Bprop) )[, 1]
        }else {
          linear_predictor_proposed <- rep(0,n)
        }
      }else{

        linear_predictor_proposed<-(as.matrix(totalW[,gammastar]) %*% (Bprop) )[, 1]


      }







      if(responseType[k]=="Gaussian"){
        if(length(gammaind)>0){
          setone <-  (start-1) + which(gamma_current[start:end]>0)
          beta_current[setone]<-Bprop
        }
        setzero <- (start -1) + which(gamma_current[start:end]==0)
        beta_current[setzero]<-0
        WBgamma <- linear_predictor_proposed
        linear_predictor_current[,k] <- WBgamma

        latent_z[,k] <- (Y[,k] - WBgamma)/d[k]


        g <- .bnpmccr_update_g_scales(
          g = g,
          gamma_block = gamma_current[start:end],
          beta_block = beta_current[start:end],
          varying_design_matrix = W,
          num_basis = num_basis,
          p = p,
          response_index = k,
          n = n
        )


        SS <- 0
        zeta2 <- compute_zeta(latent_z, invRCurr, component_membership, k)

        SS <- t(Y[,k] - linear_predictor_current[,k]) %*% (Rh*Y[,k] - Rh*(linear_predictor_current[,k]))
        SS_prop <-  t(Y[,k] - linear_predictor_current[,k] ) %*% (Y[,k] - linear_predictor_current[,k] )
        d2_prop <- rinvgamma(
          1,
          n / 2 + gaussian_variance_prior_eps,
          SS_prop / 2 + gaussian_variance_prior_eps
        )

        LogLikeProp <- (-n/2)*log(d2_prop) -(0.5/d2_prop)*SS -(1/sqrt(d2_prop))*t(Y[,k] - linear_predictor_current[,k]) %*%zeta2
        LogLike <- (-n/2)*log(d2[k]) -(0.5/d2[k])*SS -(1/sqrt(d2[k]))*t(Y[,k] - linear_predictor_current[,k] ) %*%zeta2
        LogPrior <- invgamma::dinvgamma(
          d2_prop,
          gaussian_variance_prior_eps,
          gaussian_variance_prior_eps,
          log = TRUE
        ) - invgamma::dinvgamma(
          d2[k],
          gaussian_variance_prior_eps,
          gaussian_variance_prior_eps,
          log = TRUE
        )
        LogProp <- invgamma::dinvgamma(
          d2[k],
          n / 2 + gaussian_variance_prior_eps,
          SS_prop / 2 + gaussian_variance_prior_eps,
          log = TRUE
        ) - invgamma::dinvgamma(
          d2_prop,
          n / 2 + gaussian_variance_prior_eps,
          SS_prop / 2 + gaussian_variance_prior_eps,
          log = TRUE
        )
        LogAccRatio <- LogLikeProp - LogLike  + LogPrior + LogProp


        if(log(runif(1))<LogAccRatio){
          d2[k] <- d2_prop
          d[k] <- sqrt(d2[k])
        }

        latent_z[,k] <- (Y[,k] - WBgamma)/d[k]


      }else if(responseType[k]=="gamma"){


        latent_z <- updateZCurrGamma(Y, linear_predictor_current,  latent_z, gammaAlpha, countgamma, k)

        linear_predictor_proposed_matrix <- linear_predictor_current
        linear_predictor_proposed_matrix[,k] <- linear_predictor_proposed

        ZCurrprop <-  updateZCurrGamma(Y, linear_predictor_proposed_matrix, latent_z, gammaAlpha, countgamma, k)

        zetaprop <- computeZetaGamma(ZCurrprop, component_membership, invRCurr, k)


        zeta2prop <- computeZetaGamma2(ZCurrprop, component_membership, invRCurr, k)


        zeta <- computeZetaGamma(latent_z, component_membership, invRCurr, k)

        zeta2 <- computeZetaGamma2(latent_z, component_membership, invRCurr, k)

        logLikeProp <- sum(dgamma(Y[,k], rate=gammaAlpha[countgamma]/exp(linear_predictor_proposed), shape= gammaAlpha[countgamma], log=TRUE))+
          0.5*sum(zetaprop) - sum(zeta2prop)
        logLike <- sum(dgamma(Y[,k], rate=gammaAlpha[countgamma]/exp(linear_predictor_current[,k]), shape= gammaAlpha[countgamma], log=TRUE))+
          0.5*sum(zeta) - sum(zeta2)

        V_OLD_prior <- inv_bdiagmat_OLD
        V_OLD_prior <- (V_OLD_prior + t(V_OLD_prior))/2
        V_prior <-inv_bdiagmat
        V_prior <- (V_prior + t(V_prior))/2

        logPrior<- dmvnorm(as.vector(Bprop), rep(0,length(as.vector(Bprop))), V_prior,log=TRUE) - dmvnorm(as.vector(B_kgamma),rep(0,length(as.vector(B_kgamma))), V_OLD_prior,log=TRUE)

        logAccRatio1 <- logPrior - dmvnorm( as.vector(Bprop) , muprop, inv_V,log=TRUE) + dmvnorm( as.vector(B_kgamma) ,muprop_OLD,inv_V_OLD,log=TRUE) +
          logpriorgamma
        logAccRatio <- logAccRatio1 + logLikeProp - logLike
        coefficient_proposal_total[k] <- coefficient_proposal_total[k] + 1

        if( log(runif(1)) < logAccRatio ){
          setone <- which(newGamma>0) + (k-1)*(length(newGamma))
          setzero <- which(newGamma==0) + (k-1)*(length(newGamma))
          gamma_current[setone]<-1
          beta_current[setone]<-Bprop

          gamma_current[setzero] <- 0
          beta_current[setzero]<-0
          WBgamma<-linear_predictor_proposed
          linear_predictor_current[,k]<-linear_predictor_proposed
          acceptratio[k] <- acceptratio[k] + 1
          latent_z <- ZCurrprop
          coefficient_accept_batch[k] <- coefficient_accept_batch[k] + 1
          coefficient_accept_total[k] <- coefficient_accept_total[k] + 1

        }

        if (isTRUE(adapt_coefficient_proposals) && coefficient_adapt_delay > 0L && s == coefficient_adapt_delay) {
          coefficient_accept_batch[k] <- 0
        }
        gamma_batch_size <- if (isTRUE(adapt_coefficient_proposals)) {
          adaptation_batch_size(
            iteration = s,
            burnin = burnin,
            batch_size = coefficient_adapt_batch,
            start_iteration = coefficient_adapt_delay + 1L
          )
        } else {
          0L
        }
        if (gamma_batch_size > 0L) {
          acc_rate <- coefficient_accept_batch[k] / gamma_batch_size
          coefficient_proposal_var[k] <- adapt_positive_tuning(
            current_value = coefficient_proposal_var[k],
            acceptance_rate = acc_rate,
            target_rate = coefficient_target_accept[k],
            lower_tolerance = coefficient_accept_tol_lower[k],
            upper_tolerance = coefficient_accept_tol_upper[k],
            log_step = coefficient_adapt_log_step,
            min_value = coefficient_proposal_var_min,
            max_value = coefficient_proposal_var_max
          )
          var_z[, k] <- coefficient_proposal_var[k]
          coefficient_accept_batch[k] <- 0
        }






        g <- .bnpmccr_update_g_scales(
          g = g,
          gamma_block = gamma_current[start:end],
          beta_block = beta_current[start:end],
          varying_design_matrix = W,
          num_basis = num_basis,
          p = p,
          response_index = k,
          n = n
        )





        gammaind <- which(gamma_current[start:end]>0)
        Wgamma <- totalW[,gammaind]
        B_k <- beta_current[start:end]
        B_kgamma <- B_k[gammaind]
        WBgamma <-  Wgamma %*% as.matrix(B_kgamma)

        additional_param_proposal_total[k] <- additional_param_proposal_total[k] + 1
        gammaalphaProp <- rnorm(1, mean= log(gammaAlpha[countgamma]), sd = additional_param_proposal_sd[k])

        latent_z <- updateZCurrGamma(Y, linear_predictor_current, latent_z, gammaAlpha, countgamma, k)

        ZCurrprop <-  updateZCurrGamma(Y, linear_predictor_current, latent_z, exp(gammaalphaProp), 1, k)



        zetaprop <- computeZetaGamma(ZCurrprop, component_membership, invRCurr, k)

        zeta2prop <- computeZetaGamma2(ZCurrprop, component_membership, invRCurr, k)

        zeta <- computeZetaGamma(latent_z, component_membership, invRCurr, k)

        zeta2 <- computeZetaGamma2(latent_z, component_membership, invRCurr, k)

        logPriorRatio <- log_scale_prior_ratio(
          proposed_log_value = gammaalphaProp,
          current_log_value = log(gammaAlpha[countgamma]),
          sd = log_additional_param_prior_sd
        )

        logLikeProp <- sum(dgamma(Y[,k],rate=exp(gammaalphaProp)/exp(linear_predictor_current[,k]),shape=exp(gammaalphaProp),log = TRUE))+
          0.5*sum(zetaprop) - sum(zeta2prop)
        logLike <- sum(dgamma(Y[,k],rate=gammaAlpha[countgamma]/exp(linear_predictor_current[,k]),shape=gammaAlpha[countgamma],log=TRUE))+
          0.5*sum(zeta) - sum(zeta2)

        logPropRatio <- (gammaalphaProp) - log(gammaAlpha[countgamma])
        logAccRatio <- logPriorRatio + logLikeProp - logLike + logPropRatio

        gamma_additional_param_accepted <- FALSE
        if(log(runif(1)) < logAccRatio){

          gammaAlpha[countgamma]<-exp(gammaalphaProp)
          latent_z <- ZCurrprop
          gamma_additional_param_accepted <- TRUE
        }
        if (gamma_additional_param_accepted) {
          additional_param_accept_count[k] <- additional_param_accept_count[k] + 1
          additional_param_accept_total[k] <- additional_param_accept_total[k] + 1
        }
        if (isTRUE(adapt_additional_params) && additional_param_adapt_delay > 0L && s == additional_param_adapt_delay) {
          additional_param_accept_count[k] <- 0
        }
        current_additional_param_batch <- if (isTRUE(adapt_additional_params)) {
          adaptation_batch_size(
            iteration = s,
            burnin = burnin,
            batch_size = additional_param_adapt_batch,
            start_iteration = additional_param_adapt_delay + 1L
          )
        } else {
          0L
        }
        if (current_additional_param_batch > 0L) {
          observed_accept_rate <- additional_param_accept_count[k] / current_additional_param_batch
          additional_param_proposal_sd[k] <- adapt_rw_proposal_sd(
            current_sd = additional_param_proposal_sd[k],
            acceptance_rate = observed_accept_rate,
            target_rate = additional_param_target_accept[k],
            tolerance = additional_param_accept_tol,
            log_step = additional_param_adapt_log_step,
            min_sd = additional_param_min_proposal_sd,
            max_sd = additional_param_max_proposal_sd
          )
          additional_param_accept_count[k] <- 0
        }


        latent_z <- updateZCurrGamma(Y, linear_predictor_current, latent_z, gammaAlpha, countgamma, k)



      }else{

        if(responseType[k]=="binary"){

          if(length(gammaind)>0){
            setone <- ( start-1) + which(gamma_current[start:end]>0)
            beta_current[setone]<-Bprop
          }
          setzero <- ( start- 1)+ which(gamma_current[start:end]==0)
          beta_current[setzero]<-0
          WBgamma <- linear_predictor_proposed
          linear_predictor_current[,k] <- WBgamma
          var_z[,k] <- rep(1,n)

          g <- .bnpmccr_update_g_scales(
            g = g,
            gamma_block = gamma_current[start:end],
            beta_block = beta_current[start:end],
            varying_design_matrix = W,
            num_basis = num_basis,
            p = p,
            response_index = k,
            n = n
          )


          latent_z[,k] <- latent_linear_predictor[,k]-linear_predictor_current[,k]
          TU[, k] <- cutPoint[[k]][Y[, k] + 2] - WBgamma
          TL[, k] <- cutPoint[[k]][Y[, k] + 1] - WBgamma


        }

        if(responseType[k]=="ordinal"){

          if(length(gammaind)>0){
            setone <-  (start-1 )+ which(gamma_current[start:end]>0)
            beta_current[setone]<-Bprop
          }
          setzero <- (start- 1)+ which(gamma_current[start:end]==0)
          beta_current[setzero]<-0
          WBgamma <- linear_predictor_proposed
          linear_predictor_current[,k] <- WBgamma

          g <- .bnpmccr_update_g_scales(
            g = g,
            gamma_block = gamma_current[start:end],
            beta_block = beta_current[start:end],
            varying_design_matrix = W,
            num_basis = num_basis,
            p = p,
            response_index = k,
            n = n
          )


          latent_z[,k] <- latent_linear_predictor[,k]-linear_predictor_current[,k]
          TU[, k] <- cutPoint[[k]][Y[, k] + 2] - WBgamma
          TL[, k] <- cutPoint[[k]][Y[, k] + 1] - WBgamma

        }
        if(responseType[k]=='binomial'){
          TU[,k]<- qnorm(pbinom(Y[,k],size=nTrial[countbinom],prob = 1/(1+exp(-linear_predictor_current[,k]  ))))
          TUProp<- qnorm(pbinom(Y[,k],size=nTrial[countbinom],prob = 1/(1+exp(-linear_predictor_proposed  ))))
          TL[,k]<- qnorm(pbinom(Y[,k]-1,size=nTrial[countbinom],prob = 1/(1+exp(-linear_predictor_current[,k]  ))))
          TLProp<- qnorm(pbinom(Y[,k]-1,size=nTrial[countbinom],prob = 1/(1+exp(-linear_predictor_proposed    ))))

          current_bounds <- .bnpmccr_stabilize_truncation_bounds(TL[, k], TU[, k])
          TL[, k] <- current_bounds$lower
          TU[, k] <- current_bounds$upper
          proposed_bounds <- .bnpmccr_stabilize_truncation_bounds(TLProp, TUProp)
          TLProp <- proposed_bounds$lower
          TUProp <- proposed_bounds$upper


        }
        if(responseType[k]=='negative binomial'){
          TU[,k]<- qnorm(pnbinom(Y[,k],size=negBinomPar[countnegbinom],mu = exp(linear_predictor_current[,k]) ))
          TUProp<- qnorm(pnbinom(Y[,k],size=negBinomPar[countnegbinom],mu = exp(linear_predictor_proposed)  ))
          TL[,k]<- qnorm(pnbinom(Y[,k]-1,size=negBinomPar[countnegbinom],mu = exp(linear_predictor_current[,k]) ))
          TLProp<- qnorm(pnbinom(Y[,k]-1,size=negBinomPar[countnegbinom],mu = exp(linear_predictor_proposed)  ))

          current_bounds <- .bnpmccr_stabilize_truncation_bounds(TL[, k], TU[, k])
          TL[, k] <- current_bounds$lower
          TU[, k] <- current_bounds$upper
          proposed_bounds <- .bnpmccr_stabilize_truncation_bounds(TLProp, TUProp)
          TLProp <- proposed_bounds$lower
          TUProp <- proposed_bounds$upper

        }

        if( responseType[k]=='binomial' || responseType[k]=='negative binomial' ){

          diffLikeProp <- rep(0,n); diffLike <- rep(0,n)

          diffLike_results <- updateDiffLike(RCurr, component_membership, latent_z,k - 1,  TUProp,TLProp,TU,TL)

          diffLikeProp <- diffLike_results[[1]]
          diffLike <- diffLike_results[[2]]

          diffLikeProp[which(diffLikeProp <= 0)] <- pnorm(Inf) - pnorm(8.2)
          logLikeProp <- sum(log(diffLikeProp))
          diffLike[which(diffLike<= 0)] <- pnorm(Inf) - pnorm(8.2)
          logLike <- sum(log(diffLike))


          V_OLD_prior <- inv_bdiagmat_OLD
          V_OLD_prior <- (V_OLD_prior + t(V_OLD_prior))/2
          V_prior <-inv_bdiagmat
          V_prior <- (V_prior + t(V_prior))/2


          logPrior<- dmvnorm(as.vector(Bprop), rep(0,length(as.vector(Bprop))), V_prior,log=TRUE) - dmvnorm(as.vector(B_kgamma),rep(0,length(as.vector(B_kgamma))), V_OLD_prior,log=TRUE)

          logAccRatio1 <- logPrior - dmvnorm( as.vector(Bprop) , muprop, inv_V,log=TRUE) + dmvnorm( as.vector(B_kgamma) ,muprop_OLD,inv_V_OLD,log=TRUE) +
            logpriorgamma
          logAccRatio <- logAccRatio1 + logLikeProp - logLike
          coefficient_proposal_total[k] <- coefficient_proposal_total[k] + 1

          if( log(runif(1)) < logAccRatio ){
            setone <- which(newGamma>0) + (k-1)*(length(newGamma))
            setzero <- which(newGamma==0) + (k-1)*(length(newGamma))
            gamma_current[setone]<-1
            beta_current[setone]<-Bprop

            gamma_current[setzero] <- 0
            beta_current[setzero]<-0
            WBgamma<-linear_predictor_proposed
            linear_predictor_current[,k]<-linear_predictor_proposed
            TU[,k] <-TUProp
            TL[,k] <- TLProp
            acceptratio[k] <- acceptratio[k] + 1
            coefficient_accept_batch[k] <- coefficient_accept_batch[k] + 1
            coefficient_accept_total[k] <- coefficient_accept_total[k] + 1
          }

          if (isTRUE(adapt_coefficient_proposals) && coefficient_adapt_delay > 0L && s == coefficient_adapt_delay) {
            coefficient_accept_batch[k] <- 0
            if (responseType[k] == "negative binomial" && !is.null(negbin_coefficient_proposal_var_reset)) {
              coefficient_proposal_var[k] <- negbin_coefficient_proposal_var_reset
              var_z[, k] <- coefficient_proposal_var[k]
            }
          }
          coefficient_batch_size <- if (isTRUE(adapt_coefficient_proposals)) {
            adaptation_batch_size(
              iteration = s,
              burnin = burnin,
              batch_size = coefficient_adapt_batch,
              start_iteration = coefficient_adapt_delay + 1L
            )
          } else {
            0L
          }
          if (coefficient_batch_size > 0L) {
            acc_rate <- coefficient_accept_batch[k] / coefficient_batch_size
            coefficient_proposal_var[k] <- adapt_positive_tuning(
              current_value = coefficient_proposal_var[k],
              acceptance_rate = acc_rate,
              target_rate = coefficient_target_accept[k],
              lower_tolerance = coefficient_accept_tol_lower[k],
              upper_tolerance = coefficient_accept_tol_upper[k],
              log_step = coefficient_adapt_log_step,
              min_value = coefficient_proposal_var_min,
              max_value = coefficient_proposal_var_max
            )
            var_z[, k] <- coefficient_proposal_var[k]
            coefficient_accept_batch[k] <- 0
          }

          g <- .bnpmccr_update_g_scales(
            g = g,
            gamma_block = gamma_current[start:end],
            beta_block = beta_current[start:end],
            varying_design_matrix = W,
            num_basis = num_basis,
            p = p,
            response_index = k,
            n = n
          )

        }

        if(responseType[k]=="ordinal"){
          gammaind <- which(gamma_current[start:end] > 0)

          if (length(gammaind) > 0L) {
            Wgamma <- totalW[, gammaind, drop = FALSE]
            B_k <- beta_current[start:end]
            B_kgamma <- B_k[gammaind]
            WBgamma <- Wgamma %*% as.matrix(B_kgamma)
          } else {
            WBgamma <- rep(0, n)
          }

          var_z[, k] <- rep(1, n)
          cutPoint_TMP <- rep(NA_real_, length(cutPoint[[k]]) - 2L)
          cutPoint_TMP[1] <- 0
          cutPoint_TMP[2:length(cutPoint_TMP)] <- log(cutPoint[[k]][3:(length(cutPoint[[k]]) - 1L)] - cutPoint[[k]][2:(length(cutPoint[[k]]) - 2L)])
          ordinal_cutpoint_proposal_total[k] <- ordinal_cutpoint_proposal_total[k] + 1L
          cutpoint_proptmp <- cutPoint_TMP[2:length(cutPoint_TMP)] + stats::rnorm(
            length(cutPoint_TMP) - 1L,
            mean = 0,
            sd = ordinal_cutpoint_proposal_sd[k]
          )
          cutPoint_prop <- cutPoint[[k]]
          cutPoint_prop[3] <- exp(cutpoint_proptmp)[1]

          if (nCat[k] >= 4L) {
            for (j in 4:nCat[k]) {
              cutPoint_prop[j] <- sum(exp(cutpoint_proptmp)[1:(j - 2L)])
            }
          }

          TUProp <- cutPoint_prop[Y[, k] + 2L] - WBgamma
          TU[, k] <- cutPoint[[k]][Y[, k] + 2L] - WBgamma
          TLProp <- cutPoint_prop[Y[, k] + 1L] - WBgamma
          TL[, k] <- cutPoint[[k]][Y[, k] + 1L] - WBgamma

          logpriorRatio <- normal_prior_ratio(
            proposed_value = cutpoint_proptmp,
            current_value = cutPoint_TMP[2:length(cutPoint_TMP)],
            sd = ordinal_threshold_prior_sd
          )
          diffLike_results <- updateDiffLike(RCurr, component_membership, latent_z, k - 1L, TUProp, TLProp, TU, TL)
          diffLikeProp <- diffLike_results[[1]]
          diffLike <- diffLike_results[[2]]

          diffLikeProp[which(diffLikeProp <= 0)] <- pnorm(Inf) - pnorm(8.2)
          logLikeProp <- sum(log(diffLikeProp))

          diffLike[which(diffLike <= 0)] <- pnorm(Inf) - pnorm(8.2)
          logLike <- sum(log(diffLike))

          logAccRatio <- logpriorRatio + logLikeProp - logLike

          if (log(runif(1)) < logAccRatio) {
            cutPoint[[k]][3] <- cutPoint_prop[3]
            if (nCat[k] >= 4L) {
              for (j in 4:nCat[k]) {
                cutPoint[[k]][j] <- cutPoint_prop[j]
              }
            }
            TU[, k] <- TUProp
            TL[, k] <- TLProp
            ordinal_cutpoint_accept_batch[k] <- ordinal_cutpoint_accept_batch[k] + 1L
            ordinal_cutpoint_accept_total[k] <- ordinal_cutpoint_accept_total[k] + 1L
          }

          if (isTRUE(adapt_ordinal_cutpoint_proposals) && ordinal_cutpoint_adapt_delay > 0L && s == ordinal_cutpoint_adapt_delay) {
            ordinal_cutpoint_accept_batch[k] <- 0L
          }
          ordinal_cutpoint_batch_size <- if (isTRUE(adapt_ordinal_cutpoint_proposals)) {
            adaptation_batch_size(
              iteration = s,
              burnin = burnin,
              batch_size = ordinal_cutpoint_adapt_batch,
              start_iteration = ordinal_cutpoint_adapt_delay + 1L
            )
          } else {
            0L
          }
          if (ordinal_cutpoint_batch_size > 0L) {
            acc_rate <- ordinal_cutpoint_accept_batch[k] / ordinal_cutpoint_batch_size
            ordinal_cutpoint_proposal_sd[k] <- adapt_rw_proposal_sd(
              current_sd = ordinal_cutpoint_proposal_sd[k],
              acceptance_rate = acc_rate,
              target_rate = ordinal_cutpoint_target_accept,
              tolerance = ordinal_cutpoint_accept_tol,
              log_step = ordinal_cutpoint_adapt_log_step,
              min_sd = ordinal_cutpoint_min_proposal_sd,
              max_sd = ordinal_cutpoint_max_proposal_sd
            )
            ordinal_cutpoint_accept_batch[k] <- 0L
          }
        }


        if(responseType[k]=='negative binomial'){
          gammaind <- which(gamma_current[start:end]>0)
          Wgamma <- totalW[,gammaind]
          B_k <- beta_current[start:end]
          B_kgamma <- B_k[gammaind]
          WBgamma <-  Wgamma %*% as.matrix(B_kgamma)


          additional_param_proposal_total[k] <- additional_param_proposal_total[k] + 1
          rhoProp <- rnorm(1, mean= log(negBinomPar[countnegbinom]), sd = additional_param_proposal_sd[k])
          TUProp <- qnorm(pnbinom(Y[,k],size=exp(rhoProp),mu = exp(WBgamma) ))
          TLProp <- qnorm(pnbinom(Y[,k]-1,size=exp(rhoProp),mu = exp(WBgamma) ))

          proposed_bounds <- .bnpmccr_stabilize_truncation_bounds(TLProp, TUProp)
          TLProp <- proposed_bounds$lower
          TUProp <- proposed_bounds$upper

          logPriorRatio <- log_scale_prior_ratio(
            proposed_log_value = rhoProp,
            current_log_value = log(negBinomPar[countnegbinom]),
            sd = log_additional_param_prior_sd
          )

          diffLikeProp <- rep(0,n); diffLike <- rep(0,n)

          diffLike_results <- updateDiffLike(RCurr, component_membership,latent_z,k - 1,  TUProp,TLProp,TU,TL)

          diffLikeProp <- diffLike_results[[1]]
          diffLike <- diffLike_results[[2]]

          diffLikeProp[which(diffLikeProp==0)] <- 10^(-16)
          logLikeProp <- sum(log(diffLikeProp))


          diffLike[which(diffLike==0)] <- 10^(-16)
          logLike <- sum(log(diffLike))

          logPropRatio <- rhoProp - log(negBinomPar[countnegbinom])
          logAccRatio <- logPriorRatio + logLikeProp - logLike + logPropRatio

          negbin_additional_param_accepted <- FALSE
          if(log(runif(1)) < logAccRatio){
            TU[,k]<-TUProp
            TL[,k] <-TLProp
            negBinomPar[countnegbinom]<-exp(rhoProp)
            negbin_additional_param_accepted <- TRUE
          }
          if (negbin_additional_param_accepted) {
            additional_param_accept_count[k] <- additional_param_accept_count[k] + 1
            additional_param_accept_total[k] <- additional_param_accept_total[k] + 1
          }
          if (isTRUE(adapt_additional_params) &&
              additional_param_adapt_delay > 0L &&
              s == additional_param_adapt_delay) {
            additional_param_accept_count[k] <- 0L
            additional_param_accept_total[k] <- 0L
            additional_param_proposal_total[k] <- 0L
          }
          current_additional_param_batch <- if (isTRUE(adapt_additional_params)) {
            adaptation_batch_size(
              iteration = s,
              burnin = burnin,
              batch_size = additional_param_adapt_batch,
              start_iteration = additional_param_adapt_delay + 1L
            )
          } else {
            0L
          }
          if (current_additional_param_batch > 0L) {
            observed_accept_rate <- additional_param_accept_count[k] / current_additional_param_batch
            additional_param_proposal_sd[k] <- adapt_rw_proposal_sd(
              current_sd = additional_param_proposal_sd[k],
              acceptance_rate = observed_accept_rate,
              target_rate = additional_param_target_accept[k],
              tolerance = additional_param_accept_tol,
              log_step = additional_param_adapt_log_step,
              min_sd = additional_param_min_proposal_sd,
              max_sd = additional_param_max_proposal_sd
            )
            additional_param_accept_count[k] <- 0
          }


        }


        latent_z <- updateZCurr( RCurr, component_membership, latent_z,  k - 1,  TL,  TU )

      }
    }





    alpha_mat <- matrix(alpha_h, nrow = n, ncol = H, byrow = TRUE)
    P <- pnorm(alpha_mat + psbp_linear_predictor)
    pi_prob <- matrix(NA, n, H)

    remainLen <- rep(1, n)

    for(h in seq_len(H - 1)) {
      pi_prob[, h] <- remainLen * P[, h]
      remainLen <- remainLen - pi_prob[, h]
    }

    pi_prob[, H] <- remainLen

    copula_prob <- matrix(NA, n, H)

    epsilon <- 1e-12
    pi_prob <- pmax(pi_prob, epsilon)

    for(h in seq_len(H)) {
      log_det_Rh <- as.numeric(determinant(RCurr[[h]], logarithm = TRUE)$modulus)
      const_term <- -0.5 * log_det_Rh

      M_h <- invRCurr[[h]] - diag(1, m)

      quadform_vec <- rowSums((latent_z %*% M_h) * latent_z)

      copula_prob[, h] <- const_term + (-0.5)*quadform_vec + log(pi_prob[, h])
    }

    copula_prob <- .bnpmccr_softmax_rows(copula_prob)

    component_membership <- apply(
      copula_prob, 1,
      function(prob_j) sample.int(H, size = 1, prob = prob_j)
    )


    for(h in 1:H){
      num_pi[h] <- sum(component_membership==h)
    }

    z_star <- matrix(NA_real_, n, max(H - 1L, 0L))
    z_star <- update_Z_star(component_membership, psbp_linear_predictor, z_star, alpha_h, H)

    updateDP_results <- updateDPParameters( H, p2, component_membership,
                                            psbp_gamma_current, psbp_basis_matrix, z_star,
                                            gh, omegaC, num_pi2,
                                            alpha_h, psbp_beta_current, psbp_linear_predictor )

    psbp_gamma_current <- c(updateDP_results[[1]])
    gh <- updateDP_results[[2]]
    alpha_h <- c(updateDP_results[[3]])
    psbp_beta_current <- c(updateDP_results[[4]])
    psbp_linear_predictor <- updateDP_results[[5]]
    num_pi2 <- c(updateDP_results[[6]])


    respind <- which(responseType == "gamma")

    for(h in 1:H){
      if(num_pi[h] <= 0 ){

        Wmat[[h]] <- 0
        invSigma[[h]] <- rWishart(1,m+1, Sigma = diag(1,m))
        Sigma[[h]] <- solve(matrix(invSigma[[h]],m,m)+diag(rep(10 ^(-10), m)))
        RCurr[[h]] <- diag(1/sqrt(diag(Sigma[[h]]))) %*% Sigma[[h]] %*% diag(1/sqrt(diag(Sigma[[h]])))
        invRCurr[[h]] <- solve(RCurr[[h]])
        D[[h]] <- diag(sqrt(diag(Sigma[[h]])))

      }else{
        diag(D[[h]]) <- sqrt(1/rgamma(m,(m+1)/2, scale= 2/diag(invRCurr[[h]]) ))

        Wmat[[h]] <- latent_z[which(component_membership==h),] %*% D[[h]]
        if (length(respind) > 0L) {
          dpropor[which(component_membership == h), respind] <- t(
            matrix(rep(diag(D[[h]])[respind], num_pi[h]), length(respind), num_pi[h])
          )
        }
        num <- num_pi[h]+ m + 1

        SigmaWmat <- solve(t(Wmat[[h]])%*% Wmat[[h]]+diag(1, m),tol=1e-30)

        invSigma[[h]] <- rWishart(1,num, Sigma = SigmaWmat)
        Sigma[[h]] <- solve(matrix(invSigma[[h]],m,m)+diag(rep(10 ^(-10), m)),tol=1e-30)

        RCurr[[h]] <- diag(1/sqrt(diag(Sigma[[h]]))) %*% Sigma[[h]] %*% diag(1/sqrt(diag(Sigma[[h]])))
        invRCurr[[h]] <- solve(RCurr[[h]])
        D[[h]] <- diag(sqrt(diag(Sigma[[h]])))

        if (length(respind) > 0L) {
          dpropor[which(component_membership == h), respind] <- dpropor[which(component_membership == h), respind] * t(
            matrix(rep(1 / diag(D[[h]])[respind], num_pi[h]), length(respind), num_pi[h])
          )
        }
        WD <- Wmat[[h]] %*% diag(1 / diag(D[[h]]))
        latent_z[which(component_membership==h),] <- WD
      }
    }



    Rsave[[s]] <- matrix(unlist(RCurr), m*m, H)

    posterior_samples$beta[s,]<- beta_current
    posterior_samples$beta_x[s,] <- psbp_beta_current
    posterior_samples$gamma[s,]<- gamma_current
    posterior_samples$g[s,] <- g
    cutpointsave[[s]] <- cutPoint
    negparsave[[s]] <- negBinomPar
    gammaparsave[[s]] <- gammaAlpha
    posterior_samples$numPI[s,]<-num_pi
    posterior_samples$d[s,] <- d
    posterior_samples$g_h[s,] <- gh
    posterior_samples$neg[s,] <- unlist(negBinomPar)
    posterior_samples$alpha_h[s,] <- alpha_h



    if (isTRUE(verbose) && (s %% progress_every == 0L)) {
      message(sprintf("fit_bnpmccr: iteration %d / %d completed.", s, S))
    }
  }

  margin_labels <- paste0("margin", seq_len(m), "_", responseType)
  additional_param_acceptance_rate <- rep(NA_real_, m)
  additional_param_updated <- additional_param_proposal_total > 0
  additional_param_acceptance_rate[additional_param_updated] <-
    additional_param_accept_total[additional_param_updated] /
    additional_param_proposal_total[additional_param_updated]
  coefficient_acceptance_rate <- rep(NA_real_, m)
  coefficient_updated <- coefficient_proposal_total > 0
  coefficient_acceptance_rate[coefficient_updated] <-
    coefficient_accept_total[coefficient_updated] /
    coefficient_proposal_total[coefficient_updated]
  ordinal_cutpoint_acceptance_rate <- rep(NA_real_, m)
  ordinal_cutpoint_updated <- ordinal_cutpoint_proposal_total > 0
  ordinal_cutpoint_acceptance_rate[ordinal_cutpoint_updated] <-
    ordinal_cutpoint_accept_total[ordinal_cutpoint_updated] /
    ordinal_cutpoint_proposal_total[ordinal_cutpoint_updated]
  names(additional_param_proposal_sd) <- margin_labels
  names(additional_param_accept_total) <- margin_labels
  names(additional_param_proposal_total) <- margin_labels
  names(additional_param_acceptance_rate) <- margin_labels
  names(ordinal_cutpoint_proposal_sd) <- margin_labels
  names(ordinal_cutpoint_accept_total) <- margin_labels
  names(ordinal_cutpoint_proposal_total) <- margin_labels
  names(ordinal_cutpoint_acceptance_rate) <- margin_labels
  names(coefficient_proposal_var) <- margin_labels
  names(coefficient_accept_total) <- margin_labels
  names(coefficient_proposal_total) <- margin_labels
  names(coefficient_acceptance_rate) <- margin_labels

  output <- list(
    B = posterior_samples$beta,
    G = posterior_samples$gamma,
    cutPoint = cutpointsave,
    negPar = negparsave,
    W = W,
    g = posterior_samples$g,
    numPi = posterior_samples$numPI,
    R = RCurr,
    D = D,
    d = posterior_samples$d,
    gh = posterior_samples$g_h,
    alpha_h = posterior_samples$alpha_h,
    beta2 = posterior_samples$beta_x,
    neg = posterior_samples$neg,
    Rmat = Rsave,
    gammaPar = gammaparsave,
    var_z = var_z,
    dpropormat = dpropormat,
    acceptratio = acceptratio,
    covT = X,
    predX = original_predictor_matrix,
    Y = original_response_matrix,
    burnin = burnin,
    niter = niter,
    responseType = responseType,
    prior_settings = list(
      gaussian_variance_prior_eps = gaussian_variance_prior_eps,
      log_additional_param_prior_sd = log_additional_param_prior_sd,
      ordinal_threshold_prior_sd = ordinal_threshold_prior_sd
    ),
    proposal_settings = list(
      additional_parameter_families = c("gamma", "negative binomial"),
      binomial_has_additional_parameter = FALSE,
      adapt_additional_params = adapt_additional_params,
      additional_param_adapt_batch = additional_param_adapt_batch,
      additional_param_adapt_delay = additional_param_adapt_delay,
      gamma_additional_proposal_sd_init = gamma_additional_proposal_sd_init,
      negbin_additional_proposal_sd_init = negbin_additional_proposal_sd_init,
      gamma_additional_target_accept = gamma_additional_target_accept,
      negbin_additional_target_accept = negbin_additional_target_accept,
      additional_param_accept_tol = additional_param_accept_tol,
      additional_param_adapt_log_step = additional_param_adapt_log_step,
      additional_param_min_proposal_sd = additional_param_min_proposal_sd,
      additional_param_max_proposal_sd = additional_param_max_proposal_sd,
      ordinal_parameter_family = "ordinal cutpoints",
      adapt_ordinal_cutpoint_proposals = adapt_ordinal_cutpoint_proposals,
      ordinal_cutpoint_adapt_batch = ordinal_cutpoint_adapt_batch,
      ordinal_cutpoint_adapt_delay = ordinal_cutpoint_adapt_delay,
      ordinal_cutpoint_proposal_sd_init = ordinal_cutpoint_proposal_sd_init,
      ordinal_cutpoint_target_accept = ordinal_cutpoint_target_accept,
      ordinal_cutpoint_accept_tol = ordinal_cutpoint_accept_tol,
      ordinal_cutpoint_adapt_log_step = ordinal_cutpoint_adapt_log_step,
      ordinal_cutpoint_min_proposal_sd = ordinal_cutpoint_min_proposal_sd,
      ordinal_cutpoint_max_proposal_sd = ordinal_cutpoint_max_proposal_sd,
      coefficient_parameter_families = c("gamma", "binomial", "negative binomial"),
      adapt_coefficient_proposals = adapt_coefficient_proposals,
      coefficient_adapt_batch = coefficient_adapt_batch,
      coefficient_adapt_delay = coefficient_adapt_delay,
      coefficient_adapt_log_step = coefficient_adapt_log_step,
      coefficient_proposal_var_min = coefficient_proposal_var_min,
      coefficient_proposal_var_max = coefficient_proposal_var_max,
      gamma_coefficient_proposal_var_init = gamma_coefficient_proposal_var_init,
      binomial_coefficient_proposal_var_init = binomial_coefficient_proposal_var_init,
      negbin_coefficient_proposal_var_init = negbin_coefficient_proposal_var_init,
      gamma_coefficient_target_accept = gamma_coefficient_target_accept,
      gamma_coefficient_accept_tol = gamma_coefficient_accept_tol,
      binomial_coefficient_target_accept = binomial_coefficient_target_accept,
      binomial_coefficient_accept_tol = binomial_coefficient_accept_tol,
      negbin_coefficient_target_accept = negbin_coefficient_target_accept,
      negbin_coefficient_accept_tol = negbin_coefficient_accept_tol,
      negbin_coefficient_proposal_var_reset = negbin_coefficient_proposal_var_reset
    ),
    additional_param_proposal_sd = additional_param_proposal_sd,
    additional_param_acceptance_rate = additional_param_acceptance_rate,
    additional_param_accept_total = additional_param_accept_total,
    additional_param_proposal_total = additional_param_proposal_total,
    ordinal_cutpoint_proposal_sd = ordinal_cutpoint_proposal_sd,
    ordinal_cutpoint_acceptance_rate = ordinal_cutpoint_acceptance_rate,
    ordinal_cutpoint_accept_total = ordinal_cutpoint_accept_total,
    ordinal_cutpoint_proposal_total = ordinal_cutpoint_proposal_total,
    coefficient_proposal_var = coefficient_proposal_var,
    coefficient_acceptance_rate = coefficient_acceptance_rate,
    coefficient_accept_total = coefficient_accept_total,
    coefficient_proposal_total = coefficient_proposal_total
  )

  if (!keep_burnin && burnin > 0L) {
    keep_idx <- seq.int(burnin + 1L, niter)
    output$B <- output$B[keep_idx, , drop = FALSE]
    output$G <- output$G[keep_idx, , drop = FALSE]
    output$g <- output$g[keep_idx, , drop = FALSE]
    output$numPi <- output$numPi[keep_idx, , drop = FALSE]
    output$d <- output$d[keep_idx, , drop = FALSE]
    output$gh <- output$gh[keep_idx, , drop = FALSE]
    output$alpha_h <- output$alpha_h[keep_idx, , drop = FALSE]
    if (length(output$neg) > 0L) {
      output$neg <- output$neg[keep_idx, , drop = FALSE]
    }
    output$beta2 <- output$beta2[keep_idx, , drop = FALSE]
    output$cutPoint <- output$cutPoint[keep_idx]
    output$negPar <- output$negPar[keep_idx]
    output$gammaPar <- output$gammaPar[keep_idx]
    output$Rmat <- output$Rmat[keep_idx]
  }

  class(output) <- "bnpmccr_fit"
  return(output)



}
