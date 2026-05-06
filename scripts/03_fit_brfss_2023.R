# -----------------------------------------------------------------------------
# 03_fit_brfss_2023.R
#
# Fit the BNPMCCR model to the preprocessed BRFSS 2023 real-data application.
# Run scripts/05_download_brfss_2023.R and scripts/06_prepare_brfss_2023.R first.
# -----------------------------------------------------------------------------

# Load package ----------------------------------------------------------------

if (!requireNamespace("bnpmccr", quietly = TRUE)) {
  stop(
    "The 'bnpmccr' package is not installed. ",
    "Install it before running this script, or use devtools::load_all('.') ",
    "from the package root during development."
  )
}

library(bnpmccr)

# Input / output paths ---------------------------------------------------------

processed_rds_file <- file.path(
  "data", "processed", "brfss2023_preprocessed_for_bnpmccr.rds"
)

fit_dir <- file.path("results", "brfss_2023")
dir.create(fit_dir, recursive = TRUE, showWarnings = FALSE)

fit_rds_file <- file.path(fit_dir, "bnpmccr_brfss2023_fit.rds")
settings_rds_file <- file.path(fit_dir, "bnpmccr_brfss2023_fit_settings.rds")
session_info_file <- file.path(fit_dir, "sessionInfo_brfss2023_fit.txt")

if (!file.exists(processed_rds_file)) {
  stop(
    "Preprocessed BRFSS file was not found: ", processed_rds_file, "\n",
    "Run scripts/06_prepare_brfss_2023.R first."
  )
}

# Load preprocessed data -------------------------------------------------------

brfss_data <- readRDS(processed_rds_file)

Y <- brfss_data$Y
covT <- brfss_data$X
predX <- brfss_data$X3
responseType <- brfss_data$responseType

stopifnot(
  nrow(Y) == length(covT),
  nrow(Y) == nrow(predX),
  ncol(Y) == length(responseType)
)

# MCMC and model settings ------------------------------------------------------
# These settings match the BRFSS 2023 real-data fitting setup, translated to
# the new fit_bnpmccr() wrapper.

set.seed(2025)

niter <- 5000L
burnin <- 0L
H <- 10L
num_basis <- 12L       # 10 candidate knots + constant/linear basis terms
num_basis2 <- 12L      # 10 candidate knots + constant/linear basis terms
omega <- 0.4

fit_settings <- list(
  niter = niter,
  burnin = burnin,
  initH = H,
  num_basis = num_basis,
  num_basis2 = num_basis2,
  omegaM = omega,
  omegaC = omega,
  responseType = responseType,
  seed = 2025L,
  n = nrow(Y),
  m = ncol(Y),
  p = ncol(predX)
)

saveRDS(fit_settings, settings_rds_file)

cat("Fitting BNPMCCR to BRFSS 2023 data...\n")
cat("n =", nrow(Y), "\n")
cat("m =", ncol(Y), "\n")
cat("p =", ncol(predX), "\n")

# Fit model -------------------------------------------------------------------

output <- fit_bnpmccr(
  Y = Y,
  covT = covT,
  predX = predX,
  responseType = responseType,
  niter = niter,
  burnin = burnin,
  num_basis = num_basis,
  tau = 10,
  std = FALSE,
  initH = H,
  nTrial = c(10L),
  negBinomParInit = 1,
  AlphaInit = 1,
  num_basis2 = num_basis2,
  tau2 = 1,
  sdbeta1 = sqrt(0.1),
  initAlpha = 1e10,
  omegaM = omega,
  omegaC = omega,
  mean_alpha = 0,
  verbose = TRUE,
  progress_every = 100L,
  additional_param_max_proposal_sd = 10,
  gamma_additional_target_accept = 0.35,
  gamma_additional_proposal_sd_init = 0.2,
  additional_param_adapt_log_step = 0.03,
  additional_param_adapt_delay = 500L
)

# Save outputs ----------------------------------------------------------------

saveRDS(output, fit_rds_file)

sink(session_info_file)
print(sessionInfo())
sink()

cat("Saved fit object:", fit_rds_file, "\n")
cat("Saved settings:", settings_rds_file, "\n")
cat("Saved session info:", session_info_file, "\n")
