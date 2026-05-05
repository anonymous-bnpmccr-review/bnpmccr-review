# Helper functions for spline basis construction and simple prior calculations

log_prior <- function(omega, M, num_knots) {
  num_knots * log(1 - omega) + log(omega) - log(choose(M, num_knots))
}

select_knots <- function(X, M) {
  val <- unique(X)
  quantile(val, probs = seq(0, 1, length.out = M + 2))
}

NaturalCubicBasis <- function(X, knots) {
  # boundary knots
  tL <- knots[1]
  tU <- knots[length(knots)]

  # interior knots
  intKnots <- knots[-c(1, length(knots))]

  X <- as.matrix(X)
  W <- matrix(NA_real_, nrow(X), length(intKnots) + 1)
  W[, 1] <- X

  for (j in seq_along(intKnots)) {
    W[, j + 1] <-
      ((X - intKnots[j])^3 * ((X - intKnots[j]) >= 0) -
         (X - tU)^3 * ((X - tU) >= 0)) / (tU - intKnots[j]) -
      ((X - tL)^3 * ((X - tL) >= 0) -
         (X - tU)^3 * ((X - tU) >= 0)) / (tU - tL)
  }

  W
}
