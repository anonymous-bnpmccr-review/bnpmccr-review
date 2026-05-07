// updateZCurr.cpp
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]
#include "updateZCurr.h"
#include <cmath>

using namespace Rcpp;
using namespace arma;

//-----------------------------------------------------------------
// Truncated normal sampler on [a, b].
static inline double rtruncnorm(double a, double b, double mean, double sd) {
  double alpha = (a - mean) / sd;
  double beta  = (b - mean) / sd;

  double Phi_alpha, Phi_beta;
  if(std::isinf(alpha) && alpha < 0)
    Phi_alpha = 0.0;
  else
    Phi_alpha = R::pnorm(alpha, 0.0, 1.0, 1, 0);

  if(std::isinf(beta) && beta > 0)
    Phi_beta = 1.0;
  else
    Phi_beta = R::pnorm(beta, 0.0, 1.0, 1, 0);

  double u = R::runif(0.0, 1.0);
  double p = Phi_alpha + u * (Phi_beta - Phi_alpha);

  // Clamp p to avoid numerical issues at the boundaries.
  double eps = 1e-12;
  if (p < eps) p = eps;
  if (p > 1.0 - eps) p = 1.0 - eps;

  double z = R::qnorm(p, 0.0, 1.0, 1, 0);
  return mean + sd * z;
}

//-----------------------------------------------------------------
// Update the kth latent variable for each observation using truncated normal sampling.
// k must be passed as a 0-indexed column index.
// [[Rcpp::export]]
arma::mat updateZCurr_binom(const Rcpp::List& RCurr,
                            const arma::ivec& Z_ind,
                            arma::mat ZCurr,
                            const int k,        // 0-indexed column index
                            const arma::mat &TL, // lower-bound matrix
                            const arma::mat &TU  // upper-bound matrix
) {
  int n = ZCurr.n_rows;

  for (int i = 0; i < n; ++i) {
    // R uses 1-based indexing, so subtract 1 for C++.
    int listIndex = Z_ind[i] - 1;
    arma::mat RMat = Rcpp::as<arma::mat>(RCurr[listIndex]);

    arma::rowvec row_k = RMat.row(k);
    arma::rowvec row_k_excl = row_k;
    row_k_excl.shed_col(k);

    arma::mat RMat_sub = RMat;
    RMat_sub.shed_row(k);
    RMat_sub.shed_col(k);

    arma::rowvec RR = row_k_excl * inv(RMat_sub);

    arma::rowvec Zrow = ZCurr.row(i);
    arma::rowvec Zrow_excl = Zrow;
    Zrow_excl.shed_col(k);
    double mean_z = as_scalar(RR * Zrow_excl.t());

    arma::vec col_k = RMat.col(k);
    col_k.shed_row(k);
    double sd_z = std::sqrt(1 - as_scalar(RR * col_k));

    double newVal = rtruncnorm(TL(i,0), TU(i,0), mean_z, sd_z);
    ZCurr(i, k) = newVal;
  }

  return ZCurr;
}
