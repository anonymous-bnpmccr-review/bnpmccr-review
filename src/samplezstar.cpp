// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <limits>
#include <cmath>
using namespace Rcpp;

//--------------------------------------------------------------------------
// rtruncnorm
//--------------------------------------------------------------------------
static inline double rtruncnorm(double a, double b, double mean, double sd) {
  // Z = (X - mean) / sd
  double alpha = (a - mean) / sd;
  double beta  = (b - mean) / sd;

  double Phi_alpha, Phi_beta;
  if (std::isinf(alpha) && alpha < 0) {
    Phi_alpha = 0.0;
  } else {
    Phi_alpha = R::pnorm(alpha, 0.0, 1.0, 1, 0);
  }
  if (std::isinf(beta) && beta > 0) {
    Phi_beta = 1.0;
  } else {
    Phi_beta = R::pnorm(beta, 0.0, 1.0, 1, 0);
  }

  // Uniform(0,1)
  double u = R::runif(0.0, 1.0);
  // p = Phi_alpha + u*(Phi_beta-Phi_alpha)
  double p = Phi_alpha + u * (Phi_beta - Phi_alpha);
  // original scale
  double z = R::qnorm(p, 0.0, 1.0, 1, 0);
  return mean + sd * z;
}

//--------------------------------------------------------------------------
// update_Z_star: Z_star
//--------------------------------------------------------------------------

// [[Rcpp::export]]
NumericMatrix update_Z_star(IntegerVector Z_ind,
                            NumericMatrix WB_x,
                            NumericMatrix Z_star,
                            NumericVector alpha_h,
                            int H) {
  int n = WB_x.nrow();

  for (int j = 0; j < n; j++) {
    int i = Z_ind[j] - 1;

    if (i >= 1) {
      for (int h = 0; h < i; h++) {
        double mean_val = WB_x(j, h) + alpha_h[h];
        Z_star(j, h) = rtruncnorm(-std::numeric_limits<double>::infinity(),
               0.0,
               mean_val,
               1.0);
      }
    }

    if (i != (H - 1)) {
      double mean_val = WB_x(j, i) + alpha_h[i];
      Z_star(j, i) = rtruncnorm(0.0,
             std::numeric_limits<double>::infinity(),
             mean_val,
             1.0);
    }
  }

  return Z_star;
}
