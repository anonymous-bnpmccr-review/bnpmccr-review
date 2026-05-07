// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include "updateZCurrGamma2.h"

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
arma::mat updateZCurrGamma2(
    const arma::mat& Y,
    const arma::mat& WB,
    arma::mat&        ZCurr,
    const arma::vec& gammaAlpha,
    int              countgamma,
    int              k
) {
  int n = Y.n_rows;
  int col = k - 1;                // 1-based → 0-based
  int gamma_idx = countgamma - 1; // 1-based → 0-based
  double shape = gammaAlpha(gamma_idx);

  for (int i = 0; i < n; ++i) {
    double x      = Y(i, col);
    double wb_val = WB(i);
    double scale  = std::exp(wb_val) / shape;

    double cdf = R::pgamma(x, shape, scale, /*lower_tail=*/1, /*log_p=*/0);
        cdf = std::max(1e-12, std::min(cdf, 0.99999999999999));
    double z = R::qnorm(cdf, 0.0, 1.0, /*lower_tail=*/1, /*log_p=*/0);
    ZCurr(i, col) = z;
  }

  return ZCurr;
}
