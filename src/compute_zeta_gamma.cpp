#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector computeZetaGamma(const NumericMatrix& ZCurr,
                          const IntegerVector& Z_ind,
                          const List& invRCurr,
                          int k) {

  int n = ZCurr.nrow();

  // Initialize
  NumericVector zeta(n, 0.0);

  // for(j in 1:n)
    // Rcpp: 0-based index
  for (int j = 0; j < n; j++) {

    NumericMatrix invR = invRCurr[ Z_ind[j] - 1 ];

    double val = std::pow(ZCurr(j, k - 1), 2.0) * (1.0 - invR(k - 1, k - 1));
    zeta[j] += val;
  }

  return zeta;
}
