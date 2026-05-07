#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector computeZetaGamma2(const NumericMatrix& ZCurr,
                           const IntegerVector& Z_ind,
                           const List& invRCurr,
                           int k) {
  // n: # of rows
  int n = ZCurr.nrow();
  // m: # of cols
  int m = ZCurr.ncol();

  NumericVector zeta2(n, 0.0);

  for (int kk = 0; kk < m; kk++) {
    if (kk == (k - 1)) {
      continue;
    }

      for (int j = 0; j < n; j++) {
        NumericMatrix invR = invRCurr[ Z_ind[j] - 1 ];
        zeta2[j] += ZCurr(j, k - 1) * ZCurr(j, kk) * invR(k - 1, kk);
      }
  }

  return zeta2;
}
