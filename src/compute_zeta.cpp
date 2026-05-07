#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector compute_zeta(NumericMatrix ZCurr, List invRCurr, IntegerVector Z_ind, int k) {
  int n = ZCurr.nrow();
  int m = ZCurr.ncol();
  NumericVector zeta(n, 0.0);

  //
  for (int kk = 0; kk < m; kk++) {
    if (kk == k - 1) continue;
    for (int j = 0; j < n; j++) {
      //
      int listIndex = Z_ind[j] - 1;
      NumericMatrix invR = invRCurr[listIndex];
      zeta[j] += ZCurr(j, kk) * invR(k - 1, kk);
    }
  }
  return zeta;
}
