#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
List computeZetaAndZeta2(const NumericMatrix& ZCurr,
                         const List& invRCurr,
                         const IntegerVector& Z_ind,
                         int k) {
  int n = ZCurr.nrow();
  int m = ZCurr.ncol();


  NumericVector zeta(n);
  NumericVector zeta2(n);


    for (int j = 0; j < n; j++) {
      int idx = Z_ind[j] - 1;
      NumericMatrix invR = invRCurr[idx];
      double val = ZCurr(j, k - 1);
      double invR_kk = invR(k - 1, k - 1);
      zeta[j] = val * val * (1.0 - invR_kk);
    }


    for (int kk = 0; kk < m; kk++) {
      if (kk == (k - 1)) continue;
      for (int j = 0; j < n; j++) {
        int idx = Z_ind[j] - 1;
        NumericMatrix invR = invRCurr[idx];
        zeta2[j] += ZCurr(j, k - 1) * ZCurr(j, kk) * invR(k - 1, kk);
      }
    }

  return List::create(Named("zeta") = zeta, Named("zeta2") = zeta2);
}
