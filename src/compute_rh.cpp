#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector compute_Rh(List invRCurr, IntegerVector Z_ind, int k) {
  int n = Z_ind.size();
  NumericVector Rh(n);
  //
  int kk = k - 1;

  for (int j = 0; j < n; j++) {
    //
    int listIndex = Z_ind[j] - 1;
    NumericMatrix mat = invRCurr[listIndex];
    Rh[j] = mat(kk, kk);
  }

  return Rh;
}
