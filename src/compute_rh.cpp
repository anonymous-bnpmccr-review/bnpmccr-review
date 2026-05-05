#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector compute_Rh(List invRCurr, IntegerVector Z_ind, int k) {
  int n = Z_ind.size();
  NumericVector Rh(n);
  // R의 k는 1-indexed이므로, C++에서는 k-1을 사용
  int kk = k - 1;
  
  for (int j = 0; j < n; j++) {
    // R의 Z_ind는 1-indexed이므로, C++에서는 Z_ind[j]-1 사용
    int listIndex = Z_ind[j] - 1;
    NumericMatrix mat = invRCurr[listIndex];
    Rh[j] = mat(kk, kk);
  }
  
  return Rh;
}
