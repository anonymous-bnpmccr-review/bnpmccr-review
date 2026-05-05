#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector compute_zeta(NumericMatrix ZCurr, List invRCurr, IntegerVector Z_ind, int k) {
  int n = ZCurr.nrow();
  int m = ZCurr.ncol();
  NumericVector zeta(n, 0.0);
  
  // C++는 0-indexed이므로, R의 k를 k-1로 사용합니다.
  for (int kk = 0; kk < m; kk++) {
    if (kk == k - 1) continue;  // k번째 열은 건너뜁니다.
    for (int j = 0; j < n; j++) {
      // Z_ind는 R에서 1-indexed이므로, C++에서는 0-indexed로 조정
      int listIndex = Z_ind[j] - 1;
      NumericMatrix invR = invRCurr[listIndex];
      zeta[j] += ZCurr(j, kk) * invR(k - 1, kk);
    }
  }
  return zeta;
}
