#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
List computeZetaAndZeta2(const NumericMatrix& ZCurr, 
                         const List& invRCurr, 
                         const IntegerVector& Z_ind, 
                         int k) {
  int n = ZCurr.nrow();
  int m = ZCurr.ncol();
  
  // 초기화: 각 벡터의 길이는 n (행의 수)
  NumericVector zeta(n);
  NumericVector zeta2(n);
  
  // zeta 계산: for (j in 1:n)
    for (int j = 0; j < n; j++) {
      // R은 1-indexed이므로, C++에서는 Z_ind[j] 값에서 1을 빼야 함
      int idx = Z_ind[j] - 1;
      NumericMatrix invR = invRCurr[idx];
      // k도 R에서는 1-indexed이므로, C++에서는 (k-1)를 사용
      double val = ZCurr(j, k - 1);
      double invR_kk = invR(k - 1, k - 1);
      zeta[j] = val * val * (1.0 - invR_kk);
    }
  
  // zeta2 계산: for (kk in 1:m, kk != k)
    for (int kk = 0; kk < m; kk++) {
      if (kk == (k - 1)) continue;  // k번째 컬럼은 건너뜀
      for (int j = 0; j < n; j++) {
        int idx = Z_ind[j] - 1;
        NumericMatrix invR = invRCurr[idx];
        zeta2[j] += ZCurr(j, k - 1) * ZCurr(j, kk) * invR(k - 1, kk);
      }
    }
  
  return List::create(Named("zeta") = zeta, Named("zeta2") = zeta2);
}
