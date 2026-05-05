#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector computeZetaGamma2(const NumericMatrix& ZCurr,
                           const IntegerVector& Z_ind,
                           const List& invRCurr,
                           int k) {
  // n: 관측치 개수 (행 수)
  int n = ZCurr.nrow();
  // m: ZCurr의 열 개수
  int m = ZCurr.ncol();
  
  // 결과 벡터 zeta2 초기화
  NumericVector zeta2(n, 0.0);
  
  // R 코드: ind <- 1:m
  //         for(kk in ind[-k]) { ... }
  //   -> C++에서는 0-based 이므로 kk = 0..m-1 중에서 (k-1)은 제외
  for (int kk = 0; kk < m; kk++) {
    if (kk == (k - 1)) {
      continue;  
    }
    // R 코드: for(j in 1:n)
      //   -> C++는 0-based이므로 j=0..(n-1)
      for (int j = 0; j < n; j++) {
        // R 코드: invRCurr[[Z_ind[j]]]
        //   -> Z_ind[j]가 R에서 1-based라면 C++에서 -1 처리
        NumericMatrix invR = invRCurr[ Z_ind[j] - 1 ];
        
        // R 코드:
          //   zeta2[j] <- zeta2[j] + ZCurr[j,k] * ZCurr[j,kk] * invRCurr[[Z_ind[j]]][k,kk]
        //   여기서 k, kk 모두 R에서 1-based -> C++에선 k-1, kk 사용
        zeta2[j] += ZCurr(j, k - 1) * ZCurr(j, kk) * invR(k - 1, kk);
      }
  }
  
  return zeta2;
}
