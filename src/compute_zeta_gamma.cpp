#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector computeZetaGamma(const NumericMatrix& ZCurr,
                          const IntegerVector& Z_ind,
                          const List& invRCurr,
                          int k) {
  // n: 관측치 개수 (행 수)
  int n = ZCurr.nrow();
  
  // 결과 벡터 zeta 초기화
  NumericVector zeta(n, 0.0);
  
  // R 코드: for(j in 1:n)
    // Rcpp 는 0-based 이므로 j=0부터 시작
  for (int j = 0; j < n; j++) {
    // R 코드: invRCurr[[Z_ind[j]]] 에 해당하므로
    // Z_ind[j]가 R에서 1-based라면 C++에선 -1 해줘야 함
    NumericMatrix invR = invRCurr[ Z_ind[j] - 1 ];
    
    // R 코드: zeta[j] <- zeta[j] + ZCurr[j,k]^2 * (1 - invRCurr[[Z_ind[j]]][k,k])
    // k 역시 R에서 1-based라면 C++에선 (k-1)로 인덱싱
    double val = std::pow(ZCurr(j, k - 1), 2.0) * (1.0 - invR(k - 1, k - 1));
    zeta[j] += val;
  }
  
  return zeta;
}
