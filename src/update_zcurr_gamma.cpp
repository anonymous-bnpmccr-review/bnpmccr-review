#include <Rcpp.h>
#include <cmath>
using namespace Rcpp;


// Y, WB, ZCurr는 모두 n x m 행렬이고,
// gammaAlpha는 벡터, countgamma와 k는 정수입니다.
// 주의: R에서는 1-indexing이므로, C++에서는 k와 countgamma에 대해 -1 해줍니다.

// [[Rcpp::export]]
NumericVector updateZCurrGamma(NumericMatrix Y, NumericMatrix WB, NumericMatrix& ZCurr,
                 NumericVector gammaAlpha, int countgamma, int k) {
  int n = Y.nrow();
  int col = k - 1;         // R의 k번째 열 -> C++에서는 col = k-1
  int gamma_idx = countgamma - 1;  // R의 countgamma -> C++에서는 countgamma-1
  double shape = gammaAlpha[gamma_idx];

  for (int i = 0; i < n; i++) {
    double x = Y(i, col);
    double wb_val = WB(i, col);
    // R의 pgamma(..., rate = gammaAlpha[countgamma]/exp(WB[i,k]))와 동일하게 scale 계산:
      double scale = std::exp(wb_val) / shape;

      // lower.tail = TRUE, log.p = FALSE
      double gammacdf = R::pgamma(x, shape, scale, 1, 0);

      if (gammacdf == 0.0) {
        gammacdf = 1e-12;  // 0.000000000001
      } else if (gammacdf == 1.0) {
        gammacdf = 0.99999999999999;
      }

      // qnorm(gammacdf) : 평균 0, 표준편차 1, lower.tail = TRUE, log.p = FALSE
      double qnorm_val = R::qnorm(gammacdf, 0.0, 1.0, 1, 0);
      ZCurr(i, col) = qnorm_val;
  }

  return ZCurr;
}
