// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <limits>
#include <cmath>
using namespace Rcpp;

//--------------------------------------------------------------------------
// rtruncnorm: 역변환법을 이용하여 a ~ b로 잘린 정규분포 (mean, sd)를 샘플링
//--------------------------------------------------------------------------
static inline double rtruncnorm(double a, double b, double mean, double sd) {
  // 표준화: Z = (X - mean) / sd
  double alpha = (a - mean) / sd;
  double beta  = (b - mean) / sd;

  double Phi_alpha, Phi_beta;
  // a가 -무한대인 경우
  if (std::isinf(alpha) && alpha < 0) {
    Phi_alpha = 0.0;
  } else {
    Phi_alpha = R::pnorm(alpha, 0.0, 1.0, 1, 0);
  }
  // b가 +무한대인 경우
  if (std::isinf(beta) && beta > 0) {
    Phi_beta = 1.0;
  } else {
    Phi_beta = R::pnorm(beta, 0.0, 1.0, 1, 0);
  }

  // Uniform(0,1) 표본
  double u = R::runif(0.0, 1.0);
  // 역변환: p = Phi_alpha + u*(Phi_beta-Phi_alpha)
  double p = Phi_alpha + u * (Phi_beta - Phi_alpha);
  // 표준 정규분포의 분위수를 구한 후 원래의 scale로 변환
  double z = R::qnorm(p, 0.0, 1.0, 1, 0);
  return mean + sd * z;
}

//--------------------------------------------------------------------------
// update_Z_star: Z_star 행렬을 업데이트하는 함수
//
// 인자:
//   Z_ind   : 각 행에 대해 선택된 인덱스 (길이 n, R의 1-indexed)
//   WB_x    : n x H 행렬
//   Z_star  : n x H 행렬 (업데이트할 대상)
//   alpha_h : 길이 H인 벡터
//   H       : 범주의 총 개수 (R 기준 1-indexed)
// 반환값:
//   업데이트된 Z_star 행렬
//--------------------------------------------------------------------------

// [[Rcpp::export]]
NumericMatrix update_Z_star(IntegerVector Z_ind,
                            NumericMatrix WB_x,
                            NumericMatrix Z_star,
                            NumericVector alpha_h,
                            int H) {
  int n = WB_x.nrow();

  // 각 행에 대해 업데이트
  for (int j = 0; j < n; j++) {
    // R의 1-indexed 값을 0-indexed로 변환
    int i = Z_ind[j] - 1;

    // 만약 i > 1 (즉, 원래 R에서는 i>1)인 경우, 1:(i-1) 에 해당하는 열 업데이트
    // C++에서는 0부터 (i-1)-까지
    if (i >= 1) {
      for (int h = 0; h < i; h++) {
        double mean_val = WB_x(j, h) + alpha_h[h];
        Z_star(j, h) = rtruncnorm(-std::numeric_limits<double>::infinity(),
               0.0,
               mean_val,
               1.0);
      }
    }

    // 만약 i != H 인 경우 (R에서는 if(i != H))
    // C++에서는 H가 1-indexed이므로 비교할 때 (H - 1)과 비교
    if (i != (H - 1)) {
      double mean_val = WB_x(j, i) + alpha_h[i];
      Z_star(j, i) = rtruncnorm(0.0,
             std::numeric_limits<double>::infinity(),
             mean_val,
             1.0);
    }
  }

  return Z_star;
}
