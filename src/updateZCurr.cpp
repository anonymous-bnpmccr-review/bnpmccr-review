// updateZCurr.cpp
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]
#include "updateZCurr.h"
#include <cmath>

using namespace Rcpp;
using namespace arma;

//-----------------------------------------------------------------
// truncated normal sampler: [a, b] 구간, mean, sd
// 입력: a (하한), b (상한), mean, sd
// 반환: truncated normal 분포에서 샘플링한 값
static inline double rtruncnorm(double a, double b, double mean, double sd) {
  // 표준화
  double alpha = (a - mean) / sd;
  double beta  = (b - mean) / sd;

  double Phi_alpha, Phi_beta;
  if(std::isinf(alpha) && alpha < 0)
    Phi_alpha = 0.0;
  else
    Phi_alpha = R::pnorm(alpha, 0.0, 1.0, 1, 0);

  if(std::isinf(beta) && beta > 0)
    Phi_beta = 1.0;
  else
    Phi_beta = R::pnorm(beta, 0.0, 1.0, 1, 0);

  double u = R::runif(0.0, 1.0);
  double p = Phi_alpha + u * (Phi_beta - Phi_alpha);

  // p 값이 극단값으로 인한 문제를 방지: [eps, 1-eps] 구간으로 클램핑
  double eps = 1e-12;
  if (p < eps) p = eps;
  if (p > 1.0 - eps) p = 1.0 - eps;

  double z = R::qnorm(p, 0.0, 1.0, 1, 0);
  return mean + sd * z;
}

//-----------------------------------------------------------------
// updateZCurr: RCurr, Z_ind, ZCurr, TL, TU를 사용해
// 각 관측 i에 대해 k번째 열의 latent variable Z를 truncated normal에서 샘플링하여 업데이트합니다.
// k는 C++에서는 0-indexed로 전달되어야 합니다.
// [[Rcpp::export]]
arma::mat updateZCurr_binom(const Rcpp::List& RCurr,
                      const arma::ivec& Z_ind,
                      arma::mat ZCurr,
                      const int k,        // k: 0-indexed (예: R에서 k를 전달할 때는 k-1로 변환)
                      const arma::mat &TL, // 하한값 행렬 (n x m)
                      const arma::mat &TU  // 상한값 행렬 (n x m)
) {
  int n = ZCurr.n_rows;

  // 각 관측 i에 대해 latent variable Z의 k번째 열 업데이트
  for (int i = 0; i < n; ++i) {
    // RCurr[[Z_ind[i]]] : R의 인덱스는 1-indexed이므로 C++에서는 1을 빼줍니다.
    int listIndex = Z_ind[i] - 1;
    arma::mat RMat = Rcpp::as<arma::mat>(RCurr[listIndex]);

    // 1. RMat[k, -k] : k번째 행에서 k번째 열을 제거한 row vector
    arma::rowvec row_k = RMat.row(k);
    arma::rowvec row_k_excl = row_k; // 복사 후
    row_k_excl.shed_col(k);          // k번째 열 제거

    // 2. RMat[-k, -k] : k번째 행과 k번째 열을 제거한 하위 행렬
    arma::mat RMat_sub = RMat;
    RMat_sub.shed_row(k);
    RMat_sub.shed_col(k);

    // 3. RR = row_k_excl * inv(RMat_sub)
    arma::rowvec RR = row_k_excl * inv(RMat_sub);

    // 4. mean_z = RR %*% ZCurr[i, -k]
    arma::rowvec Zrow = ZCurr.row(i);
    arma::rowvec Zrow_excl = Zrow;
    Zrow_excl.shed_col(k);
    double mean_z = as_scalar(RR * Zrow_excl.t());

    // 5. sd_z = sqrt( 1 - RR %*% (RMat[-k, k]) )
    // RMat[-k, k] : k번째 열에서 k번째 행을 제거한 column vector
    arma::vec col_k = RMat.col(k);
    col_k.shed_row(k);
    double sd_z = std::sqrt(1 - as_scalar(RR * col_k));

    // 6. 업데이트: ZCurr[i, k] <- rtruncnorm(a= TL(i,k), b= TU(i,k), mean=mean_z, sd=sd_z)
    double newVal = rtruncnorm(TL(i,0), TU(i,0), mean_z, sd_z);
    ZCurr(i, k) = newVal;
  }

  return ZCurr;
}
