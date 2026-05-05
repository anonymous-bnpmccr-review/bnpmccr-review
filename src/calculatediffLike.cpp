// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]
#include <RcppArmadillo.h>
#include <cmath>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
List updateDiffLike(const Rcpp::List &RCurr,
                    const IntegerVector &Z_ind,
                    const arma::mat &ZCurr,
                    const int k,  // k: 0-indexed (R에서 k를 전달할 때는 k-1로 변환)
                    const arma::vec &TUProp,
                    const arma::vec &TLProp,
                    const arma::mat &TU,
                    const arma::mat &TL) {
  int n = ZCurr.n_rows;
  
  // 결과 벡터: diffLikeProp와 diffLike (길이 n)
  NumericVector diffLikeProp(n);
  NumericVector diffLike(n);
  
  for (int j = 0; j < n; ++j) {
    // RCurr[[Z_ind[j]]] : R에서 인덱스는 1-indexed이므로, C++에서는 1을 빼줍니다.
    int listIndex = Z_ind[j] - 1;
    arma::mat RMat = Rcpp::as<arma::mat>(RCurr[listIndex]);
    
    // RMat[k, -k] : k번째 행에서 k번째 열을 제외한 벡터.
    // (Armadillo는 0-indexed입니다.)
    arma::rowvec row_k = RMat.row(k);
    arma::rowvec row_k_excl = row_k;
    row_k_excl.shed_col(k);
    
    // RMat[-k, -k] : k번째 행과 열을 제거한 하위행렬.
    arma::mat RMat_sub = RMat;
    RMat_sub.shed_row(k);
    RMat_sub.shed_col(k);
    
    // RR = row_k_excl %*% inv(RMat_sub)
    arma::rowvec RR = row_k_excl * inv(RMat_sub);
    
    // mean_z = RR %*% ZCurr[j, -k]
    arma::rowvec Zrow = ZCurr.row(j);
    arma::rowvec Zrow_excl = Zrow;
    Zrow_excl.shed_col(k);
    double mean_z = as_scalar(RR * Zrow_excl.t());
    
    // sd_z = 1 - RR %*% (RMat[-k, k])
    // RMat[-k, k] : k번째 열에서 k번째 행을 제외한 벡터.
    arma::vec col_k = RMat.col(k);
    col_k.shed_row(k);
    double sd_z = 1 - as_scalar(RR * col_k);
    
    double sqrt_sd_z = std::sqrt(sd_z);
    
    // diffLikeProp[j] = pnorm((TUProp[j]-mean_z)/sqrt_sd_z) - pnorm((TLProp[j]-mean_z)/sqrt_sd_z)
    double p1 = R::pnorm((TUProp[j] - mean_z)/sqrt_sd_z, 0.0, 1.0, 1, 0);
    double p2 = R::pnorm((TLProp[j] - mean_z)/sqrt_sd_z, 0.0, 1.0, 1, 0);
    diffLikeProp[j] = p1 - p2;
    
    // diffLike[j] = pnorm((TU[j,k]-mean_z)/sqrt_sd_z) - pnorm((TL[j,k]-mean_z)/sqrt_sd_z)
    double p3 = R::pnorm((TU(j, k) - mean_z)/sqrt_sd_z, 0.0, 1.0, 1, 0);
    double p4 = R::pnorm((TL(j, k) - mean_z)/sqrt_sd_z, 0.0, 1.0, 1, 0);
    diffLike[j] = p3 - p4;
  }
  
  return List::create(Named("diffLikeProp") = diffLikeProp,
                      Named("diffLike") = diffLike);
}
