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
                    const int k,  // k is 0-indexed; pass k-1 from R.
                    const arma::vec &TUProp,
                    const arma::vec &TLProp,
                    const arma::mat &TU,
                    const arma::mat &TL) {
  int n = ZCurr.n_rows;

  // Output vectors of length n.
  NumericVector diffLikeProp(n);
  NumericVector diffLike(n);

  for (int j = 0; j < n; ++j) {
    // RCurr[[Z_ind[j]]] in R is 1-indexed, so subtract 1 in C++.
    int listIndex = Z_ind[j] - 1;
    arma::mat RMat = Rcpp::as<arma::mat>(RCurr[listIndex]);

    // RMat[k, -k]: kth row excluding the kth column.
    arma::rowvec row_k = RMat.row(k);
    arma::rowvec row_k_excl = row_k;
    row_k_excl.shed_col(k);

    // RMat[-k, -k]: submatrix after removing the kth row and kth column.
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

    // sd_z = 1 - RR %*% RMat[-k, k]
    arma::vec col_k = RMat.col(k);
    col_k.shed_row(k);
    double sd_z = 1 - as_scalar(RR * col_k);

    double sqrt_sd_z = std::sqrt(sd_z);

    double p1 = R::pnorm((TUProp[j] - mean_z)/sqrt_sd_z, 0.0, 1.0, 1, 0);
    double p2 = R::pnorm((TLProp[j] - mean_z)/sqrt_sd_z, 0.0, 1.0, 1, 0);
    diffLikeProp[j] = p1 - p2;

    double p3 = R::pnorm((TU(j, k) - mean_z)/sqrt_sd_z, 0.0, 1.0, 1, 0);
    double p4 = R::pnorm((TL(j, k) - mean_z)/sqrt_sd_z, 0.0, 1.0, 1, 0);
    diffLike[j] = p3 - p4;
  }

  return List::create(Named("diffLikeProp") = diffLikeProp,
                      Named("diffLike") = diffLike);
}
