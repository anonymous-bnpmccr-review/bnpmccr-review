// updateZCurr2.h
#ifndef UPDATEZCURR_H
#define UPDATEZCURR_H

#include <RcppArmadillo.h>

// k는 0-based 인덱스
// RCurr : List of n R-matrices (each m×m)
// Z_ind : IntegerVector of length n
// ZCurr : n×m latent matrix (값 복사 또는 참조)
// TL, TU : n×m truncated normal limits
arma::mat updateZCurr_binom(
    const Rcpp::List& RCurr,
    const arma::ivec& Z_ind,
    arma::mat          ZCurr,
    int                k,
    const arma::mat   &TL,
    const arma::mat   &TU
);
#endif // UPDATEZCURR_H
