// updateZCurr2.h
#ifndef UPDATEZCURR_H
#define UPDATEZCURR_H

#include <RcppArmadillo.h>

arma::mat updateZCurr_binom(
    const Rcpp::List& RCurr,
    const arma::ivec& Z_ind,
    arma::mat          ZCurr,
    int                k,
    const arma::mat   &TL,
    const arma::mat   &TU
);
#endif // UPDATEZCURR_H
