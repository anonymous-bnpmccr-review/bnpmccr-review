#ifndef UPDATEZCURRGAMMA2_H
#define UPDATEZCURRGAMMA2_H

#include <RcppArmadillo.h>

arma::mat updateZCurrGamma2(
  const arma::mat& Y,
  const arma::mat& WB,
  arma::mat&       ZCurr,
  const arma::vec& gammaAlpha,
  int              countgamma,
  int              k
);

#endif // UPDATEZCURRGAMMA2_H
