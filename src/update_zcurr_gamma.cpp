#include <Rcpp.h>
#include <cmath>
using namespace Rcpp;


// Y, WB, and ZCurr are n x m matrices.
// gammaAlpha is a vector, and countgamma and k are integers.
// Since R uses 1-based indexing, subtract 1 from k and countgamma in C++.

// [[Rcpp::export]]
NumericVector updateZCurrGamma(NumericMatrix Y, NumericMatrix WB, NumericMatrix& ZCurr,
                               NumericVector gammaAlpha, int countgamma, int k) {
  int n = Y.nrow();
  int col = k - 1;
  int gamma_idx = countgamma - 1;
  double shape = gammaAlpha[gamma_idx];

  for (int i = 0; i < n; i++) {
    double x = Y(i, col);
    double wb_val = WB(i, col);
    double scale = std::exp(wb_val) / shape;

    double gammacdf = R::pgamma(x, shape, scale, 1, 0);

    if (gammacdf == 0.0) {
      gammacdf = 1e-12;
    } else if (gammacdf == 1.0) {
      gammacdf = 0.99999999999999;
    }

    double qnorm_val = R::qnorm(gammacdf, 0.0, 1.0, 1, 0);
    ZCurr(i, col) = qnorm_val;
  }

  return ZCurr;
}
