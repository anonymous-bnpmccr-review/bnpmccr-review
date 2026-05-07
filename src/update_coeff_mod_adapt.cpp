// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include "updateZCurrGamma2.h"   // header only
#include "updateZCurr.h"
#include <numeric>
#include <vector>
#include <algorithm>

using namespace Rcpp;
using namespace arma;

// Combine two matrices in block-diagonal form.
arma::mat bdiag(const arma::mat &A, const arma::mat &B) {
  arma::mat C = arma::zeros<arma::mat>(A.n_rows + B.n_rows, A.n_cols + B.n_cols);
  C(arma::span(0, A.n_rows - 1), arma::span(0, A.n_cols - 1)) = A;
  C(arma::span(A.n_rows, A.n_rows + B.n_rows - 1),
    arma::span(A.n_cols, A.n_cols + B.n_cols - 1)) = B;
  return C;
}


// Draw from a multivariate normal distribution using Cholesky decomposition.
arma::vec mvnrnd(const arma::vec &mu, const arma::mat &Sigma) {
  int n = mu.n_elem;
  arma::vec Y = arma::randn(n);
  arma::mat U = arma::chol(Sigma);
  return mu + U.t() * Y;
}

// [[Rcpp::export]]
Rcpp::List updateCoefficients2(
    const int k,                         // response index (0-indexed)
    const std::string responseType,      // "Gaussian", "binary", "ordinal", "binomial", "negative binomial", "gamma"
    const arma::vec &d,                  // vector d of length m
    arma::mat& ZCurr,                    // n x m current latent responses
    const arma::mat &WB,                 // n x m working response or offset
    const Rcpp::List &invRCurr,          // list of inverse matrices or related objects
    const arma::ivec &Z_ind,             // index for accessing invRCurr
    const arma::vec &Gammacurr,          // current Gamma for variable selection
    const int start, const int end,      // subvector range of Gammacurr
    const arma::mat &totalW,             // full design matrix
    const arma::mat &Xmat,               // fixed-effect design matrix
    const int p,                         // number of regression coefficients per response
    const int num_basis,                 // number of basis functions
    const arma::vec &g,                  // vector g
    const arma::mat &blockPriorInv,      // inverse prior matrix for the first block
    const arma::vec &Rh,                 // weight vector of length n
    const arma::mat &Y,                  // n x m observed responses for Gaussian cases
    arma::mat var_z,                     // n x m variance matrix, passed by value
    const arma::mat &W,                  // design matrix used for block construction
    const arma::vec &newGamma,           // newGamma vector for the MH step
    const arma::mat &Wgammastar,         // proposal design matrix for the MH step
    int countgamma,                      // counter for gamma responses
    const arma::vec& gammaAlpha,         // gammaAlpha vector for gamma responses
    int countnegbinom,                   // counter for negative binomial responses
    const arma::vec& negBinomPar,
    int countbinom,                      // counter for binomial responses
    const arma::vec& nTrial,
    const Rcpp::List& RCurr,
    const arma::mat& dpropormat,
    const int s
) {
  const int n = ZCurr.n_rows;
  const int m = ZCurr.n_cols;

  arma::vec ZTilde_k;
  arma::mat inv_bdiagmat_OLD;
  arma::mat inv_V_OLD;
  arma::vec muprop_OLD;
  arma::vec ZTilde_k_new;
  arma::vec TUProp(n), TLProp(n);
  arma::mat ZCurr_new;

  if(responseType == "Gaussian") {
    ZTilde_k = d[k] * ZCurr.col(k) + WB.col(k);

    // Compute zeta from dependence with the other responses.
    arma::vec zeta(n, fill::zeros);
    for (int j = 0; j < n; ++j) {
      arma::mat invR = Rcpp::as<arma::mat>(invRCurr[ Z_ind[j] - 1 ]);
      for (int kk = 0; kk < m; ++kk) {
        if (kk == k) continue;
        zeta(j) += ZCurr(j, kk) * invR(k, kk) / d[k];
      }
    }

    // Select relative indices with positive entries in Gammacurr[start:end].
    std::vector<int> gammaind;
    for (int i = start; i <= end; ++i) {
      if (Gammacurr(i) > 0)
        gammaind.push_back(i - start);
    }
    arma::uvec gammaind_arma(gammaind.size());
    for (size_t i = 0; i < gammaind.size(); ++i)
      gammaind_arma(i) = gammaind[i];
    arma::mat Wgamma = totalW.cols(gammaind_arma);

    std::vector<int> gammaind1;
    for (size_t i = p; i < gammaind.size(); ++i)
      gammaind1.push_back(gammaind[i] - p);

    arma::vec sub_g = g.subvec(k * p, (k + 1) * p - 1);

    std::vector<int> subgammind;
    for (size_t i = 0; i < gammaind1.size(); ++i) {
      if (gammaind1[i] <= (num_basis - 1) - 1)
        subgammind.push_back(gammaind1[i]);
    }
    arma::uvec subgammind_arma(subgammind.size());
    for (size_t i = 0; i < subgammind.size(); ++i)
      subgammind_arma(i) = subgammind[i];

    arma::mat W_sub = W.cols(subgammind_arma);
    arma::mat bdiagmat = (1.0 / sub_g(0)) * (W_sub.t() * W_sub);
    for (int q = 2; q <= p; ++q) {
      int lowerBound = (q - 1) * (num_basis - 1);
      int upperBound = q * (num_basis - 1) - 1;
      std::vector<int> temp;
      for (size_t i = 0; i < gammaind1.size(); ++i) {
        if (gammaind1[i] >= lowerBound && gammaind1[i] <= upperBound)
          temp.push_back(gammaind1[i]);
      }
      if (!temp.empty()){
        arma::uvec temp_arma(temp.size());
        for (size_t i = 0; i < temp.size(); ++i)
          temp_arma(i) = temp[i];
        arma::mat W_temp = W.cols(temp_arma);
        arma::mat block = (1.0 / sub_g(q - 1)) * (W_temp.t() * W_temp);
        bdiagmat = bdiag(bdiagmat, block);
      }
    }
    bdiagmat = bdiag(blockPriorInv, bdiagmat);

    arma::mat V = (1.0 / (d[k] * d[k])) * (Wgamma.t() * (Wgamma.each_col() % Rh)) + bdiagmat;
    arma::mat inv_V = arma::inv(V);
    inv_V = 0.5 * (inv_V + inv_V.t());

    arma::vec muprop = inv_V * (Wgamma.t() * ((1.0 / (d[k] * d[k])) * (Rh % Y.col(k)) + zeta));
    arma::vec coefprop = mvnrnd(muprop, inv_V);
    arma::vec Bprop = coefprop;

    return List::create(Named("Bprop") = Bprop,
                        Named("ZTilde") = ZTilde_k);

  } else if(responseType == "binary" || responseType == "ordinal") {

    var_z.col(k).fill(1.0);
    ZTilde_k = ZCurr.col(k) + WB.col(k);

    arma::vec zeta(n, fill::zeros);
    for (int j = 0; j < n; ++j) {
      arma::mat invR = Rcpp::as<arma::mat>(invRCurr[ Z_ind[j] - 1 ]);
      for (int kk = 0; kk < m; ++kk) {
        if (kk == k) continue;
        zeta(j) += ZCurr(j, kk) * invR(k, kk);
      }
    }
    std::vector<int> gammaind;
    for (int i = start; i <= end; ++i) {
      if (Gammacurr(i) > 0)
        gammaind.push_back(i - start);
    }
    arma::uvec gammaind_arma(gammaind.size());
    for (size_t i = 0; i < gammaind.size(); ++i)
      gammaind_arma(i) = gammaind[i];
    arma::mat Wgamma = totalW.cols(gammaind_arma);

    std::vector<int> gammaind1;
    for (size_t i = p; i < gammaind.size(); ++i)
      gammaind1.push_back(gammaind[i] - p);

    arma::vec sub_g = g.subvec(k * p, (k + 1) * p - 1);
    std::vector<int> subgammind;
    for (size_t i = 0; i < gammaind1.size(); ++i) {
      if (gammaind1[i] <= (num_basis - 1) - 1)
        subgammind.push_back(gammaind1[i]);
    }
    arma::uvec subgammind_arma(subgammind.size());
    for (size_t i = 0; i < subgammind.size(); ++i)
      subgammind_arma(i) = subgammind[i];

    arma::mat W_sub = W.cols(subgammind_arma);
    arma::mat bdiagmat = (1.0 / sub_g(0)) * (W_sub.t() * W_sub);
    for (int q = 2; q <= p; ++q) {
      int lowerBound = (q - 1) * (num_basis - 1);
      int upperBound = q * (num_basis - 1) - 1;
      std::vector<int> temp;
      for (size_t i = 0; i < gammaind1.size(); ++i) {
        if (gammaind1[i] >= lowerBound && gammaind1[i] <= upperBound)
          temp.push_back(gammaind1[i]);
      }
      if (!temp.empty()){
        arma::uvec temp_arma(temp.size());
        for (size_t i = 0; i < temp.size(); ++i)
          temp_arma(i) = temp[i];
        arma::mat W_temp = W.cols(temp_arma);
        arma::mat block = (1.0 / sub_g(q - 1)) * (W_temp.t() * W_temp);
        bdiagmat = bdiag(bdiagmat, block);
      }
    }
    bdiagmat = bdiag(blockPriorInv, bdiagmat);

    arma::mat V = (Wgamma.t() * (Wgamma.each_col() % Rh)) + bdiagmat;
    arma::mat inv_V = arma::inv(V);
    inv_V = 0.5 * (inv_V + inv_V.t());
    arma::vec muprop = inv_V * (Wgamma.t() * ((Rh % ZTilde_k) + zeta));
    arma::vec coefprop = mvnrnd(muprop, inv_V);
    arma::vec Bprop = coefprop;

    return List::create(Named("Bprop") = Bprop,
                        Named("ZTilde") = ZTilde_k);

  } else {
    if(responseType == "binomial") {
      countbinom++;
      ZTilde_k = var_z.col(k) % ZCurr.col(k) + WB.col(k);
    } else if(responseType == "negative binomial"){
      countnegbinom++;
      ZTilde_k = var_z.col(k) % ZCurr.col(k) + WB.col(k);
    } else if(responseType == "gamma") {
      countgamma++;
      ZTilde_k = var_z.col(k) % ZCurr.col(k) + WB.col(k);
    }

    arma::vec zeta(n, fill::zeros);
    for (int j = 0; j < n; ++j) {
      arma::mat invR = Rcpp::as<arma::mat>(invRCurr[ Z_ind[j] - 1 ]);
      for (int kk = 0; kk < m; ++kk) {
        if (kk == k) continue;
        zeta(j) += ZCurr(j, kk) * invR(k, kk) / var_z(j, k);
      }
    }


    arma::mat Wk = Wgammastar;
    std::vector<int> gammastarind;
    for (size_t i = p; i < newGamma.n_elem; ++i) {
      if (newGamma(i) > 0)
        gammastarind.push_back(i - p);
    }
    arma::vec sub_g = g.subvec(k * p, (k + 1) * p - 1);
    std::vector<int> subgammind;
    for (size_t i = 0; i < gammastarind.size(); ++i) {
      if (gammastarind[i] <= (num_basis - 1) - 1)
        subgammind.push_back(gammastarind[i]);
    }
    arma::uvec subgammind_arma(subgammind.size());
    for (size_t i = 0; i < subgammind.size(); ++i)
      subgammind_arma(i) = subgammind[i];
    arma::mat W_sub = W.cols(subgammind_arma);
    arma::mat bdiagmat = (1.0 / sub_g(0)) * (W_sub.t() * W_sub);
    for (int q = 2; q <= p; ++q) {
      int lowerBound = (q - 1) * (num_basis - 1);
      int upperBound = q * (num_basis - 1) - 1;
      std::vector<int> temp;
      for (size_t i = 0; i < gammastarind.size(); ++i) {
        if (gammastarind[i] >= lowerBound && gammastarind[i] <= upperBound)
          temp.push_back(gammastarind[i]);
      }
      if (!temp.empty()){
        arma::uvec temp_arma(temp.size());
        for (size_t i = 0; i < temp.size(); ++i)
          temp_arma(i) = temp[i];
        arma::mat W_temp = W.cols(temp_arma);
        arma::mat block = (1.0 / sub_g(q - 1)) * (W_temp.t() * W_temp);
        bdiagmat = bdiag(bdiagmat, block);
      }
    }

    bdiagmat = bdiag(blockPriorInv, bdiagmat);
    arma::mat inv_bdiagmat = arma::inv(bdiagmat);

    arma::mat V = (Wk.t() * arma::diagmat(1.0 / arma::square(var_z.col(k))) * (Wk.each_col() % Rh)) + bdiagmat;
    arma::mat inv_V = arma::inv(V);
    inv_V = 0.5 * (inv_V + inv_V.t());
    arma::vec muprop = inv_V * Wk.t() * (arma::diagmat(1.0 / arma::square(var_z.col(k))) * (Rh % ZTilde_k) + zeta);
    arma::vec coefprop = mvnrnd(muprop, inv_V);
    arma::vec Bprop = coefprop;

    std::vector<int> gammaind;
    for (int i = start; i <= end; ++i) {
      if (Gammacurr(i) > 0)
        gammaind.push_back(i - start);
    }
    arma::uvec gammaind_arma(gammaind.size());
    for (size_t i = 0; i < gammaind.size(); ++i)
      gammaind_arma(i) = gammaind[i];
    arma::mat Wgamma = totalW.cols(gammaind_arma);
    std::vector<int> gammaind1;
    for (size_t i = p; i < gammaind.size(); ++i)
      gammaind1.push_back(gammaind[i] - p);
    sub_g = g.subvec(k * p, (k + 1) * p - 1);
    std::vector<int> subgammind_current;
    for (size_t i = 0; i < gammaind1.size(); ++i) {
      if (gammaind1[i] <= (num_basis - 1) - 1)
        subgammind_current.push_back(gammaind1[i]);
    }
    arma::uvec subgammind_arma_current(subgammind_current.size());
    for (size_t i = 0; i < subgammind_current.size(); ++i)
      subgammind_arma_current(i) = subgammind_current[i];
    arma::mat W_sub_current = W.cols(subgammind_arma_current);
    arma::mat bdiagmat_current = (1.0 / sub_g(0)) * (W_sub_current.t() * W_sub_current);
    for (int q = 2; q <= p; ++q) {
      int lowerBound = (q - 1) * (num_basis - 1);
      int upperBound = q * (num_basis - 1) - 1;
      std::vector<int> temp;
      for (size_t i = 0; i < gammaind1.size(); ++i) {
        if (gammaind1[i] >= lowerBound && gammaind1[i] <= upperBound)
          temp.push_back(gammaind1[i]);
      }
      if (!temp.empty()){
        arma::uvec temp_arma(temp.size());
        for (size_t i = 0; i < temp.size(); ++i)
          temp_arma(i) = temp[i];
        arma::mat W_temp = W.cols(temp_arma);
        arma::mat block = (1.0 / sub_g(q - 1)) * (W_temp.t() * W_temp);
        bdiagmat_current = bdiag(bdiagmat_current, block);
      }
    }

    if(responseType == "gamma") {
      arma::mat WB_new = Wk * Bprop;
      ZCurr_new = updateZCurrGamma2(Y, WB_new, ZCurr, gammaAlpha, countgamma, k + 1);
      ZTilde_k_new = var_z.col(k) % dpropormat.col(k) %  ZCurr_new.col(k) + WB_new;
      bdiagmat_current = bdiag(blockPriorInv, bdiagmat_current);
      inv_bdiagmat_OLD = arma::inv(bdiagmat_current);
      arma::mat Wk_OLD = Wgamma;
      arma::mat V_OLD = (Wk_OLD.t() * arma::diagmat(1.0 / arma::square(var_z.col(k))) * (Wk_OLD.each_col() % Rh)) + bdiagmat_current;
      inv_V_OLD = arma::inv(V_OLD);
      inv_V_OLD = 0.5 * (inv_V_OLD + inv_V_OLD.t());
      muprop_OLD = inv_V_OLD * Wk_OLD.t() * (arma::diagmat(1.0 / arma::square(var_z.col(k))) * (Rh % ZTilde_k_new) + zeta);

    } else{
      arma::mat WB_new = Wk * Bprop;

      if(responseType == "binomial"){
        for(int i = 0; i < n; ++i) {
          double eta   = WB_new(i);
          double p     = 1.0 / (1.0 + std::exp(-eta));

          // Upper bound: P(Y[i,k] <= y)
          double c_up  = R::pbinom( Y(i,k),
                                    nTrial(countbinom - 1),
                                    p,
                                    /*lower.tail=*/1, /*log.p=*/0 );
          TUProp(i)    = R::qnorm(c_up, 0.0, 1.0, /*lower.tail=*/1, /*log.p=*/0);

          // Lower bound: P(Y[i,k] <= y - 1)
          double c_lo  = R::pbinom( Y(i,k) - 1,
                                    nTrial(countbinom - 1),
                                    p,
                                    1, 0 );
          TLProp(i)    = R::qnorm(c_lo, 0.0, 1.0, 1, 0);
        }

        ZCurr_new = updateZCurr_binom(RCurr, Z_ind, ZCurr, k, TLProp, TUProp);
        ZTilde_k_new = var_z.col(k) % ZCurr_new.col(k) + WB_new;


      }else if(responseType == "negative binomial"){ // negative binomial

        for(int i = 0; i < n; ++i) {
          double eta = WB_new(i);
          double mu  = std::exp(eta);

          // Upper bound: P(Y[i,k] <= y)
          double c_up = R::pnbinom_mu(
            Y(i,k),
            negBinomPar(countnegbinom - 1),
            mu,
            /*lower.tail=*/1,
            /*log.p=*/0
          );
          TUProp(i) = R::qnorm(
            c_up, 0.0, 1.0,
            /*lower.tail=*/1,
            /*log.p=*/0
          );

          // Lower bound: P(Y[i,k] <= y - 1)
          double c_lo = R::pnbinom_mu(
            Y(i,k) - 1,
            negBinomPar(countnegbinom - 1),
            mu,
            1, 0
          );
          TLProp(i) = R::qnorm(c_lo, 0.0, 1.0, 1, 0);
        }

        ZCurr_new = updateZCurr_binom(RCurr, Z_ind, ZCurr, k, TLProp, TUProp);
        ZTilde_k_new = var_z.col(k) % ZCurr_new.col(k) + WB_new;
      }

      bdiagmat_current = bdiag(blockPriorInv, bdiagmat_current);
      inv_bdiagmat_OLD = arma::inv(bdiagmat_current);
      arma::mat Wk_OLD = Wgamma;
      arma::mat V_OLD = (Wk_OLD.t() * arma::diagmat(1.0 / arma::square(var_z.col(k))) * (Wk_OLD.each_col() % Rh)) + bdiagmat_current;
      inv_V_OLD = arma::inv(V_OLD);
      inv_V_OLD = 0.5 * (inv_V_OLD + inv_V_OLD.t());
      muprop_OLD = inv_V_OLD * Wk_OLD.t() * (arma::diagmat(1.0 / arma::square(var_z.col(k))) * (Rh % ZTilde_k_new) + zeta);

    }


    return List::create(Named("Bprop") = Bprop,
                        Named("muprop_OLD") = muprop_OLD,
                        Named("inv_V_OLD") = inv_V_OLD,
                        Named("inv_bdiagmat_OLD") = inv_bdiagmat_OLD,
                        Named("muprop") = muprop,
                        Named("inv_V") = inv_V,
                        Named("inv_bdiagmat") = inv_bdiagmat,
                        Named("ZTilde") = ZTilde_k,
                        Named("var_z") = var_z,
                        Named("countgamma") = countgamma,
                        Named("ZTilde_new") = ZTilde_k_new,
                        Named("countnegbinom") = countnegbinom,
                        Named("negBinomPar") = negBinomPar,
                        Named("countbinom") = countbinom,
                        Named("nTrial") = nTrial,
                        Named("TLProp") = TLProp,
                        Named("TUProp") = TUProp);
  }
}
