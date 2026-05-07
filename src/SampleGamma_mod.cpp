// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <numeric>
#include <vector>
#include <algorithm>
using namespace Rcpp;
using namespace arma;

// ===== Cholesky-based calculation =====
struct CholStats {
  double logdet;  // log|Σ|
  double quad;    // μᵀ Σ^{-1} μ
};

inline CholStats chol_stats(arma::mat Sigma, const arma::vec& mu) {
  // stability
  Sigma = 0.5 * (Sigma + Sigma.t());
  arma::mat L;
  double eps = 1e-10;
  bool ok = arma::chol(L, Sigma, "lower");
  for (int t = 0; !ok && t < 6; ++t) {
    Sigma.diag() += eps;
    ok = arma::chol(L, Sigma, "lower");
    eps *= 10.0;
  }
  if (!ok) Rcpp::stop("chol failed (Sigma not SPD)");

  const double logdet = 2.0 * arma::sum(arma::log(L.diag()));
  arma::vec y = arma::solve(arma::trimatl(L), mu); // L y = μ → y = L^{-1} μ
  const double quad = arma::dot(y, y);             // μᵀ Σ^{-1} μ
  return {logdet, quad};
}

//------------------------------------------------------------------------------
// bdiag_rcpp
//------------------------------------------------------------------------------
arma::mat bdiag_rcpp(const Rcpp::List& mats) {
  int total_rows = 0, total_cols = 0;
  int n = mats.size();

  for (int i = 0; i < n; i++) {
    arma::mat current = as<arma::mat>(mats[i]);
    total_rows += current.n_rows;
    total_cols += current.n_cols;
  }

  arma::mat ans(total_rows, total_cols, arma::fill::zeros);

  int i1 = 0, j1 = 0;
  for (int i = 0; i < n; i++) {
    arma::mat current = as<arma::mat>(mats[i]);
    int nrows = current.n_rows, ncols = current.n_cols;
    ans.submat(i1, j1, i1 + nrows - 1, j1 + ncols - 1) = current;
    i1 += nrows; j1 += ncols;
  }
  return ans;
}

//------------------------------------------------------------------------------
// sampleGammaCpp
//------------------------------------------------------------------------------
// [[Rcpp::export]]
Rcpp::IntegerVector sampleGammaCpp2(Rcpp::IntegerVector Gammacurr,
                                   const arma::mat& totalW,
                                   const arma::mat& Xmat,
                                   const arma::mat& W,
                                   const arma::mat& blockPriorInv,
                                   const arma::vec& d,
                                   const arma::mat& Y,
                                   const arma::mat& ZTilde,
                                   Rcpp::CharacterVector responseType,
                                   const arma::vec& g,
                                   int num_basis,
                                   int p,
                                   int k,
                                   int start,
                                   int end,
                                   double omegaM,
                                   const arma::vec& Rh,
                                   const arma::vec& zeta) {

  int totalW_cols = totalW.n_cols;
  int Xmat_cols   = Xmat.n_cols;

  // always included
  std::vector<int> alwaysIncluded;
  for (int j = 0; j < p; j++) {
    alwaysIncluded.push_back(p + (num_basis - 1) * j);
  }

  // cols
  for (int j = p; j < totalW_cols; j++) {
    if (std::find(alwaysIncluded.begin(), alwaysIncluded.end(), j) != alwaysIncluded.end())
      continue;

    // --- PART 1: idx, idx0  ---
    std::vector<int> idx1;
    for (int i = start; i <= start + j - 1; i++) {
      if (Gammacurr[i] > 0) idx1.push_back(i - start);
    }

    std::vector<int> idx, idx0;
    if (j == totalW_cols - 1) {
      idx = idx1; idx.push_back(j);
      idx0 = idx1;
    } else {
      std::vector<int> idx2;
      for (int i = start + j + 1; i <= end; i++) {
        if (Gammacurr[i] > 0) idx2.push_back(i - start);
      }
      idx  = idx1; idx.push_back(j); idx.insert(idx.end(), idx2.begin(), idx2.end());
      idx0 = idx1; idx0.insert(idx0.end(), idx2.begin(), idx2.end());
    }

    // --- PART 2: prgamma1 ---
    int jp   = (j + 1) - Xmat_cols; // (R: j - ncol(Xmat)) 1-based
    int numj = (jp % (num_basis - 1) == 0) ? (jp / (num_basis - 1))
      : (jp / (num_basis - 1) + 1);

    std::vector<int> gamsub;
    for (int i = start + p; i <= end; i++) gamsub.push_back(Gammacurr[i]);

    int block_start = (num_basis - 1) * (numj - 1);
    int block_end   = (num_basis - 1) * numj - 1;

    std::vector<int> gamsub_block;
    for (int i = block_start; i <= block_end; i++) {
      if (i == j - p) continue;
      gamsub_block.push_back(gamsub[i]);
    }
    int gamsum = std::accumulate(gamsub_block.begin(), gamsub_block.end(), 0);
    double prgamma1 = ((1 - omegaM) * (gamsum + 1.0)) /
      ((num_basis - 1 - gamsum) + (1 - omegaM) * (gamsum + 1.0));

    // --- PART 3: W_11, W_00 ---
    arma::uvec idx_u(idx.size());
    for (unsigned int t = 0; t < idx.size(); t++) idx_u[t] = idx[t];
    arma::mat W_11 = totalW.cols(idx_u);

    arma::uvec idx0_u(idx0.size());
    for (unsigned int t = 0; t < idx0.size(); t++) idx0_u[t] = idx0[t];
    arma::mat W_00 = totalW.cols(idx0_u);

    std::vector<int> new_idx, new_idx0;
    for (unsigned int t = p; t < idx.size();  t++) new_idx.push_back(idx[t]  - p);
    for (unsigned int t = p; t < idx0.size(); t++) new_idx0.push_back(idx0[t] - p);

    std::vector<int> subidx, subidx0;
    for (auto v : new_idx)  if (v  <= (num_basis - 1) - 1) subidx.push_back(v);
    for (auto v : new_idx0) if (v <= (num_basis - 1) - 1) subidx0.push_back(v);

    // --- PART 4: bdiag  ---
    int start_g = (k - 1) * p;
    int end_g   = k * p - 1;
    arma::vec sub_g = g.subvec(start_g, end_g);

    arma::mat bdiagmat1, bdiagmat0;
    if (!subidx.empty()) {
      arma::uvec subidx_u(subidx.size());
      for (unsigned int t = 0; t < subidx.size(); t++) subidx_u[t] = subidx[t];
      arma::mat W_sub = W.cols(subidx_u);
      bdiagmat1 = (1.0 / sub_g[0]) * (W_sub.t() * W_sub);
    }
    if (!subidx0.empty()) {
      arma::uvec subidx0_u(subidx0.size());
      for (unsigned int t = 0; t < subidx0.size(); t++) subidx0_u[t] = subidx0[t];
      arma::mat W_sub0 = W.cols(subidx0_u);
      bdiagmat0 = (1.0 / sub_g[0]) * (W_sub0.t() * W_sub0);
    }

    for (int q = 2; q <= p; q++) {
      int lower_bound = (q - 1) * (num_basis - 1);
      int upper_bound =  q      * (num_basis - 1) - 1;

      std::vector<int> subidx_q, subidx0_q;
      for (auto v : new_idx)  if (v  >= lower_bound && v  <= upper_bound) subidx_q.push_back(v);
      for (auto v : new_idx0) if (v >= lower_bound && v <= upper_bound)  subidx0_q.push_back(v);

      if (!subidx_q.empty()) {
        arma::uvec subidx_q_u(subidx_q.size());
        for (unsigned int t = 0; t < subidx_q.size(); t++) subidx_q_u[t] = subidx_q[t];
        arma::mat W_sub_q = W.cols(subidx_q_u);
        arma::mat block1  = (1.0 / sub_g[q - 1]) * (W_sub_q.t() * W_sub_q);
        if (bdiagmat1.n_elem == 0) bdiagmat1 = block1;
        else {
          Rcpp::List temp; temp.push_back(bdiagmat1); temp.push_back(block1);
          bdiagmat1 = bdiag_rcpp(temp);
        }
      }
      if (!subidx0_q.empty()) {
        arma::uvec subidx0_q_u(subidx0_q.size());
        for (unsigned int t = 0; t < subidx0_q.size(); t++) subidx0_q_u[t] = subidx0_q[t];
        arma::mat W_sub0_q = W.cols(subidx0_q_u);
        arma::mat block0   = (1.0 / sub_g[q - 1]) * (W_sub0_q.t() * W_sub0_q);
        if (bdiagmat0.n_elem == 0) bdiagmat0 = block0;
        else {
          Rcpp::List temp; temp.push_back(bdiagmat0); temp.push_back(block0);
          bdiagmat0 = bdiag_rcpp(temp);
        }
      }
    }

    { Rcpp::List temp; temp.push_back(blockPriorInv); temp.push_back(bdiagmat1); bdiagmat1 = bdiag_rcpp(temp); }
    { Rcpp::List temp; temp.push_back(blockPriorInv); temp.push_back(bdiagmat0); bdiagmat0 = bdiag_rcpp(temp); }

    // --- PART 5: mu, Sigma (Gaussian + Binary/Ordinal Probit ) ---
    arma::vec mu1, mu0;
    arma::mat Sigma1, Sigma0;

    const bool isGaussian = (as<std::string>(responseType[k - 1]) == "Gaussian");

    double sigma_scale;
    arma::vec term;

    if (isGaussian) {
      // Gaussian: term = (1/d_k)*zeta + (1/d_k^2)*(Rh .* Y_k)
      double d_k = d[k - 1];
      sigma_scale = 1.0 / (d_k * d_k);
      term = (1.0 / d_k) * zeta + sigma_scale * (Rh % Y.col(k - 1));
    } else {
      // Binary/Ordinal Probit(DA): term = zeta + (Rh .* ZTilde_k), scale = 1
      sigma_scale = 1.0;
      term = zeta + (Rh % ZTilde.col(k - 1));
    }

    mu1 = W_11.t() * term;
    mu0 = W_00.t() * term;


    arma::mat W11w = W_11; W11w.each_col() %= Rh;
    arma::mat W00w = W_00; W00w.each_col() %= Rh;

    Sigma1 = sigma_scale * (W_11.t() * W11w) + bdiagmat1;
    Sigma0 = sigma_scale * (W_00.t() * W00w) + bdiagmat0;


    bdiagmat1 = 0.5 * (bdiagmat1 + bdiagmat1.t());
    bdiagmat0 = 0.5 * (bdiagmat0 + bdiagmat0.t());
    double logdet_bdiag1 = arma::log_det_sympd(bdiagmat1);
    double logdet_bdiag0 = arma::log_det_sympd(bdiagmat0);


    CholStats s1 = chol_stats(Sigma1, mu1);
    CholStats s0 = chol_stats(Sigma0, mu0);

    double logSr1 = -0.5 * s1.logdet + 0.5 * s1.quad + 0.5 * logdet_bdiag1 + std::log(prgamma1);
    double logSr0 = -0.5 * s0.logdet + 0.5 * s0.quad + 0.5 * logdet_bdiag0 + std::log(1.0 - prgamma1);

    // stable (log-sum-exp)
    double m = std::max(logSr0, logSr1);
    double prob_gamma = std::exp(logSr1 - m) / (std::exp(logSr0 - m) + std::exp(logSr1 - m));
    int gamma_sample = (R::runif(0.0, 1.0) < prob_gamma) ? 1 : 0;

    // update
    Gammacurr[start + j] = gamma_sample;
  }

  return Gammacurr;
}
