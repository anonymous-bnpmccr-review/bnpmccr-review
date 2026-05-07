// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <cmath>
using namespace Rcpp;
using namespace arma;

// Draw from a multivariate normal distribution using Cholesky decomposition.
// mean: mean vector, sigma: covariance matrix.
// [[Rcpp::export]]
arma::vec rmvnorm_sample(const arma::vec & mean, const arma::mat & sigma) {
  int n = mean.n_elem;
  arma::vec z = arma::randn(n);
  arma::mat L = arma::chol(sigma, "lower");
  return mean + L * z;
}

// [[Rcpp::export]]
Rcpp::List updateDPParameters(
    int H,                // total number of categories; updates 1:(H-1)
    int p2,               // number of parameters per category
    Rcpp::IntegerVector Z_ind, // category index for each observation, using R's 1-based indexing
    arma::vec Gammacurr_x,     // DP sparsity indicator vector of length H*p2
    arma::mat W_x,             // design matrix
    arma::mat Z_star,          // latent-variable matrix
    arma::vec gh,              // scale parameter vector
    double omegaC,             // sparsity prior constant
    arma::ivec num_pi2,        // number of selected observations for each category
    arma::vec alpha_h,         // category-specific intercepts
    arma::vec BCurr_x,         // DP regression coefficient vector
    arma::mat WB_x             // linear predictor matrix
) {
  int n = W_x.n_rows;
  int m = W_x.n_cols;

  // In R, h = 1,...,H-1; in C++, h = 0,...,H-2.
  for (int h = 0; h < H - 1; h++) {

    int h_index = h + 1;
    int start_x = h * p2;
    int end_x = (h + 1) * p2 - 1;

    // Observation indices satisfying Z_ind >= h_index.
    std::vector<int> ind_x_vec;
    for (int i = 0; i < Z_ind.size(); i++) {
      if (Z_ind[i] >= h_index)
        ind_x_vec.push_back(i);
    }
    arma::uvec ind_x;
    if (ind_x_vec.size() > 0)
      ind_x = conv_to<uvec>::from(ind_x_vec);
    else
      ind_x = uvec();  // empty

    // Update gamma_h.
    if (ind_x.n_elem <= 1) {

      // If there are zero or one observations, fix the block.
      Gammacurr_x[start_x] = 1;
      for (int k = start_x + 1; k <= end_x; k++) {
        Gammacurr_x[k] = 0;
      }

    } else {

      // In R, j = 2,...,m; in C++, j = 1,...,m-1.
      for (int j = 1; j < m; j++) {
        std::vector<int> idx1;
        for (int k = start_x; k < start_x + j; k++) {
          if (Gammacurr_x[k] > 0)
            idx1.push_back(k - start_x);
        }
        std::vector<int> idx, idx0;
        if (j == m - 1) {
          idx = idx1;
          idx.push_back(j);
          idx0 = idx1;
        } else {
          std::vector<int> idx2;
          for (int k = start_x + j + 1; k <= end_x; k++) {
            if (Gammacurr_x[k] > 0)
              idx2.push_back(k - start_x);
          }
          idx = idx1;
          idx.push_back(j);
          idx.insert(idx.end(), idx2.begin(), idx2.end());
          idx0 = idx1;
          idx0.insert(idx0.end(), idx2.begin(), idx2.end());
        }

        arma::mat W_block = W_x;
        arma::uvec idx_u = conv_to<uvec>::from(idx);
        arma::uvec idx0_u = conv_to<uvec>::from(idx0);
        arma::mat W_1 = W_block.submat(ind_x, idx_u);
        arma::mat W_0 = W_block.submat(ind_x, idx0_u);
        arma::mat ones_W1 = ones<mat>(W_1.n_rows, 1);
        arma::mat W_11 = join_horiz(ones_W1, W_1);
        arma::mat ones_W0 = ones<mat>(W_0.n_rows, 1);
        arma::mat W_00 = join_horiz(ones_W0, W_0);

        arma::uvec col_index(1);
        col_index(0) = h;
        arma::mat Z_sub_mat = Z_star.submat(ind_x, col_index);
        arma::vec Zs_sub = Z_sub_mat.col(0);

        arma::vec mu1 = W_11.t() * Zs_sub;
        arma::vec mu0 = W_00.t() * Zs_sub;

        // If W_1 or W_0 is singular, set the corresponding gamma to zero.
        double det1 = arma::det(W_1.t() * W_1);
        double det0 = arma::det(W_0.t() * W_0);
        if (det1 == 0 || det0 == 0) {
          Gammacurr_x[start_x + j] = 0;
        } else {
          arma::mat W_idx = W_block.cols(idx_u);
          arma::mat A1 = (1.0 / gh[h]) * (W_idx.t() * W_idx);
          int rA1 = A1.n_rows;
          arma::mat block1 = zeros<mat>(1 + rA1, 1 + rA1);
          block1(0, 0) = 1;
          block1.submat(1, 1, rA1, rA1) = A1;

          arma::mat W_idx0 = W_block.cols(idx0_u);
          arma::mat A0 = (1.0 / gh[h]) * (W_idx0.t() * W_idx0);
          int rA0 = A0.n_rows;
          arma::mat block0 = zeros<mat>(1 + rA0, 1 + rA0);
          block0(0, 0) = 1;
          block0.submat(1, 1, rA0, rA0) = A0;

          arma::mat M1 = W_11.t() * W_11 + block1;
          arma::mat inv_WTW1;
          bool ok1 = inv(inv_WTW1, M1);

          inv_WTW1 = 0.5 * (inv_WTW1 + inv_WTW1.t());

          arma::mat M0 = W_00.t() * W_00 + block0;
          arma::mat inv_WTW0;
          bool ok0 = inv(inv_WTW0, M0);
          inv_WTW0 = 0.5 * (inv_WTW0 + inv_WTW0.t());
          if (!ok1 || !ok0) {
            Gammacurr_x[start_x + j] = 0;
          } else {
            double muinvmu1 = as_scalar(mu1.t() * inv_WTW1 * mu1);
            double muinvmu0 = as_scalar(mu0.t() * inv_WTW0 * mu0);

            double gamsum_x = 0;
            for (int k = start_x; k <= end_x; k++) {
              if (k == start_x + j) continue;
              gamsum_x += Gammacurr_x[k];
            }
            double prgamma1_x = (1 - omegaC) * (gamsum_x + 1) /
              ((p2 - gamsum_x) + (1 - omegaC) * (gamsum_x + 1));

            arma::mat mat_det1 = (1.0 / gh[h]) * (W_idx.t() * W_idx);
            mat_det1 = 0.5*(mat_det1 + mat_det1.t());
            double logdet_mat1 = log_det_sympd(mat_det1);
            double logdet_invWTW1 = log_det_sympd(inv_WTW1);
            double logSr1 = 0.5 * logdet_mat1 + 0.5 * muinvmu1 + 0.5 * logdet_invWTW1 + std::log(prgamma1_x);

            arma::mat mat_det0 = (1.0 / gh[h]) * (W_idx0.t() * W_idx0);
            mat_det0 = 0.5 * (mat_det0 + mat_det0.t());
            double logdet_mat0 = log_det_sympd(mat_det0);
            double logdet_invWTW0 = log_det_sympd(inv_WTW0);
            double logSr0 = 0.5 * logdet_mat0 + 0.5 * muinvmu0 + 0.5 * logdet_invWTW0 + std::log(1 - prgamma1_x);

            double prob_gamma = 1.0 / (1.0 + std::exp(logSr0 - logSr1));
            double u = R::runif(0.0, 1.0);
            int gamma_sample = (u < prob_gamma) ? 1 : 0;
            Gammacurr_x[start_x + j] = gamma_sample;
          }
        }
      } // end for j (gamma update)
    } // end else (ind_x.n_elem > 1)

    // Update alpha_h and beta_h.
    num_pi2[h] = ind_x.n_elem;

    std::vector<int> gammaind_vec;
    for (int k = start_x; k <= end_x; k++) {
      if (Gammacurr_x[k] > 0)
        gammaind_vec.push_back(k - start_x);
    }
    if (gammaind_vec.size() > 0 && ind_x.n_elem > 0) {
      std::vector<int> setone_vec, setzero_vec;
      for (int k = start_x; k <= end_x; k++) {
        if (Gammacurr_x[k] > 0)
          setone_vec.push_back(k);
        else
          setzero_vec.push_back(k);
      }

      arma::mat W_block = W_x;

      arma::uvec gammaind = conv_to<uvec>::from(gammaind_vec);
      arma::mat W_idx = W_block.cols(gammaind);
      arma::mat W_gamma = W_block.submat(ind_x, gammaind);
      arma::mat ones_Wh = ones<mat>(ind_x.n_elem, 1);
      arma::mat Wh = join_horiz(ones_Wh, W_gamma);

      arma::mat A_alpha = (1.0 / gh[h]) * (W_idx.t() * W_idx);
      int rA_alpha = A_alpha.n_rows;
      arma::mat block_alpha = zeros<mat>(1 + rA_alpha, 1 + rA_alpha);
      block_alpha(0, 0) = 1;
      block_alpha.submat(1, 1, rA_alpha, rA_alpha) = A_alpha;

      arma::mat M_alpha = Wh.t() * Wh + block_alpha;
      arma::mat inv_WtW;
      bool ok_alpha = inv(inv_WtW, M_alpha);
      if (!ok_alpha) {
        alpha_h[h] = 0;
        for (size_t i = 0; i < setone_vec.size(); i++) {
          BCurr_x[ setone_vec[i] ] = 0;
        }
        for (size_t i = 0; i < setzero_vec.size(); i++) {
          BCurr_x[ setzero_vec[i] ] = 0;
        }
        for (size_t i = 0; i < ind_x.n_elem; i++) {
          WB_x(ind_x[i], h) = 0;
        }
      } else {

        arma::uvec col_index(1);
        col_index(0) = h;
        arma::mat Z_sub_mat = Z_star.submat(ind_x, col_index);

        arma::vec Zs_sub = Z_sub_mat.col(0);

        arma::vec muh = inv_WtW * (Wh.t() * Zs_sub);
        arma::mat Vh = 0.5 * (inv_WtW + inv_WtW.t()) ;

        arma::vec new_coef = rmvnorm_sample(muh, Vh);

        alpha_h[h] = new_coef(0);
        int coef_idx = 1;
        for (size_t i = 0; i < setone_vec.size(); i++) {
          BCurr_x[ setone_vec[i] ] = new_coef(coef_idx);
          coef_idx++;
        }
        for (size_t i = 0; i < setzero_vec.size(); i++) {
          BCurr_x[ setzero_vec[i] ] = 0;
        }
        arma::vec BCurr_seg = BCurr_x.subvec(start_x, end_x);
        WB_x.col(h) = W_x * BCurr_seg;
      }
    } else {
      std::vector<int> setzero_vec;
      for (int k = start_x; k <= end_x; k++) {
        if (Gammacurr_x[k] == 0)
          setzero_vec.push_back(k);
      }

      for (size_t i = 0; i < setzero_vec.size(); i++) {
        BCurr_x[ setzero_vec[i] ] = 0;
      }
      for (size_t i = 0; i < ind_x.n_elem; i++) {
        WB_x(ind_x[i], h) = 0;
      }
      double sum_val = 0;
      for (size_t i = 0; i < ind_x.n_elem; i++) {
        sum_val += Z_star(ind_x[i], h);
      }
      double denom = std::sqrt(ind_x.n_elem + 1.0);
      double mu_alpha = sum_val / denom;
      double sigma_alpha = 1.0 / denom;
      alpha_h[h] = R::rnorm(mu_alpha, sigma_alpha);
    }

    // Update g_h.
    if (ind_x.n_elem == 0) {
      gh[h] = 1;
    } else {
      double sumWB = dot(WB_x.col(h), WB_x.col(h));
      double sum_gamma = 0;
      for (int k = start_x; k <= end_x; k++) {
        sum_gamma += Gammacurr_x[k];
      }
      double shape = 0.5 + sum_gamma / 2.0;
      double rate = (n / 2.0) + (sumWB / 2.0);
      double scale = 1.0 / rate; // R::rgamma(shape, scale) returns gamma variate with given scale
      double gamma_sample = R::rgamma(shape, scale);
      gh[h] = 1.0 / gamma_sample;
    }

  } // end for h

  return Rcpp::List::create(
    Named("Gammacurr_x") = Gammacurr_x,
    Named("gh") = gh,
    Named("alpha_h") = alpha_h,
    Named("BCurr_x") = BCurr_x,
    Named("WB_x") = WB_x,
    Named("num_pi2") = num_pi2
  );
}
