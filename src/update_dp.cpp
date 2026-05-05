// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <cmath>
using namespace Rcpp;
using namespace arma;

// helper: 다변량 정규분포 표본 추출 (Cholesky 이용)
// mean: 벡터, sigma: 공분산행렬 (대칭, 양의정부호)
// [[Rcpp::export]]
arma::vec rmvnorm_sample(const arma::vec & mean, const arma::mat & sigma) {
  int n = mean.n_elem;
  arma::vec z = arma::randn(n);
  // lower-Cholesky 분해 (sigma = L * L.t())
  arma::mat L = arma::chol(sigma, "lower");
  return mean + L * z;
}

// 정석 방식의 로그 행렬식 계산 함수
// A: 입력 행렬
// reg: 정규화를 위한 상수 (기본값 1e-6)

// [[Rcpp::export]]
Rcpp::List updateDPParameters(
    int H,                // 전체 범주 수 (DP 파라미터 업데이트에서는 1:(H-1)를 업데이트)
    int p2,               // 각 범주에 해당하는 파라미터 수 (block 길이)
    Rcpp::IntegerVector Z_ind, // 각 관측에 대한 범주 지표 (R의 1-indexed)
    arma::vec Gammacurr_x,     // DP 스파스 여부 벡터 (길이 H*p2; 각 블록마다 업데이트)
    arma::mat W_x,             // 디자인행렬 (관측 x p2)
    arma::mat Z_star,          // 잠재변수 행렬 (관측 x H)
    arma::vec gh,              // scale 파라미터 벡터 (길이 H)
    double omegaC,             // 상수 (예, 스파스 priors 관련)
    arma::ivec num_pi2,        // (출력용) 각 범주별 ind_x의 개수 (길이 H)
    arma::vec alpha_h,         // 범주별 intercept (길이 H)
    arma::vec BCurr_x,         // DP 회귀계수 벡터 (길이 H*p2; 각 블록별 업데이트)
    arma::mat WB_x             // 선형예측자 행렬 (관측 x H)
) {
  // n: 관측 수, m: W_x의 열 수 (일반적으로 m == p2)
  int n = W_x.n_rows;
  int m = W_x.n_cols;
  
  // h = 1,2,...,H-1 (R에서는 1부터 H-1; C++에서는 0부터 H-2)
  for (int h = 0; h < H - 1; h++) {
    
    //Rcpp::Rcout << h << std::endl;
    // h_index: R의 h (1-indexed)
    int h_index = h + 1;
    // 블록 내 인덱스: R에서는 start_x <- (h-1)*p2 + 1, end_x <- h*p2  
    // C++에서는 0-indexed로:
    int start_x = h * p2;
    int end_x = (h + 1) * p2 - 1;
    
    //Rcpp::Rcout << end_x << std::endl;
    // ---------------------------
    // 1. 인덱스 ind_x: Z_ind >= h_index 인 관측의 행 번호 (R에서는 which(Z_ind >= h))
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
    
    //Rcpp::Rcout << h << std::endl;
    // ---------------------------
    //Rcpp::Rcout << ind_x.n_elem << " " << h<< std::endl;
    // 2. [gamma_h 업데이트]
    if (ind_x.n_elem <= 1) {
      
      // 관측 수가 0 또는 1이면 해당 블록은 고정: Gammacurr_x[start_x] = 1, 나머지 0
      Gammacurr_x[start_x] = 1;
      for (int k = start_x + 1; k <= end_x; k++) {
        Gammacurr_x[k] = 0;
      }
      
      //Rcpp::Rcout << h << std::endl;
    } else {
      
      //Rcpp::Rcout << h << std::endl;
      // W_x의 열 수(블록내 열수 m; 일반적으로 m == p2)에서 j = 2,...,m (R: j=2:m)
      // → C++에서는 j = 1,...,m-1 (0-indexed: j+1가 R의 j)
      for (int j = 1; j < m; j++) {
        //Rcpp::Rcout << h << std::endl;
        // idx1: 블록 내에서 [start_x, start_x+j-1] (즉, 첫 j개 중)에서 Gammacurr_x > 0인 상대적 인덱스
        std::vector<int> idx1;
        for (int k = start_x; k < start_x + j; k++) {
          if (Gammacurr_x[k] > 0)
            idx1.push_back(k - start_x);  // 블록내 0-indexed
        }
        std::vector<int> idx, idx0;
        if (j == m - 1) {
          idx = idx1;
          idx.push_back(j);  // 현재 j (이미 블록내 0-indexed; R의 j+1)
          idx0 = idx1;
        } else {
          std::vector<int> idx2;
          for (int k = start_x + j + 1; k <= end_x; k++) {
            if (Gammacurr_x[k] > 0)
              idx2.push_back(k - start_x);
          }
          idx = idx1;
          idx.push_back(j);  // 가운데 j
          idx.insert(idx.end(), idx2.begin(), idx2.end());
          idx0 = idx1;
          idx0.insert(idx0.end(), idx2.begin(), idx2.end());
        }
        // (q_x = idx.size() 는 여기서 사용하지 않음)
        
        //Rcpp::Rcout << j << std::endl;
        // — 추출: 현재 블록에 해당하는 W_x의 부분행렬
        arma::mat W_block = W_x;
        //Rcpp::Rcout << j << std::endl;
        // idx, idx0를 uvec로 변환 (이들은 블록내 열 번호, 0-indexed)
        arma::uvec idx_u = conv_to<uvec>::from(idx);
        arma::uvec idx0_u = conv_to<uvec>::from(idx0);
        // W_1 = W_block(ind_x, idx); W_0 = W_block(ind_x, idx0)
        arma::mat W_1 = W_block.submat(ind_x, idx_u);
        arma::mat W_0 = W_block.submat(ind_x, idx0_u);
        // W_11 = [1 | W_1] (열 벡터 1을 앞에 붙임)
        arma::mat ones_W1 = ones<mat>(W_1.n_rows, 1);
        arma::mat W_11 = join_horiz(ones_W1, W_1);
        arma::mat ones_W0 = ones<mat>(W_0.n_rows, 1);
        arma::mat W_00 = join_horiz(ones_W0, W_0);
        
        //Rcpp::Rcout << j << std::endl;
        // Z_star의 해당 행(ind_x)와 h번째 열 (C++에서는 h가 0-indexed)
        //arma::vec Zs_sub = Z_star(ind_x, h);
        // ind_x는 arma::uvec 타입의 행 인덱스라고 가정합니다.
        // h는 0-indexed인 열 번호입니다.
        arma::uvec col_index(1);
        col_index(0) = h;
        // 하위행렬(submatrix) 추출 (하나의 열을 포함한 matrix)
        arma::mat Z_sub_mat = Z_star.submat(ind_x, col_index);
        // matrix의 첫번째(유일한) 열을 vector로 변환
        arma::vec Zs_sub = Z_sub_mat.col(0);
        
        arma::vec mu1 = W_11.t() * Zs_sub;
        arma::vec mu0 = W_00.t() * Zs_sub;
        
        //Rcpp::Rcout <<" i" << std::endl;
        
        // 만약 W_1 또는 W_0의 (자기곱) 행렬식이 0이면 해당 gamma를 0으로 설정
        double det1 = arma::det(W_1.t() * W_1);
        double det0 = arma::det(W_0.t() * W_0);
        if (det1 == 0 || det0 == 0) {
          Gammacurr_x[start_x + j] = 0;
          //Rcpp::Rcout << start_x + j << std::endl;
        } else {
          // --- block diagonal 행렬 구성 ---
          // (1) idx에 대해: A1 = (1/gh[h]) * t(W_block.cols(idx)) * W_block.cols(idx)
          arma::mat W_idx = W_block.cols(idx_u);
          arma::mat A1 = (1.0 / gh[h]) * (W_idx.t() * W_idx);
          // block1 = bdiag(1, A1)
          int rA1 = A1.n_rows;
          arma::mat block1 = zeros<mat>(1 + rA1, 1 + rA1);
          block1(0, 0) = 1;
          block1.submat(1, 1, rA1, rA1) = A1;
          
          // (2) idx0에 대해
          arma::mat W_idx0 = W_block.cols(idx0_u);
          arma::mat A0 = (1.0 / gh[h]) * (W_idx0.t() * W_idx0);
          int rA0 = A0.n_rows;
          arma::mat block0 = zeros<mat>(1 + rA0, 1 + rA0);
          block0(0, 0) = 1;
          block0.submat(1, 1, rA0, rA0) = A0;
          
          // inv_WTW1 = solve( t(W_11)%*%W_11 + block1 )
          arma::mat M1 = W_11.t() * W_11 + block1;
          arma::mat inv_WTW1;
          bool ok1 = inv(inv_WTW1, M1);
          
          inv_WTW1 = 0.5 * (inv_WTW1 + inv_WTW1.t());
          
          // ok1 =  0.5 * (ok1 + ok1.t());
          // inv_WTW0 = solve( t(W_00)%*%W_00 + block0 )
          arma::mat M0 = W_00.t() * W_00 + block0;
          arma::mat inv_WTW0;
          bool ok0 = inv(inv_WTW0, M0);
          // ok0 =  0.5 * (ok0 + ok0.t());
          inv_WTW0 = 0.5 * (inv_WTW0 + inv_WTW0.t());
          if (!ok1 || !ok0) {
            //Rcpp::Rcout <<" i" << std::endl;
            Gammacurr_x[start_x + j] = 0;
          } else {
            double muinvmu1 = as_scalar(mu1.t() * inv_WTW1 * mu1);
            double muinvmu0 = as_scalar(mu0.t() * inv_WTW0 * mu0);
            
            //Rcpp::Rcout << m << std::endl;
            // gamsum_x = sum( Gammacurr_x[start_x:end_x] except the (j)th element
            double gamsum_x = 0;
            for (int k = start_x; k <= end_x; k++) {
              if (k == start_x + j) continue;
              gamsum_x += Gammacurr_x[k];
            }
            double prgamma1_x = (1 - omegaC) * (gamsum_x + 1) /
              ((p2 - gamsum_x) + (1 - omegaC) * (gamsum_x + 1));
            
            // log determinant 계산
            arma::mat mat_det1 = (1.0 / gh[h]) * (W_idx.t() * W_idx);
            mat_det1 = 0.5*(mat_det1 + mat_det1.t());
            //mat_det1 = 0.5 * (mat_det1 + mat_det1.t());
            double logdet_mat1 = log_det_sympd(mat_det1);
            double logdet_invWTW1 = log_det_sympd(inv_WTW1);
            double logSr1 = 0.5 * logdet_mat1 + 0.5 * muinvmu1 + 0.5 * logdet_invWTW1 + std::log(prgamma1_x);
            
            arma::mat mat_det0 = (1.0 / gh[h]) * (W_idx0.t() * W_idx0);
            mat_det0 = 0.5 * (mat_det0 + mat_det0.t());
            //inv_WTW0 = 0.5 * (inv_WTW0 + inv_WTW0.t());
            double logdet_mat0 = log_det_sympd(mat_det0);
            double logdet_invWTW0 = log_det_sympd(inv_WTW0);
            double logSr0 = 0.5 * logdet_mat0 + 0.5 * muinvmu0 + 0.5 * logdet_invWTW0 + std::log(1 - prgamma1_x);
            
            double prob_gamma = 1.0 / (1.0 + std::exp(logSr0 - logSr1));
            // rbinom(1,1,prob=prob_gamma): 간단히 uniform과 비교
            double u = R::runif(0.0, 1.0);
            int gamma_sample = (u < prob_gamma) ? 1 : 0;
            Gammacurr_x[start_x + j] = gamma_sample;
          }
        }
      } // end for j (gamma update)
    } // end else (ind_x.n_elem > 1)
    
    //Rcpp::Rcout << h << std::endl;
    //Rcpp::Rcout <<" i" << std::endl;
    
    // ---------------------------
    // 3. [alpha_h 및 beta_h 업데이트]
    // ind_x (행 인덱스)는 위에서 구한 것을 그대로 사용.
    num_pi2[h] = ind_x.n_elem;  // (R: num_pi2[h] <- length(ind_x))
    
    //Rcpp::Rcout << num_pi2[h] << std::endl;
    // gammaind_x: 해당 블록 내에서 Gammacurr_x > 0인 상대 인덱스 (0-indexed)
    std::vector<int> gammaind_vec;
    for (int k = start_x; k <= end_x; k++) {
      if (Gammacurr_x[k] > 0)
        gammaind_vec.push_back(k - start_x);
    }
    // if (length(gammaind_x)>0 & length(ind_x)>0)
    if (gammaind_vec.size() > 0 && ind_x.n_elem > 0) {
      // setone_x: 글로벌 인덱스 (in BCurr_x) where Gammacurr_x in block > 0
      // setzero_x: where == 0
      std::vector<int> setone_vec, setzero_vec;
      for (int k = start_x; k <= end_x; k++) {
        if (Gammacurr_x[k] > 0)
          setone_vec.push_back(k);
        else
          setzero_vec.push_back(k);
      }
      
      //Rcpp::Rcout <<" i" << std::endl;
      
      // Wh = cbind( rep(1, length(ind_x)), W_x[ind_x, gammaind_x] )
      arma::mat W_block = W_x;
      
      // Rcpp::Rcout << W_x.nco << std::endl;
      
      arma::uvec gammaind = conv_to<uvec>::from(gammaind_vec);
      //arma::uvec idx_u = conv_to<uvec>::from(idx);
      arma::mat W_idx = W_block.cols(gammaind);
      arma::mat W_gamma = W_block.submat(ind_x, gammaind);
      arma::mat ones_Wh = ones<mat>(ind_x.n_elem, 1);
      arma::mat Wh = join_horiz(ones_Wh, W_gamma);
      
      //Rcpp::Rcout << W_gamma.n_cols << ' ' << h << std::endl;
      
      //Rcpp::Rcout <<" i" << std::endl;
      
      // block diagonal: bdiag(1, 1/gh[h]*t(W_x[, gammaind])%*%W_x[, gammaind])
      arma::mat A_alpha = (1.0 / gh[h]) * (W_idx.t() * W_idx);
      int rA_alpha = A_alpha.n_rows;
      arma::mat block_alpha = zeros<mat>(1 + rA_alpha, 1 + rA_alpha);
      block_alpha(0, 0) = 1;
      block_alpha.submat(1, 1, rA_alpha, rA_alpha) = A_alpha;
      
      arma::mat M_alpha = Wh.t() * Wh + block_alpha;
      arma::mat inv_WtW;
      bool ok_alpha = inv(inv_WtW, M_alpha);
      //Rcpp::Rcout << ok_alpha << std::endl;
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
        
        //Rcpp::Rcout <<" i" << std::endl;
        // ind_x는 arma::uvec 타입의 행 인덱스라고 가정합니다.
        // h는 0-indexed인 열 번호입니다.
        arma::uvec col_index(1);
        col_index(0) = h;
        // 하위행렬(submatrix) 추출 (하나의 열을 포함한 matrix)
        arma::mat Z_sub_mat = Z_star.submat(ind_x, col_index);
        
        //arma::vec col3 = Z_sub_mat.col(1);
        
        // matrix의 첫번째(유일한) 열을 vector로 변환
        arma::vec Zs_sub = Z_sub_mat.col(0);
        
        arma::vec muh = inv_WtW * (Wh.t() * Zs_sub);
        arma::mat Vh = 0.5 * (inv_WtW + inv_WtW.t()) ;
        // new_coef ~ N(muh, Vh); 차원: 1 + (number of gammaind)
        
        //Rcpp::Rcout << h << std::endl;
        
        arma::vec new_coef = rmvnorm_sample(muh, Vh);
        
        //Rcpp::Rcout <<h << std::endl;
        alpha_h[h] = new_coef(0); 
        // setone_x에 해당하는 BCurr_x는 new_coef의 나머지 요소로 업데이트
        int coef_idx = 1;
        for (size_t i = 0; i < setone_vec.size(); i++) {
          BCurr_x[ setone_vec[i] ] = new_coef(coef_idx);
          // Rcpp::Rcout << coef_idx << std::endl;
          coef_idx++;
        }
        // setzero_x는 0으로
        for (size_t i = 0; i < setzero_vec.size(); i++) {
          BCurr_x[ setzero_vec[i] ] = 0;
        }
        // WB_x[, h] = W_x %*% BCurr_x[start_x:end_x]
        arma::vec BCurr_seg = BCurr_x.subvec(start_x, end_x);
        WB_x.col(h) = W_x * BCurr_seg;
        
        
        // Rcpp::Rcout <<" i" << std::endl;
      }
    } else {
      // 만약 gammaind가 없거나 ind_x가 없다면
      std::vector<int> setzero_vec;
      for (int k = start_x; k <= end_x; k++) {
        //Rcpp::Rcout <<" i" << std::endl;
        if (Gammacurr_x[k] == 0)
          setzero_vec.push_back(k);
      }
      
      //Rcpp::Rcout <<" i" << std::endl;
      for (size_t i = 0; i < setzero_vec.size(); i++) {
        BCurr_x[ setzero_vec[i] ] = 0;
      }
      for (size_t i = 0; i < ind_x.n_elem; i++) {
        WB_x(ind_x[i], h) = 0;
      }
      // alpha_h[h] = rnorm(1, mean = 1/sqrt(length(ind_x)+1)*sum(Z_star[ind_x,h]), sd = 1/sqrt(length(ind_x)+1))
      double sum_val = 0;
      for (size_t i = 0; i < ind_x.n_elem; i++) {
        sum_val += Z_star(ind_x[i], h);
      }
      double denom = std::sqrt(ind_x.n_elem + 1.0);
      double mu_alpha = sum_val / denom;
      double sigma_alpha = 1.0 / denom;
      alpha_h[h] = R::rnorm(mu_alpha, sigma_alpha);
      
      //Rcpp::Rcout <<" i" << std::endl;
    }
    
    // ---------------------------
    // 4. [g_h 업데이트]
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
    
    //Rcpp::Rcout << h << std::endl;
  } // end for h
  
  // 업데이트된 파라미터들을 List로 반환
  return Rcpp::List::create(
    Named("Gammacurr_x") = Gammacurr_x,
    Named("gh") = gh,
    Named("alpha_h") = alpha_h,
    Named("BCurr_x") = BCurr_x,
    Named("WB_x") = WB_x,
    Named("num_pi2") = num_pi2
  );
}