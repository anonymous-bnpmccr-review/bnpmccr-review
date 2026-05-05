// 파일명도 ASCII로: update_coeff_mod.cpp

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include "updateZCurrGamma2.h"   // 헤더만!
#include "updateZCurr.h"
// #include "update_zcurr_gamma_mod.cpp"  ← 삭제!
#include <numeric>
#include <vector>
#include <algorithm>

using namespace Rcpp;
using namespace arma;
// helper: 두 행렬을 block-diagonal 형태로 결합하는 함수
arma::mat bdiag(const arma::mat &A, const arma::mat &B) {
  arma::mat C = arma::zeros<arma::mat>(A.n_rows + B.n_rows, A.n_cols + B.n_cols);
  C(arma::span(0, A.n_rows - 1), arma::span(0, A.n_cols - 1)) = A;
  C(arma::span(A.n_rows, A.n_rows + B.n_rows - 1),
    arma::span(A.n_cols, A.n_cols + B.n_cols - 1)) = B;
  return C;
}


// helper: 다변량 정규분포 표본 추출 (Cholesky 분해 이용)
// 주의: arma::chol()는 상삼각행렬을 반환하므로 U.t()를 곱해줍니다.
arma::vec mvnrnd(const arma::vec &mu, const arma::mat &Sigma) {
  int n = mu.n_elem;
  arma::vec Y = arma::randn(n);
  arma::mat U = arma::chol(Sigma); // 상삼각행렬
  return mu + U.t() * Y;
}

// [[Rcpp::export]]
Rcpp::List updateCoefficients2(
    const int k,                         // response index (0-indexed)
    const std::string responseType,      // "Gaussian", "binary", "ordinal", "binomial", "negative binomial", "gamma"
    const arma::vec &d,                  // 벡터 d (길이는 m)
    arma::mat& ZCurr,              // n x m, 현재 잠재 반응값들
    const arma::mat &WB,                 // n x m, working response 또는 offset
    const Rcpp::List &invRCurr,          // 각 관측치에 대한 역행렬(또는 관련 객체) 리스트
    const arma::ivec &Z_ind,             // invRCurr 접근용 인덱스 (0-indexed)
    const arma::vec &Gammacurr,          // 현재 Gamma (변수 선택 관련)
    const int start, const int end,      // Gammacurr의 하위 벡터 범위 (0-indexed)
    const arma::mat &totalW,             // 전체 design 행렬 (회귀계수 업데이트에 사용)
    const arma::mat &Xmat,               // 고정효과 디자인 행렬
    const int p,                       // 각 반응별 회귀계수 개수
    const int num_basis,               // basis 함수 개수
    const arma::vec &g,                // 벡터 g (길이는 적어도 m*p 이상)
    const arma::mat &blockPriorInv,    // 사전 분포의 역행렬 (첫 블록)
    const arma::vec &Rh,               // 가중치 벡터 (길이 n)
    const arma::mat &Y,                // n x m, (Gaussian의 경우) 실제 반응값
    arma::mat var_z,                   // n x m, 이후 업데이트될 분산 행렬 (by value)
    const arma::mat &W,                // 블록 결합에 사용될 design 행렬
    const arma::vec &newGamma,         // MH 단계에서 사용하는 newGamma 벡터
    const arma::mat &Wgammastar,       // MH 제안 행렬
    int countgamma,                    // gamma 응답의 경우 카운터
    const arma::vec& gammaAlpha,       // gamma 응답용 gammaAlpha 벡터
    int countnegbinom,                 // neg binom 응답의 경우 카운터
    const arma::vec& negBinomPar,
    int countbinom,                     // binom 응답의 경우 카운터
    const arma::vec& nTrial,
    const Rcpp::List& RCurr,
    const arma::mat& dpropormat,
    const int s
) {
  const int n = ZCurr.n_rows;
  const int m = ZCurr.n_cols;
  
  // k번째 반응에 대한 업데이트된 latent 변수 (ZTilde)를 저장할 벡터
  arma::vec ZTilde_k;
  arma::mat inv_bdiagmat_OLD;
  arma::mat inv_V_OLD;
  arma::vec muprop_OLD;
  arma::vec ZTilde_k_new;    // gamma 분기에서만 계산하든, 아니든
  arma::vec TUProp(n), TLProp(n);
  arma::mat ZCurr_new;
  
  if(responseType == "Gaussian") {
    // Gaussian인 경우: ZTilde[,k] = d[k] * ZCurr[,k] + WB[,k]
    ZTilde_k = d[k] * ZCurr.col(k) + WB.col(k);
    
    // 다른 반응들과의 상호작용에 의한 zeta 계산
    arma::vec zeta(n, fill::zeros);
    for (int j = 0; j < n; ++j) {
      arma::mat invR = Rcpp::as<arma::mat>(invRCurr[ Z_ind[j] - 1 ]);
      for (int kk = 0; kk < m; ++kk) {
        if (kk == k) continue;
        zeta(j) += ZCurr(j, kk) * invR(k, kk) / d[k];
      }
    }
    
    // Gammacurr의 start~end 구간 중 양수인 원소의 상대 인덱스 선택
    std::vector<int> gammaind;
    for (int i = start; i <= end; ++i) {
      if (Gammacurr(i) > 0)
        gammaind.push_back(i - start);
    }
    arma::uvec gammaind_arma(gammaind.size());
    for (size_t i = 0; i < gammaind.size(); ++i)
      gammaind_arma(i) = gammaind[i];
    arma::mat Wgamma = totalW.cols(gammaind_arma);
    
    // fixed effects 처리를 위해 Xmat에 해당하는 첫 요소들을 제외한 후 p를 빼줌
    //const int ncolX = Xmat.n_cols;
    std::vector<int> gammaind1;
    for (size_t i = p; i < gammaind.size(); ++i)
      gammaind1.push_back(gammaind[i] - p);
    
    // k번째 반응에 해당하는 g의 하위 벡터: 인덱스 [k*p, (k+1)*p - 1]
    arma::vec sub_g = g.subvec(k * p, (k + 1) * p - 1);
    
    // gammaind1 중 (num_basis-1) 미만인 원소 선택
    std::vector<int> subgammind;
    for (size_t i = 0; i < gammaind1.size(); ++i) {
      if (gammaind1[i] <= (num_basis - 1) - 1)
        subgammind.push_back(gammaind1[i]);
    }
    arma::uvec subgammind_arma(subgammind.size());
    for (size_t i = 0; i < subgammind.size(); ++i)
      subgammind_arma(i) = subgammind[i];
    
    // W의 하위 행렬을 이용하여 block-diagonal 행렬 구성
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
    // blockPriorInv와 결합
    bdiagmat = bdiag(blockPriorInv, bdiagmat);
    
    // V 및 그 역행렬 계산
    arma::mat V = (1.0 / (d[k] * d[k])) * (Wgamma.t() * (Wgamma.each_col() % Rh)) + bdiagmat;
    arma::mat inv_V = arma::inv(V);
    inv_V = 0.5 * (inv_V + inv_V.t());
    
    // 제안(mean) 계산 후 다변량 정규분포 표본 추출
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
    
    //const int ncolX = Xmat.n_cols;
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
    //Rcpp:: Rcout << countnegbinom <<"\n";
    if(responseType == "binomial") {
      countbinom++;
      ZTilde_k = var_z.col(k) % ZCurr.col(k) + WB.col(k);
    } else if(responseType == "negative binomial"){
      //Rcpp:: Rcout << "A\n";
      countnegbinom++;
      ZTilde_k = var_z.col(k) % ZCurr.col(k) + WB.col(k);
    } else if(responseType == "gamma") {
      countgamma++; // 카운터 증가
     // double trig_val = R::trigamma(gammaAlpha[countgamma - 1])/2 ;
      ZTilde_k = var_z.col(k) % ZCurr.col(k) + WB.col(k);
    }
    
    arma::vec zeta(n, fill::zeros);
    for (int j = 0; j < n; ++j) {
      arma::mat invR = Rcpp::as<arma::mat>(invRCurr[ Z_ind[j] - 1 ]);
      for (int kk = 0; kk < m; ++kk) {
        if (kk == k) continue;
        zeta(j) += ZCurr(j, kk) * invR(k, kk) / var_z(j, k);
        //Rcpp::Rcout<< var_z(j,k) << std::endl;
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
    //const int ncolX = Xmat.n_cols;
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
   // Rcpp::Rcout << "A" ;
    if(responseType == "gamma") {
      arma::mat WB_new = Wk * Bprop;
      ZCurr_new = updateZCurrGamma2(Y, WB_new, ZCurr, gammaAlpha, countgamma, k + 1);
     // Rcpp::Rcout << "A \n" ;
      ZTilde_k_new = var_z.col(k) % dpropormat.col(k) %  ZCurr_new.col(k) + WB_new;
     // Rcpp::Rcout << "A \n" ;
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
        // 1) Allocate
       
        // 2) Compute the raw truncation limits
        for(int i = 0; i < n; ++i) {
          double eta   = WB_new(i);
          double p     = 1.0 / (1.0 + std::exp(-eta)); 
          
          // P(Y[i,k] <= y)
          double c_up  = R::pbinom( Y(i,k),
                                    nTrial(countbinom - 1),
                                    p,
                                    /*lower.tail=*/1, /*log.p=*/0 );
          TUProp(i)    = R::qnorm(c_up, 0.0, 1.0, /*lower.tail=*/1, /*log.p=*/0);
          
          // P(Y[i,k] <= y-1)
          double c_lo  = R::pbinom( Y(i,k) - 1,
                                    nTrial(countbinom - 1),
                                    p,
                                    1, 0 );
          TLProp(i)    = R::qnorm(c_lo, 0.0, 1.0, 1, 0);
        }
        
        //Rcpp:: Rcout << TLProp.n_rows;
        ZCurr_new = updateZCurr_binom(RCurr, Z_ind, ZCurr, k, TLProp, TUProp);
        // Rcpp::Rcout << "A \n" ;
        ZTilde_k_new = var_z.col(k) % ZCurr_new.col(k) + WB_new;

        
      }else if(responseType == "negative binomial"){ // negative binomial
        // 내부에 이미 int n = ZCurr.n_rows, int k, arma::vec WBprop, arma::uvec negBinomPar, int countnegbinom 등이 정의되어 있다고 가정
        
        // 1) TL/TU 벡터 확보
       
        // 2) truncation limits 계산
        for(int i = 0; i < n; ++i) {
          double eta = WB_new(i);            // linear predictor
          double mu  = std::exp(eta);        // NB 평균 파라미터
          
          // 상한: P(Y[i,k] ≤ y)
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
          
          // 하한: P(Y[i,k] ≤ y-1)
          double c_lo = R::pnbinom_mu(
            Y(i,k) - 1,
            negBinomPar(countnegbinom - 1),
            mu,
            1, 0
          );
          TLProp(i) = R::qnorm(c_lo, 0.0, 1.0, 1, 0);
        }
        //Rcpp:: Rcout << TLProp.n_rows;
        ZCurr_new = updateZCurr_binom(RCurr, Z_ind, ZCurr, k, TLProp, TUProp);
        // Rcpp::Rcout << "A \n" ;
        ZTilde_k_new = var_z.col(k) % ZCurr_new.col(k) + WB_new;
        
        
       
        // 이제 TLProp, TUProp을 rtruncnorm 등에 바로 사용할 수 있습니다.
        // 예: ZCurr(i,k) = rtruncnorm(TLProp(i), TUProp(i), mean_z, sd_z);
        
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
