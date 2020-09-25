#include <RcppArmadillo.h>
#include <truncnorm.h>
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo, RcppDist)]]

/* *********** NOTES *********** 
* X = solve( A, B ); Solve a dense system of linear equations, A*X = B, where X is unknown
* i.e., solve( A, B ) gives A^{-1}B
* 
*/

/* *********** GENERAL/HELPER FUNCTIONS *********** 
* 
* Used for creating matrices and such needed by the samplers below.
* Includes a univariate normal sampler, etc.
* 
* arma::mvnrnd(mu, Sigma, n) samples from multivariate normal.
* rnormArma(double mu, double sig_sq) samples from univariate normal.
* 
*/

// [[Rcpp::export]]
double rnormArma(double mu, double sig_sq) {
  /* A single draw of a univariate normal distribution.
  */ 
  return mu + pow(sig_sq, 0.5) * arma::randn( );
}

//[[Rcpp::export]]
bool is_inf( double x){
  /* 
   * General function for determining whether a value is infinite.
   */
  return traits::is_infinite<REALSXP>(x);
}

/* *********** ERROR SAMPLING FUNCTIONS *********** 
* 
* Used for sampling error terms for A_i = B_i + E_i,
* with A_i and B_i known in Gibbs step and E~N(0,diag(\sigma_1^2,...,\sigma_P^2)).
* 
*/

// [[Rcpp::export]]
double sample_sigsq_p(double a_sig, double b_sig, int N, double RSS) {
  /* 
  * General function for sampling \sigma_p^2, given 1/\sigma_p^2 ~ Ga(a_sig/2,b_sig/2)
  * and you observe N samples with residual sum of squares RSS.
  * Note that arma::randg uses scale parameterization, i.e. 1/rate. 
  */
  arma::vec samp = arma::randg(1, arma::distr_param((a_sig+N)/2,2/(b_sig+RSS)));
  return 1/samp(0);
}

// [[Rcpp::export]]
arma::vec sample_sigsq_x(double a_sig, double b_sig, arma::mat D_min_mu, 
                         std::vector< std::string > type) {
  /* A draw of diag(\sigma_1^2, \sigma_P^2), given 1/\sigma_p^2 ~ Ga(a_sig/2,b_sig/2),
  * D_min_mu is the P x N matrix of data minus its mean (P is no of vars, N is no of obs).
  * 
  */ 
  int P = D_min_mu.n_rows;
  int N = D_min_mu.n_cols;
  arma::vec sigsq_all(P);
  sigsq_all.fill(1.0); // This will be refilled with something else if variable is not binary.
  std::string bin ("binary"); 
  for( int p=0; p<P; p++ ){
    if( bin.compare(type[p]) != 0 ){ // If variable p is NOT binary, sample sigsq.
      double RSS = sum( square( D_min_mu.row(p) ) );
      sigsq_all(p) = sample_sigsq_p(a_sig, b_sig, N, RSS);
    }
  }
  return sigsq_all;
}

// [[Rcpp::export]]
arma::vec sample_sigsq_y(double a_sig, double b_sig, arma::mat D_min_mu, 
                         arma::mat obs_Y, bool homo_var) {
  /* A draw of diag(\sigma_1^2, \sigma_P^2), given 1/\sigma_p^2 ~ Ga(a_sig/2,b_sig/2),
  * D_min_mu is the P x N matrix of data minus its mean (P is no of vars, N is no of obs).
  * obs_Y is the P x N matrix where (p,i) is 1 if Y is observed and 0 else.
  * homo_var says whether or not variance is homoscedastic (T) or hetero (F).
  * 
  */ 
  int P = D_min_mu.n_rows;
  int N = D_min_mu.n_cols;
  arma::vec sigsq_all(P);
  double RSS;
  double N_obs;
  if( homo_var ){
    RSS = 0.0;
    N_obs = 0.0;
    for( int p=0; p<P; p++ ){
      arma::rowvec D_min_mu_p = D_min_mu.row(p);
      arma::rowvec obs_Y_p = obs_Y.row(p);
      for(int i=0; i<N; i++){
        double obs_Y_p_i = obs_Y_p(i);
        if( obs_Y_p_i>0 ){ // If you observe 1 or more response values for chem i with dose d.
          RSS += obs_Y_p_i * pow(D_min_mu_p(i), 2.0);
          N_obs += 1.0;
        }
      }
    }
    sigsq_all.fill( sample_sigsq_p(a_sig, b_sig, N_obs, RSS) );
  } else{
    for( int p=0; p<P; p++ ){
      arma::rowvec D_min_mu_p = D_min_mu.row(p);
      arma::rowvec obs_Y_p = obs_Y.row(p);
      RSS = 0.0;
      N_obs = 0;
      for(int i=0; i<N; i++){
        if( obs_Y_p(i)==1 ){
          RSS += pow(D_min_mu_p(i), 2.0);
          N_obs += 1;
        }
      }
      sigsq_all(p) = sample_sigsq_p(a_sig, b_sig, N_obs, RSS);
    }
  }

  return sigsq_all;
}

// [[Rcpp::export]]
arma::vec sample_sigsq_longy(double a_sig, double b_sig, arma::vec D_min_mu_long, 
                             arma::vec dind_long, bool homo_var, int P) {
  /* A draw of diag(\sigma_1^2, \sigma_P^2), given 1/\sigma_p^2 ~ Ga(a_sig/2,b_sig/2),
  * D_min_mu is the long vector of data minus its mean (each of the N_obs rows are one observation).
  * dind_long is the N_obs vector giving the dose index at each observation of D_min_mu.
  * homo_var says whether or not variance is homoscedastic (T) or hetero (F).
  * 
  */ 
  int N_obs = D_min_mu_long.n_rows;
  arma::vec sigsq_all(P);
  double RSS;
  int dind;
  arma::vec RSS_vec(P);
  arma::vec nobs_vec(P);
  if( homo_var ){
    RSS = 0.0;
    for( int i=0; i<N_obs; i++ ){
      RSS += pow(D_min_mu_long(i), 2.0);
    }
    sigsq_all.fill( sample_sigsq_p(a_sig, b_sig, N_obs, RSS) );
  } else{
    RSS_vec.fill(0.0);
    nobs_vec.fill(0.0);
    for( int i=0; i<N_obs; i++ ){
      dind = dind_long(i);
      RSS_vec(dind) += pow(D_min_mu_long(i), 2.0);
      nobs_vec(dind) += 1.0;
    }
    for( int p=0; p<P; p++ ){
      sigsq_all(p) = sample_sigsq_p(a_sig, b_sig, nobs_vec(p), RSS_vec(p));
    }
  }
  
  return sigsq_all;
}

// [[Rcpp::export]]
arma::mat get_X_min_mu(arma::mat X, arma::mat Theta, arma::mat eta,
                       arma::mat xi, arma::mat nu, arma::vec Zmean) {
  /* Get the value of X minus its mean (given all mean params Gibbs).
  * 
  */
  int N = X.n_cols;
  arma::mat E_X = Theta * eta + xi * nu; // S x N (doesn't include Zmean, this gets added next line)
  for( int i=0; i<N; i++){ E_X.col(i) = E_X.col(i) + Zmean; }
  return X - E_X;
}

// [[Rcpp::export]]
arma::vec get_Y_min_mu_long(arma::vec Y_long, arma::mat Lambda, arma::mat eta,
                            arma::ivec IDs_long, arma::ivec dind_long, arma::vec Ymean) {
  /* Get the value of Y minus its mean (given all mean params Gibbs).
   * 
   */ 
  int N_long = Y_long.n_rows;
  arma::mat mu = Lambda*eta; // D times N
  arma::vec Y_min_mu_long(N_long);
  for(int i=0; i<N_long; i++){ Y_min_mu_long(i) = Y_long(i) - mu( dind_long(i), IDs_long(i) ) - Ymean( dind_long(i) ); }
  return Y_min_mu_long;
}

// [[Rcpp::export]]
arma::mat get_Y_min_mu(arma::mat Y, arma::mat Lambda, arma::mat eta) {
  /* Get the value of Y minus its mean (given all mean params Gibbs).
  * 
  */ 
  return Y - Lambda*eta;
}

// [[Rcpp::export]]
arma::vec get_tau(arma::vec delta) {
  /* Get tau = cumulative product of delta.
  * 
  */ 
  return arma::cumprod(delta);
}

/* *********** SCORE SAMPLING FUNCTIONS *********** 
 * 
 * Used for sampling the joint score components \eta = [\eta_1, ..., \eta_N]
 * and X-specific score components \nu = [\nu_1, ..., \nu_N].
 * 
 */

// [[Rcpp::export]]
arma::mat sample_eta_all(arma::mat Y, arma::mat X, arma::mat xi, arma::mat nu, arma::mat Lambda, 
                         arma::mat Theta, arma::vec sigsq_y, arma::vec sigsq_x, arma::mat obs_Y) {
  /* A draw of the K x N matrix [\eta_1,...,\eta_N], the joint scores associated with both X and Y.
  * Here D = [Y; X - \xi \nu] is (D+S) x N vector (cols are obs-specific),
  * Omega = [\Lambda; \Theta] is (D+S) x K matrix (common across obs),
  * and sigsq_D_inv is (D+S) x 1 vector of inv diag entries of independent noise cov matrix for D.
  * 
  */ 
  arma::mat Omega = join_cols(Lambda, Theta); // (D + S) x K
  arma::mat Omega_t = Omega.t(); // K x (D + S)
  arma::mat D = join_cols(Y, X - xi*nu); // (D + S) x N
  arma::mat D_t = D.t(); // N x (D + S)
  int num_doses = Y.n_rows; // Equivalent to D
  int num_feats = X.n_rows; // Equivalent to S
  int P = D.n_rows; // P = D + S
  int K = Omega.n_cols;
  int N = D.n_cols;
  arma::mat eta_all(K, N);
  arma::mat eta_all_t(N, K);
  // Set up what will be the covariance matrices used to sample \eta_i.
  arma::mat SigInv_Omeg(P, K);
  arma::mat SigStar(K, K);
  for( int i=0; i<N; i++ ){
    arma::vec obs_Y_i = obs_Y.col(i);
    int P_obs = sum(obs_Y_i>0) + num_feats; // Number of doses at which Y_i is observed.
    arma::uvec ind_tmp(P_obs);
    int d_tmp = 0;
    arma::vec sigsq_y_tmp = sigsq_y;
    for( int d=0; d<num_doses; d++ ){ // Not all Y are observed.
      if( obs_Y_i(d)>0 ){
        ind_tmp(d_tmp) = d;
        d_tmp += 1;
        sigsq_y_tmp(d) = sigsq_y(d)/obs_Y_i(d);
      }
    }
    for( int ss=num_doses; ss<P; ss++ ){ // All X are observed.
      ind_tmp(d_tmp) = ss;
      d_tmp += 1;
    }
    // Get covariance matrix used when sampling each individual \eta_i.
    arma::vec sigsq_D_inv = join_cols(pow(sigsq_y_tmp, -1.0), pow(sigsq_x, -1.0));
    for( int p=0; p<P; p++ ){ SigInv_Omeg.row(p) = sigsq_D_inv(p) * Omega.row(p); }
    SigStar = Omega_t.cols(ind_tmp) * SigInv_Omeg.rows(ind_tmp);
    for( int k=0; k<K; k++ ){ SigStar(k,k) = SigStar(k,k) + 1; }
    SigStar = SigStar.i();
    // Get (common) mean multiplication matrix used when sampling each individual \eta_i.
    arma::mat mu_mult(K, P_obs);
    mu_mult = SigStar * SigInv_Omeg.rows(ind_tmp).t(); 
    arma::vec D_t_i = D_t.row(i).t();
    // Sample \eta_i for i=1,...,N
    eta_all_t.row(i) = arma::mvnrnd(mu_mult * D_t_i.rows(ind_tmp), SigStar, 1).t();
  }
  eta_all = eta_all_t.t();
  return eta_all; 
}

// [[Rcpp::export]]
arma::mat sample_nu_all(arma::mat X, arma::mat xi, arma::mat eta,
                        arma::mat Theta, arma::vec sigsq_x) {
  /* A draw of the J x N matrix [\nu_1,...,\nu_N], the scores associated with just X.
  * Here D = X - \Theta \eta is S x N vector (cols are obs-specific),
  * and sigsq_x_inv is S x 1 vector of inv diag entries of independent noise cov matrix for X.
  * 
  */ 
  arma::mat xi_t = xi.t();
  arma::mat D = X - Theta*eta;
  arma::mat D_t = D.t(); // N x S
  arma::vec sigsq_x_inv = pow(sigsq_x, -1.0);
  int S = D.n_rows;
  int J = xi.n_cols;
  int N = D.n_cols;
  arma::mat nu_all(J, N);
  arma::mat nu_all_t(N, J);
  // Get (common) covariance matrix used when sampling each individual \nu_i.
  arma::mat SigInv_xi(S, J);
  arma::mat SigStar(J, J);
  for( int s=0; s<S; s++ ){ SigInv_xi.row(s) = sigsq_x_inv(s) * xi.row(s); }
  SigStar = xi_t * SigInv_xi;
  for( int j=0; j<J; j++ ){ SigStar(j,j) = SigStar(j,j) + 1; }
  SigStar = SigStar.i();
  // Get (common) mean multiplication matrix used when sampling each individual \nu_i.
  arma::mat mu_mult(J, S);
  mu_mult = SigStar * SigInv_xi.t(); 
  // Sample \nu_i for i=1,...,N
  for( int i=0; i<N; i++ ){
    nu_all_t.row(i) = arma::mvnrnd(mu_mult * D_t.row(i).t(), SigStar, 1).t();
  }
  nu_all = nu_all_t.t();
  return nu_all; 
}

/* *********** XI SAMPLING FUNCTIONS *********** 
 * 
 * Used for sampling the factor loading matrix \xi
 * and associated hyper-parameters phi_xi, delta_xi (and thus tau_xi).
 * 
 */

// [[Rcpp::export]]
arma::mat sample_xi(arma::mat X, arma::mat nu, arma::mat eta, arma::mat Theta, 
                    arma::vec sigsq_x, arma::mat phi_xi, arma::vec tau_xi) {
  /* A draw of the S x J matrix \xi, the loadings associated with just X.
  * Here D = X - \Theta \eta is S x N vector (cols are obs-specific),
  * and sigsq_x_inv is S x 1 vector of inv diag entries of independent noise cov matrix for X.
  * 
  */ 
  arma::mat D = X - Theta*eta;
  arma::vec sigsq_x_inv = pow(sigsq_x, -1.0);
  int S = D.n_rows;
  int J = nu.n_rows;
  arma::mat xi(S,J);
  arma::mat nu_nuT = nu * nu.t();
  arma::mat SigInv(J,J);
  for( int s=0; s<S; s++ ){
    SigInv = sigsq_x_inv(s) * nu_nuT;
    for( int j=0; j<J; j++ ){
      SigInv(j,j) = SigInv(j,j) + phi_xi(s,j)*tau_xi(j);
    }
    SigInv = SigInv.i();
    xi.row(s) = arma::mvnrnd(SigInv*nu*(D.row(s).t())*sigsq_x_inv(s), SigInv, 1).t();
  }
  return xi; 
}

// [[Rcpp::export]]
arma::mat sample_phi_xi(double gamma_xi, arma::vec tau_xi, arma::mat xi) {
  /* A draw of the S x J matrix \phi_{xi}.
  * The prior on each entry \phi_{\xi, s,j} is Ga(gamma_xi/2, gamma_xi/2)
  * 
  */ 
  int S = xi.n_rows;
  int J = xi.n_cols;
  arma::mat phi_xi(S,J);
  for( int s=0; s<S; s++ ){
    for( int j=0; j<J; j++ ){
      arma::vec samp = arma::randg(1, arma::distr_param((gamma_xi+1)/2,
                                                        2/(gamma_xi+pow(xi(s,j),2.0)*tau_xi(j))));
      phi_xi(s,j) = 1/samp(0);
    }
  }
  return phi_xi;
}

// [[Rcpp::export]]
arma::vec sample_delta_xi(double a1, double a2, arma::mat xi, 
                          arma::mat phi_xi, arma::vec delta_xi) {
  /* 
  * Function for sampling delta = [\delta_{1},...,\delta_{J}]', 
  * given \delta_{1} ~ Ga(a1,1) a priori
  * and \delta_{h} ~ Ga(a2,1) h \in {2,...,J} a priori. 
  * xi and phi_xi are S x J, delta_xi is J x 1 (this is old delta sample).
  * Note that arma::randg uses scale parameterization, i.e. 1/rate.
  * 
  */
  int S = xi.n_rows;
  int J = xi.n_cols;
  double shape;
  double rate;
  
  for( int h=0; h<J; h++ ){
  
    // Specify shape and initial rate parameters
    if(h==0){
      shape = a1 + S*J/2;
      rate = 1.0;
    } else{ 
      shape = a2 + S*(J-h)/2; 
      rate = 1.0/delta_xi(h-1); // Evan's modification
    }
    
    // Calculate the components of the Gamma posterior
    double tau_minus = 1.0;
    for( int j=0; j<J; j++ ){
      if( j != h ){
        tau_minus = tau_minus * delta_xi(j);
      }
      if( j >= h ){
        double tmp_sum = 0.0;
        for( int s=0; s<S; s++ ){
          tmp_sum += phi_xi(s,j) * xi(s,j) * xi(s,j);
        }
        rate += 0.5*tau_minus*tmp_sum;
      }
    }
    // Sample a new value for \delta_h
    delta_xi(h) = arma::randg( 1,arma::distr_param(shape,1/rate) )(0);
  
  }
  
  return delta_xi;
}

/* *********** THETA SAMPLING FUNCTIONS *********** 
 * 
 * Used for sampling the factor loading matrix \Theta
 * and associated hyper-parameters betasq_th (global), 
 * gammasq_th (local), and tau_th (column).
 * 
 */

// [[Rcpp::export]]
arma::mat sample_Theta(arma::mat X, arma::mat nu, arma::mat eta, arma::mat xi, arma::vec sigsq_x, 
                       double betasq_th, arma::mat gammasq_th, arma::vec tau_th) {
  /* A draw of the S x K matrix \Theta, the X loadings associated with the shared scored.
  * Here D = X - \xi \nu is S x N vector (cols are obs-specific),
  * and sigsq_x_inv is S x 1 vector of inv diag entries of independent noise cov matrix for X.
  * Note tau_th is precision term (column-specific shrinkage),
  * and betasq_th (global shrink) and gammasq_th (local shrink) are variance terms.
  * 
  */ 
  arma::mat D = X - xi*nu;
  arma::vec sigsq_x_inv = pow(sigsq_x, -1.0);
  arma::mat eta_t = eta.t(); // N x J
  int S = D.n_rows;
  int K = eta.n_rows;
  arma::mat Theta(S,K);
  arma::mat eta_etaT = eta * eta_t;
  arma::mat SigInv(K,K);
  double prec_bet = pow(betasq_th, -1.0);
  arma::mat prec_gamma = pow(gammasq_th, -1.0);
  for( int s=0; s<S; s++ ){
    SigInv = sigsq_x_inv(s) * eta_etaT;
    for( int k=0; k<K; k++ ){
      SigInv(k,k) = SigInv(k,k) + prec_bet*prec_gamma(s,k)*tau_th(k);
    }
    SigInv = SigInv.i();
    Theta.row(s) = arma::mvnrnd(SigInv*eta*(D.row(s).t())*sigsq_x_inv(s), SigInv, 1).t();
  }
  return Theta; 
}

// [[Rcpp::export]]
double sample_betasq_th(double t, arma::mat Theta, arma::mat gammasq_th, 
                        arma::vec tau_th) {
  /* 
  * Function for sampling global variance term \beta^2, from Makalic paper.
  * A priori \beta^2|t ~ IG(1/2,1/t).
  * Theta and gamma_th (variance term) are S x K, tau_th is K x 1 (precision term).
  * Note that arma::randg uses scale parameterization, i.e. 1/rate.
  * 
  */
  int S = Theta.n_rows;
  int K = Theta.n_cols;
  double shape;
  double rate;
  double samp;
  arma::mat prec_gamma = pow(gammasq_th, -1.0);
  
  // Specify shape and initial rate parameters
  shape = (S*K+1)/2;
  rate = 1.0/t;
  
  // Calculate the components of the Gamma posterior
  for( int k=0; k<K; k++ ){
    double tmp_sum = 0.0;
    for( int s=0; s<S; s++ ){
      tmp_sum += prec_gamma(s,k) * Theta(s,k) * Theta(s,k);
    }
    rate += 0.5*tau_th(k)*tmp_sum;
  }
  
  samp = arma::randg( 1,arma::distr_param(shape,1/rate) )(0);
  
  return 1/samp; // INVERSE gamma distributed bc variance term
}

// [[Rcpp::export]]
arma::mat sample_gammasq_th(arma::mat s_mat, arma::mat Theta, double betasq_th, 
                            arma::vec tau_th) {
  /* 
  * Function for sampling local variance term gammasq_th, from Makalic paper.
  * A priori gammasq_th_{sk}|s_{sk} ~ IG(1/2,1/s_{sk}).
  * Theta and gamma_th (variance term) are S x K, tau_th is K x 1 (precision term).
  * Note that arma::randg uses scale parameterization, i.e. 1/rate.
  * 
  */
  int S = Theta.n_rows;
  int K = Theta.n_cols;
  arma::mat gammasq_th(S,K);
  double shape;
  double rate;
  
  // Specify shape and initial rate parameters
  shape = 1.0;
  
  // Calculate the components of the Gamma posterior
  for( int s=0; s<S; s++ ){
    for( int k=0; k<K; k++ ){
      rate = 1.0/s_mat(s,k) + Theta(s,k)*Theta(s,k)*tau_th(k)/(2*betasq_th);
      gammasq_th(s,k) = arma::randg( 1,arma::distr_param(shape,1/rate) )(0);
    }
  }
  
  return 1/gammasq_th; // INVERSE gamma distributed bc variance term
}

// [[Rcpp::export]]
arma::mat sample_s_mat(arma::mat gammasq_th) {
  /* 
  * Function for sampling hyper-prior for gammasq_th, from Makalic paper.
  * A priori gammasq_th_{sk}|s_{sk} ~ IG(1/2,1/s_{sk}).
  * 
  */
  int S = gammasq_th.n_rows;
  int K = gammasq_th.n_cols;
  arma::mat s_mat(S,K);
  double shape;
  double rate;
  
  // Specify shape and initial rate parameters
  shape = 1.0;
  
  // Calculate the components of the Gamma posterior
  for( int s=0; s<S; s++ ){
    for( int k=0; k<K; k++ ){
      rate = 1.0 + 1.0/gammasq_th(s,k);
      s_mat(s,k) = arma::randg( 1,arma::distr_param(shape,1/rate) )(0);
    }
  }
  
  return 1/s_mat; // INVERSE gamma distributed bc variance term
}

// [[Rcpp::export]]
double sample_t(double betasq_th) {
  /* 
  * Function for sampling hyper-prior for betasq_th, from Makalic paper.
  * A priori betasq_th|t ~ IG(1/2,1/t), t ~ IG(1/2,1/2)
  * 
  */
  double shape = 1.0;
  double rate = 1.0 + 1.0/betasq_th;
  betasq_th = arma::randg( 1,arma::distr_param(shape,1/rate) )(0);
  
  return 1/betasq_th; // INVERSE gamma distributed bc variance term
}

/* *********** MEAN SAMPLING FUNCTIONS *********** 
 * 
 * Used for sampling the mean of Y and X, getting covariance matrix
 * for GP for Y, and sampling associated hyper-parameters.
 * Specifically, YplusMean = Ymean + Lambda eta + epsilon and 
 * XplusMean = Xmean + Theta eta + Xi nu + e.
 * 
 */

// [[Rcpp::export]]
arma::vec sample_meanY(arma::mat YplusMean, arma::mat Lambda, arma::mat eta, 
                       arma::vec sigsq_y, arma::mat alphCovDD, arma::mat obs_Y){
  /* 
   * Function for sampling D-dimensional mean vector for Y assuming a GP prior.
   * YplusMean is D x N matrix of data.
   * Lambda is D x K factor loadings matrix.
   * eta is K x N matrix where column N is the length-K vec \eta_i associated with obs i.
   * alphCovDD is the D x D matrix with entry (d, d') being sig*exp(-r||d-d'||^2).
   * sigsq_y is the D vector of variance terms for Y.
   * obs_Y is the D x N matrix where (d,i) is 1 if Y is observed and 0 else.
   * 
   */
  int D = YplusMean.n_rows;
  int N = YplusMean.n_cols;
  arma::mat Yst = YplusMean - Lambda*eta;
  arma::vec S_st(D);
  arma::vec w_sum(D);
  arma::vec gp_mu(D);
  arma::mat gp_cov(D,D);
  S_st.fill(0.0);
  w_sum.fill(0.0);
  
  // Only count elements that are observed 
  for( int i=0; i<N; i++ ){ 
    for( int d=0; d<D; d++ ){
      double obs_Y_tmp = obs_Y(d,i);
      if( obs_Y_tmp>0 ){
        double w_tmp = (obs_Y_tmp / sigsq_y(d)); // inverse variance weights
        S_st(d) += w_tmp * Yst(d,i);
        w_sum(d) += w_tmp;
      }
    }
  }
  arma::vec w_sum_inv = 1/w_sum; // D-dimensional
  S_st = S_st % w_sum_inv; // % is element-wise multiplication
  // Calculate GP mean and covariance.
  arma::mat tmp_mult = alphCovDD;
  for( int d=0; d<D; d++ ){ tmp_mult(d,d) += w_sum_inv(d); } // + nugget(d)
  tmp_mult = alphCovDD * tmp_mult.i();
  gp_mu = tmp_mult * S_st;
  gp_cov = alphCovDD - tmp_mult * alphCovDD;
  
  return arma::mvnrnd(gp_mu, gp_cov, 1); // Ymean
}

// [[Rcpp::export]]
double sample_psi_Ymn(double g, arma::vec Ymean, 
                      arma::mat covDD, double nugget){
  /* 
   * Function for sampling global precision term psi_ymn, 
   * for Ymean ~ GP(0_D, c_k()) where c_k(d,d') = (1/psi_ymn) exp(-r||d-d'||^2).
   * Here covDD is the matrix with entry (i,j) being exp(-r*||d_i-d_j||^2),
   * A priori, psi_lam ~ Ga(g/2, g/2).
   * nugget is the value added to the diagonal for computational stability.
   * 
   */
  
  int D = Ymean.n_rows;
  covDD.diag() += nugget;
  arma::mat covDD_inv = covDD.i();
  double shape;
  double rate;
  double samp;
  
  // Specify shape parameter and rate parameter.
  shape = (g+D)/2.0;
  rate = g/2.0;
  rate += ((Ymean.t() * covDD_inv * Ymean )/2.0).eval()(0);
  
  samp = arma::randg( 1,arma::distr_param(shape,1/rate) )(0);
  return samp;
}

// [[Rcpp::export]]
arma::vec sample_meanZ(arma::mat ZplusMean, arma::mat Theta, arma::mat eta, 
                       arma::mat xi, arma::mat nu, arma::vec sigsq_z, arma::vec tau_Zmn){
  /* 
   * Function for sampling S-dimensional mean vector for Z assuming a N(0,1/tau_s) prior.
   * ZplusMean is S x N matrix of data.
   * Theta is S x K factor loadings matrix.
   * eta is K x N matrix where column N is the length-K vec \eta_i associated with obs i.
   * xi is S x J factor loadings matrix.
   * nu is J x N matrix where column N is the length-J vec \nu_i associated with obs i.
   * sigsq_z is the D vector of variance terms for Z.
   * tau_Zmn is the S vector of precision terms for Zmean.
   * 
   */
  int S = ZplusMean.n_rows;
  int N = ZplusMean.n_cols;
  arma::mat Zst = ZplusMean - (Theta * eta + xi * nu); // S x N
  arma::vec Zmean(S);
  arma::vec Ztmp(N);
  double tau_s;
  double sigsq_s;
  double post_mean;
  double post_sigsq;
  
  // Only count elements that are observed 
  for( int s=0; s<S; s++ ){
    Ztmp = Zst.row(s);
    tau_s = tau_Zmn(s);
    sigsq_s = sigsq_z(s);
    post_sigsq = 1.0 / ( 1.0/tau_s + N/sigsq_s );
    post_mean = post_sigsq * arma::sum(Ztmp) / sigsq_s;
    Zmean(s) = rnormArma(post_mean, post_sigsq);
  }
  
  return Zmean;
}

// [[Rcpp::export]]
arma::vec sample_tau_Zmn(arma::vec Zmean) {
  /* 
   * Function for sampling local variance term gammasq_th, from Makalic paper.
   * A priori gammasq_th_{sk}|s_{sk} ~ IG(1/2,1/s_{sk}).
   * Theta and gamma_th (variance term) are S x K, tau_th is K x 1 (precision term).
   * Note that arma::randg uses scale parameterization, i.e. 1/rate.
   * 
   */
  int S = Zmean.n_rows;
  arma::vec tau_Zmn(S);
  double shape;
  double rate;
  
  // Specify shape and initial rate parameters
  shape = 1.0;
  
  // Calculate the components of the Gamma posterior
  for( int s=0; s<S; s++ ){
    rate = 0.5 + Zmean(s)*Zmean(s)/2.0;
    tau_Zmn(s) = arma::randg( 1,arma::distr_param(shape,1/rate) )(0);
  }
  
  return tau_Zmn; // Precision term
}

    

/* *********** LAMBDA SAMPLING FUNCTIONS *********** 
 * 
 * Used for sampling the factor loading matrix \Lambda, getting covariance matrix
 * for GP, and smapling associated precision hyper-parameter psi_lam (global);
 * alpha_lam (K x 1) is \alpha_lam_{k} = (psi_lam * prod{h=1}^k delta_ome_{h})^{-1}.
 * 
 */

// [[Rcpp::export]]
arma::mat get_covDD(arma::vec d_vec, double l){
  /* 
  * Function for getting D x D matrix where entry (d, d') is exp(-(1/l^2)||d-d'||^2).
  * The lengthscale l determines the length of the 'wiggles' in your function. 
  * In general, you won't be able to extrapolate more than l units away from your data.
  * See https://www.cs.toronto.edu/~duvenaud/cookbook/
  * 
  */
  int D = d_vec.n_elem;
  double inv_lsq = pow(l, -2.0);
  arma::mat covDD(D,D);
  for( int i=0; i<D; i++ ){
    for( int j=(i+1); j<D; j++ ){
      double d1 = d_vec(i);
      double d2 = d_vec(j);
      double tmp_cov = exp(-0.5*inv_lsq*pow(fabs(d1-d2),1.9999)); // vs. 2.0 for numeric stability
      covDD(i, j) = tmp_cov;
      covDD(j, i) = tmp_cov;
    }
  }
  for( int i=0; i<D; i++ ){ covDD(i, i) = 1.0; }
  
  return covDD;
}

// [[Rcpp::export]]
arma::mat sample_Lambda(arma::mat Y, arma::mat Lambda, arma::mat eta, 
                        arma::vec alpha_lam, arma::vec sigsq_y, 
                        arma::mat covDD, arma::mat obs_Y){
  /* 
  * Function for sampling each column of Lambda = [\lambda_{:,1},...,\lambda_{:,K}].
  * Each column is D-dimensional vector and is given a GP prior, i.e., 
  * \lambda_{:,k} ~ GP(0_D, c_k()) where c_k(d,d') = \alpha_k exp(-r||d-d'||^2).
  * eta is K x N matrix where column N is the length-K vec \eta_i associated with obs i.
  * covDD is the D x D matrix with entry (d, d') being exp(-r||d-d'||^2).
  * sigsq_y is the D vector of variance terms for Y.
  * obs_Y is the D x N matrix where (d,i) is 1 if Y is observed and 0 else.
  * 
  */
  int D = Y.n_rows;
  int N = Y.n_cols;
  int K = eta.n_rows;
  arma::rowvec eta_k_save(N); 
  arma::mat Y_st(D,N);
  arma::mat covDD_k(D,D);
  arma::vec S_st(D);
  arma::vec w_sum(D);
  arma::vec gp_mu_k(D);
  arma::mat gp_cov_k(D,D);
  
  // Sample each length-D column vector \lambda_{1:D,k} for k=1,...,K
  for( int k=0; k<K; k++ ){
    eta_k_save = eta.row(k);
    // Set kth row of eta to zero so kth col of Lambda doesn't contribute to Y_st.
    eta.row(k).fill(0.0);
    // Get Y* (col i is Y*_i = (Y_i - \sum_{h!=k} \lambda_{col h}\eta_{h,i})/\eta_{k,i} )
    Y_st = Y - Lambda*eta;
    S_st.fill(0.0);
    w_sum.fill(0.0);
    for( int i=0; i<N; i++ ){ 
      double eta_k_i = eta_k_save(i);
      Y_st.col(i) = (Y_st.col(i)) / eta_k_i;
      for( int d=0; d<D; d++ ){
        double obs_Y_tmp = obs_Y(d,i);
        if( obs_Y_tmp>0 ){
          double w_tmp = eta_k_i * eta_k_i * (obs_Y_tmp / sigsq_y(d)); // inverse variance weights
          S_st(d) += w_tmp * Y_st(d,i);
          w_sum(d) += w_tmp;
        }
      }
    }
    arma::vec w_sum_inv = 1/w_sum; // D-dimensional
    S_st = S_st % w_sum_inv; // % is element-wise multiplication
    // Get column-specific GP covariance matrix covDD_k.
    covDD_k = alpha_lam(k)*covDD;
    // Calculate GP mean and covariance.
    arma::mat tmp_mult = covDD_k;
    for( int d=0; d<D; d++ ){ tmp_mult(d,d) += w_sum_inv(d); } // + nugget(d)
    tmp_mult = covDD_k * tmp_mult.i();
    gp_mu_k = tmp_mult * S_st;
    gp_cov_k = covDD_k - tmp_mult * covDD_k;
    // Sample kth column of Lambda.
    Lambda.col(k) = arma::mvnrnd(gp_mu_k, gp_cov_k, 1);
    // Set kth row of eta back to original value.
    eta.row(k) = eta_k_save;
  }
  
  return Lambda;
}

// [[Rcpp::export]]
Rcpp::List sample_Lambda_err(arma::mat Y, arma::mat Lambda, arma::mat eta, 
                            arma::vec alpha_lam, arma::vec sigsq_y, 
                            arma::mat covDD, arma::mat obs_Y){
  /* 
  * Function for sampling each column of Lambda = [\lambda_{:,1},...,\lambda_{:,K}].
  * Each column is D-dimensional vector and is given a GP prior, i.e., 
  * \lambda_{:,k} ~ GP(0_D, c_k()) where c_k(d,d') = \alpha_k exp(-r||d-d'||^2).
  * eta is K x N matrix where column N is the length-K vec \eta_i associated with obs i.
  * covDD is the D x D matrix with entry (d, d') being exp(-r||d-d'||^2).
  * sigsq_y is the D vector of variance terms for Y.
  * obs_Y is the D x N matrix where (d,i) is 1 if Y is observed and 0 else.
  * 
  */
  int D = Y.n_rows;
  int N = Y.n_cols;
  int K = eta.n_rows;
  arma::rowvec eta_k_save(N); 
  arma::mat Y_st(D,N);
  arma::mat covDD_k(D,D);
  arma::vec S_st(D);
  arma::vec w_sum(D);
  arma::vec gp_mu_k(D);
  arma::mat gp_cov_k(D,D);
  arma::mat Lambda_old = Lambda;
  bool bad = false;
  
  // Sample each length-D column vector \lambda_{1:D,k} for k=1,...,K
  for( int k=0; k<K; k++ ){
    eta_k_save = eta.row(k);
    // Set kth row of eta to zero so kth col of Lambda doesn't contribute to Y_st.
    eta.row(k).fill(0.0);
    // Get Y* (col i is Y*_i = (Y_i - \sum_{h!=k} \lambda_{col h}\eta_{h,i})/\eta_{k,i} )
    Y_st = Y - Lambda*eta;
    S_st.fill(0.0);
    w_sum.fill(0.0);
    for( int i=0; i<N; i++ ){ 
      double eta_k_i = eta_k_save(i);
      Y_st.col(i) = (Y_st.col(i)) / eta_k_i;
      for( int d=0; d<D; d++ ){
        double obs_Y_tmp = obs_Y(d,i);
        if( obs_Y_tmp > 0 ){
          double w_tmp = eta_k_i * eta_k_i * (obs_Y_tmp / sigsq_y(d)); // inverse variance weights
          S_st(d) += w_tmp * Y_st(d,i);
          w_sum(d) += w_tmp;
        }
      }
    }
    arma::vec w_sum_inv = 1/w_sum; // D-dimensional
    S_st = S_st % w_sum_inv; // % is element-wise multiplication
    // Get column-specific GP covariance matrix covDD_k.
    covDD_k = alpha_lam(k)*covDD;
    // Calculate GP mean and covariance.
    arma::mat tmp_mult = covDD_k;
    for( int d=0; d<D; d++ ){ tmp_mult(d,d) += w_sum_inv(d); } // + nugget(d)
    tmp_mult = covDD_k * tmp_mult.i();
    gp_mu_k = tmp_mult * S_st;
    gp_cov_k = covDD_k - tmp_mult * covDD_k;
    // Sample kth column of Lambda.
    try{
      Lambda.col(k) = arma::mvnrnd(gp_mu_k, gp_cov_k, 1);
    } catch (...) {
      bad = true;
      break;
    }
    // Set kth row of eta back to original value.
    eta.row(k) = eta_k_save;
  }
  
  if( bad ){
    return Rcpp::List::create(Rcpp::Named("Lambda") = Lambda_old,
                              Rcpp::Named("bad") = bad);
  } else{
    return Rcpp::List::create(Rcpp::Named("Lambda") = Lambda,
                              Rcpp::Named("bad") = bad);
  }
}

// [[Rcpp::export]]
double sample_psi_lam(double g, arma::mat Lambda, arma::vec delta_lam, 
                      arma::mat covDD, double nugget){
  /* 
  * Function for sampling global precision term psi_lam, 
  * where \alpha_lam{k} = (psi_lam \prod_{h=1}^k delta_lam_{k})^(-1)
  * for \lambda_{:,k} ~ GP(0_D, c_k()) where c_k(d,d') = \alpha_k exp(-r||d-d'||^2).
  * Here covDD is the matrix with entry (i,j) being exp(-r*||d_i-d_j||^2),
  * and the kth covariance matrix without psi is (1/prod(delta_lam(1:k)))*covDD.
  * A priori, psi_lam ~ Ga(g/2, g/2).
  * nugget is the value added to the diagonal for computational stability.
  * 
  */
  
  int D = Lambda.n_rows;
  int K = Lambda.n_cols;
  covDD.diag() += nugget;
  arma::mat covDD_inv = covDD.i();
  arma::vec lam_k(D);
  double thet_tmp = 1.0;
  double shape;
  double rate;
  double samp;
  
  // Specify shape parameter, initialize rate parameter.
  shape = (g+D*K)/2.0;
  rate = g/2.0;
  
// Calculate the components of the Gamma posterior
  for( int k=0; k<K; k++ ){
    thet_tmp = thet_tmp/delta_lam(k);
    lam_k = Lambda.col(k);
    rate += ((lam_k.t() * thet_tmp * covDD_inv * lam_k )/2.0).eval()(0);
  }
  
  samp = arma::randg( 1,arma::distr_param(shape,1/rate) )(0);
  return samp;
}

/* *********** SHARED COLUMN SHRINKAGE SAMPLING FUNCTIONS *********** 
 * 
 * Used for sampling the multiplicative shrinkage term for the precision
 * element shared by both Lambda and Theta, tau_ome.
 * 
 */

// [[Rcpp::export]]
arma::vec sample_delta_ome(double a1, double a2, arma::vec delta_ome, 
                           double psi_lam, arma::mat Lambda, arma::mat covDD, double nugget,
                           arma::mat Theta, double betasq_th, arma::mat gammasq_th){
  /* 
  * Function for sampling multiplicative shrinkage terms delta_lam = [\delta_1,...,\delta_K], 
  * shared between BOTH Lambda and Theta, and Omega = [Lambda; Theta] (row-bound),
  * where \alpha_lam{k} = (psi_lam \prod_{h=1}^k delta_ome_{k})^(-1)
  * for \lambda_{:,k} ~ GP(0_D, c_k()) where c_k(d,d') = \alpha_k exp(-r||d-d'||^2)
  * and \theta_{jk} ~ N(0, \beta^2 \gamma_{jk}^2 \tau_ome_{k}^(-1)) independently.
  * Here covDD is the matrix with entry (i,j) being exp(-r*||d_i-d_j||^2), and the
  * hth cov matrix w/out \delta_k (h>=k) is (psi*prod(delta_lam(1:h,h!=k)))^(-1)*covDD.
  * Note the cov matrix for column k is (psi*prod(delta_lam(1:h)))^(-1)*covDD,
  * so the prec matrix can be written delta_lam(k)*psi*prod(delta_lam(1:h,h!=k))*covDD^(-1).
  * A priori, delta_ome_{1} ~ Ga(a1, 1), and delta_ome_{h} ~ Ga(a2, 1), h>1.
  * 
  */
  
  int D = Lambda.n_rows;
  int S = Theta.n_rows;
  int K = Lambda.n_cols;
  double prec_bet = pow(betasq_th, -1.0);
  arma::mat prec_gamma = pow(gammasq_th, -1.0);
  covDD.diag() += nugget;
  arma::mat psi_covDD_inv = psi_lam*covDD.i();
  double shape;
  double rate;
  
  for( int h=0; h<K; h++ ){
    
    // Specify shape and initial rate parameters
    if(h==0){
      shape = a1 + K*(D+S)/2;
      rate = 1.0;
    } else{ 
      shape = a2 + (K-h)*(D+S)/2; 
      rate = 1.0;
    }
    
    // Calculate the components of the Gamma posterior
    double tau_minus = 1.0;
    for( int k=0; k<K; k++ ){
      if( k != h ){
        tau_minus = tau_minus * delta_ome(k); // \prod_{l=1,l!=h}^{k} \delta_l
      }
      if( k >= h ){
        // Calculate the components of the Gamma posterior
        // First add in Lambda contribution
        arma::vec lam_k = Lambda.col(k);
        rate += (0.5*(lam_k.t() * tau_minus * psi_covDD_inv * lam_k )).eval()(0);
        // Then add in Theta contribution
        double tmp_sum = 0.0;
        for( int s=0; s<S; s++ ){
          tmp_sum += prec_gamma(s,k) * Theta(s,k) * Theta(s,k);
        }
        rate += 0.5*prec_bet*tau_minus*tmp_sum;
      }
    }
    
    // Sample a new value for \delta_h
    delta_ome(h) = arma::randg( 1,arma::distr_param(shape,1/rate) )(0);
    
  } // for( int h=0; h<K; h++ ){
  
  return delta_ome;
}

/* *********** MISSING DATA SAMPLING FUNCTIONS *********** 
 * 
 * Used for sampling new values of Y when Y itself isn't observed.
 * 
 */


// [[Rcpp::export]]
arma::mat sample_Y_miss(arma::mat Lambda, arma::mat eta, arma::vec sigsq_y,
                        arma::mat Y, arma::mat obs_Y, arma::vec Ymean){
  /* 
  * Function for sampling Y when Y isn't observed at all locations.
  * Actual relationship: Y = Lambda_true %*% eta_true + e_y.
  * Lambda is the D x K matrix of factor loadings.
  * Eta is the K x N matrix of factor scores.
  * sigsq_y is the D x 1 vector of error variance terms.
  * obs_Y is the D x N matrix where (d,i) is n (n is the no of times chem i is observed at dose d).
  * Ymean is the D x 1 vector with the mean of Y, i.e. Y = Ymean + Lambda eta + epsilon.
  * 
  */
  
  int D = Lambda.n_rows;
  int N = eta.n_cols;
  arma::mat E_Y = Lambda * eta; // Expectated value of Y (D x N).
  arma::mat Y_samp = Y;
  
  for( int d=0; d<D; d++ ){
    arma::rowvec obs_Y_d = obs_Y.row(d);
    double sigsq_y_d = sigsq_y(d);
    for( int n=0; n<N; n++ ){
      if( obs_Y_d(n)==0 ){ // If Y(d,n) is not observed...
        Y_samp(d,n) = rnormArma(E_Y(d,n) + Ymean(d), sigsq_y_d);
      }
    } // for( int n=0; n<N; n++ ){
  } // for( int d=0; d<D; d++ ){
  
  return Y_samp;
}

/* *********** SAMPLING LATENT VARS FOR NON-NORMAL DATA *********** 
 * 
 * Currently, count and binary data are supported.
 * 
 */

// [[Rcpp::export]]
Rcpp::List sample_X(std::vector< std::string > type, arma::mat X_original, arma::vec sigsq_x,
                    arma::mat Theta, arma::mat eta, arma::mat xi, arma::mat nu, arma::vec Zmean){
  /* 
  * Function for sampling the latent variable x underlying x_{original}.
  * For binary x_{original}, x|{all},x_o=0 ~ N_(E[x_o], 1) and x|{all},x_o=1 ~ N+(E[x_o], 1).
  * For count x_o, x|{all},x_o=0 ~ N_(E[x_o], \sig_x) and x|{all},x_o=j ~ N[j-1,j)(E[x_o], \sig_x).
  * type is a length-S vector specifying {"continuous", "binary", "count"} for each variable in x_o.
  * X_original is the S x N matrix of the observed X values.
  * sigsq_x is the S x 1 noise variance terms for the S variables in X.
  * E(X) = Theta (S x K) eta (K x N) + xi (S x J) nu (J x N).
  * See https://cran.r-project.org/web/packages/RcppDist/vignettes/RcppDist.pdf
  * 
  */
  int S = type.size();
  int N = X_original.n_cols;
  arma::mat X = X_original;
  double inf = std::numeric_limits<double>::infinity();
  std::string conti ("continuous");
  std::string binar ("binary");
  std::string count ("count");
  int inf_s=0; // Number of times an infinite sample is drawn (breaks sampler).
  arma::mat E_X = Theta * eta + xi * nu; // S x N (doesn't include Zmean, this gets added next line)
  for( int i=0; i<N; i++){ E_X.col(i) = E_X.col(i) + Zmean;}
  
  for( int s=0; s<S; s++ ){
    if( conti.compare(type[s]) != 0 ){ // If variable s is NOT continuous.
      arma::vec E_Xs = E_X.row(s).t();
      arma::vec Xs_original = X_original.row(s).t();
      double sig_xs = pow(sigsq_x(s), 0.5);
      arma::vec X_samp(N);
      for( int i=0; i<N; i++ ){
        double Xsi_original = Xs_original(i);
        if( binar.compare(type[s]) == 0 ){ // If variable s is binary.
          if( Xsi_original==0 ){
            X_samp(i) = r_truncnorm(E_Xs[i], sig_xs, -inf, 0);
            if( is_inf(X_samp(i)) ){ X_samp(i)=-50; inf_s=inf_s+1; }
          } else{
            X_samp(i) = r_truncnorm(E_Xs[i], sig_xs, 0, inf);
            if( is_inf(X_samp(i)) ){ X_samp(i)=50; inf_s=inf_s+1; }
          }
        } else{ // If variable s is count.
          if( Xsi_original==0 ){
            X_samp(i) = r_truncnorm(E_Xs[i], sig_xs, -inf, 0);
            if( is_inf(X_samp(i)) ){ X_samp(i)=-50; inf_s=inf_s+1; }
          } else{
            X_samp(i) = r_truncnorm(E_Xs[i], sig_xs, Xsi_original-1, Xsi_original);
          }
        }
      }
      X.row(s) = X_samp.t(); // Fill in row s with sampled X_{i,s} for all i.
    }
  }
  
  return Rcpp::List::create(Rcpp::Named("Z") = X,
                            Rcpp::Named("inf_samples") = inf_s);
}

// [[Rcpp::export]]
arma::mat sample_X_init(std::vector< std::string > type, arma::mat X_original, arma::vec sigsq_x){
  /* 
  * Function for sampling the latent variable x underlying x_{original}.
  * For binary x_{original}, x|{all},x_o=0 ~ N_(E[x_o], 1) and x|{all},x_o=1 ~ N+(E[x_o], 1).
  * For count x_o, x|{all},x_o=0 ~ N_(E[x_o], \sig_x) and x|{all},x_o=j ~ N[j-1,j)(E[x_o], \sig_x).
  * type is a length-S vector specifying {"continuous", "binary", "count"} for each variable in x_o.
  * X_original is the S x N matrix of the observed X values.
  * See https://cran.r-project.org/web/packages/RcppDist/vignettes/RcppDist.pdf
  * 
  */
  int S = type.size();
  int N = X_original.n_cols;
  arma::mat X = X_original;
  double inf = std::numeric_limits<double>::infinity();
  std::string conti ("continuous");
  std::string binar ("binary");
  std::string count ("count");
  
  for( int s=0; s<S; s++ ){
    if( conti.compare(type[s]) != 0 ){ // If variable s is NOT continuous.
      arma::vec Xs_original = X_original.row(s).t();
      double sig_xs = pow(sigsq_x(s), 0.5);
      arma::vec X_samp(N);
      for( int i=0; i<N; i++ ){
        double Xsi_original = Xs_original(i);
        if( conti.compare(type[s]) != 0 ){ // If variable s is binary.
          if( Xsi_original==0 ){
            X_samp(i) = r_truncnorm(0.0, sig_xs, -inf, 0);
          } else{
            X_samp(i) = r_truncnorm(0.0, sig_xs, 0, inf);
          }
        } else{ // If variable s is count.
          if( Xsi_original==0 ){
            X_samp(i) = r_truncnorm(0.0, sig_xs, -inf, 0);
          } else{
            X_samp(i) = r_truncnorm(Xsi_original+0.5, sig_xs, Xsi_original-1, Xsi_original);
          }
        }
      }
      X.row(s) = X_samp.t(); // Fill in row s with sampled X_{i,s} for all i.
    }
  }
  
  return X;
}

/* *********** MISC ADDITIONAL FUNCTIONS *********** 
 * 
 * For use in examples and display purposes, or in the FPCA (Y-only) sampler.
 * 
 */

// [[Rcpp::export]]
arma::mat get_sqexp_kernel(arma::vec d_vec, double l, double sig, double nugget){
  /* 
  * Function for getting D x D matrix where entry (d, d') is sig^2*exp(-(1/l^2)||d-d'||^2),
  * plus nugget along the diagonal (i.e., if d=d').
  * The lengthscale l determines the length of the 'wiggles' in your function;
  * in general, you won't be able to extrapolate more than l units away from your data.
  * The function variance sig determines how close to the mean the function stays.
  * See https://www.cs.toronto.edu/~duvenaud/cookbook/
  * 
  */
  int D = d_vec.n_elem;
  double inv_lsq = pow(l, -2.0);
  double sigsq = pow(sig, 2.0);
  arma::mat covDD(D,D);
  for( int i=0; i<D; i++ ){
    for( int j=(i+1); j<D; j++ ){
      double d1 = d_vec(i);
      double d2 = d_vec(j);
      double tmp_cov = sigsq*exp(-0.5*inv_lsq*pow(d1-d2,2.0));
      covDD(i, j) = tmp_cov;
      covDD(j, i) = tmp_cov;
    }
  }
  for( int i=0; i<D; i++ ){ covDD(i, i) = 1.0 + nugget; }
  
  return covDD;
}

// [[Rcpp::export]]
arma::mat sample_eta_fpca(arma::mat Y, arma::mat Lambda, arma::vec sigsq_y, arma::mat obs_Y) {
  /* A draw of the K x N matrix [\eta_1,...,\eta_N], the scores associated with Y.
  * Here sigsq_D_inv is (D+S) x 1 vector of inv diag entries of independent noise cov matrix for D.
  * 
  */ 
  arma::mat Lambda_t = Lambda.t(); // K x D
  arma::mat Y_t = Y.t(); // N x D
  int num_doses = Y.n_rows; // Equivalent to D
  int K = Lambda.n_cols;
  int N = Y.n_cols;
  arma::mat eta_all(K, N);
  arma::mat eta_all_t(N, K);
  // Set up what will be the covariance matrices used to sample \eta_i.
  arma::mat SigInv_Lam(num_doses, K);
  arma::mat SigStar(K, K);
  for( int i=0; i<N; i++ ){
    arma::vec obs_Y_i = obs_Y.col(i);
    int P_obs = sum(obs_Y_i>0); // Number of doses at which Y_i is observed.
    arma::uvec ind_tmp(P_obs);
    int d_tmp = 0;
    arma::vec sigsq_y_tmp = sigsq_y;
    for( int d=0; d<num_doses; d++ ){ // Not all Y are observed.
      if( obs_Y_i(d)>0 ){
        ind_tmp(d_tmp) = d;
        d_tmp += 1;
        sigsq_y_tmp(d) = sigsq_y(d)/obs_Y_i(d);
      }
    }
    for( int ss=num_doses; ss<num_doses; ss++ ){ // All X are observed.
      ind_tmp(d_tmp) = ss;
      d_tmp += 1;
    }
    // Get covariance matrix used when sampling each individual \eta_i.
    arma::vec sigsq_D_inv = pow(sigsq_y_tmp, -1.0);
    for( int p=0; p<num_doses; p++ ){ SigInv_Lam.row(p) = sigsq_D_inv(p) * Lambda.row(p); }
    SigStar = Lambda_t.cols(ind_tmp) * SigInv_Lam.rows(ind_tmp);
    for( int k=0; k<K; k++ ){ SigStar(k,k) = SigStar(k,k) + 1; }
    SigStar = SigStar.i();
    // Get (common) mean multiplication matrix used when sampling each individual \eta_i.
    arma::mat mu_mult(K, P_obs);
    mu_mult = SigStar * SigInv_Lam.rows(ind_tmp).t(); 
    arma::vec Y_t_i = Y_t.row(i).t();
    // Sample \eta_i for i=1,...,N
    eta_all_t.row(i) = arma::mvnrnd(mu_mult * Y_t_i.rows(ind_tmp), SigStar, 1).t();
  }
  eta_all = eta_all_t.t();
  return eta_all; 
}

// [[Rcpp::export]]
arma::vec sample_delta_lam(double a1, double a2, arma::vec delta_ome, 
                           double psi_lam, arma::mat Lambda, arma::mat covDD, double nugget){
  /* 
  * Function for sampling multiplicative shrinkage terms delta_lam = [\delta_1,...,\delta_K], 
  * for column shrinkage of Lambda,
  * where \alpha_lam{k} = (psi_lam \prod_{h=1}^k delta_ome_{k})^(-1)
  * for \lambda_{:,k} ~ GP(0_D, c_k()) where c_k(d,d') = \alpha_k exp(-r||d-d'||^2).
  * Here covDD is the matrix with entry (i,j) being exp(-r*||d_i-d_j||^2), and the
  * hth cov matrix w/out \delta_k (h>=k) is (psi*prod(delta_lam(1:h,h!=k)))^(-1)*covDD.
  * Note the cov matrix for column k is (psi*prod(delta_lam(1:h)))^(-1)*covDD,
  * so the prec matrix can be written delta_lam(k)*psi*prod(delta_lam(1:h,h!=k))*covDD^(-1).
  * A priori, delta_ome_{1} ~ Ga(a1, 1), and delta_ome_{h} ~ Ga(a2, 1), h>1.
  * 
  */
  
  int D = Lambda.n_rows;
  int K = Lambda.n_cols;
  covDD.diag() += nugget;
  arma::mat psi_covDD_inv = psi_lam*covDD.i();
  double shape;
  double rate;
  
  for( int h=0; h<K; h++ ){
    
    // Specify shape and initial rate parameters
    if(h==0){
      shape = a1 + K*D/2;
      rate = 1.0;
    } else{ 
      shape = a2 + (K-h)*D/2; 
      rate = 1.0;
    }
    
    // Calculate the components of the Gamma posterior
    double tau_minus = 1.0;
    for( int k=0; k<K; k++ ){
      if( k != h ){
        tau_minus = tau_minus * delta_ome(k); // \prod_{l=1,l!=h}^{k} \delta_l
      }
      if( k >= h ){
        // Calculate the components of the Gamma posterior
        // Add in Lambda contribution
        arma::vec lam_k = Lambda.col(k);
        rate += (0.5*(lam_k.t() * tau_minus * psi_covDD_inv * lam_k )).eval()(0);
      }
    }
    
    // Sample a new value for \delta_h
    delta_ome(h) = arma::randg( 1,arma::distr_param(shape,1/rate) )(0);
    
  } // for( int h=0; h<K; h++ ){
  
  return delta_ome;
}
