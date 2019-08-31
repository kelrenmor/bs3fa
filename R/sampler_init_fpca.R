sampler_init_fpca <- function(random_init, N, D, K){

  if(random_init){ # INITIALIZE RANDOMLY
    small_sd = 0.001
    med_sd = 1
    
    # Loadings components
    Lambda = matrix(rnorm(D*K, 0, small_sd), nrow=D, ncol=K)
    # Scores components
    eta = matrix(rnorm(K*N, 0, med_sd), nrow=K, ncol=N)
    # Error components
    sigsq_y_vec = matrix(med_sd^2, nrow=D, ncol=1)
    # Hyperparams for Lambda
    psi_lam = 1
    alpha_lam = matrix(1, nrow=K, ncol=1)
    # Hyperparameters for shared column shrinkage
    tau_ome = matrix(1, nrow=K, ncol=1)
    delta_ome = matrix(1, nrow=K, ncol=1)
    
  } else{ # INITIALIZE TO SVD SOLUTIONS
    svd_y = svd(Y)
    
    # Loadings components
    Lambda = norm_bycol(as.matrix(svd_y$u[,1:K]))
    # Joint scores components
    eta = t(Lambda) %*% Y
    # Error components
    sigsq_y_vec = matrix(rep(mean((Y-Lambda%*%eta)^2),D))
    # Hyperparams for Lambda
    psi_lam = 1
    alpha_lam = matrix(1, nrow=K, ncol=1)
    # Hyperparameters for shared column shrinkage
    tau_ome = matrix(1, nrow=K, ncol=1)
    delta_ome = matrix(1, nrow=K, ncol=1)
  }

  init_list = list("Lambda"=Lambda, "eta"=eta, "sigsq_y_vec"=sigsq_y_vec, 
                   "psi_lam"=psi_lam, "alpha_lam"=alpha_lam, "tau_ome"=tau_ome, "delta_ome"=delta_ome)
  return(init_list)
}
