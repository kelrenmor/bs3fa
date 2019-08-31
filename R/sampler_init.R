sampler_init <- function(random_init, N, D, S, K, J, X_type, X){

  if(random_init){ # INITIALIZE RANDOMLY
    small_sd = 0.001
    med_sd = 1
    
    # Loadings components
    Lambda = matrix(rnorm(D*K, 0, small_sd), nrow=D, ncol=K)
    Theta = matrix(rnorm(S*K, 0, small_sd), nrow=S, ncol=K)
    xi = matrix(rnorm(S*J, 0, small_sd), nrow=S, ncol=J)
    # Scores components
    eta = matrix(rnorm(K*N, 0, med_sd), nrow=K, ncol=N)
    nu = matrix(rnorm(J*N, 0, small_sd), nrow=J, ncol=N)
    # Error components
    sigsq_y_vec = matrix(med_sd^2, nrow=D, ncol=1)
    sigsq_x_vec = matrix(small_sd^2, nrow=S, ncol=1)
    # Hyperparams for xi
    phi_xi = matrix(1, nrow=S, ncol=J)
    tau_xi = matrix(1, nrow=J, ncol=1)
    delta_xi = matrix(1, nrow=J, ncol=1)
    # Hyperparams for Theta
    betasq_th = 1
    gammasq_th = matrix(1, nrow=S, ncol=K)
    # Hyper-hyper params for horseshoe prior on Theta 
    s_mat = matrix(1, nrow=S, ncol=K)
    t = 1
    # Hyperparams for Lambda
    psi_lam = 1
    alpha_lam = matrix(1, nrow=K, ncol=1)
    # Hyperparameters for shared column shrinkage
    tau_ome = matrix(1, nrow=K, ncol=1)
    delta_ome = matrix(1, nrow=K, ncol=1)
    # Latent variable Z for non-continuous X
    Z = sample_X_init(X_type, X, sigsq_x_vec)
    
  } else{ # INITIALIZE TO SVD SOLUTIONS
    svd_xy = svd(rbind(X,Y))
    
    # Joint loadings components
    Theta = norm_bycol(as.matrix(svd_xy$u[1:S,1:K]))
    Lambda = norm_bycol(as.matrix(svd_xy$u[(S+1):(D+S),1:K]))
    # Joint scores components
    eta = t(Lambda) %*% Y
    # Individual components
    svd_xresid = svd(X-Theta%*%eta)
    xi = norm_bycol(as.matrix(svd_xresid$u[,1:J]))
    nu = t(xi) %*% (X-Theta%*%eta)
    # Error components
    sigsq_y_vec = matrix(rep(mean((Y-Lambda%*%eta)^2),D))
    sigsq_x_vec = matrix(rep(mean((X-Theta%*%eta-xi%*%nu)^2),S))
    # Hyperparams for xi
    phi_xi = matrix(1, nrow=S, ncol=J)
    tau_xi = matrix(1, nrow=J, ncol=1)
    delta_xi = matrix(1, nrow=J, ncol=1)
    # Hyperparams for Theta
    betasq_th = 1
    gammasq_th = matrix(1, nrow=S, ncol=K)
    # Hyper-hyper params for horseshoe prior on Theta 
    s_mat = matrix(1, nrow=S, ncol=K)
    t = 1
    # Hyperparams for Lambda
    psi_lam = 1
    alpha_lam = matrix(1, nrow=K, ncol=1)
    # Hyperparameters for shared column shrinkage
    tau_ome = matrix(1, nrow=K, ncol=1)
    delta_ome = matrix(1, nrow=K, ncol=1)
    Z = sample_X_init(X_type, X, sigsq_x_vec)
  }

  init_list = list("Lambda"=Lambda, "Theta"=Theta, "xi"=xi, "eta"=eta, "nu"=nu, "sigsq_y_vec"=sigsq_y_vec, 
                   "sigsq_x_vec"=sigsq_x_vec, "phi_xi"=phi_xi, "tau_xi"=tau_xi, "delta_xi"=delta_xi, 
                   "betasq_th"=betasq_th, "gammasq_th"=gammasq_th, "s_mat"=s_mat, "t"=t, "psi_lam"=psi_lam,
                   "alpha_lam"=alpha_lam, "tau_ome"=tau_ome, "delta_ome"=delta_ome)
  return(init_list)
}
