############# Bs^3FA data simulation #############

normalize_and_rotate <- function(mat){
  library(pracma)
  mat = pracma::gramSchmidt(mat)$Q # Orthogonalize Theta (if not already orthogonal)
  if(ncol(mat)!=1){
    mat = mat %*% varimax(mat)$rotmat # Rotate and normalize
  }
  return(mat)
}
simulate_lambda <- function(k, p){
  sapply(1:k, function(j) simulate_lambda_column(p, j))
}
simulate_lambda_column <- function(p, j){ # p is dimension of data, j is column number to simulate
  value = runif(n = p, min = .5, max = 1) * sample(c(-1, 1), size = p, replace=TRUE)
  nonzero = rbinom(n = p, size = 1, p = .4 + .2/j) # higher column number means more chance of 0s
  value[!nonzero] = 0
  value
}

simulate_data <- function(N,D,S,K,J,std_error_y,std_error_x,prob_miss=0,real_Y=NULL,prob_rep=0,prop_zero=0,
                          gp_ls=0.3, gp_nv=1, gp_nug=1e-6, sample_real=NULL, wt_xfacs=1, taper_lam=NULL,
                          X_type=rep("continuous",S)){
  # Load packages
  library(mvtnorm) # Use for multivariate normal sampling
  library(sparseEigen) # Use to make sparse orthogonal matrix.
  # Set up necessary matrices for Y
  if( is.null(real_Y) ){ # If no real dose response data are provided, simulate from a GP.
    doses_long=dvec_unique=(1:D)/D
    avg_dose_resp=rep(0,D)
    #Lambda_true=t( mvtnorm::rmvnorm(n=K, mean=avg_dose_resp, sigma=get_sqexp_kernel(dvec_unique, gp_ls, gp_nv, gp_nug)) )
    obs_tmp=t( mvtnorm::rmvnorm(n=K, mean=avg_dose_resp, sigma=get_sqexp_kernel(dvec_unique, gp_ls, gp_nv, gp_nug)) )
    svd_tmp = svd(t(obs_tmp))
    Lambda_true = as.matrix(svd_tmp$v[,1:K])
    #plot(Lambda_true, type="l")
    if( !is.null(sample_real) ){
      print('Warning: sample_real set to a number, but no real data provided. Need to provide real_Y to use.')
    }
  } else{ # If real dose response data are provided, smooth over those curves to generate simulated 'real' data.
    if( !is.null(sample_real) ){ # make real_Y a resampled version of itself
      real_Y = real_Y[sample(1:nrow(real_Y), sample_real, replace=T),]
    }
    doses = as.numeric(colnames(real_Y))
    doses_long = seq(min(doses), max(doses), length=D)
    Y_smooth = matrix(NA,nrow=nrow(real_Y),ncol=D)
    for(i in 1:nrow(real_Y)){
      lo = loess(real_Y[i,]~doses) # , span=0.7
      Y_smooth[i,] = predict(lo, doses_long)
      #plot(doses,real_Y[i,]); lines(doses_long, predict(lo, doses_long), col='red', lwd=2)
    }
    #plot(doses_long, apply(Y_smooth,2,mean), type="l", xlab="dose", ylab="response")
    avg_dose_resp = apply(Y_smooth,2,mean)
    Y_smooth = scale(Y_smooth, scale=FALSE) # subtract average curve from data
    svd_y = svd(Y_smooth)
    Lambda_true = as.matrix(svd_y$v[,1:K])
  }
  Lambda_true=normalize_and_rotate(Lambda_true)
  if( !is.null(taper_lam) ){
    # Make max 1 and scale others relative to max, sort decreasing
    taper_lam = sort(taper_lam/max(taper_lam), decreasing=T)
    for(k in 1:K){
      Lambda_true[,k] = Lambda_true[,k]*taper_lam[k]
    }
  }
  eta_true=matrix(rnorm(K*N), nrow=K, ncol=N)
  if(prop_zero>0){ # Make some number of chemicals just have constant dose-response curves
    if( !((prop_zero<=1) & (prop_zero>=0)) ){print("Warning: Need prop_zero in range of [0,1]")}
    eta_true[,sample(1:N,round(N*prop_zero))] = 0
  }
  # Get Y itself and make some values unobserved (if prob_miss > 0)
  true_curve = Lambda_true %*% eta_true
  if(prob_rep==0){ # Singly observed values per dose.
    e_y=matrix(rnorm(D*N, mean=0, sd=std_error_y),nrow=D,ncol=N)
    Y = true_curve + e_y
    for(ii in 1:nrow(Y)){for(jj in 1:ncol(Y)){if( rbinom(1,1,prob=prob_miss) ){ Y[ii,jj] = NA }}}
  } else{ # Multiply observed values per dose (prob_rep probability of multiple obs at a given dose).
    Y = matrix(NA, nrow=0, ncol=3)
    e_y = c()
    for(ii in 1:N){
      rep_tmp = rbinom(1,1,prob=prob_rep)
      if(rep_tmp){
        nreps = sample(2:5, 1)
      } else{
        nreps = 1
      }
      true_curve_tmp = true_curve[,ii]
      dose_samps = sort(sample(1:D, round((1-prob_miss)*length(1:D))))
      for(dd in dose_samps){ # If you observe values at a given dose.
        e_y_tmp = rnorm(nreps,mean=0,sd=std_error_y)
        e_y = c(e_y, e_y_tmp)
        obs_tmp = cbind(ii, doses_long[dd], true_curve_tmp[dd] + e_y_tmp)
        Y = rbind(Y, obs_tmp)
      }
    }
    Y = as.data.frame(Y)
    colnames(Y) = c("id","dose","resp")
    
  }
  
  # Set up necessary matrices for X
  tmp=matrix(rnorm(K*S),nrow=K,ncol=S)
  Theta_true=sparseEigen::spEigen(t(tmp) %*% tmp, q=K, rho=0.2)$vectors
  if( J!=0 ){
    xi_true=norm_bycol(simulate_lambda(J, S))
  } else{
    xi_true=matrix(0,nrow=S,ncol=J)
  }
  nu_true=matrix(rnorm(J*N), nrow=J, ncol=N)
  e_x=matrix(rnorm(S*N, mean=0, sd=std_error_x),nrow=S,ncol=N)
  # Get X itself
  X = Theta_true %*% eta_true + wt_xfacs * xi_true %*% nu_true + e_x
  colnames(X) = 1:N
  
  # Round depending on X_type variable
  # "continuous", "binary", "count"
  if( !(length(X_type)==S) ){print('Warning: Need length(X_type) equal to S.')}
  for(s in 1:S){
    type_tmp = X_type[s]
    if(type_tmp=="continuous"){
      # Do nothing
    } else if(type_tmp=="binary"){
      X[s,] = 1*(X[s,]>0)
    } else if(type_tmp=="count"){
      X[s,(X[s,]==0)] = 0
      X[s,] = round(5*X[s,])
    } else{
      print('Warning: Need X_type to only contain "continuous", "binary", "count".')
    }
  }

  dat_list = list("X"=X,"Y"=Y,"eta_true"=eta_true,"Lambda_true"=Lambda_true,"e_y"=e_y,
                  "Theta_true"=Theta_true,"xi_true"=xi_true,"nu_true"=nu_true,"e_x"=e_x,
                  "avg_dose_resp"=avg_dose_resp,"K"=K,"J"=J,"doses"=doses_long, 
                  "true_curve"=true_curve, "X_type"=X_type)
  return(dat_list)
}
