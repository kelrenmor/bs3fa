### Functions used when evaluating performance of each method ###

# Get pairwise distances between observations.
get_dists = function(eta, inds=1:ncol(eta), Lambda=NULL){
  # eta should be K x N (i.e., obs are along columns)
  # inds can either be the test inds or left blank for all indices.
  # If Lambda is provided, distance will be weighted by column
  # norms of Lambda for each draw of eta / Lambda.
  if(is.null(Lambda)){
    dists = as.matrix(dist(t(eta)))[inds,inds]
    dists = c( dists[lower.tri(dists, diag=FALSE)] )
  } else{
    wts = apply(Lambda, 2, function(x) sum(x^2))
    wtdeta = sweep(t(eta), 2, wts, function(x,y) x * sqrt(y))
    dists = as.matrix(dist(wtdeta))[inds,inds]
    dists = c( dists[lower.tri(dists, diag=FALSE)] )
  }
  return(dists)
}

# Split an array into lists along dimension n,
# e.g. turn an N x P x S array a into a length-S list of 
# N x P matrices with split.along.dim(a, 3).
split.along.dim <- function(a, n){
  setNames(lapply(split(a, arrayInd(seq_along(a), dim(a))[, n]),
                  array, dim = dim(a)[-n], dimnames(a)[-n]),
           dimnames(a)[[n]])
}

# Get norm of vector
norm_vec = function(x) sqrt(sum(x^2))

# Get norm of each row of matrix
norm_byrow = function(mat){
  if( is.matrix(mat) ){ # When K>1.
    n_rows = nrow(mat)
    for(i in 1:n_rows){
      mat[i,] = mat[i,]/norm_vec(mat[i,])
    }
  } else{ # When K=1.
    mat = mat/norm_vec(mat)
  }
  return(mat)
}

# Get norm of each col of matrix
norm_bycol = function(mat){
  if( is.matrix(mat) ){ # When K>1.
    n_cols = ncol(mat)
    for(j in 1:n_cols){
      mat[,j] = mat[,j]/norm_vec(mat[,j])
    }
  } else{ # When K=1.
    mat = mat/norm_vec(mat)
  }
  return(mat)
}

# Function to get the MSE of two vectors
get_mse = function(true_vec, pred_vec){
  return( mean( (true_vec-pred_vec)^2 ) )
}

# Function to produce summary statistics (mean and +/- sd)
data_summary = function(x) {
  m <- mean(x)
  ymin <- m-sd(x)
  ymax <- m+sd(x)
  return(c(y=m,ymin=ymin,ymax=ymax))
}

# Function for whether to flip sign of vector
get_sign = function(eta_true, eta){
  if(cor(c(eta_true), c(eta))<0){
    sign = -1
  } else{ sign = 1 }
  return(sign)
}

# Function to match up eta, Lambda, and Theta to their true indices
reorder_entries = function(eta_true, eta, Lambda_true, Lambda, Theta_true, Theta, 
                           is_bsssfa=T, alpha=0.05, cred_band="simultaneous"){
  # Figure out the correlation between pairwise correlations of eta_true and eta
  K_tru = nrow(eta_true)
  K_mod = nrow(eta)
  cor_mat = matrix(NA,nrow=K_tru,ncol=K_mod)
  rownames(cor_mat) = 1:K_tru
  for(k_true in 1:K_tru){
    et_tru = c(eta_true[k_true,])
    for(k_mod in 1:K_mod){
      if(is_bsssfa){
        et_mod = apply(eta[k_mod,,],1,mean)
      } else{ et_mod = c(eta[k_mod,]) }
      cor_mat[k_true, k_mod] = cor(et_tru,et_mod)
    }
  }
  # Get estimates for BSSSFA model if necessary
  if(is_bsssfa){
    eta_tmp = eta_low = eta_upp = matrix(NA,nrow=K_mod,ncol=ncol(eta_true))
    Lambda_tmp = Lambda_low = Lambda_upp = matrix(NA,nrow=nrow(Lambda_true),ncol=K_mod)
    Theta_tmp = Theta_low = Theta_upp = matrix(NA,nrow=nrow(Theta_true),ncol=K_mod)
    for(k_mod in 1:K_mod){
      # Get mean
      eta_tmp[k_mod,] = apply(eta[k_mod,,],1,mean)
      Lambda_tmp[,k_mod] = apply(Lambda[,k_mod,],1,mean)
      Theta_tmp[,k_mod] = apply(Theta[,k_mod,],1,mean)
      if(cred_band=="simultaneous"){
        credBands_eta = get_credBands(sampFuns=t(eta[k_mod,,]), alpha=alpha)
        credBands_Lam = get_credBands(sampFuns=t(Lambda[,k_mod,]), alpha=alpha)
        credBands_The = get_credBands(sampFuns=t(Theta[,k_mod,]), alpha=alpha)
        eta_low[k_mod,] = credBands_eta[,1]
        eta_upp[k_mod,] = credBands_eta[,2]
        Lambda_low[,k_mod] = credBands_Lam[,1]
        Lambda_upp[,k_mod] = credBands_Lam[,2]
        Theta_low[,k_mod] = credBands_The[,1]
        Theta_upp[,k_mod] = credBands_The[,2]
      } else{ # cred_band="pointwise"
        # Get lower 2.5% of 95% credible interval
        eta_low[k_mod,] = apply(eta[k_mod,,],1,function(x) unname(quantile(x,alpha/2)))
        Lambda_low[,k_mod] = apply(Lambda[,k_mod,],1,function(x) unname(quantile(x,alpha/2)))
        Theta_low[,k_mod] = apply(Theta[,k_mod,],1,function(x) unname(quantile(x,alpha/2)))
        # Get upper 97.5% of 95% credible interval
        eta_upp[k_mod,] = apply(eta[k_mod,,],1,function(x) unname(quantile(x,1-alpha/2)))
        Lambda_upp[,k_mod] = apply(Lambda[,k_mod,],1,function(x) unname(quantile(x,1-alpha/2)))
        Theta_upp[,k_mod] = apply(Theta[,k_mod,],1,function(x) unname(quantile(x,1-alpha/2)))
      }
    }
    eta = eta_tmp; Lambda = Lambda_tmp; Theta = Theta_tmp
  } else{eta_low=eta_upp=Lambda_low=Lambda_upp=Theta_low=Theta_upp=NULL}
  # Now get indices for re-indexing
  new_mod_ind = rep(NA,K_tru)
  for(k in 1:K_tru){
    max_ind = which(abs(cor_mat) == max(abs(cor_mat)), arr.ind=TRUE) 
    # max_ind[1] is the TRUE k, max_ind[2] is the best fitting model k
    new_mod_ind[max_ind[1]] = max_ind[2]
    # Change sign of each vector if necessary
    eta[max_ind[2],] = sign(cor_mat[max_ind]) * eta[max_ind[2],]
    Lambda[,max_ind[2]] = sign(cor_mat[max_ind]) * Lambda[,max_ind[2]]
    Theta[,max_ind[2]] = sign(cor_mat[max_ind]) * Theta[,max_ind[2]]
    eta_low[max_ind[2],] = sign(cor_mat[max_ind]) * eta_low[max_ind[2],]
    Lambda_low[,max_ind[2]] = sign(cor_mat[max_ind]) * Lambda_low[,max_ind[2]]
    Theta_low[,max_ind[2]] = sign(cor_mat[max_ind]) * Theta_low[,max_ind[2]]
    eta_upp[max_ind[2],] = sign(cor_mat[max_ind]) * eta_upp[max_ind[2],]
    Lambda_upp[,max_ind[2]] = sign(cor_mat[max_ind]) * Lambda_upp[,max_ind[2]]
    Theta_upp[,max_ind[2]] = sign(cor_mat[max_ind]) * Theta_upp[,max_ind[2]]
    # Set to 0 because now the kth row of eta_true is "done"
    cor_mat[max_ind[1],] = 0
    cor_mat[,max_ind[2]] = 0
  }
  # Finally, re-index
  new_mod_ind = c(new_mod_ind,setdiff(1:K_mod,new_mod_ind))
  eta = eta[new_mod_ind,]
  Lambda = Lambda[,new_mod_ind] 
  Theta = Theta[,new_mod_ind]
  eta_low = eta_low[new_mod_ind,]
  Lambda_low = Lambda_low[,new_mod_ind] 
  Theta_low = Theta_low[,new_mod_ind]
  eta_upp = eta_upp[new_mod_ind,]
  Lambda_upp = Lambda_upp[,new_mod_ind] 
  Theta_upp = Theta_upp[,new_mod_ind]
  # Save everything to a list, and voila!
  res = list("eta" = eta, "Lambda"=Lambda, "Theta"=Theta,
             "eta_low" = eta_low, "Lambda_low"=Lambda_low, "Theta_low"=Theta_low,
             "eta_upp" = eta_upp, "Lambda_upp"=Lambda_upp, "Theta_upp"=Theta_upp)
  return(res)
}

pred_drcurve <- function(Lambda_mod, eta_mod, rescale=1, alpha=0.05, cred_band="simultaneous"){
  nsims = dim(Lambda_mod)[3]
  D = dim(Lambda_mod)[1]
  N = dim(eta_mod)[2]
  # Get predicted dose response curve for each Gibbs iter
  dr_mod = array(NA, dim=c(D,N,nsims))
  for(s in 1:nsims){ dr_mod[,,s] = Lambda_mod[,,s] %*% eta_mod[,,s] }
  dr_mod = dr_mod / rescale
  # Get estimated dose response curve and 95% credible interval
  dr_est = dr_low = dr_upp = matrix(NA, nrow=D, ncol=N)
  for(i in 1:N){
    # Get mean and lower 2.5% and upper 97.5% of 95% credible interval
    dr_est[,i] = apply(dr_mod[,i,],1,mean)
    if(cred_band=="simultaneous"){
      credBands = get_credBands(sampFuns=t(dr_mod[,i,]), alpha=alpha)
      dr_low[,i] = credBandsY[,1]
      dr_upp[,i] = credBandsY[,2]
    } else{ # cred_band="pointwise"
      dr_low[,i] = apply(dr_mod[,i,],1,function(x) unname(quantile(x,alpha/2)))
      dr_upp[,i] = apply(dr_mod[,i,],1,function(x) unname(quantile(x,1-alpha/2)))
    }
  }
  dr_list = list("dr_est"=dr_est, "dr_low"=dr_low, "dr_upp"=dr_upp)
  return(dr_list)
}

