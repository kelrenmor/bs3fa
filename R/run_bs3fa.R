run_bs3fa <- function(X, Y, K, J, X_type=rep("continuous", nrow(X)), post_process=T,
                      D=ifelse(ncol(X)==ncol(Y), nrow(Y), length(unique(Y[,2]))),
                      dvec_unique= if(ncol(X)==ncol(Y)) 1:nrow(Y) else sort(unique(Y[,2])),
                      nsamps_save=500, thin=10, burnin=5000, nugget=1e-8, l=D*0.0008, 
                      update_ls=list("type"="auto", "niter_max"=500, "l_diff"=1/(10*D), 
                                     "reset_ls"=round(3*burnin/4), "l_new"=NULL),
                      homo_Y=T, print_progress=T, scale_X=T,
                      a_sig_y=NULL, b_sig_y=NULL, a_sig_x=1, b_sig_x=1, 
                      num_ls_opts=50, ls_opts='auto', save_original_data=F) # change these to sample length-scale!
{
  # X - S x N chemical feature matrix, where S is the number of features and N is the no of obs.
  #     If Y is provided in 'long' format colnames(X) must give IDs used for Y.
  # Y - D x N dose response curve matrix, where S is the number of doses and N is the no of obs.
  #     (if Y provided in this format, and no explicit labels are included in X/Y, col alignment assumed)
  #     Here Missing data should be coded as NA.
  #     Alternatively, T x 3 matrix where column 1 is ID (must also be provided in first row of X), 
  #     column 2 is dose, and column 3 is response. 
  # X_type - length-S vector giving 'type' of each variable in X ("continuous", "binary", "count" supported).
  # dvec_unique - 1:D by default, else vector or D x 1 matrix of doses corresponding to rows of Y.
  # K - The maximum dimension for the common latent space.
  # J - The maximum dimension for the feature-specific latent space.
  # post_process - If T, correct for rotational ambihuity, label/sign switching.
  # nsamps_save - The number of samples of Lambda, Theta, eta, missing Y to be saved.
  # thin - Every thin-th sample will be kept, the rest discarded.
  # burnin - The number of initial samples to be thrown out.
  #          Note that the total samples overall are burnin + thin*nsamps_save.
  # nugget - Add for numerical stability in inversion of CovDD.
  # l - GP length-scale; IMPORTANT parameter, set conservatively before initialization.
  # update_ls - A list with entries type (gives type of updating, either "auto", "manual", "sample", or "none"),
  #             niter_max (for auto type, max times to try new l to see if it works),
  #             l_diff (for auto type, difference by which to bump up in l at each step of tuner),
  #             l_new (for manual type, new l to switch to after some burn-in period),
  #             reset_ls (for manual/auto type, at what ss to reset length-scale param),
  #             NOTE if update_ls=="sample" a grid of num_ls_opts length-scale values will be tried
  #             OR set type to "none" to keep the same l throughout burnin and sampling.
  # random_init - 
  # homo_Y - Set to T for homoscedastic variance, F for hetero.
  # print_progress - Set to T to print sample number every iteration.
  # scale_X - Scale X (set variable means to 0 and SDs to 1).
  # a1_delta_xi/a1_delta_om - Hyperparameters for B&D shrinkage prior.
  # a2_delta_xi/a2_delta_om - Hyperparameters for B&D shrinkage prior.

  # Load libraries and Cpp functions
  library(abind)
  
  # Set some default settings
  fr_normalize=T
  random_init=T # Set to T to initialize with random numbers, F to initialize to SVD solution.
  a1_delta_xi=2.1; a1_delta_om=2.1; a2_delta_xi=3.1; a2_delta_om=3.1 # Hyperparameters for B&D shrinkage prior.
  bad_samp_tol=nsamps_save # Total number of bad_samps before killing sampler.
  return_original_scale=T # Get things back to their original scale before returning output
  
  ##### Do some checks:
  types = unique(X_type)
  cond = (sum(sapply(types, function(x) !(x %in% c("continuous","binary","count"))))==0)
  if( !cond ){ stop('X_type must be length-S vector containing only {"continuous","binary","count"}') }
  
  ##### Do preliminary data manipulation and normalization
  N = ncol(X); N1 = ncol(Y)
  if( N==N1 ){ # Then number of obs per chemical/dose combo is 1 (when observed).
    obs_Y = 1*(!is.na(Y))
    D = nrow(Y)
    longY = FALSE
  }else if(N1==3){ # Manipulate
    if( is.null(colnames(X)) ){ stop("If ncol(Y)=3 (ID, dose, response), then colnames(X) must give ID.") }
    IDs = as.character(colnames(X))
    IDs_y = as.character(Y[,1])
    dvec_unique = sort(unique(Y[,2]))
    D = length(dvec_unique) # Number of uniquely observed dose values.
    longY = TRUE
    Y_long = matrix(Y[,3]) # save long format Y for sampling noise variance term for Y
    dind_long = matrix(sapply(Y[,2], function(x) which(dvec_unique==x) - 1)) # indexed for C++
    IDs_long = rep(NA, nrow(Y))
    obs_Y = Y_new = matrix(NA, nrow=D, ncol=N)
    for(i in 1:length(IDs)){
      id = IDs[i]
      ind_tmp = which(IDs_y == id)
      IDs_long[ind_tmp] = (as.numeric(as.character(id))-1) # for C++ indexing
      doses_tmp = Y[ind_tmp,2]
      resps_tmp = Y[ind_tmp,3]
      for(d in 1:D){
        dind_tmp = which(doses_tmp == dvec_unique[d])
        tmp_sum = length(dind_tmp)
        obs_Y[d, i] = tmp_sum
        if(tmp_sum>0){
          Y_new[d, i] = mean(resps_tmp[dind_tmp])
        }
      }
    }
    Y = Y_new
  } else{
    stop("Must have either ncol(X)=ncol(Y) or ncol(Y)=3 (ID, dose, response).")
  }
  # Need to randomly init for non-observed sol'n
  all_obs=(sum(is.na(Y))==0); if( !all_obs ){random_init=T}
  all_nobs_mat = matrix(0, nrow=D, ncol=N) # for sampling Y_save
  ##### Do more data manipulation and normalization
  # Scale dvec_unique to be between 0 and 1 if not already so
  dvec_unique_original=dvec_unique
  dvec_unique=(dvec_unique-min(dvec_unique))/(max(dvec_unique)-min(dvec_unique))
  # Scale columns of X to have unit variance and 0 mean for continuous variables
  if( scale_X ){
    scaled_X = t(scale(t(X)))
    X[X_type=="continuous", ] = scaled_X[X_type=="continuous", ]
  }
  # Get number of unique values by row (so S total) for X
  num_un = apply(X, 1, function(vec) length(unique(vec)))
  # Automatically make X binary (0/1) if only two values
  for(s in 1:nrow(X)){
    if(num_un[s]==2){
      un_vals = unique(X[s,])
      X[s,] = 1*(X[s,] == un_vals[1]) # recode as 0/1
      X_type[s] = "binary"
    }
  }
  cond = !(num_un==1)
  X = X[cond,]
  if(sum(cond)==1){ X=matrix(X,nrow=1) }
  X_type = X_type[cond]
  S = nrow(X) # Now save row dimension of X (i.e., no of features)
  if( sum(!cond) > 0 ){print(paste(sum(!cond),"X variables have no variation, removed."))}
  not_cont = !(X_type=="continuous")
  # Normalize data by Frobenius norm (make Y on same relative scale as X)
  if(fr_normalize){
    Y_to_norm = Y
    for(j in 1:nrow(Y_to_norm)){
      tmp = Y_to_norm[j,]
      tmp[is.na(tmp)] = mean(tmp, na.rm=T)
      Y_to_norm[j,] = tmp
    }
    norm_Y = norm(Y_to_norm, type="F")
    norm_X = norm(X, type="F")
  } else{
    norm_Y = 1
    norm_X = 1
  }
  norm_rescale=norm_X/norm_Y
  Y = norm_rescale*Y # scale Y so relative weight is the same as that of X
  if(longY){Y_long = norm_rescale*Y_long}
  a_y = ifelse( is.null(a_sig_y), 1, a_sig_y)
  b_y = ifelse( is.null(b_sig_y), 1, norm_rescale^2 * b_sig_y)
  
  # Initialize parameters and define hyperparameter values
  init_list = sampler_init(random_init, N, D, S, K, J, X_type, X)
  list2env(init_list, environment()) # puts list elements in environment
  g_xi = g_psi = 1; 
  covDD = get_covDD(matrix(dvec_unique), l);
  
  # Set up framework to sample the length-scale, if user sets num_ls_opts > 1
  if(num_ls_opts>1){
    update_ls_bool = TRUE
    er_to_l = function(er){sqrt(er/6)} # function to go from effective range to length scale l
    # get_covDD() is parameterized as sig^2 exp(-0.5 ||d-d'||^2 / l^2)
    # For sig^2 exp(-phi ||d-d'||^2), back of envelope is effective range is 3/phi
    # phi = 0.5 l^(-2) ----> 3/phi = 6 l^2 (small ranges, i.e. smaller than range of data, cause issues)
    # so l = (effective_range / 6) ^(0.5)
    sample_ls = TRUE
    er_min = 1/D + 1e-2 # corresponds to minimum effective range
    er_max = 1 - 1e-2 # corresponds roughly to eff range spanning all data
    l_opts = seq(er_to_l(er_min), er_to_l(er_max), length.out=num_ls_opts)
    l_new = median(l_opts)
    # Pre-compute the covDD matrices, and the log determinants and inverses of each.
    covDD_all = covDDinv_all = array(NA, dim=c(D,D,num_ls_opts))
    logdetCovDD_all = rep(NA, num_ls_opts)
    for(j in 1:num_ls_opts){
      covDD_all[,,j] = get_covDD(matrix(dvec_unique), l_opts[j])
      covDDinv_all[,,j] = solve(covDD_all[,,j])
      logdetCovDD_all[j] = determinant(covDD_all[,,j], logarithm=TRUE)$modulus[1]
    }
    reset_ls = round(3*burnin/4)
  } else{
    # Handle l updating
    if( update_ls[["type"]]=="none" ){
      update_ls_bool = FALSE
    } else if( update_ls[["type"]]=="auto" ){
      l_diff = update_ls[["l_diff"]]
      niter_max = update_ls[["niter_max"]]
      reset_ls = update_ls[["reset_ls"]]
      update_ls_bool = TRUE
    } else if( update_ls[["type"]]=="manual" ){
      l_new = update_ls[["l_new"]]
      reset_ls = update_ls[["reset_ls"]]
      update_ls_bool = TRUE
    } else{
      stop("update_ls[['type']] must be one of 'auto', 'manual', 'none'")
    }
  }
  
  # Make matrices to save the samples of Lambda, Theta, eta, nu, and xi
  Theta_save = array(NA, dim=c(S,K,nsamps_save))
  Lambda_save = array(NA, dim=c(D,K,nsamps_save))
  eta_save = array(NA, dim=c(K,N,nsamps_save))
  xi_save = array(NA, dim=c(S,J,nsamps_save))
  nu_save = array(NA, dim=c(J,N,nsamps_save))
  if(homo_Y){ sigsq_y_save = rep(NA, nsamps_save) }else{ sigsq_y_save = matrix(NA, nrow=D, ncol=nsamps_save) }
  sigsq_x_save = matrix(NA, nrow=S, ncol=nsamps_save)
  Y_save = DRcurve_save = array(NA, dim=c(D,N,nsamps_save))
  X_save = array(NA, dim=c(S,N,nsamps_save))
  if(num_ls_opts>1){l_save = rep(NA, nsamps_save)}else{l_save=NULL}
  tau_save = matrix(NA, nrow=K, ncol=nsamps_save)
  tauxi_save = matrix(NA, nrow=J, ncol=nsamps_save)
  
  ##### Run sampler
  init=T # Whether or not intialization is needed (will be changed to F upon initialization in sampler)
  ind=1 # Starting index for saving values.
  bad_samps=0 # Number of samples for which cov matrix is not symmetric PD
  inf_samps=0 # Number of samples for which the non-continuous X samps come up inf
  nsamps = nsamps_save*thin + burnin # Total number of times to loop through sampler.
  update_samps = seq(1, nsamps, round(nsamps/20))
  psi_lam_min = Inf # Initialize to infinity so anything is smaller.
  for(ss in 1:nsamps){
    if( print_progress & ss%in%update_samps ){
      print(paste(sep="",round(100*ss/nsamps),"% done sampling"))
      print(paste(sep="",bad_samps," bad samples"))
    }
    
    ##### Update length-scale hyperparameter to be as 'smooth' as possible.
    if( init & update_ls_bool ){ 
      if( ss>reset_ls ){
        if( (update_ls[["type"]]=="auto") & (!(num_ls_opts>1)) ){
          l = update_l(l, l_diff, dvec_unique, niter_max, Y, Lambda, eta, 
                       alpha_lam_tmp, sigsq_y_vec, obs_Y)
        } else{ l = l_new }
        covDD = get_covDD(matrix(dvec_unique), l)
        init = F
      }
    }
    
    ##### Sample Y-specific factor loading matrix \Lambda and shrinkage params  #####
    
    # Loadings matrix
    Lam_samp = sample_Lambda_err(Y, Lambda, eta, alpha_lam, sigsq_y_vec, covDD, obs_Y)
    Lambda = Lam_samp$Lambda
    bad_samps = bad_samps + Lam_samp$bad
    if( bad_samps>bad_samp_tol){
      print(paste(sep="","Error: bad_samps=",bad_samps,"with l=",l,", try smaller l"))
      return(-1)
    }
    # Hyper-params
    psi_lam = sample_psi_lam(g_psi, Lambda, delta_ome, covDD, nugget)
    delta_ome = sample_delta_ome(a1_delta_om, a2_delta_om, delta_ome, psi_lam, 
                                 Lambda, covDD, nugget, Theta, betasq_th, gammasq_th)
    tau_ome = get_tau(delta_ome)
    alpha_lam = 1/(psi_lam*tau_ome)
    # Save min psi_lam value for checking smoothness parameter during burn-in
    if( ss<burnin & ss>min(round(burnin/2), round(reset_ls/2)) ){
      psi_lam_min = min(psi_lam_min, psi_lam)
      alpha_lam_tmp = 1/(psi_lam_min*get_tau(delta_ome))
    }
    
    # Length-scale (if that's the user-specified choice)
    if( (num_ls_opts>1) & (ss>reset_ls) ){ 
      lind = sample_lind(Lambda, alpha_lam, covDDinv_all, logdetCovDD_all)
      l = l_opts[lind]
      covDD = covDD_all[,,lind]
    }
    
    ##### Sample latent variable Z corresponding to non-continuous X  #####
    
    Z_samp = sample_X(X_type, X, sigsq_x_vec, Theta, eta, xi, nu)
    Z = Z_samp$Z
    inf_samps = inf_samps + 1*(sum(Z_samp$inf_samples)>1)
    if( inf_samps>bad_samp_tol){
      print(paste(sep="","Error: inf_samps=",inf_samps,", try increasing J"))
      return(-1)
    }
    
    ##### Sample X-specific factor loading matrix \xi, scores \nu, and shrinkage params  #####
    
    # Score vectors
    nu = sample_nu_all(Z, xi, eta, Theta, sigsq_x_vec)
    # Loadings matrix
    xi = sample_xi(Z, nu, eta, Theta, sigsq_x_vec, phi_xi, tau_xi)
    # Hyper-params
    phi_xi = sample_phi_xi(g_xi, tau_xi, xi)
    delta_xi = sample_delta_xi(a1_delta_xi, a2_delta_xi, xi, phi_xi, delta_xi)
    tau_xi = get_tau(delta_xi)
    
    ##### Sample X-specific factor loading matrix \Theta and shrinkage params  #####
    
    # Loadings matrix
    Theta = sample_Theta(Z, nu, eta, xi, sigsq_x_vec, betasq_th, gammasq_th, tau_ome)
    # Hyper-params (note delta_ome (and thus tau_ome) sampled in Lambda region later)
    betasq_th = sample_betasq_th(t, Theta, gammasq_th, tau_ome)
    gammasq_th = sample_gammasq_th(s_mat, Theta, betasq_th, tau_ome)
    # Hyper-hyper-params
    s_mat = sample_s_mat(gammasq_th)
    t = sample_t(betasq_th)
    
    ##### Sample shared factor score eta = [\eta_1, ..., \eta_N] #####
    
    eta = sample_eta_all(Y, Z, xi, nu, Lambda, Theta, sigsq_y_vec, sigsq_x_vec, obs_Y)
    
    ##### Sample error terms #####

    # Error terms for Y
    if(longY){
      Y_min_mu = get_Y_min_mu_long(Y_long, Lambda, eta, IDs_long, dind_long)
      sigsq_y_vec = sample_sigsq_longy(a_y, b_y, Y_min_mu, dind_long, homo_Y, D)
    } else{
      Y_min_mu = get_Y_min_mu(Y, Lambda, eta)
      sigsq_y_vec = sample_sigsq_y(a_y, b_y, Y_min_mu, obs_Y, homo_Y)
    }
    
    # Error terms for X
    X_min_mu = get_X_min_mu(Z, Theta, eta, xi, nu)
    sigsq_x_vec = sample_sigsq_x(a_sig_x, b_sig_x, X_min_mu, X_type)
    
    ##### Save samples #####
    if( ss>burnin & ss%%thin==0 ){
      Theta_save[,,ind] = Theta
      Lambda_save[,,ind] = Lambda
      eta_save[,,ind] = eta
      xi_save[,,ind] = xi
      nu_save[,,ind] = nu
      if(homo_Y){ sigsq_y_save[ind] = sigsq_y_vec[1] }else{ sigsq_y_save[,ind] = sigsq_y_vec }
      sigsq_x_save[,ind] = sigsq_x_vec
      Y_save[,,ind] = sample_Y_miss(Lambda, eta, sigsq_y_vec, Y, all_nobs_mat)
      DRcurve_save[,,ind] = Lambda %*% eta
      X_save[not_cont,,ind] = Z[not_cont,] # only save sampled X vals
      if(num_ls_opts>1){l_save[ind] = l}
      tau_save[,ind] = tau_ome;
      tauxi_save[,ind] = tau_xi;
      ind = ind + 1
    }

  }
  
  # Reset the Y-specific things back to their original scales!
  Lambda_save = Lambda_save / ifelse(return_original_scale,norm_rescale,1)
  sigsq_y_save = sigsq_y_save / ifelse(return_original_scale,norm_rescale^2,1)
  Y_save = Y_save / ifelse(return_original_scale,norm_rescale,1)
  DRcurve_save = DRcurve_save / ifelse(return_original_scale,norm_rescale,1)
  
  if( post_process & (K>1) ){
    
    if( print_progress ){print("Done sampling, beginning post-processing.")}
    
    ### Correct ambiguity for joint components
    # Correct rotational ambiguity.
    Omega_save = abind(Lambda_save, Theta_save, along=1)
    Omega_rotated = mcrotfact(split.along.dim(Omega_save,3), file=FALSE, ret='both')
    Omega_save = abind(Omega_rotated$samples, along=3)
    Lambda_save = Omega_save[1:D,,]; Theta_save = Omega_save[(D+1):(D+S),,]
    if(S==1){ Theta_save=array(Theta_save, dim=c(S,K,nsamps_save)) }
    # Apply inverse rotation to eta so Lam * eta remains unchanged by Lam rotation
    eta_save = invrot(split.along.dim(eta_save,3), Omega_rotated$rots)
    # Fix possible label/sign switching ambiguity.
    pivot = Omega_rotated$samples[[ round( length(Omega_rotated$samples)/2 ) ]]
    permsignList = lapply(Omega_rotated$samples, msfOUT, pivot)
    Theta_save = abind(applier(split.along.dim(Theta_save,3), permsignList), along=3)
    Lambda_save = abind(applier(split.along.dim(Lambda_save,3), permsignList), along=3)
    # Double t in two lapply calls since eta transposed relative to Lambda and Theta
    eta_save = abind(lapply(applier(lapply(eta_save, t), permsignList), t), along=3) 
    # Get out predicted mean for Lambda (on original scale if specified), eta, and Theta
    Theta_mean = apply(Theta_save,c(1,2),mean)
    Lambda_mean = apply(Lambda_save,c(1,2),mean)
    eta_mean = apply(eta_save,c(1,2),mean)
    
    ### Correct ambiguity for X-specific components
    # Correct rotational ambiguity.
    Xi_rotated = mcrotfact(split.along.dim(xi_save,3), file=FALSE, ret='both')
    xi_save = abind(Xi_rotated$samples, along=3)
    if(S==1){ xi_save=array(xi_save, dim=c(S,J,nsamps_save)) }
    # Apply inverse rotation to nu so Xi * nu remains unchanged by Xi rotation
    nu_save = invrot(split.along.dim(nu_save,3), Xi_rotated$rots)
    # Fix possible label/sign switching ambiguity.
    pivot = Xi_rotated$samples[[ round( length(Xi_rotated$samples)/2 ) ]]
    permsignList = lapply(Xi_rotated$samples, msfOUT, pivot)
    xi_save = abind(applier(split.along.dim(xi_save,3), permsignList), along=3)
    # Double t in two lapply calls since eta transposed relative to Lambda and Theta
    nu_save = abind(lapply(applier(lapply(nu_save, t), permsignList), t), along=3) 
    # Get out predicted mean for Lambda (on original scale if specified), eta, and Theta
    Xi_mean = apply(xi_save,c(1,2),mean)
    nu_mean = apply(nu_save,c(1,2),mean)
  } else{
    if( print_progress ){print("Done sampling, no post-processing performed.")}
    Theta_mean = 'Theta not identifiable, posterior mean not calculated'
    Lambda_mean = 'Lambda not identifiable, posterior mean not calculated'
    eta_mean = 'eta not identifiable, posterior mean not calculated'
    Xi_mean = 'Xi not identifiable, posterior mean not calculated'
    nu_mean = 'nu not identifiable, posterior mean not calculated'
  }
  
  ##### Get out predicted mean and 95% credible interval for Y (on original scale if specified)
  Y_mean = apply(DRcurve_save,c(1,2),mean)
  Y_ll = apply(Y_save,c(1,2),function(x) quantile(x, 0.025))
  Y_ul = apply(Y_save,c(1,2),function(x) quantile(x, 0.975))
  DR_ll = apply(DRcurve_save,c(1,2),function(x) quantile(x, 0.025))
  DR_ul = apply(DRcurve_save,c(1,2),function(x) quantile(x, 0.975))
  
  
  ##### Save everything in a list and return said list.
  res = list("Theta_save"=Theta_save, "Lambda_save"=Lambda_save, "eta_save"=eta_save, 
             "Xi_save" = xi_save, "nu_save" = nu_save, 
             "l_save" = l_save, "tau_save" = tau_save, "tauxi_save" = tauxi_save,
             "sigsq_y_save"=sigsq_y_save, "sigsq_x_save"=sigsq_x_save,
             "DRcurve_save"=DRcurve_save, "Y_save"=Y_save, "Z"=X_save, 
             # (ABOVE) All saves for parameters
             "DRcurve_mean"=Y_mean, "DR_ll"=DR_ll, "DR_ul"=DR_ul,
             "Y_mean"=Y_mean, "Y_ll"=Y_ll, "Y_ul"=Y_ul, 
             "Lambda_mean"=Lambda_mean, "Theta_mean"=Theta_mean, "eta_mean"=eta_mean, 
             "Xi_mean"=Xi_mean, "nu_mean"=nu_mean)
             # (ABOVE) Summary stats for parameters
  if(save_original_data){ # Input data (save to output for reference)
    res_dat = list("dvec_unique"=dvec_unique_original, "l"=l, "covDD"=covDD, 
                   "Y"=Y, "X"=X, "not_cont_X_vars"=not_cont,
                   "kept_X_vars_original"=cond, "S"=S, "D"=D, "K"=K, "J"=J)
    res = c(res, res_dat)
  }
  return(res)
  
}
