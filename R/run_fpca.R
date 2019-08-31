run_fpca <- function(Y, K, dvec_unique=1:nrow(Y), post_process=T, Y_format='long',
                     nsamps_save=500, thin=10, burnin=5000, nugget=1e-8, l=nrow(Y)*0.0008, 
                     update_ls=list("type"="auto", "niter_max"=500, "l_diff"=1/(10*nrow(Y)), 
                                       "reset_ls"=round(3*burnin/4), "l_new"=NULL),
                     random_init=T, homo_Y=T, print_progress=T, a1_delta_om=2.1, a2_delta_om=3.1,
                     a_sig_y=1, b_sig_y=1, bad_samp_tol=nsamps_save, debug=F, sigsq_Y_fixed=NULL){
  
  # Y - D x N dose response curve matrix, where S is the number of doses and N is the no of obs.
  #     (if Y provided in this format, and no explicit labels are included in Y, col alignment assumed)
  #     Here Missing data should be coded as NA. If this is format, specify Y_format='wide'.
  #     Alternatively, (Y_format='long') T x 3 matrix where column 1 is ID, column 2 is dose, and column 3 is response. 
  # dvec_unique - 1:D by default, else vector or D x 1 matrix of doses corresponding to rows of Y.
  # K - The maximum dimension for the common latent space.
  # post_process - If T, correct for rotational ambihuity, label/sign switching.
  # nsamps_save - The number of samples of Lambda, Theta, eta, missing Y to be saved.
  # thin - Every thin-th sample will be kept, the rest discarded.
  # burnin - The number of initial samples to be thrown out.
  #          Note that the total samples overall are burnin + thin*nsamps_save.
  # nugget - Add for numerical stability in inversion of CovDD.
  # l - GP length-scale; SUPER IMPORTANT parameter, set conservatively before initialization.
  # update_ls - A list with entries type (gives type of updating, either "auto", "manual", or "none"),
  #             niter_max (for auto type, max times to try new l to see if it works),
  #             l_diff (for auto type, difference by which to bump up in l at each step of tuner),
  #             l_new (for manual type, new l to switch to after some burn-in period),
  #             reset_ls (for manual/auto type, at what ss to reset length-scale param),
  #             OR set type to "none" to keep the same l throughout burnin and sampling.
  # 
  # random_init - Set to T to initialize with random numbers, F to initialize to SVD solution.
  # homo_Y - Set to T for homoscedastic variance, F for hetero.
  # print_progress - Set to T to print sample number every iteration.
  # update_ls - Set to T to update length-scale hyperparam partway through sampling.
  # a1_delta_om/a2_delta_om - Hyperparameters for B&D shrinkage prior.
  # bad_samp_tol - Total number of bad_samps before killing sampler.
  # debug - For internal debugging, will print when Lambda sample is bad.
  # sigsq_Y_fixed - Default NULL, but to fix sigsq_Y_vec set to numeric value.
  
  # Load libraries and Cpp functions
  library(abind)
  
  ##### Do preliminary data manipulation and normalization
  if( Y_format=='wide' ){ # Then number of obs per chemical/dose combo is 1 (when observed).
    obs_Y = 1*(!is.na(Y))
    D = nrow(Y)
  }else if(ncol(Y)==3 & Y_format=='long'){ # Manipulate
    IDs = unique(as.character(Y[,1]))
    N = length(IDs)
    IDs_y = as.character(Y[,1])
    dvec_unique = sort(unique(Y[,2]))
    D = length(dvec_unique) # Number of uniquely observed dose values.
    obs_Y = Y_new = matrix(NA, nrow=D, ncol=N)
    for(i in 1:length(IDs)){
      id = IDs[i]
      ind_tmp = which(IDs_y == id)
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
    stop("Must have either ncol(Y)=[no. of unique doses] and specify Y_format='wide', or ncol(Y)=3 (ID, dose, response) and specify Y_format='long'.")
  }
  # Need to randomly init for non-observed sol'n (FIXME?)
  all_obs=(sum(is.na(Y))==0); if( !all_obs ){random_init=T}
  ##### Do more data manipulation and normalization
  # Scale dvec_unique to be between 0 and 1 if not already so
  dvec_unique_original=dvec_unique
  dvec_unique=(dvec_unique-min(dvec_unique))/(max(dvec_unique)-min(dvec_unique))
  norm_rescale=1 # legacy from other code with X involved
  
  # Initialize parameters and define hyperparameter values
  init_list = sampler_init_fpca(random_init, N, D, K)
  g_psi = 1
  list2env(init_list, environment()) # puts list elements in environment
  covDD = get_covDD(matrix(dvec_unique), l);
  # Handle l updating
  update_ls_bool = T
  if( update_ls[["type"]]=="none" ){
    update_ls_bool = F
  } else if( update_ls[["type"]]=="auto" ){
    l_diff = update_ls[["l_diff"]]
    niter_max = update_ls[["niter_max"]]
    reset_ls = update_ls[["reset_ls"]]
  } else if( update_ls[["type"]]=="manual" ){
    l_new = update_ls[["l_new"]]
    reset_ls = update_ls[["reset_ls"]]
  } else{
    stop("update_ls[['type']] must be one of 'auto', 'manual', 'none'")
  }
  
  # Make matrices to save the samples of Lambda, Theta, and eta
  Lambda_save = array(NA, dim=c(D,K,nsamps_save))
  eta_save = array(NA, dim=c(K,N,nsamps_save))
  sigsq_y_save = matrix(NA, nrow=D, ncol=nsamps_save)
  Y_save = array(NA, dim=c(D,N,nsamps_save))

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
        if( update_ls[["type"]]=="auto" ){
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
    if( debug ){
      print(ss)
      if( Lam_samp$bad==1 ){
        print("BAD LAMBDA SAMPLE!")
      }
    }
    if( bad_samps>bad_samp_tol){
      print(paste(sep="","Error: bad_samps=",bad_samps,"with l=",l,", try smaller l"))
      return(-1)
    }
    # Hyper-params
    psi_lam = sample_psi_lam(g_psi, Lambda, delta_ome, covDD, nugget)
    delta_ome = sample_delta_lam(a1_delta_om, a2_delta_om, delta_ome, psi_lam, 
                                 Lambda, covDD, nugget)
    tau_ome = get_tau(delta_ome)
    alpha_lam = 1/(psi_lam*tau_ome)
    # Save min psi_lam value for checking smoothness parameter during burn-in
    if( ss<burnin & ss>min(round(burnin/2), round(reset_ls/2)) ){
      psi_lam_min = min(psi_lam_min, psi_lam)
      alpha_lam_tmp = 1/(psi_lam_min*get_tau(delta_ome))
    }
    
    ##### Sample shared factor score eta = [\eta_1, ..., \eta_N] #####
    
    eta = sample_eta_fpca(Y, Lambda, sigsq_y_vec, obs_Y)
    
    ##### Sample error terms #####
    
    # Error terms for Y
    Y_min_mu = get_Y_min_mu(Y, Lambda, eta)
    if(is.null(sigsq_Y_fixed)){
      sigsq_y_vec = sample_sigsq_y(a_sig_y, b_sig_y, Y_min_mu, obs_Y, homo_Y)
    } else{
      sigsq_y_vec = matrix(rep(norm_rescale^2 * sigsq_Y_fixed, D))
    }

    ##### Save samples #####
    if( ss>burnin & ss%%thin==0 ){
      Lambda_save[,,ind] = Lambda
      eta_save[,,ind] = eta
      sigsq_y_save[,ind] = sigsq_y_vec
      Y_save[,,ind] = sample_Y_miss(Lambda, eta, sigsq_y_vec, Y, obs_Y)
      ind = ind + 1
    }
    
  }
  
  if( post_process & (K>1) ){
    if( print_progress ){print("Done sampling, beginning post-processing.")}
    # Correct rotational ambiguity.
    Omega_rotated = mcrotfact(split.along.dim(Lambda_save,3), file=FALSE, ret='both')
    Lambda_save = abind(Omega_rotated$samples, along=3)
    # Apply inverse rotation to eta so Lam * eta remains unchanged by Lam rotation
    eta_save = invrot(split.along.dim(eta_save,3), Omega_rotated$rots)
    # Fix possible label/sign switching ambiguity.
    pivot = Omega_rotated$samples[[ round( length(Omega_rotated$samples)/2 ) ]]
    permsignList = lapply(Omega_rotated$samples, msfOUT, pivot)
    Lambda_save = abind(applier(split.along.dim(Lambda_save,3), permsignList), along=3)
    # Double t in two lapply calls since eta transposed relative to Lambda and Theta
    eta_save = abind(lapply(applier(lapply(eta_save, t), permsignList), t), along=3) 
    # Get out predicted mean for Lambda, and eta
    Lambda_mean = apply(Lambda_save,c(1,2),mean)
    eta_mean = apply(eta_save,c(1,2),mean)
  } else{
    if( print_progress ){print("Done sampling, no post-processing performed.")}
    Lambda_mean = 'Lambda not identifiable, posterior mean not calculated'
    eta_mean = 'eta not identifiable, posterior mean not calculated'
  }
  
  ##### Get out predicted mean for Y
  Y_mean = matrix(0, nrow=dim(Lambda_save)[1], ncol=dim(eta_save)[2])
  for(i in 1:nsamps_save){ Y_mean = Y_mean + Lambda_save[,,i] %*% eta_save[,,i] }
  Y_mean = Y_mean/(nsamps_save)
  
  ##### Save everything in a list and return said list.
  res = list("Lambda_save"=Lambda_save, "eta_save"=eta_save, 
             "sigsq_y_save"=sigsq_y_save, "Y_save"=Y_save, "dvec_unique"=dvec_unique, 
             "dvec_unique_original"=dvec_unique_original, "Y_mean"=Y_mean, "Lambda_mean"=Lambda_mean, "eta_mean"=eta_mean,
             "l"=l, "covDD"=covDD, "Y"=Y, "norm_rescale"=norm_rescale, "D"=D, "K"=K)
  return(res)
  
}