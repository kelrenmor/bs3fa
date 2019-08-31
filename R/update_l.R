update_l = function(l, l_diff, dvec_unique, niter_max, Y, Lambda, eta, 
                    alpha_lam_tmp, sigsq_y_vec, obs_Y){
  l_new = l + l_diff
  still_adding = T
  dvec_mat = matrix(dvec_unique)
  while( still_adding ){
    covDD = get_covDD(dvec_mat, l_new);
    niter=1; success = T
    while( success & niter<niter_max ){
      success = F
      err = tryCatch({
        sample_Lambda(Y, Lambda, eta, alpha_lam_tmp, sigsq_y_vec, covDD, obs_Y)
        success = T
      }, error = function(e) {}
      )
      niter = niter + 1
    }
    if( niter==niter_max ){
      l_new = l_new + l_diff
    } else{
      l_new = l_new - l_diff
      l_new = (l + l_new)/2
      still_adding = F
    }
  }
  return(l_new)
}
