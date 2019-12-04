sample_lind = function(Lambda, sigsqLam, covDDinv_all, logdetCovDD_all){
  # Sample from discrete distribution, return index.
  num_ls = length(logdetCovDD_all)
  log_prob = rep(NA, num_ls)
  K = ncol(Lambda)
  precLam = 1/sigsqLam
  for(i in 1:num_ls){
    covDD_inv = covDDinv_all[,,i]
    log_prob[i] = 0.0
    for(k in 1:K){
      Lamk = Lambda[,k,drop=F]
      quad_tmp = precLam[k] * t(Lamk) %*% covDD_inv %*% Lamk
      log_prob[i] = log_prob[i] - 0.5*(logdetCovDD_all[i] + quad_tmp)
    }
  }
  maxlg = max(log_prob)
  # https://stats.stackexchange.com/questions/66616/converting-normalizing-very-small-likelihood-values-to-probability
  probls = exp(log_prob - maxlg)
  probls = probls/sum(probls)
  which_ls = sample(x=(1:num_ls), 1, replace=T, prob=probls)
  return(which_ls)
}