#----------------------------------------------------------------------------
# Compute Simultaneous Credible Bands
#
# Compute (1-alpha)\% credible BANDS for a function based on MCMC samples using Crainiceanu et al. (2007). 
# Function from fosr package (https://github.com/drkowal/fosr).
#
# @param sampFuns \code{Nsims x m} matrix of \code{Nsims} MCMC samples and \code{m} points along the curve
# @param alpha confidence level
#
# @return \code{m x 2} matrix of credible bands; the first column is the lower band, the second is the upper band
#
# @note The input needs not be curves: the simultaneous credible "bands" may be computed
# for vectors. The resulting credible intervals will provide joint coverage at the (1-alpha)%
# level across all components of the vector.
#
get_credBands = function(sampFuns, alpha = .05){
  
  N = nrow(sampFuns); m = ncol(sampFuns)
  
  # Compute pointwise mean and SD of f(x):
  Efx = colMeans(sampFuns); SDfx = apply(sampFuns, 2, sd)
  
  # Compute standardized absolute deviation:
  Standfx = abs(sampFuns - tcrossprod(rep(1, N), Efx))/tcrossprod(rep(1, N), SDfx)
  
  # And the maximum:
  Maxfx = apply(Standfx, 1, max)
  
  # Compute the (1-alpha) sample quantile:
  Malpha = quantile(Maxfx, 1-alpha)
  
  # Finally, store the bands in a (m x 2) matrix of (lower, upper)
  cbind(Efx - Malpha*SDfx, Efx + Malpha*SDfx)
}
#----------------------------------------------------------------------------
# Compute Simultaneous Band Scores (SimBaS)
#
# Compute simultaneous band scores (SimBaS) from Meyer et al. (2015, Biometrics).
# SimBaS uses MC(MC) simulations of a function of interest to compute the minimum
# alpha such that the joint credible bands at the alpha level do not include zero.
# This quantity is computed for each grid point (or observation point) in the domain
# of the function. Function from fosr package (https://github.com/drkowal/fosr).
#
# @param sampFuns \code{Nsims x m} matrix of \code{Nsims} MCMC samples and \code{m} points along the curve
#
# @return \code{m x 1} vector of simBaS
#
# @note The input needs not be curves: the simBaS may be computed
# for vectors to achieve a multiplicity adjustment.
#
# @note The minimum of the returned value, \code{PsimBaS_t},
# over the domain \code{t} is the Global Bayesian P-Value (GBPV) for testing
# whether the function is zero everywhere.
#
get_simBaSc = function(sampFuns){
  
  N = nrow(sampFuns); m = ncol(sampFuns)
  
  # Compute pointwise mean and SD of f(x):
  Efx = colMeans(sampFuns); SDfx = apply(sampFuns, 2, sd)
  
  # Compute standardized absolute deviation:
  Standfx = abs(sampFuns - tcrossprod(rep(1, N), Efx))/tcrossprod(rep(1, N), SDfx)
  
  # And the maximum:
  Maxfx = apply(Standfx, 1, max)
  
  # And now compute the SimBaS scores:
  PsimBaS_t = rowMeans(sapply(Maxfx, function(x) abs(Efx)/SDfx <= x))
  
  PsimBaS_t
}