# Apply CAPout output to a list

# ARGUMENTS: l: list of matrices to apply to permutations to
#            p: permutation / sign output from  CAPout

applier = function(l, p){
  return(mapply(aplr, l, p, SIMPLIFY = FALSE))
}

# rotate factors to enforce identifiability
# based on varimax rotation 

# ARGUMENTS: lambda: file path to factor matrix sample list (or just the list);
#            tolerance: rotation algorithm stopping tolerance;
#            ncores: number of cores
#            normalize: logical. whether to normalize lambda samples
#            file: logical. whether lambda was passed directly or as an Rds
#            ret: character, whether to return rotated loadings themselves ('loadings'),
#                 the rotation matrix ('rotmat'), or both ('both')

mcrotfact = function(lambda, 
                     tolerance = 1e-5, ncores = 1,
                     normalize = FALSE, file = TRUE, ret='loadings'){
  library(parallel)
  if(file) lambda = readRDS(lambdafile) # read in factor samples
  n = length(lambda)           # initialize shared attributes
  if(normalize) lambda = mclapply(lambda, scale,
                                  mc.cores = ncores)
  
  # Varimax rotations
  
  if(ret=='loadings'){
    Lr = mclapply(lambda,
                  function(lam) as(varimax(lam, eps = tolerance)[["loadings"]],
                                   "matrix"), 
                  mc.cores = ncores)
    LrMean = Reduce("+", Lr) / n
    class(LrMean) = "matrix"
    return(list(samples = Lr, mean = LrMean))
  }
  
  if(ret=='rotmat'){
    rots = mclapply(lambda,
                    function(lam) as(varimax(lam, eps = tolerance)[["rotmat"]],
                                     "matrix"), 
                    mc.cores = ncores)
    return(list(rots = rots))
  }
  
  if(ret=='both'){
    both = mclapply(lambda,
                    function(lam) varimax(lam, eps = tolerance), 
                    mc.cores = ncores)
    Lr = mclapply(lapply(both, '[[', 1), function(lam) as(lam, 'matrix'))
    rot = mclapply(lapply(both, '[[', 2), function(ro) as(ro, 'matrix'))
    LrMean = Reduce("+", Lr) / n
    class(LrMean) = "matrix"
    return(list(samples = Lr, mean = LrMean, rots=rot))
  }
  
}


# inverse rotate scores based on loading rotation
# because lam eta = lam rot rot^T eta
# To see how this works, run the code below!
# n=5; k=2; p=3
# Lam = matrix(rnorm(n*k),nrow=n,ncol=k)
# eta = matrix(rnorm(k*p),nrow=k,ncol=p)
# # The following should be equal:
# Lam %*% eta
# as(varimax(Lam, eps=1e-5)[["loadings"]],"matrix") %*% t(varimax(Lam, eps=1e-5)$rotmat) %*% eta


# ARGUMENTS: eta: score matrix sample list;
#            rots: rotation matrices from varimax on Lambda;

invrot = function(eta, rots){
  library(parallel)
  n=length(eta)
  trots = lapply(rots,t)
  eta_rot = Map('%*%',trots,eta)
  return(eta_rot)
}

# function to permute in minimal-switch order
# will output repeated orders for i > k!
