get_perfectCorInds = function(Xmat, tol=1e-12){
  rm_inds = "none"
  Xmat = as.data.frame(Xmat)
  colnames(Xmat) = 1:ncol(Xmat)
  z = cor(Xmat)
  z[lower.tri(z,diag=TRUE)]=NA  # Prepare to drop duplicates and meaningless information
  z=as.data.frame(as.table(z))  # Turn into a 3-column table
  z=na.omit(z)                  # Get rid of the junk we flagged above
  z=z[order(-abs(z$Freq)),]     # Sort by highest abs correlation (whether +ve or -ve)
  bigz_inds = which( abs(z$Freq) > max(1-tol, 0.99) )
  if( length(bigz_inds)>0 ){
    bigz_inds_new = bigz_inds[ sapply( bigz_inds, function(i) all.equal(1, z$Freq[i], tolerance=tol)==TRUE ) ]
    if( length(bigz_inds_new)>0 ){
      rm_inds = sapply( z$Var1[bigz_inds_new], function(x) which(colnames(Xmat)==x) )
    }
  }
  return(rm_inds)
}
