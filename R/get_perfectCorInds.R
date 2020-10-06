get_perfectCorInds = function(Xmat, tol=1e-12){
  z = cor(Xmat)
  z[lower.tri(z,diag=TRUE)]=NA  # Prepare to drop duplicates and meaningless information
  z=as.data.frame(as.table(z))  # Turn into a 3-column table
  z=na.omit(z)                  # Get rid of the junk we flagged above
  z=z[order(-abs(z$Freq)),]     # Sort by highest abs correlation (whether +ve or -ve)
  bigz_inds = which( abs(z$Freq) > max(1-tol, 0.99) )
  bigz_inds = bigz_inds[ sapply( bigz_inds, function(i) all.equal(1, z$Freq[i], tolerance=tol)==TRUE ) ]
  rm_inds = sapply( z$Var1[bigz_inds], function(x) which(colnames(Xmat)==x) )
  return(rm_inds)
}
