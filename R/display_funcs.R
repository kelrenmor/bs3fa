plot_matrix <- function(A, type="other", tit=NULL, include_legend=T){
  library(reshape2)
  library(ggplot2)
  longA = melt(A)
  if( type=="Lambda" ){
    xLab="k"; yLab="d"
    if(is.null(tit)){tit = expression(Lambda~entries)}
  } else if( type=="Theta" ){
    xLab="k"; yLab="s"
    if(is.null(tit)){tit = expression(Theta~entries)}
  } else{
    xLab="k"; yLab="p"
    if(is.null(tit)){tit = "Matrix entries"}
  }
  if(include_legend){
  ggplot(longA, aes(x = Var2, y = Var1)) + 
    geom_tile(aes(fill=value), colour="grey20") + 
    scale_fill_gradient2(low = "#3d52bf", high = "#33961b", mid = "white") + 
    theme(panel.grid.major = element_blank(), panel.border = element_blank(), 
          axis.text = element_blank(), panel.background = element_blank(), 
          axis.ticks = element_blank(), legend.title = element_text(), 
          plot.title = element_text(hjust = 0.5)) + 
    xlab(xLab) + ylab(yLab) + scale_y_reverse() + 
    ggtitle(tit)
  } else{
    ggplot(longA, aes(x = Var2, y = Var1)) + 
      geom_tile(aes(fill=value), colour="grey20") + 
      scale_fill_gradient2(low = "#3d52bf", high = "#33961b", mid = "white") + 
      theme(panel.grid.major = element_blank(), panel.border = element_blank(), 
            axis.text = element_blank(), panel.background = element_blank(), 
            axis.ticks = element_blank(), legend.position = "none", 
            plot.title = element_text(hjust = 0.5)) + 
      xlab(xLab) + ylab(yLab) + scale_y_reverse() + 
      ggtitle(tit)
  }
}

plot_Lambda_true <- function(Lambda, doses=1:nrow(Lambda)/nrow(Lambda), inds=1:3){
  library(gridExtra)
  library(ggplot2)
  Lambda = as.data.frame(Lambda)
  xLim = range(doses)
  yLim = range(Lambda[,inds])
  if(length(inds)>3){stop("inds must be of length <= 3")}
  if(length(inds)>=1){
    p1 = qplot(doses, Lambda[,inds[1]], geom="line") + xlim(xLim) + ylim(yLim) + 
      theme_minimal() + theme_bw() + xlab("dose") + ylab("loading") + 
      ggtitle(paste("Loading",inds[1]))
  }
  if(length(inds)>=2){
    p2 = qplot(doses, Lambda[,inds[2]], geom="line") + xlim(xLim) + ylim(yLim) +
      theme_minimal() + theme_bw() + xlab("dose") + ylab("loading") + 
      ggtitle(paste("Loading",inds[2]))
  }
  if(length(inds)>=3){
    p3 = qplot(doses, Lambda[,inds[3]], geom="line") + xlim(xLim) + ylim(yLim) +
      theme_minimal() + theme_bw() + xlab("dose") + ylab("loading") + 
      ggtitle(paste("Loading",inds[3]))
  }
  if(length(inds)==1){
    grid.arrange(p1, nrow = 1)
  }else if(length(inds)==2){
    grid.arrange(p1, p2, nrow = 1)
  }else{
    grid.arrange(p1, p2, p3, nrow = 1)
  }
}

plot_Lambda <- function(Lambda_low, Lambda, Lambda_upp, Lambda_true, k, dose=1:nrow(Lambda)/nrow(Lambda)){
  library(reshape2)
  library(ggplot2)
  df = as.data.frame(cbind(dose, Lambda_low[,k], Lambda[,k], Lambda_upp[,k], Lambda_true[,k]))
  colnames(df) = c("dose", "ll", "est", "ul", "truth")
  
  ggplot(data = df, aes(dose, est)) +
    geom_ribbon(aes(x=dose, ymax=ll, ymin=ul), fill="grey", alpha=.5) +
    geom_line(aes(y = ul), lty = 3, colour = 'black') +
    geom_line(aes(y = ll), lty = 3, colour = 'black')+
    geom_line(lty = 2) + geom_line(aes(y = truth), colour = 'red') +
    ylab("loading") + ggtitle(paste("Column",k,"of loadings matrix for Y")) + 
    theme_minimal() + theme_bw()
}

plot_Lambda_mod_tru_k <- function(df, k, ylab_head, title_head, xLim, yLim){
  if(!is.null(title_head)){
    tmp_title = paste(title_head,k)
  } else{tmp_title=''}
  ggplot(data = df, aes(dose, est)) +
    geom_ribbon(aes(x=dose, ymax=ll, ymin=ul), fill="grey", alpha=.5) +
    geom_line(aes(y = ul), lty = 3, colour = 'black') +
    geom_line(aes(y = ll), lty = 3, colour = 'black')+
    geom_line(lty = 2) + geom_line(aes(y = truth), colour = 'red') +
    ylab(ylab_head) + ggtitle(tmp_title) + ylim(yLim) + xlim(xLim) + 
    theme_minimal() + theme_bw()
}

plot_Lambda_mod_tru <- function(Lambda_low, Lambda, Lambda_upp, Lambda_true, 
                                doses=1:nrow(Lambda)/nrow(Lambda), inds=1:3,
                                ylab_head="loading", title_head="Loading",
                                mean_curve=rep(0,length(doses))){
  library(gridExtra)
  library(ggplot2)
  xLim = range(doses)
  for(i in 1:length(inds)){
    Lambda_low[,inds[i]] = Lambda_low[,inds[i]] + mean_curve
    Lambda[,inds[i]] = Lambda[,inds[i]] + mean_curve
    Lambda_upp[,inds[i]] = Lambda_upp[,inds[i]] + mean_curve
    Lambda_true[,inds[i]] = Lambda_true[,inds[i]] + mean_curve
  }
  yLim = range(c(Lambda_low[,inds],Lambda_upp[,inds],Lambda_true[,inds]))
  if(length(inds)>3){stop("inds must be of length <= 3")}
  if(length(inds)>=1){
    df = as.data.frame(cbind(doses, Lambda_low[,inds[1]], Lambda[,inds[1]], 
                             Lambda_upp[,inds[1]], Lambda_true[,inds[1]]))
    colnames(df) = c("dose", "ll", "est", "ul", "truth")
    p1 = plot_Lambda_mod_tru_k(df, inds[1], ylab_head, title_head, xLim, yLim)
  }
  if(length(inds)>=2){
    df = as.data.frame(cbind(doses, Lambda_low[,inds[2]], Lambda[,inds[2]],
                             Lambda_upp[,inds[2]],Lambda_true[,inds[2]]))
    colnames(df) = c("dose", "ll", "est", "ul", "truth")
    p2 = plot_Lambda_mod_tru_k(df, inds[2], ylab_head, title_head, xLim, yLim)
  }
  if(length(inds)>=3){
    df = as.data.frame(cbind(doses, Lambda_low[,inds[3]], Lambda[,inds[3]], 
                             Lambda_upp[,inds[3]], Lambda_true[,inds[3]]))
    colnames(df) = c("dose", "ll", "est", "ul", "truth")
    p3 = plot_Lambda_mod_tru_k(df, inds[3], ylab_head, title_head, xLim, yLim)
  }
  if(length(inds)==1){
    grid.arrange(p1, nrow = 1)
  }else if(length(inds)==2){
    grid.arrange(p1, p2, nrow = 1)
  }else{
    grid.arrange(p1, p2, p3, nrow = 1)
  }
}

plot_data <- function(Y, true_curve, avg_dose_resp=rep(0,nrow(Y)), 
                      doses=1:nrow(Y)/nrow(Y), inds=1:3){
  library(gridExtra)
  library(ggplot2)
  Y = as.data.frame(Y)
  for(j in 1:ncol(Y)){
    Y[,j] = Y[,j]+avg_dose_resp
    true_curve[,j] = true_curve[,j]+avg_dose_resp
  }
  xLim = range(doses)
  yLim = range(c(Y[,inds],true_curve[,inds]), na.rm=T)
  p1 = qplot(doses, Y[,inds[1]]) + xlim(xLim) + ylim(yLim) + 
    theme_minimal() + theme_bw() + xlab("dose") + ylab("response") + 
    ggtitle("Chemical 1") + geom_line(aes(x=doses,y=true_curve[,inds[1]]))
  p2 = qplot(doses, Y[,inds[2]]) + xlim(xLim) + ylim(yLim) +
    theme_minimal() + theme_bw() + xlab("dose") + ylab("response") + 
    ggtitle("Chemical 2") + geom_line(aes(x=doses,y=true_curve[,inds[2]]))
  p3 = qplot(doses, Y[,inds[3]]) + xlim(xLim) + ylim(yLim) +
    theme_minimal() + theme_bw() + xlab("dose") + ylab("response") + 
    ggtitle("Chemical 3") + geom_line(aes(x=doses,y=true_curve[,inds[3]]))
  
  grid.arrange(p1, p2, p3, nrow = 1)
}
