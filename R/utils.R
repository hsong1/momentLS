#' compute XtX
#'@export
makeXtX <-
function(x, s_x = NULL){
  XtX = outer(x,x)
  XtX = (1+XtX)/(1-XtX)
  if(is.null(s_x)){
    s_x = sqrt((1+x^2)/(1-x^2))
  }
  XtX = XtX / outer(s_x, s_x)
}
 
#' compute Xtr
#'@export
computeXtr = function(x, r, s_x = NULL){
  Xtr=rep(0,length(x))
  n = length(r)
  exponents=seq(0,n-1,by=1)
  
  for (i in 1:length(x)){
    Xtr[i]=2*sum((x[i]^exponents)*r)-r[1]
  }
  
  if(is.null(s_x)){
    s_x = sqrt((1+x^2)/(1-x^2))
  }
  Xtr = Xtr / s_x
  return(Xtr)
}

#' create a grid of alpha values
#'@export
makeGrid <-
  function(upper_threshold,
           nX,
           cm = FALSE,
           scale = c("log", "equidist")) {
    
    scale = match.arg(scale, choices = c("log","equidist"))
    xTemp=NULL
    
    if (scale=="log"){
      if(nX%%2!=1){nX=nX+1; warning("provide odd nX when scale=log")}
      if(!cm){nX2 = (nX+1)/2}else{nX2=nX}
      x0=seq(0,log(1-upper_threshold),length.out=nX2)
      x0=1-exp(x0)
      
      if (cm){
        xTemp=x0
      }else{
        xTemp=c(-x0[length(x0):2],x0)
      }
      
    }else{
      if (cm){
        xTemp=seq(0,upper_threshold,length.out=nX)
      }else{
        xTemp=seq(-upper_threshold,upper_threshold,length.out=nX)
      }
      
    }
    return(xTemp)
  }



#' compute asymptotic variance of a moment sequence (weights, support)
#'@export
asympVariance<-function(weights,support){
  avar=sum(weights*(1+support)/(1-support))
  return(avar)
}



findIndices = function(support, alphaGrid){
  if(length(support)==0){ind=integer(0)}else{
    ind =sapply(1:length(support), function(j) which(abs(alphaGrid  - support[j])<1e-10) )
    if(length(ind)!= length(support))stop("not all support points are in alphaGrid")
  }
  
  return(ind)
}

#' compute first M moments of a moment sequence (weights, support)
#'@export
computeMoments = function(support,weights,M=100){
  moments=rep(0,M)
  for (i in 1:M){
    for (j in 1:length(support)){
      moments[i]=moments[i]+support[j]^(i-1)*weights[j]
    }
  }
  return(moments)
}
 



