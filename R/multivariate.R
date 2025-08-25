mLSE_ij = function(fits,a,b,method){
  len = sapply(fits, function(x) {length(x$support)})
  support = unlist(lapply(fits, function(x) {x$support}))
  weights_all = lapply(fits, function(x) {x$weights})
  if(method==1){
    weights = 
      c(weights_all[[1]],rep(0,len[2]+len[3]))-
      a^2*c(rep(0,len[1]), weights_all[[2]],rep(0,len[3])) - 
      b^2*c(rep(0,len[1]+len[2]), weights_all[[3]])
    weights = weights/(2*a*b)
  }else if(method==2){
    weights = 
      c(weights_all[[1]],rep(0,len[2]))-
      c(rep(0,len[1]), weights_all[[2]])
    
    weights = weights/(4*a*b)
  }
  
  # order
  weights = weights[order(support)]
  support = support[order(support)]
  list(support=support, weights=weights)
}

# 
# M=10000; rho = 0.4
# ch = generateAR1Chain(M = M, rho = rho)
# g = function(x){cbind(x,x^2)}
# x= matrix(g(ch$x),ncol=2)
# delta = NULL;
# estimate_delta_ij = FALSE;
# init =NULL;
# tol = 10^-8; # alg. tolerance
# maxit = 1000;
# gridControl = list(
#   nX = 1001, # number of grid points
#   cm = FALSE, # completely monotone?
#   scale = "log"
# );
# scale_chains = TRUE;
# type = "mLSE";
# method = 2

#'@export
mtvMLSE = function(
    x, # input data 
    delta = NULL,
    estimate_delta_ij = FALSE,
    init =NULL,
    tol = 10^-8, # alg. tolerance
    maxit = 1000,
    gridControl = list(
      nX = 1001, # number of grid points
      cm = FALSE, # completely monotone?
      scale = "log"
    ),
    scale_chains = TRUE,
    type = "mLSE",
    method = 2,
    psd_adjust = TRUE
){  
  
  
  d = dim(x)[2]; M = dim(x)[1]
  sigmaEst = matrix(NA,nrow=d,ncol=d)
  deltaMat = matrix(NA,nrow=d,ncol=d)
  if(!method%in%c(1,2)){stop("method has to be 1 or 2")}
  if(method==2){deltaMat2 = matrix(NA,nrow=d,ncol=d)}else{deltaMat2 = NULL}
  mLSEs = vector("list",length = d*d)
  dim(mLSEs)=c(d,d)
  
  # estimate diagonal matrices
  for(i in 1:d){
    
    if(type=="mLSE"){
      ri = autocov(x[,i])
      if(is.null(delta)){
        deltaMat[i,i] = tune_delta(x = x[,i],nSplits = 5,c_M_const = 0,fcn.summarize = "mean")$delta * 0.8
      }else{
        deltaMat[i,i] = delta
      }
      
      fit=SR1(r=ri, delta=deltaMat[i,i],init = init,tol = tol,maxit = maxit,gridControl = gridControl)
      sigmaEst[i,i]= asympVariance(weights = fit$weights, support = fit$support)
      mLSEs[i,i][[1]] = list(support=fit$support, weights=fit$weights)
      
    }else if(type=="init.pos"){
      sigmaEst[i,i] = mcmc::initseq(x[,i])$var.pos
      fit = NULL
    }else if(type=="init.conv"){
      sigmaEst[i,i] = mcmc::initseq(x[,i])$var.con
      fit = NULL
    }
    
    
  }
  
  i=1;j=2
  for(i in 1:(d-1)){
    for(j in (i+1):d){
      
      if(scale_chains){
        # a = 1/sqrt(sigmaEst[i,i]); b = 1/sqrt(sigmaEst[j,j])
        a = 1/sqrt(autocov(x[,i])[1]); b= 1/sqrt(autocov(x[,j])[1])
      }else{
        a =1; b=1
      }
      
      cur_ch = a*x[,i]+b*x[,j]
      if(method==2) {cur_ch2 = a*x[,i]-b*x[,j]}
      
      if(type=="mLSE"){
        
        rij=autocov(cur_ch)
        if(method==2){rij2 = autocov(cur_ch2)}
        
        # tune delta if not null
        if(!is.null(delta)){
          deltaMat[i,j] = delta
          if(method==2){deltaMat2[i,j] = delta}
          
        }else{
          
          if(estimate_delta_ij){
            deltaMat[i,j] =  tune_delta(x = cur_ch, nSplits = 5,c_M_const = 0,fcn.summarize = "mean")$delta * 0.8
            if(method==2){deltaMat2[i,j] =  tune_delta(x = cur_ch2, nSplits = 5,c_M_const = 0,fcn.summarize = "mean")$delta * 0.8}
            
          }else{
            deltaMat[i,j] = min(deltaMat[i,i],deltaMat[j,j])
            if(method==2){deltaMat2[i,j] = deltaMat[i,j]}
          }
        }
        
        
        # fit momentLS for cur_ch = ax[,i]+bx[,j]
        fit=SR1(r=rij, delta=deltaMat[i,j], init = init,tol = tol,maxit = maxit,gridControl = gridControl)
        vij = asympVariance(weights = fit$weights, support = fit$support)
        
        if(method==2){
          fit2=SR1(r=autocov(cur_ch2), delta=deltaMat2[i,j], init = init,tol = tol,maxit = maxit,gridControl = gridControl)
          vij2 = asympVariance(weights = fit2$weights, support = fit2$support)
          
        }
        
      }else if(type=="init.pos"){
        vij = mcmc::initseq(cur_ch)$var.pos
        if(method==2){vij2 = mcmc::initseq(cur_ch2)$var.pos}
        fit = NULL; fit2 = NULL
      }else if(type=="init.conv"){
        vij = mcmc::initseq(cur_ch)$var.con
        if(method==2){vij2 = mcmc::initseq(cur_ch2)$var.con}
        fit = NULL; fit2 = NULL
      }
      
      if(method==1){
        # Sigma_ij = (aVar(fit)-a^2sigmaEst[i,i] - b^2sigmaEst[j,j])/2ab
        sigmaEst[i,j] = (vij-a^2*sigmaEst[i,i]-b^2*sigmaEst[j,j])/(2*a*b)
        fits = list(fit,mLSEs[i,i][[1]],mLSEs[j,j][[1]])
        
        
      }else if(method==2){
        sigmaEst[i,j] = (vij-vij2)/(4*a*b)
        fits = list(fit,fit2)
      }
      
      if(type=="mLSE"){
        mLSEs[i,j][[1]] = mLSE_ij(fits,a,b,method=method)
        mLSEs[j,i][[1]] = mLSE_ij(fits,a,b,method=method)
      }
      sigmaEst[j,i] = sigmaEst[i,j]
    }
  }
  
  evd = eigen(sigmaEst)
  is.psd = min(evd$values)>=0
  
  if(type=="mLSE"& psd_adjust){
    
    if(!is.psd){
      # cov_adjust = TRUE
      mtvMLSE_psd = 
        mtvMLSE2(x = x,dir = evd$vectors,delta = min(diag(deltaMat)),
                 init=init,tol=tol,maxit=maxit, gridControl = gridControl,type=type)
      sigmaEst.f = mtvMLSE_psd$sigmaEst
    }else{
      sigmaEst.f = sigmaEst
      sigmaEst = NULL
      mtvMLSE_psd = NULL
    }
    
  }else{
    sigmaEst.f = sigmaEst
    sigmaEst = NULL
    mtvMLSE_psd = NULL
  }  
  
  
  return(list(cov = sigmaEst.f, 
              mLSEs = mLSEs, 
              deltaMat =deltaMat, 
              deltaMat2 = deltaMat2,
              is.psd = is.psd,
              cov.before.adj = sigmaEst, 
              mtvMLSE_psd = mtvMLSE_psd))
}


mtvMLSE2 = function(
    x, # input data
    dir,
    delta = NULL,
    init = NULL,
    tol = 10^-8, # alg. tolerance
    maxit = 1000,
    gridControl = list(
      nX = 1001, # number of grid points
      cm = FALSE, # completely monotone?
      scale = "log"
    ),
    type = "mLSE"
){
  d = dim(x)[2]; M = dim(x)[1]
  sigmaEst = matrix(0,nrow=d,ncol=d)
  
  xD = as.matrix(x)%*%dir
  mLSEs = vector("list",length = d)
  deltaVec = c()
  
  for(i in 1:d){
    if(type=="mLSE"){
      ri = autocov(xD[,i])
      if(is.null(delta)){
        deltaVec[i] = tune_delta(x = xD[,i],nSplits = 5,c_M_const = 0,fcn.summarize = "mean")$delta * 0.8
      }else{
        deltaVec[i] = delta
      }
      fit=SR1(r=ri, delta = deltaVec[i], init = init,tol = tol,maxit = maxit,gridControl = gridControl)
      sigmaEst[i,i]= asympVariance(weights = fit$weights, support = fit$support)
    }else if(type=="init.pos"){
      sigmaEst[i,i] = mcmc::initseq(xD[,i])$var.pos
    }else if(type=="init.conv"){
      sigmaEst[i,i] = mcmc::initseq(xD[,i])$var.con
    }
    mLSEs[[i]] = fit
  }
  
  
  return(list(sigmaEst = dir%*%diag(x = diag(sigmaEst))%*%t(dir),
              sigmaEst0 = sigmaEst,
              deltaVec=deltaVec,
              mLSEs=mLSEs))
}



mtvAsympVariance = function(coefs,support){
  if(is.null(dim(coefs))){coefs= matrix(coefs,ncol=1)}
  Reduce("+", lapply(1:nrow(coefs), function(k) outer(coefs[k,],coefs[k,])*(1+support[k])/(1-support[k]) ))
}
