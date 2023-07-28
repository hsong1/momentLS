#' Tune delta for momentLS estimators
#' 
#' @description
#' Tune based on a modification of an adaptive bandwidth selection method by Politis <https://doi.org/10.1080/10485250310001604659>. See also our paper <https://arxiv.org/abs/2207.12705>.
#'@export
tune_delta = function(x, 
                      nSplits = 1,
                      method= "ft", 
                      seq_type = "even",
                      const = 1,
                      c_M_const = 0.1){
  
  nSplits = round(nSplits)
  if(nSplits>1){
    ind = cut(1:length(x),breaks = nSplits, labels = FALSE)
  }else{
    ind = rep(1,length(x))
  }
  
  d_all = sapply(1:nSplits, function(j){
    r = autocov(x[ind==j]); m = length(r)
    # compute trunctation pt for each split
    c(m, 
      truncation_point(r = r,method = method, seq_type = seq_type,
                       c1 = 1,c2 = c_M_const*sqrt(log(m))))
  })
  
  d = apply(d_all,1,mean)
  m = d[1] # average chain length per split
  m_trunc_pt = d[2]; # average truncation point
  
  dhat = max(1-exp(-0.5*( log(m)/m_trunc_pt)), 1/m)
  d_all = data.frame(t(d_all))
  colnames(d_all) = c("m","trunc_pt")
  return(list(delta = const*dhat, d_all = d_all))
  
}

#'@export
truncation_point = function(r, 
                            method= "ft", 
                            seq_type = "even",
                            c1=1,c2=1){
  
  # compute empirical autocovariance
  M = length(r)
  method = match.arg(method, c("init","ft")) 
  seq_type = match.arg(seq_type, c("raw","even"))
  rc = r[-1]/r[1]
  
  if(seq_type=="even"){
    # check even lags 
    # trunc_pt: k such that 
    # method == "ft": rho_M(k+2)<c_M sqrt(logM/M) or
    # method == "init": rho_M(k+2)<0
    
    
    if(method=="init"){
      index = seq(from=2,to=M-1,by=2)
      k = min(which(rc[index]<0)) # rc[index[k-1]] > 0; rc[index[k]] <0
      if(!is.infinite(k)){trunc_pt = index[k-1]}else{trunc_pt = max(index)}
      
    }else{
      trunc_pt = choose_mhat(rc=rc, c1=c1,c2=c2,M=M, incr=2)
    }
    
  }else if(seq_type=="raw"){
    
    # check all lags 
    # trunc_pt: k such that 
    # method == "ft": rho_M(k+1)<c_M sqrt(logM/M) 
    
    if(method=="init"){stop("if method==init, seq_type needs to be even")}
    trunc_pt = choose_mhat(rc=rc, c1=c1,c2=c2,M=M, incr=1)
  }
  
  # dhat  = max(1-exp(-0.5*( log(M)/trunc_pt)), 1/M)
  
  return(trunc_pt)
}


#'@export
choose_mhat = function(rc,c1,c2,M, incr = 1){
  # return first m such that
  # all(abs(rc[mhat+ind])<= c2*sqrt(log(M)/M))
  # when c1=1, first time when rc[mhat+incr]<= c2*sqrt(log(M)/M)
  cond = FALSE; mhat = 0
  max_ind = max(incr, c1)
  ind = seq(from=incr, to = max_ind, by=incr)
  while(1){
    cond = all(abs(rc[mhat+ind])<= c2*sqrt(log(M)/M))
    if(cond||is.na(cond)){break}else{mhat = mhat+incr}
  }
  mhat = max(1,mhat)
  return(mhat)
}

#'@export
compute_delta_hat = function(mhat,m){
  1-exp(-0.5*( log(m)/mhat))
}