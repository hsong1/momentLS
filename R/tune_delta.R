choose_mhat = function(rc, c1,c2, incr = 1){
  cond = FALSE; mhat = 0
  M = length(rc)
  K_N=c1*round(sqrt(log(M)))
  ind = seq(from=incr, to = K_N, by=incr)
  while(1){
    cond = all(abs(rc[mhat+ind])<= c2*sqrt(log(M)/M))
    if(cond){break}else{mhat = mhat+incr}
  }
  mhat = max(1,mhat)
  return(mhat)
}

#'@export
tune_delta = function(r, 
                      method= c("init","ft"), 
                      seq_type = c("raw","even","average"),
                      const = .5,
                      c1=1,c2=1){
  
  # compute empirical autocovariance
  M = length(r)
  method = match.arg(method, c("init","ft")) 
  seq_type = match.arg(seq_type, c("raw","even","average"))
  
  if(seq_type =="average"){
    Gamma=r[seq(from=1,to=(M-1),by=2)]+r[seq(from=2,to=M,by=2)]
    
    if(method=="init"){
      ind_G = min(which(Gamma<0))-1 # last j s.t. Gamma[j] > 0, i.e., Gamma(0),Gamma(1),...,Gamma(indG-1) > 0
      # Gamma(ind_G-1) = gamma(2*(ind_G-1))+gamma(2*(ind_G-1) +1); needs to include 1,..., 2(ind_G-1)+2
    }else{
      rc = Gamma[-1] / Gamma[1]
      ind_G = choose_mhat(rc = rc,c1 = c1,c2 = c2,incr = 1)
    }
    trunc_pt = 2*(ind_G-1) +2
  }else if(seq_type=="even"){
    if(method=="init"){
      trunc_pt = min(which(r[seq(from=2,to=M,by=2)] <0 ))-2
    }else{
      rc = r[-1]/r[1]
      trunc_pt = choose_mhat(rc=rc, c1=c1,c2=c2,incr=2)
    }
    
  }else if(seq_type=="raw"){
    if(method=="init"){stop("if method==init, seq_type needs to be average or even")}
    rc = r[-1]/r[1]
    trunc_pt = choose_mhat(rc=rc, c1=c1,c2=c2,incr=1)
  }
  
  dhat  = max(1-exp(-0.5*( log(M)/trunc_pt)), 1/M)
  
  return(c(trunc_pt=trunc_pt, dhat=const*dhat))
}
