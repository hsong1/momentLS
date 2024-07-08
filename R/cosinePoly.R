cosinePoly = function(coef=c(0,1), side=2){
  a <- as.numeric(coef)
  while ((la <- length(a)) > 1 && a[la] == 0) a <- a[-la]
  structure(list(coef=a, side =side), class="cosinePoly")
}
#x = cosinePoly()

#' @export
#' @method print cosinePoly
print.cosinePoly=function(x,...,quote=FALSE){
  a = round(x$coef,3)
  la = length(a)
  if(la>6){a = a[1:6]}
  print_cosine = paste0("*cos(",0:(length(a)-1),"w)")
  
  msg = paste0(a,print_cosine,collapse = " + ")
  if(x$side==2){
    print_cosine2 = paste0("*cos(-",1:(length(a)-1),"w)")
    msg =  paste(msg, paste0(a[-1],print_cosine2,collapse = " + "),sep = " + ")
  }
  if(la>6){msg=paste(msg,"+ ...",sep = "")}
  cat(x$side,"-sided cosine polynomial:\np(w)=",msg,"\n",sep = "")
}

#print(x)


eval=function(w,x){
  # returns cosine polynomial x evaluated at w
  if(!is(x, "cosinePoly")) stop("x should be cosinePoly object")
  p = x$coef
  L = length(p)-1
  # if side = 1; sum_{k=0}^{L} p(k) cos(wk)
  # if side = 2; sum_{k=-L}^{L} p(k) cos(wk)
  ret = numeric(length = length(w))
  for(i in 1:length(ret)){
    if(x$side==2){
      out = p*cos( w[i]*(0:L))
      ret[i] = 2*sum(out)-p[1]
    }else{
      out = p*cos( w[i]*(0:L) )
      ret[i] = sum(out)
    }
  }
  
  return(ret)
}
x1=cosinePoly(side = 2)
x=cosinePoly(side=1)

convert_twosides <-function(x){
  a = x$coef
  if(length(a)>1){x$coef = c(a[1],a[-1]/2)}
  x$side = 2
  return(x)
}

multiply<-function(x1,x2){
  
  if(!is(x1, "cosinePoly") || !is(x2, "cosinePoly")) stop("x1 and x2 should be cosinePoly object")
  
  if(x1$side==1 & x2$side==1){side=1
  }else if(x1$side==2 & x2$side==2){side=2
  }else if(x1$side != x2$side){
    if(x1$side==1){x1 = convert_twosides(x1)}
    if(x2$side==1){x2 = convert_twosides(x2)}
    side= 2
  }
  p1 = x1$coef; p2 = x2$coef
  
  M1=length(p1)
  M2=length(p2)
  
  M=2*(M1+M2-1)-1
  p1aug=c(p1,rep(0,M-M1))
  p2aug=c(p2,rep(0,M-M2))
  
  if (side==1){
    
    #multiply 
    # x1 = \sum_{k=0}^{M1-1}a_k cos(kw) with 
    # x2 = \sum_{k=0}^{M2-1}b_k cos(kw)
    # x1*x2 = \sum_{k=0}^{M1+M2-2}c_k cos(kw)}
    #return c_k, k=0,...,M1+M2-2
    
    
    phi1=(fft(p1aug)+fft(p1aug,inverse=TRUE))/2
    phi2=(fft(p2aug)+fft(p2aug,inverse=TRUE))/2
    phi=phi1*phi2
    
    p3=Re(fft(phi,inverse = TRUE))[1:(M1+M2-1)]/M
    if (length(p3)>1){p3[2:length(p3)]=2*p3[2:length(p3)]}
    
    x3 = cosinePoly(coef = p3,side = side)
    return(x3)
  }else{
    
    #assuming symmetric a_k, b_k sequences such that a_k=a_{-k}, b_k=b_{-k}
    #multiply 
    # x1 = \sum_{k=-(M1-1)}^{M1-1}a_k cos(kw)
    # x2 = \sum_{k=-(M2-1)}^{M2-1}b_k cos(kw)
    # x1*x2 = \sum_{k=-(M1+M2-2)}^{M1+M2-2}c_k cos(kw)
    #return c_k, k=0,1...,M1+M2-2; c_{-k}=c_k
    
    phi1=(fft(p1aug)+fft(p1aug,inverse=TRUE))-p1aug[1]
    phi2=(fft(p2aug)+fft(p2aug,inverse=TRUE))-p2aug[1]
    phi=phi1*phi2
    p3=Re(fft(phi,inverse = TRUE))[1:(M1+M2-1)]/M
    x3 = cosinePoly(coef = p3,side = side)
    return(x3)
  }
}


cosinePoly_from_beta = function(beta){
  # return a_k s.t. 
  # prod_{i=1}^n (1+beta_i^2 - 2*beta_i cos(w)) = \sum_{k=0}^{n} a_k cow(wk)
  mat = cbind(1+beta^2, -2*beta)
  cp_list = lapply(1:nrow(mat), function(i) cosinePoly(coef = mat[i,],side = 1))
  p = cp_list[[1]]
  if(nrow(mat)>1){
    for(i in 2:nrow(mat)){
      p = multiply(p,cp_list[[i]])
    }
  }
  
  return(p)
}

