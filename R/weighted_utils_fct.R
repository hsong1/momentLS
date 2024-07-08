#'@export
poissonKernel = function(rho,w) {(1-rho^2)/(1-2*cos(w)*rho + rho^2)}


#### roots and factorization ####

find_roots_mu<-function(mu){
  
  support=mu$support
  weights=mu$weights
  
  # find ai,bi,qi
  # phi(w) = a0 + sum_{i=1}^n ai / (cosw - bi) 
  # phi(w) = c_phi \prod_{i=1}^L (cosw - qi) / \prod_{i=1}^n (cosw - bi)
  # for L = n, n-1, n-2
  
  # check inputs
  if(length(support)!=length(weights)) stop("length(support)!=length(weights)")
  #if(any( weights < .Machine$double.eps^{1/2})) stop("all weights have to be positive ")
  
  # check whether support contains 0
  ind_0 = which(abs(support) < .Machine$double.eps^{1/2})
  
  if(length(ind_0)>0){
    a0 = weights[ind_0]
    support_0 = support[-ind_0]
    weights_0 = weights[-ind_0]
  }else{
    a0 = 0
    support_0 = support; weights_0 = weights
  }
  
  # n = number of non-zero support points
  n = length(support_0)
  
  ### a few special cases ###
  if(n==0){
    # case I and n_phi = 1 (single support at 0)
    out = list(a0 = a0, ai = numeric(),bi = numeric(),qi = numeric(), c_phi = a0, n=n,  weights= weights, support = support)
    return(out)
  }
  
  # if n>0 compute ai and bi
  ai = -weights_0*(1-support_0^2)/(2*support_0)
  bi = (1+support_0^2)/(2*support_0)
  ai = ai[order(bi)]
  bi = bi[order(bi)]
  
  if((n==1) & (length(ind_0)==0) & (abs(sum(ai)) > .Machine$double.eps^{1/2})){
    # case II and n_phi = 1 (single non zero support)
    out = list(a0 = a0,ai = ai,bi = bi,qi = numeric(), c_phi = ai, n=n, weights= weights, support = support)
    return(out)
  }
  
  if ((n==2) & (length(ind_0)==0) & (abs(sum(ai))<.Machine$double.eps^{1/2})){
    # case III and n_phi = 2 (two non-zero supports with opposite signs with a_1+a_2 = 0)
    out = list(a0 = a0,ai = ai,bi = bi,qi = numeric(), c_phi = NULL,n=n, weights= weights, support = support)
    return(out)
  }
  
  
  # find q(x)
  q_fcn = function(x){
    val = c()
    for(i in 1:length(x)){
      val[i] = a0+sum(ai/(x[i]-bi))  
    }
    
    return(val)
  }
  
  ### call root finding subroutine ###
  qi = find_roots_fcn(ai = ai,bi = bi,a0=a0,q_fcn = q_fcn,tol = .Machine$double.eps)
  c_phi = ifelse(a0!=0,a0,sum(ai))
  out = list(a0=a0,ai=ai,bi=bi,qi=qi, c_phi= c_phi, n=n, weights= weights, support = support)
  
  return(out)
} 


find_roots_fcn = function(ai,bi,a0,q_fcn, 
                          min_val = -1e6, max_val = 1e6, 
                          tol = .Machine$double.eps^{1/2}){
  
  # find roots of q_fcn(x) = a0 + sum_{i=1}^n ai / (x-bi)
  if(!all.equal(bi, sort(bi))){stop("ai, bi needs to be sorted")}
  
  nroots_n =  (abs(a0) > tol) # check whether a0 != 0
  nroots_n1 = (!nroots_n) & (abs(sum(ai)) > tol) # check whether sum(ai) != 0
  nroots_n2 = (!nroots_n) & (abs(sum(ai)) <= tol) 
  
  
  interval_endpoints = 
    cbind(c(min_val, bi), c(bi, max_val))
  # whether b1<0<bn
  sign_change = ifelse(sign(bi[1])*sign(bi[length(bi)])==-1, TRUE,FALSE) 
  if(sign_change){
    idx_sign_change = which((interval_endpoints[,1]<0) & (interval_endpoints[,2] >=0))
  }
  idx_root = 1:nrow(interval_endpoints)
  
  #is there a root in (-Inf,b1)?
  leftRoot=NULL
  
  #is there a root in (bn,Inf)?
  rightRoot=NULL
  
  if(nroots_n){
    ## remove 1 from idx_root
    # nroots = length(support)
    if(sign_change){
      # (-infty, b1),(b1,b2),...,(bn,infty) -> remove idx_sign_change
      idx = setdiff(idx_root, idx_sign_change)
      leftRoot=TRUE
      rightRoot=TRUE
    }else if(!sign_change){
      if(bi[length(bi)]<0){
        # (-infty,b1),...,(bn-1,bn) -> remove last
        idx = idx_root[-length(idx_root)]
        leftRoot=TRUE
        rightRoot=FALSE
      }else if(bi[1]>0){
        # (b1,b2),...,(bn-1,bn), (bn, infty) -> remove first
        idx = idx_root[-1]
        leftRoot=FALSE
        rightRoot=TRUE
      }
    }
    
  }
  
  if(nroots_n1){
    ## remove 2 from idx_root
    # nroots = length(support)-1
    if(sign_change){
      idx = setdiff(idx_root, idx_sign_change)
      if(sum(ai)>0){
        # (b1, b2),...,(bn, infty) -> remove first, idx_sign_change
        idx = idx[-1]
        leftRoot=FALSE
        rightRoot=TRUE
      }else if(sum(ai)<0){
        # (-infty,b1),...,(bn-1,bn) -> remove last, idx_sign_change
        idx = idx[-length(idx)]
        leftRoot=TRUE
        rightRoot=FALSE
      }
      
    }else if(!sign_change){
      # (b1,b2),...,(bn-1,bn) -> remove first, last
      idx = idx_root[-c(1,length(idx_root))]
      leftRoot=FALSE
      rightRoot=FALSE
    }
  }
  
  if(nroots_n2){
    
    #(b1,b2),(b2,b3),...,(b_{J-1},b_J),(b_{J+1},b_{J+2}),...,(b_{n-1},b_{n})
    
    #there is necessarily a sign change    
    idx = setdiff(idx_root, idx_sign_change)
    
    #remove left and right
    leftRoot=FALSE
    rightRoot=FALSE
    idx = idx[-1]
    idx = idx[-length(idx)]
    
    # nroots = length(support)-2
    #stop("not yet implemented")
  }
  
  endpoints = interval_endpoints[idx,,drop=F]
  endpoints[,1] = endpoints[,1]*exp(sign(endpoints[,1])*10^-8) #endpoints[,1]+1e-10
  endpoints[,2] = endpoints[,2]*exp(-sign(endpoints[,2])*10^-8) #endpoints[,2]-1e-10
  qfcn_vals = 
    cbind(q_fcn(endpoints[,1]),
          q_fcn(endpoints[,2]))
  
  #endpoints
  #qfcn_vals
  if(!all(apply(sign(qfcn_vals),1,prod)== -1)){
    warning("sign not opposite; potentially increase min/max boundary values")
  }
  
  
  qi = sapply(1:nrow(endpoints), 
              function(j) {
                lower=endpoints[j,1]
                upper=endpoints[j,2]
                extendInt="no"
                if ((leftRoot) & (j==1)){
                  extendInt="downX"
                }
                if ((rightRoot) & (j==nrow(endpoints))){
                  extendInt="upX"
                }
                
                return(uniroot(q_fcn, lower = endpoints[j,1],
                               upper = endpoints[j,2],
                               extendInt = extendInt,
                               tol = tol)$root
                )})
  return(qi)
  
}


q_alpha = function(ai){
  if(any(abs(ai)<.Machine$double.eps^{1/2})){stop("a!=0 in q_alpha(a)")}
  val = (1+ai^2)/(2*ai)
  return(val)
} 

q_alpha_inverse = function(qi){
  if(!all((qi>1) | (qi<(-1)))) stop("wrong qi range")
  return(ifelse(qi>1, qi-sqrt(qi^2-1),qi+sqrt(qi^2-1)))
}


eval_fct_roots = function(w, fct){
  out = sapply(w, function(wi) with(fct, c_phi*prod((cos(wi)-qi))/prod((cos(wi)-bi))) )
  return(out)
}

phi_fct = function(mu){
  
  # root finding
  fct_root = find_roots_mu(mu)
  
  # expose outputs of findRoots_mu
  a0 = fct_root$a0
  ai = fct_root$ai
  bi = fct_root$bi
  qi = fct_root$qi
  c_phi = fct_root$c_phi
  n = fct_root$n
  weights= fct_root$weights
  support = fct_root$support
  
  # wheter support contains 0 
  ind_0 = which(abs(support) < .Machine$double.eps^{1/2})
  
  if(length(ind_0)>0){
    a0 = weights[ind_0]
    support_0 = support[-ind_0]
    weights_0 = weights[-ind_0]
  }else{
    a0 = 0
    support_0 = support; weights_0 = weights
  }
  
  if(length(ind_0) > 0){
    # case I
    if(n==0){beta=0; xi = numeric(); C_phi = a0 }
    if(n> 0){
      beta = sort(support); 
      xi = sort(q_alpha_inverse(qi))
      C_phi = a0* prod(-2*support_0) / prod(-2*xi)
    }
  }else if( (length(ind_0)==0) & (abs(sum(ai)) > .Machine$double.eps^{1/2})){
    # case II
    if(n==1){beta = sort(support); xi = numeric(); C_phi = sum(ai)*(-2*support_0)}
    if(n> 1){
      beta = sort(support); xi = sort(q_alpha_inverse(qi));
      C_phi = sum(ai) * prod(-2*support_0) / prod(-2*xi)
    }
  }else if( (length(ind_0)==0) & (abs(sum(ai)) < .Machine$double.eps^{1/2})){
    # case III
    stop("not implemented")
    if(n==1){stop("cannot happen; error")}
    if(n==2){
      beta = sort(support); xi = numeric()
    }
    if(n>2){
      beta = sort(support); xi = sort(c(0,q_alpha_inverse(qi)))
    }
  }
    
  out = list(beta=beta,xi=xi,C_phi = C_phi,fct_root = fct_root)
  return(out)
}

eval_fct=function(w,fct){
  sapply(w, function(wi) with(fct, C_phi*prod(1-2*xi*cos(wi)+xi^2)/prod(1-2*beta*cos(wi)+beta^2) ))
}


phi_inverse_pf=function(mu){
 
  fct = phi_fct(mu = mu)
  beta = sort(fct$beta)
  xi = sort(fct$xi)
  
  # xi==0 case (Case III) not implemented yet
  stopifnot(all(abs(xi) >= .Machine$double.eps^{1/2})) 
  
  n_phi = length(beta)
  
  if(n_phi==1){
    out = list(
      C_phi_inv = 1/fct$C_phi,
      phi_inv_cp = cosinePoly_from_beta(beta),
      pf_coef = NULL,
      xi = xi
    )
  }else if(n_phi==2){
    out=list(
      C_phi_inv = 1/fct$C_phi,
      phi_inv_cp = cosinePoly_from_beta(beta),
      pf_coef = 1/(1-xi^2),
      xi = xi
    )
  }else{
    # when n_phi > 2
    beta_r1 = beta[-c(1,length(beta))]
    beta_r2 = beta[c(1,length(beta))]
    
    # check whether length(beta) +1 == length(xi)
    stopifnot(length(beta_r1)+1 == length(xi))
    
    # check whether 0 is contained in beta_r1
    ind_0 = which(abs(beta_r1) < .Machine$double.eps^{1/2})
    
    
    # solve system of linear equations
    RHS = c(rep(0, length(beta_r1)), prod(1+beta_r1^2) / prod(1+xi^2))
    # mat%*% pf_coef == RHS
    mat = matrix(NA,nrow = length(xi), ncol = length(xi))
    fcn_pf = function(x){(1-xi^2)/(1-2*xi*x+xi^2) }
    if(length(ind_0)>0){xseq = (1+beta_r1[-ind_0]^2)/(2*beta_r1[-ind_0])}else{
      xseq=(1+beta_r1^2)/(2*beta_r1)
    }
    
    for(i in 1:length(xseq)){
      x = xseq[i]
      mat[i,] = fcn_pf(x)
    }
    if(length(ind_0)==1){
      mat[(nrow(mat)-1),] = (1-xi^2)/(-2*xi)
    }
    
    mat[nrow(mat),] = fcn_pf(0)
    
    C_phi_inv = 1/fct$C_phi
    phi_inv_cp = cosinePoly_from_beta(beta_r2)
    pf_coef = solve(mat,RHS)
    out=list(C_phi_inv = C_phi_inv, phi_inv_cp = phi_inv_cp, pf_coef=pf_coef, xi = xi)
  }
  return(out)
}



# #### fcns related to chebyshev polynomials ####
# # p(x) = sum_k a_k T_k(x) = sum_k a_k cos(kx)
# # e.g., p = c(a1,a2,a3) represents a1*T_0(x)+a2*T_1(x)+a3*T_2(x)
# 
# multiply_cp<-function(p1,p2){
#   d1=length(p1)-1
#   d2=length(p2)-1
#   
#   d3=d1+d2
#   p3=rep(0,d3+1)
#   
#   
#   for (i in 1:(d1+1)){
#     i1=i-1
#     for (j in 1:(d2+1)){
#       j1=j-1
#       
#       p3[i1+j1+1]=p3[i1+j1+1]+0.5*p1[i]*p2[j]
#       p3[abs(i1-j1)+1]=p3[abs(i1-j1)+1]+0.5*p1[i]*p2[j]
#     }
#   }
#   return(p3)
# }
# 
# 
# 
# 
# eval_cp=function(w,p, two_sided = FALSE){
#   L = length(p)-1
#   # if two_sided = FALSE; sum_{k=0}^{L} p(k) cos(wk)
#   # if two_sided = TRUE; sum_{k=-L}^{L} p(k) cos(wk)
#   ret = numeric(length = length(w))
#   for(i in 1:length(ret)){
#     if(two_sided){
#       out = p*cos( w[i]*(0:L))
#       ret[i] = 2*sum(out)-p[1]
#     }else{
#       out = p*cos( w[i]*(0:L) )
#       ret[i] = sum(out)
#     }
#   }
#   
#   return(ret)
# }
# 
# 
# 
