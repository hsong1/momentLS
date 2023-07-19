#' compute ||m||^2 for a moment sequence m
#'@export
momentNorm2<-function(support,weights){
  # ||m||^2 = int (1+alpha_1 alpha_2)/ (1-alpha_i alpha_j) dF(alpha_1,alpha_2)
  # solutionNorm2 in SupportReduction.R
  
  aat = outer(support,support)
  norm2 = sum((1+aat)/(1-aat) * outer(weights, weights))
  
  return(norm2)
}

#' compute sum_k r(k)^2 for k = -M,...,M for r=(r(0),r(1),...,r(M))
#'@export
l2Norm2 = function(r){
  return(2*sum(r^2)-r[1]^2)
}

#' inner product between r in L2(Z) and m in M_infty
#'@export
innerProduct_L2Moment = function(r,support ,weights, precomputed=list(norm2_r = NULL,X_unscl_Atr= NULL)){
  
  # <r,m> = sum_alpha <r,x_alpha> w_alpha
  if(length(weights)==0){innerProduct = 0}else{
    
    if(!is.null(precomputed$X_unscl_Atr)){
      # XAtr = [<x_alpha/s_alpha, r>]_{alpha in A}
      innerProduct = sum(precomputed$X_unscl_Atr * weights)
    }else{
      # compute sum_alpha <x_alpha,r> w_alpha
      n = length(r)
      innerProduct = sum(sapply(1:length(weights), function(i){
        x=support[i]^(0:(n-1)); return((2*sum(x*r)-r[1])*weights[i])}
      ))
    }
  
  }
  
  return(innerProduct)
  
}

#' inner product between m1,m2 in M_infty
#'@export
innerProduct_twoMoments = function(support1, weights1, support2, weights2){
  # inner product between m1,m2 in M_infty
  aat = outer(support1,support2)
  ip = sum((1+aat)/(1-aat) * outer(weights1, weights2))
  return(ip)
}
  
#' compute || r - m ||^2 for r in L2(Z) and m in M_infty
#'@export
L2diff_L2Moment = function(r,support ,weights, precomputed=list(norm2_r = NULL,X_unscl_Atr= NULL)){
  # || r - m ||^2 = ||r||^2 + ||m||^2 - 2 <r,m>
  
  if(is.null(precomputed$norm2_r)){norm2_r=2*sum(r^2)-r[1]^2}else{
    norm2_r = precomputed$norm2_r
  }
  if(length(weights)==0){
    f = norm2_r
  }else{
    norm2_m=momentNorm2(weights = weights,support = support)
    innerProduct = innerProduct_L2Moment(r = r,support = support,weights = weights,precomputed = precomputed)
    
    f = norm2_r+norm2_m-2*innerProduct
  }
  if(f<0)stop("wrong precomputed inputs")
  return(f)
  
}

#' compute || m1 - m2 ||^2 for m1,m2 in M_infty
#'@export
L2diff_twoMoments = function(support1, weights1, support2, weights2){
  
  val = momentNorm2(support = support1,weights = weights1)+
    momentNorm2(support = support2, weights = weights2)-
    2* innerProduct_twoMoments(support1 = support1,weights1 = weights1,support2 = support2,weights2 = weights2)
  return(val)

}

