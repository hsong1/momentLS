#' Tune delta for momentLS estimators
#' 
#' @description
#' Tune based on a modification of an adaptive bandwidth selection method by Politis <https://doi.org/10.1080/10485250310001604659>. 
#' See also our paper <https://arxiv.org/abs/2207.12705>.
#' 
#'@export
tune_delta = function(x, 
                      nSplits = 1,
                      c_M_const = 0,
                      fcn.summarize = "median"){
  
  # split x into nSplits sub-chains
  nSplits = round(nSplits)
  if(nSplits>1){
    ind = cut(1:length(x),breaks = nSplits, labels = FALSE)
  }else{
    ind = rep(1,length(x))
  }
  
  mhat_all = sapply(1:nSplits, function(j){
    r_j = autocov(x[ind==j])
    m_j = length(r_j)
    # find first mhat such that 
    # rc_j(mhat+2) <= c_M*sqrt(log(m_j)/m_j)
    # where c_M = c_M_const*sqrt(log(m_j))
    mhat = find_mhat(r = r_j, c_M = c_M_const*sqrt(log(m_j)))
    dhat = compute_delta_hat(M = m_j,mhat = mhat)
    c(m_j, mhat, dhat)
  })
  mhat_all = data.frame(t(mhat_all))
  colnames(mhat_all)= c("M","m_hat","d_hat")
  
  fcn.summarize = match.arg(fcn.summarize, choices=c("mean","median"))
  
  if(fcn.summarize =="mean"){
    mhat_summ = apply(mhat_all,2,mean)
  }else if(fcn.summarize=="median"){
    mhat_summ = apply(mhat_all,2,median)
  }
  
  m = mhat_summ[1] # average chain length per split
  mhat = mhat_summ[2]; # average truncation point
  
  delta_hat = compute_delta_hat(mhat = mhat,M = m)
  
  return(list(delta = delta_hat, mhat_all = mhat_all))
}


#'@export
find_mhat = function(r,c_M){
  # return first m such that
  # rc[mhat+2]<= c2*sqrt(log(M)/M))
  # where M = length of r

  cond = FALSE; mhat = 0
  M = length(r)
  rc = r[-1]/r[1] # empirical autocorrelations
  
  while(1){
    cond = (rc[mhat+2]<= c_M*sqrt(log(M)/M) )
    if(cond||is.na(cond)){break}else{mhat = mhat+2}
  }
  return(mhat)
}

#'@export
compute_delta_hat = function(M,mhat){
  if(length(mhat)!=length(M)) {stop("length(mhat)!=length(M")}
  deltas = foreach(i = 1:length(mhat),.combine="c")%do%{
    if (mhat[i] > 0) {delta = 1 - exp(-0.5 * (log(M[i])/mhat[i]))}
    else if (mhat[i] == 0) {delta = 1}
    delta = max(delta, 1/M)
  }
  return(deltas)
}
