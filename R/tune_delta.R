#' Tune delta for momentLS estimators
#' 
#' @description
#' Tune based on a modification of an adaptive bandwidth selection method by Politis <https://doi.org/10.1080/10485250310001604659>. 
#' See also our paper <https://arxiv.org/abs/2207.12705>.
#'@param x input chain
#'@param nSplits number of splits
#'@param c_M_const c_M constant in \eqn{\hat{m} = \min\{t\in 2\mathbb{N}; \hat{\rho}_M(t+2) \le c_M\sqrt{\log M/M}\}}
#'@param fcn.summarize if fcn.summarize =="mean", the average of deltas from each split will be returned. If fcn.summarize =="median", median of deltas from each split will be returned. 
#'@return A numeric value (delta) and a data frame (mhat_all)
#'\itemize{
#'   \item delta: estimated delta (mean or median from \eqn{\hat{delta}_M} from each split.)
#'   \item mhat_all: estimated \eqn{\hat{m}} and \eqn{\hat{delta}_M} from each split.
#' }
#'@export
tune_delta = function(x, 
                      nSplits = 1,
                      c_M_const = 0,
                      fcn.summarize = "mean"){
  
  r_j_all = compute_autocov_splits(x,nSplits)
  
  mhat_all = sapply(1:nSplits, function(j){
    r_j =  r_j_all[[j]]
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
    delta_hat = mean(mhat_all$d_hat)
  }else if(fcn.summarize=="median"){
    mhat_summ = median(mhat_all$d_hat)
  }
  
  
  return(list(delta = delta_hat, mhat_all = mhat_all))
}

#'@export
compute_autocov_splits = function(x,L){
  
  # return the list of length L
  # each element l contains 
  # l=0 B^{-1} sum t=0,...,B-1-k g(Xt)g(Xt+k)
  # l=2,...,L B^{-1} sum t=(l-1)*B-k to l*B-k g(Xt)g(Xt+k)
  
  # center
  xc = x-mean(x)
  B = floor(length(xc)/L)
  xc = xc[1:(L*B)] # make M to be the multiple of L
  
  cum_sums = foreach(l=1:L)%do%{
    (l*B)*autocov(xc[1:(l*B)],center = F)[1:B]
  }
  
  
  autocov_batches = foreach(l=1:L)%do%{
    if(l==1){(1/B)*cum_sums[[1]]}else{
      (1/B)*(cum_sums[[l]]-cum_sums[[l-1]])
    }
  }
  return(autocov_batches)
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

#' @export
compute_delta_hat = function(M,mhat){
  if(length(mhat)!=length(M)) {stop("length(mhat)!=length(M")}
  deltas = foreach(i = 1:length(mhat),.combine="c")%do%{
    if (mhat[i] > 0) {delta = 1 - exp(-0.5 * (log(M[i])/mhat[i]))}
    else if (mhat[i] == 0) {delta = 1}
    delta = max(delta, 1/M)
  }
  return(deltas)
}
