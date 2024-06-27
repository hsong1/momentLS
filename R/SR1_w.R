#' SR1_w
#' 
#' Fit the weighted MomentLS estimator given the input autocovariance sequence estimator
#' 
#' @param r input autocovariance estimate
#' @param delta parameter to specify appropriate moment space; use tune_delta if unknown
#' @param alphaGrid a grid of x values from [-1+delta,1-delta]; will be automatically created if left unspecified
#' @param init initialization for the SR algorithm
#' @param tol algorithm tolerance
#' @param maxit maximum number of iterations for support reduction algorithm
#' @param gradTrace trace gradients if TRUE
#' @param gridControl list which contains parameters to specify alphaGrid
#' @param precomputed an optional list which can contain some precomputed quantities for computational speed-up
#' @return a SR1fit object containing estimated support and weights
#'@useDynLib momentLS
#'@import Rcpp
#'@export
SR1_w <-
  function(r, # input response vector
           delta = NULL,
           alphaGrid = NULL,
           phi = NULL,
           init = NULL,
           tol = 10^-8, # alg. tolerance
           maxit = 1000,
           gradTrace = FALSE,
           n_alphas = 501,
           gridControl = list(
             cm = FALSE, # completely monotone?
             scale = "log"
           ),
           precomputed = list(
             s_alpha = NULL,
             XtX = NULL,
             Xtr = NULL,
             alphaGrid = NULL,
             input = NULL)
  ) {
    
    
    n = length(r)
    
    ## input checks (XtX, Xtr, s_alpha, alphaGrid) ##
    
    ## alphaGrid = points in the grid from [-threshold,threshold] (or [0,threshold] if cm = TRUE)
    if(!is.null(alphaGrid)){
      # if both alphaGrid and delta are provided, delta value is ignored
      if(!is.null(delta)){warning("provided delta not used")}
    }
    if(is.null(alphaGrid)){
      if(is.null(delta)){stop("at least one of delta or alphaGrid needs to be provided")}
      alphaGrid = makeGrid(nX = n_alphas, upper_threshold = 1-delta, cm = gridControl$cm, scale = gridControl$scale)
    }
    
    ## weight phi ##
    if(is.null(phi)){
      m_uw = SR1(r,alphaGrid = alphaGrid); phi = m_uw; weightType="momentLS"
    }else{
      
      if(class(phi)=="SRfit1"){m_uw=phi; weightType="momentLS"}else{
        # wseq and phi_wseq 
        if(is.null(phi$wseq) || is.null(phi$phi_wseq)){
          stop("wseq and phi_wseq need to be provided")}
        stopifnot(length(phi$wseq)==length(phi$phi_wseq))
        stopifnot(length(phi$wseq)>=length(r))
        weightType="others"
      }
      
    }
    
    
    ## XtX ##

    if(weightType=="momentLS"){
      if(is.null(precomputed$XtX_w)){
        XtX_w = makeXtX_w(alphaGrid = alphaGrid, m_uw = m_uw)
        s_alpha = sqrt(diag(XtX_w))
        XtX_w = XtX_w/ outer(s_alpha,s_alpha)
      }else{
        if(is.null(precomputed$alphaGrid)){
          stop("when XtX_w is not NULL, the precomputed list should contain alphaGrid")}else{
            stopifnot( abs(precomputed$alphaGrid-alphaGrid) <1e-4)}
        if(is.null(precomputed$s_alpha)){
          stop("s_alpha needs to be provided")
        }
        XtX_w = precomputed$XtX_w
      }
      
      
      ## Xtr ##
      # Xtr[i] = <x_alpha[i], r>_phi / s_alpha[i]
      if(is.null(precomputed$Xtr_w)){
        Xtr_w = computeXtr_w(alphaGrid = alphaGrid,r = r,m_uw = m_uw)
        Xtr_w = Xtr_w / s_alpha
      }else{
        # check that Xtr was computed using the correct alphaGrid and input r
        if(is.null(precomputed$input)){
          stop("When Xtr_w is not NULL, the precomputed list should contain input autocovariance estimate")}else{
            stopifnot( abs(precomputed$input-r) <1e-4)
          }
        if(is.null(precomputed$alphaGrid)){
          stop("When Xtr_w is not NULL, the precomputed list should contain alphaGrid")}else{
            stopifnot( abs(precomputed$alphaGrid-alphaGrid) <1e-4)
          }
        stopifnot(length(precomputed$Xtr_w)==length(alphaGrid))
        Xtr_w = precomputed$Xtr_w
      }
    }else{
      raug = c(r, rep(0, length(phi$wseq)-length(r)))
      FT_r_wseq = Re(2*fft(raug)) - r[1]
      message("numerical integration for XtX_w and Xtr_w")
      Xtr_w = computeXtr_w_cpp(alphaGrid = alphaGrid,FT_r_wseq = FT_r_wseq,wseq = phi$wseq,phi_wseq = phi$phi_wseq)
      XtX_w = makeXtX_w_cpp(alphaGrid = alphaGrid,wseq = phi$wseq,phi_wseq = phi$phi_wseq,diag = TRUE)
      s_alpha = sqrt(diag(XtX_w))
      XtX_w = XtX_w/ outer(s_alpha,s_alpha)
      Xtr_w = Xtr_w / s_alpha
    }    
    
    
    ##### relative tolerance:
    norm2_r = 2*sum(r^2)-r[1]^2
    tol=tol*sqrt(norm2_r)
    
    ## support reduction algorithm ##
    SR_out = supportReduction(XtX = XtX_w,Xtr = Xtr_w,s_alpha = s_alpha,init = init,gradTrace = gradTrace, tol = tol)
    
    #### output ####
    inds = SR_out$active_inds
    beta = SR_out$beta
    
    #reorder support for convenience
    support=alphaGrid[inds]
    
    #unstandardize beta
    weights=beta/s_alpha[inds]
    weights=weights[order(support)]
    inds = inds[order(support)]
    support=support[order(support)]
    
    
    ## compute objective value
    #fvals = L2diff_L2Moment(r = r, support = support, weights = weights, 
    #                        precomputed=list(norm2_r = norm2_r, X_unscl_Atr = Xtr[inds]*s_alpha[inds] ) )
    
    res = structure(
      list(weights = weights,
           support = support,
           active_inds = inds,
           alphaGrid = alphaGrid,
           SRinputList = list(XtX_w=XtX_w,Xtr_w=Xtr_w,s_alpha=s_alpha),
           opt = list(iter=SR_out$iter,fvals=NA,convergence=SR_out$convergence,gradient=SR_out$gradient,gradientExt = SR_out$gradientExt),
           phi = phi,
           M = length(r)
      ),
      class = "SRfit1")
    
    return(res)
    
  }

