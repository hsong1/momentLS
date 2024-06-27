#' SR1
#' 
#' Fit the MomentLS estimator given the input autocovariance sequence estimator
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
SR1 <-
  function(r, # input response vector
           delta = NULL,
           alphaGrid = NULL,
           init = NULL,
           tol = 10^-8, # alg. tolerance
           maxit = 1000,
           gradTrace = FALSE,
           n_alphas = 1001,
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
    
    ## XtX ##
    ## each column of X is [..., alpha^1, 1, alpha^1, alpha^2, alpha^3,...] / s_alpha
    ## where s_alpha = sqrt( sum_k alpha^{2|k|} ) = sqrt((1+alpha^2)/(1-alpha^2))
    
    if(is.null(precomputed$s_alpha)){
      s_alpha  = sqrt((1+alphaGrid^2)/(1-alphaGrid^2))
      }else{s_alpha = precomputed$s_alpha}
    
    if(is.null(precomputed$XtX)){
      XtX=makeXtX(alphaGrid, s_x = s_alpha)
    }else{
      # check that XtX was computed using the correct alphaGrid
      if(is.null(precomputed$alphaGrid)){
        stop("when XtX is not NULL, the precomputed list should contain alphaGrid")}else{
        stopifnot( abs(precomputed$alphaGrid-alphaGrid) <1e-4)
      }
      XtX = precomputed$XtX
      
    }
    
    ## Xtr ##
    # Xtr[i] = <x_alpha[i], r> / s_alpha[i]
    if(is.null(precomputed$Xtr)){
      Xtr = Xtr_cpp(x=alphaGrid, a=r)
    }else{
      # check that Xtr was computed using the correct alphaGrid and input r
      if(is.null(precomputed$input)){
        stop("When Xtr is not NULL, the precomputed list should contain input autocovariance estimate")}else{
          stopifnot( abs(precomputed$input-r) <1e-4)
        }
      if(is.null(precomputed$alphaGrid)){
        stop("When Xtr is not NULL, the precomputed list should contain alphaGrid")}else{
        stopifnot( abs(precomputed$alphaGrid-alphaGrid) <1e-4)
      }
      stopifnot(length(precomputed$Xtr)==length(alphaGrid))
      Xtr = precomputed$Xtr
    }
    
    
    ##### relative tolerance:
    norm2_r = 2*sum(r^2)-r[1]^2
    tol=tol*sqrt(norm2_r)
    
    ## support reduction algorithm ##
    SR_out = supportReduction(XtX = XtX,Xtr = Xtr,s_alpha = s_alpha,init = init,gradTrace = gradTrace, tol = tol, maxit=maxit)
    
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
    fvals = L2diff_L2Moment(r = r, support = support, weights = weights, 
                            precomputed=list(norm2_r = norm2_r, X_unscl_Atr = Xtr[inds]*s_alpha[inds] ) )
    
    res = structure(
      list(weights = weights,
           support = support,
           active_inds = inds,
           alphaGrid = alphaGrid,
           SRinputList = list(XtX=XtX,Xtr=Xtr,s_alpha=s_alpha),
           opt = list(iter=SR_out$iter,fvals=fvals,convergence=SR_out$convergence,gradient=SR_out$gradient,gradientExt = SR_out$gradientExt),
           M = length(r)
      ),
      class = "SRfit1")
  
    return(res)
    
  }

