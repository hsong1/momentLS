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
#' @useDynLib momentLS
#' @import Rcpp
#' @export
SR1_w <-
  function(r, # input response vector
           delta = NULL,
           alphaGrid = NULL,
           phi = NULL,
           n_phi = length(r),
           comp_method = c("exact","num"),
           init = NULL,
           tol = 10^-6, # alg. tolerance
           maxit = 1000,
           gradTrace = FALSE,
           n_alphas = 501,
           gridControl = list(
             cm = FALSE, # completely monotone?
             scale = "log"
           ),
           precomputed = list(
             s_alpha = NULL,
             XtX_w = NULL,
             Xtr_w = NULL,
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
      m_uw = SR1(r,alphaGrid = alphaGrid); weightType="momentLS"
      phi = m_uw
    }
    
    ## set up weightType ##
    comp_method = match.arg(comp_method,choices = c("exact","num"))
    
    if(is(phi,"SRfit1")&comp_method=="exact"){
      m_uw=phi; weightType="momentLS"
    }else if(is(phi,"SRfit1")&comp_method=="num"){
      wseq = (0:(n_phi-1))*2*pi/n_phi
      phi_wseq = phi_cpp(wseq,support = phi$support,weights = phi$weights)
      weightType = "others"
    }else{
      # wseq and phi_wseq 
      if(is.null(phi$wseq) || is.null(phi$phi_wseq)){
        stop("wseq and phi_wseq need to be provided")}
      stopifnot(length(phi$wseq)==length(phi$phi_wseq))
      stopifnot(length(phi$wseq)>=length(r))
      message("numerical integrations are based on provided wseq")
      phi_wseq = phi$phi_wseq
      wseq = phi$wseq
      weightType="others"
      }
      
    
    
    
    ## XtX ##
    if(!is.null(precomputed$XtX_w)){
      # checks
      if(is.null(precomputed$alphaGrid)){stop("when XtX_w is not NULL, the precomputed list should contain alphaGrid")}else{
        stopifnot( abs(precomputed$alphaGrid-alphaGrid) <1e-4)}
      if(is.null(precomputed$s_alpha)){
        stop("s_alpha needs to be provided")}
      
      XtX_w = precomputed$XtX_w
      s_alpha = precomputed$s_alpha
      
    }else{
      # when no precomputed XtX_w
      if(weightType=="momentLS"){
        XtX_w = makeXtX_w(alphaGrid = alphaGrid, m_uw = m_uw)
        s_alpha = sqrt(diag(XtX_w))
        XtX_w = XtX_w/ outer(s_alpha,s_alpha)
      }else{
        message("numerical integration for XtX_w")
        XtX_w = makeXtX_w_cpp(alphaGrid = alphaGrid,wseq = wseq,phi_wseq = phi_wseq,diag = TRUE)
        s_alpha = sqrt(diag(XtX_w))
        XtX_w = XtX_w/ outer(s_alpha,s_alpha)
      }
    }
    
    
     
      ## Xtr ##
      # Xtr[i] = <x_alpha[i], r>_phi / s_alpha[i]
    
    if(!is.null(precomputed$Xtr_w)){
      # check that Xtr was computed using the correct alphaGrid and input r
      if(is.null(precomputed$input)){
        stop("When Xtr_w is not NULL, the precomputed list should contain input autocovariance estimate")}else{
          stopifnot( abs(precomputed$input-r) <1e-4)}
      if(is.null(precomputed$alphaGrid)){
        stop("When Xtr_w is not NULL, the precomputed list should contain alphaGrid")}else{
          stopifnot( abs(precomputed$alphaGrid-alphaGrid) <1e-4)
        }
      stopifnot(length(precomputed$Xtr_w)==length(alphaGrid))
      
      Xtr_w = precomputed$Xtr_w
    }else{
      # when no precomputed Xtr_w
      if(weightType=="momentLS"){
        Xtr_w = computeXtr_w(alphaGrid = alphaGrid,r = r,m_uw = m_uw)
        Xtr_w = Xtr_w / s_alpha
        }else{
          message("numerical integration for Xtr_w")
          raug = c(r, rep(0, length(wseq)-length(r)))
          FT_r_wseq = Re(2*fft(raug)) - r[1]
          Xtr_w = computeXtr_w_cpp(alphaGrid = alphaGrid,FT_r_wseq = FT_r_wseq,
                                   wseq = wseq,phi_wseq = phi_wseq)
          Xtr_w = Xtr_w / s_alpha
        }
    }
    
    
    
    ##### relative tolerance:
    norm2_r = 2*sum(r^2)-r[1]^2
    tol=tol*sqrt(norm2_r)
    
    ## initialization ##
    if(is.null(init)){
      #no support points and weights to start out with (length 0)
      inds = integer(0) #current active set
      beta = numeric(0) #values of the weights corresponding to xTemp[inds]
    }else{
      #init = list(weights, support)
      if(any(init$weights<0)){stop("all weights need to be positive")}
      inds = findIndices(support = init$support,alphaGrid = alphaGrid)
      beta = init$weights/ s_alpha[inds]
    }
    
    ## support reduction algorithm ##
    SR_out = supportReduction(XtX = XtX_w,Xtr = Xtr_w,s_alpha = s_alpha,beta = beta, inds=inds, gradTrace = gradTrace, tol = tol)
    
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

