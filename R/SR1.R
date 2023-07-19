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
           delta,
           alphaGrid = NULL,
           init = NULL,
           tol = 10^-8, # alg. tolerance
           maxit = 1000,
           gradTrace = FALSE,
           gridControl = list(
             nX = 1001, # number of grid points
             cm = FALSE, # completely monotone?
             scale = "log"
           ),
           precomputed = list(
             norm2_r = NULL,
             s_alpha = NULL,
             XtX = NULL,
             Xtr = NULL,
             alphaGrid = NULL,
             input = NULL)
  ) {
    
    
    
    if(gradTrace){gradientExt=list(gradMin=NULL,gradMax=NULL,gradMinScaled=NULL,gradMaxScaled=NULL)}else{gradientExt = NULL}
    n = length(r)
    
    ## input checks ##
    
    ## create a design matrix X ##
    ## xTemp = points in the grid from [-threshold,threshold] (or [0,threshold] if cm = TRUE)
    ## each column of X is [..., alpha^1, 1, alpha^1, alpha^2, alpha^3,...] / s_alpha
    ## where s_alpha = sqrt( sum_k alpha^{2|k|} ) = sqrt((1+alpha^2)/(1-alpha^2))
    if(is.null(alphaGrid)){
      xTemp = makeGrid(nX = gridControl$nX, upper_threshold = 1-delta, cm = gridControl$cm, scale = gridControl$scale)
    }else{
      xTemp = alphaGrid
    }
    if(is.null(precomputed$norm2_r)){norm2_r = 2*sum(r^2)-r[1]^2}else{norm2_r = precomputed$norm2_r}
    
    ##### relative tolerance:
    tol=tol*sqrt(norm2_r)
    
    
    ## precompute Xtr and XtX ##
    ## Compute Xtr
    if(is.null(precomputed$Xtr)){
      Xtr = Xtr_cpp(x=xTemp, a=r)
    }else{
      # check that Xtr was computed using the correct alphaGrid, input
      if(is.null(precomputed$input)){
        stop("When Xtr is not NULL, the precomputed list should contain input autocovariance estimator")}else{
          stopifnot( abs(precomputed$input-r) <1e-4)
        }
      if(is.null(precomputed$alphaGrid)){stop("When Xtr is not NULL, the precomputed list should contain alphaGrid")}else{
        stopifnot( abs(precomputed$alphaGrid-xTemp) <1e-4)
      }
      if(length(precomputed$Xtr)!=length(xTemp)){stop("length(Xtr)!=length(xTemp)")}
      Xtr = precomputed$Xtr
    }
    
    
    ## Compute XtX
    
    if(is.null(precomputed$s_alpha)){s_xTemp  = sqrt((1+xTemp^2)/(1-xTemp^2))}else{s_xTemp = precomputed$s_alpha}
    
    if(is.null(precomputed$XtX)){
      XtX=makeXtX(xTemp, s_x = s_xTemp)
    }else{
      # check that XtX was computed using the correct alphaGrid
      if(is.null(precomputed$alphaGrid)){stop("when XtX is not NULL, the precomputed list should contain alphaGrid")}else{
        stopifnot( abs(precomputed$alphaGrid-xTemp) <1e-4)
      }
      XtX = precomputed$XtX
      
    }
    ###if(is.null(precomputed$s_alpha)){s_xTemp  = sqrt((1+xTemp^2)/(1-xTemp^2))}else{s_xTemp = precomputed$s_alpha}
    
    ## support reduction algorithm ##
    
    ### Step 1: initialization ###
    if(is.null(init)){
      #no support points and weights to start out with (length 0)
      beta = numeric(0) #values of the weights corresponding to xTemp[inds]
      inds = integer(0) #current active set
    }else{
      #init = list(weights, support)
      if(any(init$weights<0)){stop("all weights need to be positive")}
      inds = findIndices(support = init$support,alphaGrid = alphaGrid)
      beta = init$weights/ s_xTemp[inds]
    }
    
    #### Step2 - 4: Loop ####
    count = 0 #iteration count
    convergence  = FALSE
    
    while (count< maxit){
      count=count+1
      
      ### sequential unrestricted minimization and support reduction ###
      SR_conv = FALSE
      if(any(beta<0)) stop("weights should be always non-negative")
      current = beta #current feasible solution
      
      while (!SR_conv){
        ##### unrestricted minimization ####
        if(length(inds)==0){
          betaNew = numeric(0) # inds = empty; then solution = null measure
        }else{
          cholesky=FALSE
          if (!cholesky){
            betaNew=solve(XtX[inds,inds,drop=FALSE],Xtr[inds])
          }else{
            R=chol(XtX[inds,inds,drop=FALSE])
            Qty=solve(t(R),Xtr[inds])
            betaNew=solve(R,Qty)
          }
        }
        
        if(all(betaNew>0)){
          # if all betaNew > 0, we are good; go to step 2
          SR_conv  = TRUE; beta = betaNew
          
        }else{
          # current = (c1,...,cK)
          # betaNew = (d1,...,dK)
          # Mass(t) = (1-t)(currentMass) + t newMass ( t=0: old; t=1: new )
          # largest 0<=t<=1 to make Mass(t)>=0 for all K?
          
          
          t = - current / (betaNew - current)
          t = t[-length(t)]
          # t[t<0] = -Inf; t[t>1] = -Inf
          # t = t[which.max(t)]
          t[t<0] = Inf; t[t>1] = Inf
          t_star = unique(t[which.min(t)])
          current = (1-t_star)*current + t_star*betaNew
          inds = inds[-which.min(t)]
          current = current[-which.min(t)]
          
        }
      }
      
      
      ### Step 2: find new support point
      ### 2-1: compute a gradient vector ###
      XtXw = XtX[,inds,drop=FALSE] %*% beta
      gradient = XtXw - Xtr
      if (gradTrace){gradientExt=updateGradient(gradientExt,gradient,s_xTemp)}
      
      ### 2-2 decide whether to add support point ###
      gradMin=min(gradient)
      # cat("gradMin at count",count,":",gradMin,"\n")
      
      #### if all gradients > -tol at current beta, we are done
      if (gradMin > -tol){ #changed to: if (gradMin > (-tol*sqrt(norm2_r))){...} ?
        convergence = TRUE
        break
      }
      
      ### Step 3: extended current support set ###
      ##### otherwise, add the support point where the gradient is largest
      inds = c(inds,which.min(gradient)) 
      beta = c(beta,0) 
    }
    
    if(!convergence){warning("SR algorithm did not converge after iter=",maxit," iterations")}
    
    #### output ####

    
    #reorder support for convenience
    support=xTemp[inds]
    
    #unstandardize beta
    weights=beta/s_xTemp[inds]
    weights=weights[order(support)]
    inds = inds[order(support)]
    support=support[order(support)]
    
    
    
    
    
    #compute objective value
    
    fvals = L2diff_L2Moment(r = r, support = support, weights = weights, 
                            precomputed=list(norm2_r = norm2_r, X_unscl_Atr = Xtr[inds]*s_xTemp[inds] ) )
    
    
    res = structure(
      list(weights = weights,
           support = support,
           inds = inds,
           XtX = XtX,
           s_alpha = s_xTemp,
           Xtr = Xtr,
           iter = count,
           fvals = fvals,
           minGrad = min(gradient),
           convergence = convergence,
           gradient = gradient,
           M = length(r),
           alphaGrid = xTemp,
           gradientExt = gradientExt
      ),
      class = "SRfit1")
    return(res)
    
  }

updateGradient<-function(gradientExt,grad,scale){
  gradientExt$gradMin=c(gradientExt$gradMin,min(-grad))
  gradientExt$gradMax=c(gradientExt$gradMax,max(-grad))
  gradientExt$gradMinScaled=c(gradientExt$gradMinScaled,min(-grad*scale))
  gradientExt$gradMaxScaled=c(gradientExt$gradMaxScaled,max(-grad*scale))
  return(gradientExt)
}
