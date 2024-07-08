supportReduction = function(XtX,Xtr, s_alpha, beta, inds, gradTrace = FALSE, tol = 1e-6, maxit=1000){
  
  ## support reduction algorithm ##
  ### Input: scaled XtX and Xtr, scaling factor s_alpha, and initialization (beta,inds)
  
  if(gradTrace){
    gradientExt=list(gradMin=NULL,
                     gradMax=NULL,
                     gradMinScaled=NULL,
                     gradMaxScaled=NULL)
    }else{gradientExt = NULL}
  
  
  
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
    if (gradTrace){gradientExt=updateGradient(gradientExt,gradient,s_alpha)}
    
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

  
  return(list(beta=beta, active_inds = inds, iter=count, convergence=convergence,gradient = gradient, gradientExt = gradientExt))
}
    

updateGradient<-function(gradientExt,grad,scale){
  gradientExt$gradMin=c(gradientExt$gradMin,min(-grad))
  gradientExt$gradMax=c(gradientExt$gradMax,max(-grad))
  gradientExt$gradMinScaled=c(gradientExt$gradMinScaled,min(-grad*scale))
  gradientExt$gradMaxScaled=c(gradientExt$gradMaxScaled,max(-grad*scale))
  return(gradientExt)
}
