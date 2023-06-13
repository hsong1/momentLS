#'@import doParallel

ind_min1se = function(m_vals,se_vals){
  # choose the delta with the smallest m_vals (tie: return the smallest)
  indmin <- min(which(m_vals==min(m_vals))) 
  
  ind <-  intersect(which(abs(m_vals-m_vals[indmin])<=se_vals[indmin]),(1:indmin))
  if(length(ind)==0){
    ind1se <-  indmin
  } else {
    ind1se <-  min(ind)
  }
  return(list(indmin = indmin, ind1se = ind1se))
}


update_lst_precomp = function(delta_c, lst_precomp){
  # input: lst_precomp containing alphaGrid, XtX, Xtr, and s_alpha
  # based on delta_c, subset XtX, Xtr, alphaGrid, and recompute s_alpha
  # return: list containing sub matrix / vector of alphaGrid, XtX, Xtr, and s_alpha
  
  ind = abs(lst_precomp$alphaGrid)<=1-delta_c
  XtX_sub = lst_precomp$XtX[ind,ind]
  Xtr_sub = lst_precomp$Xtr[ind]
  s_alpha_sub = sqrt((1+lst_precomp$alphaGrid[ind]^2)/(1-lst_precomp$alphaGrid[ind]^2))
  return(list(alphaGrid = lst_precomp$alphaGrid[ind], XtX = XtX_sub, Xtr = Xtr_sub, s_alpha= s_alpha_sub))
}


precomp_X = function(alphaGrid){
  s_alpha = sqrt((1+alphaGrid^2)/(1-alphaGrid^2))
  XtX = makeXtX(x = alphaGrid,s_x = s_alpha)
  
  lst = list(alphaGrid = alphaGrid, XtX=XtX, s_alpha= s_alpha)
  return(lst)
}

precomp_SR = function(lst_X, r){
  Xtr = Xtr_cpp(x = lst_X$alphaGrid, a = r)
  lst = list(alphaGrid = lst_X$alphaGrid, XtX = lst_X$XtX, s_alpha = lst_X$s_alpha,
             Xtr = Xtr, input = r)
  return(lst)
}




#'@export
bootrisk = function(deltaseq, 
                    spectrf = NULL, 
                    spectr.alt = list(type = c("L2seq","momentseq"), seq=NULL),
                    B=10 , 
                    x  = NULL, 
                    M0 = length(x), 
                    boot.type  = "WB", 
                    parallel = FALSE,
                    method = NULL 
){
  
  ## output
  ## mean_b|| gamma(spectrf) - Pi_delta(r(D_b); delta)||^2 where r(D_b) ~ boot(x,spectrf)
  
  ## input
  # deltaseq = gap values
  # spectrf = input periodogram OR spectr.alt (either L2seq or momentseq)
  # method == c("numInt","exact")
  # x is needed when boot.type = "RB"
  # M0 is needed when boot.type = "WB"
  
  if(is.null(method)){
    # default method: If spectrf is provided, "numInt". If spectr.alt is provided, "exact"
    if(!is.null(spectrf)){method = "numInt"}else{method="exact"}
  }else{
    # default method can be overwritten by the user choice
    if(!method %in% c("numInt","exact")) stop("method needs to be one of numInt, exact")
  }
  
  boot.type = match.arg(boot.type, choices = c("RB","WB"))
  if(boot.type=="RB"){if(is.null(x)) stop("x needs to be provided")}
  if(boot.type=="WB"){if(M0==0) stop("M0 needs to be provided")}
  
  # Obtain spectrf based on L2 or momentseq.
  # spectrf is needed for bootstrapping
  if(is.null(spectrf)){
    if(is.null(spectr.alt$seq)){stop("when spectrf is not provided, spectr.alt needs to be provided")}
    spectr.alt$type = match.arg(spectr.alt$type, c("L2seq","momentseq"))
    ft_ = ft_seq(type = spectr.alt$type, seq = spectr.alt$seq)
    spectrf = approxfun(x = ft_$FFreq, y = ft_$spectr, rule=2)
    
  }
  
  deltaseq = sort(deltaseq) #smallest to largest
  # a grid for SR
  alphaGrid = makeGrid(upper_threshold = 1-min(deltaseq), nX = 1001, cm = F, scale = "log")
  lst_X = precomp_X(alphaGrid = alphaGrid)
  
  
  
  #### parallel environment
  if(parallel){
    cl <- makeCluster(parallel::detectCores()-2)
    registerDoSNOW(cl)
    
    pb <- txtProgressBar(max = B, style = 3)
    progress <- function(n) setTxtProgressBar(pb, n)
    opts <- list(progress = progress)
    
    clusterExport(cl,list=c("M0","boot.type","x","parallel"),envir = environment())
  }else{
    foreach::registerDoSEQ()
    opts = NULL
    pb = progress::progress_bar$new(total=B)
  }
  
  
  vals = foreach(b=1:B,.combine="cbind",.packages=c("momentLS","foreach"),.options.snow=opts) %dopar%{
    
    if(!parallel){pb$tick()}
    
    if(boot.type=="RB"){
      Fx = subsetFcoefs(fft(x))
      boot = with(Fx, TFTboot(a_k=a_k,b_k=b_k, M = M0,spectrf =  spectrf,type = "RB"))
    }else{
      boot = TFTboot(a_k=NULL,b_k=NULL,M = M0, spectrf = spectrf,type = "WB")
    }
    
    x_boot = with(boot, reconstrTS(a_k_boot,b_k_boot,M0))
    
    ### compute Pi_d(r(D_b)) for d in deltaseq
    r = autocov(x_boot)
    lst_precomp = precomp_SR(lst_X,r)
    
    mLSE = list(); lst_precomp_c = list()
    for(i in 1:length(deltaseq)){
      delta_c = deltaseq[i]
      lst_precomp_c[[i]] = update_lst_precomp(delta_c,lst_precomp = lst_precomp)
      mLSE[[i]] = SR1(r = r, delta = delta_c, alphaGrid = lst_precomp_c[[i]]$alphaGrid, precomputed = lst_precomp_c[[i]])  
    }
    
    v = L2diff_gammaSeq(mLSE = mLSE, spectrf = spectrf,spectr.alt = spectr.alt,method = method,w = 2*pi*(0:(M0-1))/M0,dw = 1/M0)
    
    return(v)
    
  }
  rownames(vals) = deltaseq
  if(parallel) {stopCluster(cl)}
  
  m_vals = apply(vals,1,mean)
  se_vals = apply(vals,1,sd)/sqrt(B)
  
  inds=ind_min1se(m_vals,se_vals)
  
  options = list(boot.type = boot.type, method=method)
  ret_ls= list(vals= vals, m_vals = m_vals, se_vals = se_vals, deltaseq = deltaseq, delta.min = deltaseq[inds$indmin], delta.1se = deltaseq[inds$ind1se], options=options)
  return(ret_ls)
  
}



