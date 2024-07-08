#' compute XtX_w
#'@export
computeXtr_w = function(alphaGrid, 
                        r,m_uw, 
                        numchck_control = list(chck_TOL = .Machine$double.eps^{1/2},
                                               rplc_TOL = .Machine$double.eps^{1/3},
                                               N_wseq = 1e6)){
  
  # factorize 1/phi(w) from m_uw
  inv_fct = phi_inverse_pf(m_uw)
  p1 = multiply(inv_fct$phi_inv_cp,inv_fct$phi_inv_cp)
  
  # compute (two-sided) chebyshev coefficients of p(w)=rhat(w)*p1(w)
  FT_r = cosinePoly(coef = r,side = 2)
  p = multiply(p1,FT_r)
  tilde_bk = p$coef
  # p2 = c(r[1],r[-1]*2)
  # p = multiply_cp(p1,p2)
  # tilde_bk = c(p[1], p[-1]/2)
  tilde_ck = compute_tilde_ck(tilde_bk)
  
  # compute Bi
  xi = inv_fct$xi
  C_phi_inv = inv_fct$C_phi_inv
  pf_coef = inv_fct$pf_coef
  cK = function(a,b){a*(1-b^2)/((1-a*b)*(a-b))}
  
  if(length(xi)==0){
    # out = C_phi^{-2}*(1/2pi) int_w K(a,w)p(w) dw
    out = C_phi_inv^2*int_Ka1p(a1 = alphaGrid,tilde_bk = tilde_bk)
  }else if(length(xi)==1){
    # out = C_phi^{-2}*(1/2pi) int_w K(a,w)p(w){d_phi^2 K(xi_1,w)^2} dw
    grid_mat = expand.grid(alphaGrid,xi)
    out = (C_phi_inv*pf_coef)^2*int_Ka1Ka2pow2p(a1 = grid_mat$Var1,a2 = grid_mat$Var2,tilde_bk = tilde_bk,tilde_ck = tilde_ck)
  }else{
    
    grid_mat = lapply(1:length(xi), function(k) expand.grid(alphaGrid,xi[k]))
    # alphaGrid = column1; xi[1] = column2; ...
    grid_mat = cbind(grid_mat[[1]]$Var1,sapply(grid_mat,function(x) x$Var2)) 
    
    term1 = matrix(nrow=length(alphaGrid),ncol = length(xi))
    for(k in 1:length(xi)){
      # int_eval1 = 1/(2pi)int p(w)K(a_i,w)K(xi_k,w)^2 dw
      int_eval1 = int_Ka1Ka2pow2p(a1 = grid_mat[,1], a2 = grid_mat[,(1+k)],tilde_bk = tilde_bk,tilde_ck = tilde_ck)
      term1[,k] = int_eval1 * pf_coef[k]^2
    }
    term1 = apply(term1,1,sum)
    
    term2 = matrix(nrow=length(alphaGrid), ncol = choose(length(xi),2))
    it = 1
    for(k in 1:(length(xi)-1)){
      for(l in (k+1):length(xi)){
        
        term2_1 = cK(xi[k],xi[l])*int_Ka1Ka2p(a1 = grid_mat[,1],a2 = grid_mat[,(1+k)],tilde_bk = tilde_bk,tilde_ck = tilde_ck)
        term2_2 = cK(xi[l],xi[k])*int_Ka1Ka2p(a1 = grid_mat[,1],a2 = grid_mat[,(1+l)],tilde_bk = tilde_bk,tilde_ck = tilde_ck)
        term2[,it] = pf_coef[k]*pf_coef[l]*(term2_1 + term2_2)
        it = it+1
      }
    }
    rm(it)
    term2 = apply(term2,1,sum)*2
    out = (term1 + term2) * C_phi_inv^2
  }
  
  # check 
  chck_TOL = numchck_control$chck_TOL
  rplc_TOL = numchck_control$rplc_TOL
  N_wseq = numchck_control$N_wseq
  
  if(length(xi)>0){
    
    dist= abs(grid_mat[,-1,drop=F] - grid_mat[,1])
    ind_chck = which(dist < chck_TOL,arr.ind = T)
    ind_chck = unique(ind_chck[,1])
    
    if(length(ind_chck)>0){
      message("numerical checks for ",length(ind_chck)," alphas")
      
      # numerically compute (1/2pi) int K(a,w)FT_r(w)/phi(w)^2 dw
      # for a in alphaGrid[ind_chck]
      wseq = 2*pi*(0:(N_wseq-1))/N_wseq
      phi_wseq = phi_cpp(wseq = wseq, support = m_uw$support,weights = m_uw$weights)
      
      # compute FT(r)(wseq)
      raug = c(r, rep(0, length(wseq)-length(r)))
      FT_r_wseq = Re(2*fft(raug)) - r[1]
      # numerical integration
      out_chck = computeXtr_w_cpp(alphaGrid = alphaGrid[ind_chck], FT_r_wseq = FT_r_wseq, wseq = wseq,phi_wseq = phi_wseq)
      
      # replace
      ind_rplc = which(abs(out_chck - out[ind_chck]) > rplc_TOL)
      if(length(ind_rplc)>0){
        warning("<x_a,r>_phi numerically computed for ", length(ind_rplc), " alphas")
        out[ind_chck[ind_rplc]] = out_chck[ind_rplc]}
    } # fi for ind_chck
    
  }
  
  return(out)
}





#' compute Xtr_w
#'@export
makeXtX_w =function(alphaGrid, m_uw, 
                    numchck_control = list(chck_TOL = .Machine$double.eps^{1/2},
                                           rplc_TOL = .Machine$double.eps^{1/3},
                                           N_wseq = 1e6)){
  # factorize 1/phi(w) from m_uw
  inv_fct = phi_inverse_pf(m_uw)
  p1 = multiply(inv_fct$phi_inv_cp,inv_fct$phi_inv_cp)
  p1 = p1$coef

  if(length(p1)<5){
    p1_pad = rep(0,5)
    p1_pad[1:length(p1)] = p1
    p1 = p1_pad; rm(p1_pad)
  }
  # compute Aij
  xi = inv_fct$xi
  C_phi_inv = inv_fct$C_phi_inv
  pf_coef = inv_fct$pf_coef

  # output
  p=length(alphaGrid)
  H=matrix(0,p,p)
  # H[i,l] = 1/(2pi) int K(alpha_i,w)K(alpha_l,w)/ phi(w)^2 dw
  
  if(length(xi)==0){
    xi = 0;
    pf_coef = 1
  }
  
  # evaluate
  # 1/(2pi*C_phi^2) int_w p1(w)^2 \{K(alpha_i,w)K(alpha_l,w)\} \{sum_j sum_k d_k d_l K(xj,w)K(xl,w)} dw
  # Find first 5 (two-sided) Fourier Coefficients of h_i(w) = K(alpha_i,w) (d_k d_l K(xj,w)K(xl,w)) 
  # for each alpha_i in alphaGrid
  cab=matrix(0,p,p)
  cab=outer(alphaGrid,alphaGrid,function(a,b){a*(1-b^2)/((a-b)*(1-a*b))})
  
  ns=length(xi)
  fourierCoefficients=matrix(0,p,5)
  for (i in 1:p){
    for (j in 1:ns){
      for (k in 1:ns){
        for (n in 1:5){
          a=alphaGrid[i]
          b=xi[j]
          c=xi[k]
          cj=pf_coef[j]
          ck=pf_coef[k]
          fourierCoefficients[i,n]=fourierCoefficients[i,n]+h_abck(a,b,c,n-1)*cj*ck
        }
      }
    }
  }
  # compute int_w p1(w)^2 h_i(w) dw
  intVec=as.numeric((fourierCoefficients%*%p1))
  
  # compute
  # H[i,j] = 1/(2pi*C_phi^2) int_w c_K(i,j) p1(w)^2 h_i(w) + c_K(j,i) p1(w)^2 h_j(w)
  
  K=intVec*cab
  H=K+t(K)
  
  # recompute diagonals
  # Find first 5 (two-sided) Fourier Coefficients of 
  # h2_i(w) = K(alpha_i,w)^2 (d_k d_l K(xj,w)K(xl,w)) 
  fourierCoefficients2=matrix(0,p,5)
  
  for (i in 1:p){
    for (j in 1:ns){
      for (k in 1:ns){
        for (n in 1:5){
          a=alphaGrid[i]
          b=xi[j]
          c=xi[k]
          cj=pf_coef[j]
          ck=pf_coef[k]
          fourierCoefficients2[i,n]=fourierCoefficients2[i,n]+h_a2bck(a,b,c,n-1)*cj*ck
        }
      }
    }
  }
  
  intVec2=(fourierCoefficients2%*%p1)
  diag(H)=intVec2
  
  H=H*C_phi_inv^2
  
  # check 
  chck_TOL = numchck_control$chck_TOL
  rplc_TOL = numchck_control$rplc_TOL
  N_wseq = numchck_control$N_wseq
  
  alphaGridDist = as.matrix(dist(alphaGrid) )
  ind_chck = which(0 < alphaGridDist & alphaGridDist < chck_TOL, arr.ind=T)
  ind_chck = ind_chck[ind_chck[,1]<ind_chck[,2], ]
  if(nrow(ind_chck)>0){
    message("numerical checks for ",nrow(ind_chck)," alphas")
    wseq = 2*pi*(0:(N_wseq-1))/N_wseq
    phi_wseq = phi_cpp(wseq = wseq, support = m_uw$support,weights = m_uw$weights)
    out_chck = vector("numeric",length = nrow(ind_chck))
    for(i in 1:nrow(ind_chck)){
      out_chck[i] = makeXtX_w_cpp(alphaGrid = c(alphaGrid[ind_chck[i,1]],alphaGrid[ind_chck[i,2]]), wseq = wseq, phi_wseq = phi_wseq,diag = F)[2,1]
    }
    
   ind_rplc= which(abs(H[ind_chck] - out_chck) > rplc_TOL)
   if(length(ind_rplc)>0){
     warning("<x_a,x_a'>_phi numerically computed for ", length(ind_rplc), " alpha pairs")
     H[ind_chck[ind_rplc,]] = out_chck[ind_rplc]}
  }
  
  return(H) 
}
