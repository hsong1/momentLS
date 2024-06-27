compute_tilde_ck = function(tilde_bk){
  
  L = length(tilde_bk)-1
  if(L==0){
    # int_Kp = tilde_b0; d/da int_Kp = 0
    tilde_ck = 0
  }else if(L==1){
    # int_Kp = tilde_b0 + 2tilde_b1a; d/da int_Kp = 2tilde_b1
    tilde_ck = tilde_bk[2]*2
  }else{
    tilde_ck = c()
    tilde_ck = (2:L)*tilde_bk[3:(L+1)]
    tilde_ck = c(tilde_bk[2]*2, tilde_ck)
  }
  
  return(tilde_ck)
}



#### I1 ####
int_Ka1p = function(a1,tilde_bk,approx = TRUE){
  # evaluate (1/2pi) int K(a1,w) p(w)dw
  # p(w) = \sum tilde_bk cos(wk)
  if(approx){
    val = Xtr_cpp(x = a1,a = tilde_bk,standardization = FALSE)
    #print("Xtr_cpp")
  }else{
    print("exact!")
    exponents = seq(0, length(tilde_bk)-1,by=1)
    val = outer(a1,exponents,"^")%*%tilde_bk*2-tilde_bk[1]
    val = as.numeric(val)
  }
  return(val)
}

int_Ka1pow2p = function(a1, tilde_bk, tilde_ck, approx= TRUE){
  # evaluate (1/2pi) int K(a1,w)^2 p(w)dw
  # tilde_ck = compute_tilde_ck(tilde_bk)
  #a1*Xtr_cpp(x = a1, a = tilde_ck, standardization = FALSE) + 
  #  (a1^2+1)/(1-a1^2)*Xtr_cpp(x = a1, a = tilde_bk, standardization = FALSE)
  a1*int_Ka1p(a1 = a1, tilde_bk = tilde_ck, approx = approx) + 
    (a1^2+1)/(1-a1^2)*int_Ka1p(a1 = a1, tilde_bk = tilde_bk, approx=approx)
}


int_Ka1pow3p = function(a1,tilde_bk, tilde_ck, tilde_dk,approx=TRUE){
  # evaluate (1/2pi) int K(a1,w)^3 p(w)dw
  #tilde_ck = compute_tilde_ck(tilde_bk)
  #tilde_dk = compute_tilde_ck(tilde_ck)
  term1 = (a1^2+1)/(1-a1^2) * int_Ka1pow2p(a1 = a1,tilde_bk = tilde_bk,tilde_ck = tilde_ck,approx=approx)
  
  int2 = int_Ka1p(a1 = a1,tilde_bk = tilde_ck,approx=approx)
  
  term2_1 = 4*a1 / (1-a1^2)^2 * int_Ka1p(a1 = a1, tilde_bk = tilde_bk, approx=approx)
  term2_2 = (a1^2+1)/(1-a1^2)*int2
  term2_3 = int2
  term2_4 = a1*int_Ka1p(a1 = a1,tilde_bk = tilde_dk,approx=approx)
  term2 = (a1/2)*(term2_1+term2_2+term2_3+term2_4)
  
  return(term1+term2)
}

#### I2 ####
int_Ka1Ka2p <- function(a1,a2,tilde_bk,tilde_ck=NULL,approx=TRUE, TOL = .Machine$double.eps^{1/3}){
  # evaluate (1/2pi) int K(a1,w)K(a2,w) p(w)dw
  
  if(length(a1)!=length(a2)) stop("length(a1)!=length(a2)")
  
  is_a1_0 = (abs(a1)<TOL)
  is_a2_0 = (abs(a2)<TOL)
  a1_equals_a2 = (a1==a2)
  val = vector(mode="numeric",length = length(a1))
  
  # when a1==0 & a2==0
  ind_1 = which(is_a1_0 & is_a2_0)
  if(length(ind_1)>0){
    val[ind_1] = tilde_bk[1]
  }
  
  # when a1!=0  & a2==0
  ind_2 = which(!is_a1_0 & is_a2_0)
  if(length(ind_1)>0){
    val[ind_2] = int_Ka1p(a1 = a1[ind_2],tilde_bk = tilde_bk,approx=approx)
  }
  
  # when a1==0  & a2!=0
  ind_3 = which(is_a1_0 & !is_a2_0)
  if(length(ind_3)>0){
    val[ind_3] = int_Ka1p(a1 = a2[ind_3],tilde_bk = tilde_bk,approx=approx)
  }
  
  # when a1, a2 !=0 and a1 !=a2
  ind_4 = which(!is_a1_0 & !is_a2_0 & !a1_equals_a2)
  # (1/2pi) int K(a1,w)K(a2,w) p(w)dw for a1 != a2
  if(length(ind_4)>0){
    sum_1 = matrix(int_Ka1p(a1 = c(a1[ind_4],a2[ind_4]),tilde_bk = tilde_bk,approx=approx),ncol=2)
    K_wei = cbind(-(1-a2[ind_4]^2)/(2*a2[ind_4]), (1-a1[ind_4]^2)/(2*a1[ind_4]))
    val[ind_4]= apply(sum_1*K_wei,1,sum)/(q_alpha(a1[ind_4])-q_alpha(a2[ind_4]))
  }
  
  ind_5 = which(!is_a1_0 & !is_a2_0 & a1_equals_a2)
  # (1/2pi) int K(a1,w)^2 p(w)dw for a1 == a2
  if(length(ind_5)>0){
    val[ind_5] = int_Ka1pow2p(a1 = a2[ind_5],tilde_bk = tilde_bk,tilde_ck,approx=approx)
  }
  return(val)
  
}

g_int_Ka1Ka2p = function(a1,a2,tilde_bk,tilde_ck, approx= TRUE, TOL=.Machine$double.eps^{1/3}){
  
  if(min(abs(a1),abs(a2)) < TOL){stop("both a1 and a2 !=0")}
  if(any(abs(a1-a2) < TOL)) stop("a1!=a2")
  
  # evaluate d/da2 int K(a1,w)K(a2,w) p(w)dw
  wei=matrix(nrow = length(a1),ncol = 3)
  wei[,1] = (1+a2^2)/(2*a2^2*(q_alpha(a1)-q_alpha(a2))) + (1-a2^2)^2 / (4*a2^3*(q_alpha(a1)-q_alpha(a2))^2)
  wei[,2] = -(1-a1^2)*(1-a2^2)/(4*a1*a2^2*(q_alpha(a1)-q_alpha(a2))^2)
  wei[,3] = (1-a1^2)/(2*a1*(q_alpha(a1)-q_alpha(a2)))
  
  # tilde_ck = compute_tilde_ck(tilde_bk)
  # sum_1 = Xtr_cpp(x = c(a1,a2),a = tilde_bk,standardization = FALSE)
  sum_1 = matrix(int_Ka1p(a1 = c(a1,a2),
                          tilde_bk = tilde_bk,approx=approx),ncol=2)
  sum_2 = int_Ka1p(a1 = a2,tilde_bk = tilde_ck,approx=approx)
  
  apply(wei*cbind(sum_1,sum_2),1,sum)
  
}


int_Ka1Ka2pow2p = function(a1,a2,tilde_bk,tilde_ck,approx=TRUE,TOL = .Machine$double.eps^{1/3}){
  # evaluate (1/2pi) int K(a1,w)K(a2,w)^2 p(w)dw
  if(any(abs(a2)< TOL)) stop("a2 has to be always non-zero")
  
  
  val = vector(mode = "numeric", length = length(a1))
  
  ind_1 = which(abs(a1)<TOL)
  if(length(ind_1)>0){
    # when a1 ==0  (1/2pi) int K(a2,w)^2 p(w)dw
    val[ind_1] = int_Ka1pow2p(a1 = a2[ind_1],tilde_bk = tilde_bk, tilde_ck = tilde_ck,approx=approx)
  }
  
  # when a1!=0 and a1!=a2  
  ind_2 = which(abs(a1)>=TOL & abs(a1-a2)>=TOL)
  if(length(ind_2)>0){
    # (1/2pi) int K(a1,w) K(a2,w)^2 p(w)dw  
    term1 = (a2[ind_2]^2+1)/(1-a2[ind_2]^2)*int_Ka1Ka2p(a1 = a1[ind_2],a2 = a2[ind_2],tilde_bk = tilde_bk, approx=approx) 
    term2 = a2[ind_2]*g_int_Ka1Ka2p(a1 = a1[ind_2],a2 = a2[ind_2],tilde_bk = tilde_bk,tilde_ck = tilde_ck, approx=approx)
    val[ind_2] = term1+term2
  }
  
  # when a1!=0 and a1==a2  
  ind_3 = which(abs(a1)>=TOL & abs(a1-a2)< TOL)
  if(length(ind_3)>0){
    # (1/2pi) int K(a2,w)^3 p(w)dw  
    tilde_dk = compute_tilde_ck(tilde_ck)
    val[ind_3] = int_Ka1pow3p(a1 = a2[ind_3],tilde_bk = tilde_bk, tilde_ck = tilde_ck, tilde_dk = tilde_dk,approx=approx)
  }
  
 
  return(val)
  
}
 
#### discrete convolutions ####
h_abck<-function(a,b,c,k){
  if (k==0){
    return((-1 - b*c + a*(b + c)*(-1 + b*c) + a^2*b*c*(1 + b*c))/
             ((-1 + a*b)*(-1 + a*c)*(-1 + b*c)))
  }
  if (k==1){
    return(-((a + b + c - a*b*c*(b*c + a*(b + c)))/((-1 + a*b)*(-1 + a*c)*
                                                      (-1 + b*c))))
  }
  if (k==2){
    return(((-b)*c - c^2 + a*(b + c)*(-1 + b*c) + b^2*(-1 + c^2) + 
              a^2*(-1 + b^2 + b*c + c^2))/((-1 + a*b)*(-1 + a*c)*
                                             (-1 + b*c)))
  }
  if (k==3){
    return((a^3*(-1 + b^2)*(-1 + b*c)*(-1 + c^2) + 
              (b + c)*(-c^2 + b^2*(-1 + c^2)) - 
              a*(-1 + b*c)*((-b)*c - c^2 + b^2*(-1 + c^2)) - 
              a^2*(b + c)*(1 - b*c - c^2 + b^2*(-1 + c^2)))/
             ((-1 + a*b)*(-1 + a*c)*(-1 + b*c)))
  }
  if (k==4){
    return((-a^4 - a^3*b - a^2*b^2 + a^4*b^2 - a*b^3 + a^3*b^3 - b^4 + 
              a^2*b^4 - a^3*c - a^2*b*c + a^4*b*c - a*b^2*c + 
              2*a^3*b^2*c - b^3*c + 2*a^2*b^3*c - a^4*b^3*c + a*b^4*c - 
              a^3*b^4*c - a^2*c^2 + a^4*c^2 - a*b*c^2 + 2*a^3*b*c^2 - 
              b^2*c^2 + 3*a^2*b^2*c^2 - a^4*b^2*c^2 + 2*a*b^3*c^2 - 
              2*a^3*b^3*c^2 + b^4*c^2 - a^2*b^4*c^2 - a*c^3 + a^3*c^3 - 
              b*c^3 + 2*a^2*b*c^3 - a^4*b*c^3 + 2*a*b^2*c^3 - 
              2*a^3*b^2*c^3 + b^3*c^3 - 2*a^2*b^3*c^3 + a^4*b^3*c^3 - 
              a*b^4*c^3 + a^3*b^4*c^3 - c^4 + a^2*c^4 + a*b*c^4 - 
              a^3*b*c^4 + b^2*c^4 - a^2*b^2*c^4 - a*b^3*c^4 + a^3*b^3*c^4)/
             ((-1 + a*b)*(-1 + a*c)*(-1 + b*c)))
  }
}

h_a2bck<-function(a,b,c,k){
  if (k==0){
    return((1 + b*c - 2*a*(b + c)*(-1 + b*c) - 2*a^5*b*c*(b + c)*
              (-1 + b*c) - a^6*b^2*c^2*(1 + b*c) + 
              2*a^3*(b + c)*(-1 + b*c)*(1 + b*c) + 
              a^2*(1 + b^3*c - c^2 + b*c*(-5 + c^2) - b^2*(1 + 2*c^2)) + 
              a^4*(-c^2 + b*c*(2 + c^2) + b^2*(-1 + 5*c^2) + 
                     b^3*(c - c^3)))/((-1 + a^2)*(-1 + a*b)^2*(-1 + a*c)^2*
                                        (-1 + b*c)))
  }
  if (k==1){
    return((b + c - a^6*b^2*c^2*(b + c) + a^2*b*c*(b + c)*(-4 + b*c) + 
              2*a^3*(b + c)^2*(-1 + b*c) + a^4*(b + c)*(-1 + 4*b*c) + 
              a*(2 - 2*b^2*c^2) + a^5*(2*b*c - 2*b^3*c^3))/
             ((-1 + a^2)*(-1 + a*b)^2*(-1 + a*c)^2*(-1 + b*c)))
  }
  if (k==2){
    return((b^2 + b*c + c^2 - b^2*c^2 - 2*a*(b + c)*(-1 + b*c) - 
              2*a^5*b*c*(b + c)*(-1 + b*c) + 2*a^3*(b + c)*(-1 + b*c)*
              (1 + b*c) - a^6*b*c*(-1 + b^2 + b*c + c^2) + 
              a^2*(3 - 3*b^2 - 4*b*c - 3*c^2 + b^3*c^3) + 
              a^4*(-1 + b*c*(4*b*c + 3*c^2 - 3*b^2*(-1 + c^2))))/
             ((-1 + a^2)*(-1 + a*b)^2*(-1 + a*c)^2*(-1 + b*c)))
  }
  if (k==3){
    return((b^3 + b^2*c + b*c^2 - b^3*c^2 + c^3 - b^2*c^3 - 
              a^6*(b + c)*(-1 + b^2 + c^2) - 2*a^5*(-1 + b*c)*
              (-1 + b^2 + b*c + c^2) + 2*a*(-1 + b*c)*
              ((-b)*c - c^2 + b^2*(-1 + c^2)) - 
              a^4*(b + c)*(4 - 4*b*c - 3*c^2 + 3*b^2*(-1 + c^2)) - 
              2*a^3*(-1 + b*c)*(2 - 2*b*c - 3*c^2 + b^2*(-3 + 2*c^2)) + 
              a^2*(b + c)*(3 - 4*b*c - 3*c^2 + b^2*(-3 + 4*c^2)))/
             ((-1 + a^2)*(-1 + a*b)^2*(-1 + a*c)^2*(-1 + b*c)))
  }
  if (k==4){
    return((5*a^4 - 3*a^6 + 4*a^3*b - 6*a^5*b + 2*a^7*b + 3*a^2*b^2 - 
              9*a^4*b^2 + 4*a^6*b^2 + 2*a*b^3 - 6*a^3*b^3 + 6*a^5*b^3 - 
              2*a^7*b^3 + b^4 - 3*a^2*b^4 + 3*a^4*b^4 - a^6*b^4 + 
              4*a^3*c - 6*a^5*c + 2*a^7*c + 3*a^2*b*c - 12*a^4*b*c + 
              8*a^6*b*c - a^8*b*c + 2*a*b^2*c - 12*a^3*b^2*c + 
              14*a^5*b^2*c - 4*a^7*b^2*c + b^3*c - 7*a^2*b^3*c + 
              15*a^4*b^3*c - 8*a^6*b^3*c + a^8*b^3*c - 2*a*b^4*c + 
              6*a^3*b^4*c - 6*a^5*b^4*c + 2*a^7*b^4*c + 3*a^2*c^2 - 
              9*a^4*c^2 + 4*a^6*c^2 + 2*a*b*c^2 - 12*a^3*b*c^2 + 
              14*a^5*b*c^2 - 4*a^7*b*c^2 + b^2*c^2 - 10*a^2*b^2*c^2 + 
              19*a^4*b^2*c^2 - 9*a^6*b^2*c^2 + a^8*b^2*c^2 - 4*a*b^3*c^2 + 
              14*a^3*b^3*c^2 - 14*a^5*b^3*c^2 + 4*a^7*b^3*c^2 - b^4*c^2 + 
              4*a^2*b^4*c^2 - 6*a^4*b^4*c^2 + 4*a^6*b^4*c^2 - 
              a^8*b^4*c^2 + 2*a*c^3 - 6*a^3*c^3 + 6*a^5*c^3 - 2*a^7*c^3 + 
              b*c^3 - 7*a^2*b*c^3 + 15*a^4*b*c^3 - 8*a^6*b*c^3 + 
              a^8*b*c^3 - 4*a*b^2*c^3 + 14*a^3*b^2*c^3 - 14*a^5*b^2*c^3 + 
              4*a^7*b^2*c^3 - b^3*c^3 + 8*a^2*b^3*c^3 - 15*a^4*b^3*c^3 + 
              7*a^6*b^3*c^3 - a^8*b^3*c^3 + 2*a*b^4*c^3 - 6*a^3*b^4*c^3 + 
              6*a^5*b^4*c^3 - 2*a^7*b^4*c^3 + c^4 - 3*a^2*c^4 + 
              3*a^4*c^4 - a^6*c^4 - 2*a*b*c^4 + 6*a^3*b*c^4 - 
              6*a^5*b*c^4 + 2*a^7*b*c^4 - b^2*c^4 + 4*a^2*b^2*c^4 - 
              6*a^4*b^2*c^4 + 4*a^6*b^2*c^4 - a^8*b^2*c^4 + 2*a*b^3*c^4 - 
              6*a^3*b^3*c^4 + 6*a^5*b^3*c^4 - 2*a^7*b^3*c^4 - 
              a^2*b^4*c^4 + 3*a^4*b^4*c^4 - 3*a^6*b^4*c^4 + a^8*b^4*c^4)/
             ((-1 + a^2)*(-1 + a*b)^2*(-1 + a*c)^2*(-1 + b*c)))
  }
}