#'@export
# subset "relevant" Fourier coefficients
# Fourier coefficients corresponding to w_1,...w_N 
# where w_j = 2*pi*j /M and N = M/2-1 if M is even or M = M/2-1/2 if M is odd
subsetFcoefs =  function(Fx){
  # input: Fx = fft([x[1],...,x[M]])
  # Fx contians FT of x at w_0,...,w_{M-1}
  M = length(Fx)
  # Re(Fx[1]) == sum(x) 
  N = floor((M-1)/2)
  # when M is even, Re(Fx[N+2]) == sum((-1)^{0:(M-1)}*x) 
  # take N coefficients starting from w_1 
  a_k = Re(Fx[(1+1):(N+1)]) 
  b_k = Im(Fx[(1+1):(N+1)])
  
  # check (when M is even)
  # max(abs(Re(Fx[(N+3):M])- a_k[(N:1)]))
  # max(abs(Im(Fx[(N+3):M])+ b_k[(N:1)]))

  # check (when M is odd)
  # max(abs(Re(Fx[(N+2):M])- a_k[(N:1)]))
  # max(abs(Im(Fx[(N+2):M])+ b_k[(N:1)]))
  
  return(list(a_k=a_k, b_k=b_k, M = M))
}

#'@export
# reconstruct original signal (other than mean)
# based on FT coefficients from w_1,...,w_N
reconstrTS = function(a_k, b_k, M){
  fullFx = vector(mode = "complex",length = M)
  N = floor((M-1)/2)
  fullFx[2:(N+1)] = complex(real = a_k,imaginary = b_k)
  if(M%%2==0){
    fullFx[(N+3):M] = complex(real = a_k[N:1], imaginary = -b_k[N:1])
  }else{
    fullFx[(N+2):M] = complex(real = a_k[N:1], imaginary = -b_k[N:1])
  }
  
  iFx = fft(fullFx,inverse = T)/M
  if(max(abs(Im(iFx))) > 1e-8) stop("iFx has to be a real-valued signal")
  return(Re(iFx))
}

# bootstrap based on Fourier coefficients
# a_k/sqrt(M), b_k/sqrt(M) a~ N(0, 1/2 Fgamma(w_k))
# Fgamma(w) = sum_k gamma(k) e^{-iwk}
#'@export
TFTboot = function(a_k, b_k, M, spectrf, type=c("RB","WB")){
  
  type = match.arg(type, choices = c("RB","WB"))
  
  # Obtain spectr = [Fgamma(w_k) ; for k = 1,...,N]
  N = floor((M-1)/2)
  freq = 2*pi*(1:N)/M; spectr = spectrf(freq)
  
  if(type=="RB"){
    stilde_vec = c(a_k/sqrt((M/2)*spectr), b_k/sqrt((M/2)*spectr))
    svec = (stilde_vec - mean(stilde_vec))/sd(stilde_vec) # standardization
    svec_resamp = sample(svec,size = length(svec),replace = T)
  }else{
    svec_resamp = rnorm(n = 2*N)
  }
  
  a_k_resamp = svec_resamp[1:N] * sqrt((M/2)*spectr)
  b_k_resamp = svec_resamp[(N+1):(2*N)]* sqrt((M/2)*spectr)
  return(list(a_k_boot= a_k_resamp, b_k_boot = b_k_resamp, M=M))
}

# poissonKernel(w,rho) =   (F rho^(|*|)) (w) = {sum_{t in Z} rho^{|t|} e^{-itw} = (1-rho^2) / (1-2cos(w)rho + rho^2)
#'@export
poissonKernel = function(w,rho) {(1-rho^2)/(1-2*cos(w)*rho + rho^2)}


# for m(k) = sum_j w_j (rho_j)^{|k|} 
# (Fm)(w) = sum_j w_j poissonKernel(rho_j,w)


#'@export
momentSpectrDens = function(momentSeq, FFreq){
  nSupports = length(momentSeq$support)
  res = matrix(nrow = length(FFreq),ncol = nSupports)
  for(i in 1:nSupports){
    res[,i] = poissonKernel(w = FFreq,rho = momentSeq$support[i])*momentSeq$weights[i]
  }
  # res = foreach (i = 1:nSupports,.combine = "cbind")%do%{
  #   poissonKernel(w = FFreq,rho = momentSeq$support[i])*momentSeq$weights[i]
  # }
  if(nSupports>1){
    res= apply(res,1,sum)
  }
  
  return(res)
}


#'@export
ft_seq = function(type, seq){

  # estimated F(seq)(FFreq) based on type, seq
  if(type == "momentseq"){
    M = seq$M
    FFreq = 2*pi*(0:(M-1))/M
    y = momentSpectrDens(momentSeq = seq,FFreq = FFreq)
    
  }else if(type == "L2seq"){
    M = length(seq)
    FFreq = 2*pi*(0:(M-1))/M
    y = Mod(fft(seq))
  }
  return(list(FFreq=FFreq, spectr = y))

}


# # check
# Fx_sub = subsetFcoefs(Fx)
# iFx= with(Fx_sub, reconstrTS(a_k,b_k,M))
# max(abs(autocov(iFx)- autocov(x)))



