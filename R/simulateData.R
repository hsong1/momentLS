#' Simulate some reversible Markov chains
#'@import doSNOW
#'@import doParallel
#'@import parallel
#'@import foreach
#'@import stats
#'@useDynLib momentLS
#'
#'@examples
#' # Simulate an AR1 chain
#' chainParams = list(type="AR",  M = 10000, rho = 0.5)
#' ch_ar1 = generateChain(chainParams)
#' 
#' # Simulate a finite state space MH chain
#' chainParams = list(type="MH",  M = 10000,  nStates = 100, g = matrix(rnorm(100*1),ncol=1), discreteMC = simulate_discreteMC(nStates = 100), d = NULL)
#' ch_mh = generateChain(chainParams)
#'@export
generateChain = function(chainParams){
  # chainParams = list(type="AR",  M = 10000, rho = rho)
  # chainParams = list(type="MH",  M = 10000,  nStates = 100, g = matrix(rnorm(100*1),ncol=1), discreteMC = simulate_discreteMC(nStates = nStates), d = NULL)
  
  M = chainParams$M
  if(chainParams$type=="AR"){
    rho = chainParams$rho
    x = generateAR1Chain(M = M, rho = rho)
    F_weights=x$F_weights; F_support = x$F_support
    varTruth = sum(F_weights*(1+F_support)/(1-F_support))
    opt=list()
    x = x$x
  }else if(chainParams$type =="MH"){
    x= generateMHChain(M = M,nStates = chainParams$nStates, discreteMC = chainParams$discreteMC, g = chainParams$g, d = chainParams$d)
    F_weights=x$F_weights; F_support = x$F_support
    
    if(is.null(chainParams$d)){chainParams$d=1}
    if(chainParams$d==1){
      varTruth = sum(F_weights*(1+F_support)/(1-F_support))
    }else if(chainParams$d>1){
      varTruth = mtvAsympVariance(coefs = x$g_coefs, support = F_support)
    }
    
    # g_coef_j= <g,Phi_j>_pi
    opt=list(g_coefs=x$g_coefs)
    x = x$x
  }
  
  return(list(type=chainParams$type, x=x, varTruth = varTruth,
              F_weights=F_weights,F_support = F_support, opt = opt))
}

mtvAsympVariance = function(coefs,support){
  if(is.null(dim(coefs))){coefs= matrix(coefs,ncol=1)}
  Reduce("+", lapply(1:nrow(coefs), function(k) outer(coefs[k,],coefs[k,])*(1+support[k])/(1-support[k]) ))
}

#'@export
generateAR1Chain = function(M,rho){
  
  ### Sample Chain ###
  x = numeric(M); x[1] = rnorm(1, mean=0, sd = sqrt(1/(1-rho^2)))
  for(t in 2:M){
    x[t] = rho*x[t-1] + rnorm(1)
  }
  
  return(list(x=x, F_weights=1/(1-rho^2),F_support=rho ))
}


#'@export
generateMHChain = function(M,nStates, discreteMC, g = NULL, d = NULL){
  
  if(is.null(discreteMC)){
    warning("invariant pi and MH Kernel Q simulated")
    discreteMC = simulate_discreteMC(nStates = nStates) # simulate pi_vec, Q
  }
  
  v = with(discreteMC, simulateChain(pi_vec,Q,M)) # simulate length n chain based on pi_vec,Q
  
  # compute g(Vt) 
  if(is.null(g)){
    if(is.null(d)){d=1}
    warning("g simulated from N(0,Id)")
    g = matrix(rnorm(nStates*d),ncol=d) # simulate g
  }else{
    if( is.null(dim(g))){
      if(length(g)!=nStates){stop("length(g) needs to be the same as nStates")}else{
        g = matrix(g,ncol=1)
      }
    }else if(nrow(g)!=nStates){stop("nrow(g) needs to be the same as nStates")}
  }
  
  trueRepMeasure = with(discreteMC, discreteRepresentingMeasure(Q,pi_vec,g))
  
  x = g[v,]
  # x=with(discreteMC,simulateChain(pi_vec,Q,M,g)$x)
  return(list(x=x, F_weights = trueRepMeasure$weights, F_support = trueRepMeasure$support, g_coefs = trueRepMeasure$coefs))
}



#'@export
simulate_discreteMC = function(nStates){
  
  # Simulate random target distribution (pi_vec) and Q matrix based on random proposal distribution
  
  # number of states
  N=nStates
  
  # target vector pi = distribution on state N
  pi_vec=runif(N)
  pi_vec=pi_vec/sum(pi_vec)
  
  # generate Q distribution
  Q=matrix(0,ncol=N,nrow=N)
  
  # random proposal distribution
  proposal=matrix(0,N,N)
  for (i in 1:N){
    proposal[i,]=runif(N)
    proposal[i,]=proposal[i,]/sum(proposal[i,])
  }
  
  # Metropolis - Hastings Kernel Q (based on target vector pi_vec, and proposal distribution)
  for (i in 1:N){
    for (j in 1:N){
      if (j!=i){
        if (proposal[i,j]>0){
          Aij=min(1,pi_vec[j]/pi_vec[i]*proposal[j,i]/proposal[i,j])
          Q[i,j]=Aij*proposal[i,j]
          Q[i,i]=Q[i,i]+(1-Aij)*proposal[i,j]
        }
      }
      if (j==i){
        Q[i,i]=Q[i,i]+proposal[i,j]
      }
    }
  }
  #return stationary pi_vec, proposal P, and kernel Q
  # g = rnorm(nStates)
  return(list(Q=Q,pi_vec=pi_vec,proposal = proposal))
}

simulateChain<-function(pi_vec,Q,n){
  v=rep(0,n)
  v[1]=which(rmultinom(1,1,prob=pi_vec)!=0)
  for (i in 2:n){
    v[i]=which(rmultinom(1,1,prob=Q[v[i-1],])!=0)
  }
  
  return(v)
}

# mcsettings = simulate_discreteMC(nStates)
# Q = mcsettings$Q
# pi_vec = mcsettings$pi_vec
# d=2
# g = matrix(rnorm(nStates*d),ncol=d)

discreteRepresentingMeasure = function(Q,pi_vec,g){
  piCrossproduct = function(u,v,pi_vec){
    sum(u*v*pi_vec)
  }
  
  eigs=eigen(Q)
  vectors=eigs$vectors
  values=eigs$values
  
  # eigenvectors Phi_1 ,..., Phi_{nStates}
  # <Phi_i,Phi_j>_pi = 1[i=j]
  PhiVectors = sapply(1:ncol(vectors),function(i){
    vectors[,i]/sqrt(piCrossproduct(vectors[,i],vectors[,i],pi_vec))
  })

  
  # For gamma(k) = <Q^k bar{g}, bar{g}>_pi, 
  # the representing measure is sum of point mass at eig$values with weights = <g,Phi_j>_pi^2
  # supports
  meanZero=which(values<(1-10^-8)) # eigenvectors corr. to values != 1.
  F_support=values[meanZero]
  
  # weights
  # g = sum_{j=1}^d <g,Phi_j>_pi Phi_j
  # coef_j = <g,Phi_j>_pi, weights = coef_j^2
  if(is.null(dim(g))){d=1}else{d = ncol(g)}
  coefs     = matrix(nrow = length(meanZero),ncol = d)
  F_weights = matrix(nrow = length(meanZero),ncol = d)
  # k=1
  for(k in 1:ncol(g)){
    coef = sapply(1:ncol(PhiVectors),function(i) piCrossproduct(g[,k],PhiVectors[,i],pi_vec))
    coefs[, k] = coef[meanZero]
    F_weights[,k] = coefs[,k]^2
  }
  if(dim(F_weights)[2]==1){dim(F_weights)=NULL}
  if(dim(coefs)[2]==1){dim(coefs)=NULL}
  return(list(support= F_support, weights = F_weights, coefs= coefs))
}



