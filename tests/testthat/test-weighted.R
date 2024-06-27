library(momentLS)
#### Examples #####
set.seed(123)
# Simulate a finite state space MH chain
#chainParams = list(type ="AR", M=10000, rho=0.9)
chainParams = list(type="MH",  M = 10000,  nStates = 100, g = matrix(rnorm(100*1),ncol=1), discreteMC = simulate_discreteMC(nStates = 100), d = NULL)
ch = generateChain(chainParams)
dhat = tune_delta(ch$x,nSplits = 5)
r = autocov(ch$x)
m_uw1 = SR1(r = r,delta = dhat$delta*0.8,n_alphas = 101)
m_uw2 = SR1(r = r,delta = dhat$delta*0.8,n_alphas = 101)
m_uw3 = SR1(r = r,delta = dhat$delta*0.8,n_alphas = 101)

m_uw2$support = 0.1
m_uw2$weights = m_uw2$weights[2]

# m_uw3$support = c(0,m_uw3$support[1])
# m_uw3$weights = c(0.1,m_uw3$weights[1])

m_uw3$support = c(0,m_uw3$support)
m_uw3$weights = c(0.1,m_uw3$weights)

w = runif(5,-pi,pi)
mu = m_uw3
m_uw = mu

############
#source("weighted.R")

wseq = seq(-pi,pi,length.out=100000)
hatrM_seq= eval_cp(w = wseq,p = r,two_sided = T)
phi_wseq = sapply(wseq, function(wi) sum(mu$weights * poissonKernel(rho = mu$support, w = wi)))
alphaGrid = mu$alphaGrid
eval1 = computeXtr_w(alphaGrid,r = r,m_uw = mu)
eval2 = sapply(alphaGrid, function(a) sum(poissonKernel(rho = a,w = wseq)*hatrM_seq / phi_wseq^2)/length(wseq))
print(max(abs(eval1 - eval2)))

plot(eval1,eval2)
abline(0,1)



Mat = matrix(ncol=length(alphaGrid), nrow=length(alphaGrid))
for(i in 1:length(alphaGrid)){
  for(j in i:length(alphaGrid)){
    alpha1 = alphaGrid[i]; alpha2 = alphaGrid[j]
    Mat[i,j] = sum(poissonKernel(alpha1,wseq)* poissonKernel(alpha2,wseq) / phi_wseq^2) / length(wseq)
    Mat[j,i]= Mat[i,j]
  }
}

eval2 = makeXtX_w(alphaGrid,m_uw)

print(max(Mat-eval2))
plot(Mat,eval2, pch="*")
abline(0,1)


# numerical checks

eval1=computeXtr_w(alphaGrid = m_uw$alphaGrid,r,m_uw)
eval2=computeXtr_w(alphaGrid = m_uw$alphaGrid,r,m_uw,
                   numchck_control = list(chck_TOL = 0.01,rplc_TOL=0,N_wseq=100000))
max(abs(eval1-eval2))



eval1=makeXtX_w(alphaGrid = m_uw$alphaGrid,m_uw = m_uw)
eval2=makeXtX_w(alphaGrid = m_uw$alphaGrid,m_uw,
                numchck_control = list(chck_TOL = 0.01,rplc_TOL=0,N_wseq=100000))
max(abs(eval1-eval2))


rm(m_uw)

fit_w=SR1_w(r,delta = 0.1,n_alphas = 101)
fit_w$phi
N = 1e5
phi=list(wseq = 2*pi*(0:(N-1))/N,
         phi_wseq = phi_cpp(wseq = 2*pi*(0:(N-1))/N,support = fit_w$phi$support,weights = fit_w$phi$weights))
debugonce(SR1_w)
fit_w2=SR1_w(r = r,delta = 0.1,phi = phi,n_alphas = 101)

fit_w
fit_w2
