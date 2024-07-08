set.seed(123)
# Simulate a finite state space MH chain
#chainParams = list(type ="AR", M=10000, rho=0.9)
chainParams = list(type="MH",  M = 10000,  nStates = 100, g = matrix(rnorm(100*1),ncol=1), discreteMC = simulate_discreteMC(nStates = 100), d = NULL)
ch = generateChain(chainParams)
dhat = tune_delta(ch$x,nSplits = 5)
r = autocov(ch$x)
m_uw = SR1(r = r,delta = dhat$delta*0.8,n_alphas = 101)
wseq = runif(min = -pi,max = pi,n=10) # sample w

phifct_out = phi_fct(m_uw)

# phi_factorization check
eval1=eval_fct(w = wseq,fct = phifct_out)
eval2 = apply(sapply(1:length(m_uw$weights), function(i) 
  m_uw$weights[i]*poissonKernel_cpp(rho = m_uw$support[i],wseq = wseq)),1,sum)

max(abs(eval1-eval2))
testthat::expect_equal(eval1, eval2)

# phi_inverse_pf check
phiInvfct_out=phi_inverse_pf(m_uw)
eval1=1/sapply(wseq,function(a) eval_fct(a,phi_fct(m_uw)))
frac_xi = t(sapply(wseq,function(w) poissonKernel(phiInvfct_out$xi,w)))
eval2=as.numeric(phiInvfct_out$C_phi_inv*eval(wseq, phiInvfct_out$phi_inv_cp)*frac_xi%*%phiInvfct_out$pf_coef)
max(abs(eval1-eval2))
testthat::expect_equal(eval1, eval2)

