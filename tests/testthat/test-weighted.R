library(momentLS)
#### Examples #####
set.seed(1234)
# Simulate a finite state space MH chain
#chainParams = list(type ="AR", M=10000, rho=0.9)
chainParams = list(type="MH",  M = 10000,  nStates = 100, g = matrix(rnorm(100*1),ncol=1), discreteMC = simulate_discreteMC(nStates = 100), d = NULL)
ch = generateChain(chainParams)
dhat = tune_delta(ch$x,nSplits = 5)
r = autocov(ch$x)

# three examples
m_uw1 = SR1(r = r,delta = dhat$delta*0.8,n_alphas = 101)
# single support
m_uw2 = SR1(r = r,delta = dhat$delta*0.8,n_alphas = 101)
m_uw2$support = m_uw2$support[2]
m_uw2$weights = m_uw2$weights[2]
# add 0
m_uw3 = m_uw1
m_uw3$support = c(0,m_uw3$support)
m_uw3$weights = c(0.1,m_uw3$weights)

# test sequence
w = runif(5,-pi,pi)

############
#mu = m_uw3
mu_list=list(m_uw1,m_uw2,m_uw3)
mu_names = c("m_uw1","m_uw2","m_uw3")
wseq = seq(-pi,pi,length.out=400000)

#hatrM_seq = eval(w = wseq,x = cosinePoly(r)) # time-consuming
#save(hatrM_seq,file="tests/testthat/testdata/hatrM_seq.Rdata")
#load("tests/testthat/testdata/hatrM_seq.Rdata")
load(testthat::test_path("testdata","hatrM_seq.Rdata"))
i=3
for(i in 1:length(mu_list)){
  mu = mu_list[[i]]
  mu_name = mu_names[i]
  
  # Xtr_w
  # numerical evaluation
  phi_wseq = sapply(wseq, function(wi) sum(mu$weights * poissonKernel(rho = mu$support, w = wi)))
  alphaGrid = mu$alphaGrid
  eval1 = computeXtr_w(alphaGrid,r = r,m_uw = mu)
  eval2 = sapply(alphaGrid, function(a) sum(poissonKernel(rho = a,w = wseq)*hatrM_seq / phi_wseq^2)/length(wseq))
  print(max(abs(eval1-eval2)))
  testthat::expect_equal(eval1, eval2,tolerance = 1e-4)
  plot(eval1,eval2)
  abline(0,1)
  
  # XtX_w
  # XtX_w_num = matrix(ncol=length(alphaGrid), nrow=length(alphaGrid))
  # for(i in 1:length(alphaGrid)){
  #   for(j in i:length(alphaGrid)){
  #     alpha1 = alphaGrid[i]; alpha2 = alphaGrid[j]
  #     XtX_w_num[i,j] = sum(poissonKernel(alpha1,wseq)* poissonKernel(alpha2,wseq) / phi_wseq^2) / length(wseq)
  #     XtX_w_num[j,i]= XtX_w_num[i,j]
  #   }
  # }
  # save(XtX_w_num, file = paste0("tests/testthat/testdata/XtX_w_",mu_name,".Rdata"))
  # load(paste0("tests/testthat/testdata/XtX_w_",mu_name,".Rdata"))
  load(testthat::test_path("testdata",paste0("XtX_w_",mu_name,".Rdata")))
  eval1 = XtX_w_num
  eval2 = makeXtX_w(alphaGrid,mu)
  
  print(max(abs(eval1-eval2)))
  testthat::expect_equal(eval1, eval2,tolerance = 1e-4)
  
  plot(eval1,eval2, pch="*")
  abline(0,1)
  
  
  # numerical checks
  eval1=computeXtr_w(alphaGrid = mu$alphaGrid,r,mu)
  suppressWarnings({eval2=computeXtr_w(alphaGrid = mu$alphaGrid,r,mu,
                     numchck_control = list(chck_TOL = 0.01,rplc_TOL=0,N_wseq=100000))})
  print(max(abs(eval1-eval2)))
  testthat::expect_equal(eval1, eval2)
  
  eval1=makeXtX_w(alphaGrid = mu$alphaGrid,m_uw = mu)
  suppressWarnings({eval2=makeXtX_w(alphaGrid = mu$alphaGrid,mu,
                  numchck_control = list(chck_TOL = 0.01,rplc_TOL=0,N_wseq=100000))})
  print(max(abs(eval1-eval2)))
  testthat::expect_equal(eval1, eval2)
}


fit_w=SR1_w(r,delta = 0.1,n_alphas = 101,comp_method = "exact")
N = 1e5
phi=list(wseq = 2*pi*(0:(N-1))/N,
         phi_wseq = phi_cpp(wseq = 2*pi*(0:(N-1))/N,
                            support = fit_w$phi$support,
                            weights = fit_w$phi$weights))

fit_w2=SR1_w(r = r,delta = 0.1,phi = phi, n_alphas = 101)
#save(fit_w2, file = paste0("tests/testthat/testdata/fit_w2.Rdata"))
#load(testthat::test_path("testdata","fit_w2.Rdata"))

eval1=fit_w
eval2=fit_w2
testthat::expect_equal(eval1$weights,eval2$weights,tolerance = 1e-6)
testthat::expect_equal(eval1$support,eval2$support,tolerance = 1e-6)

fit_w3= SR1_w(r = r,delta = 0.1,comp_method = "num",
              n_alphas = 101,phi =fit_w$phi,n_phi = 1e5)

eval1=fit_w2
eval2=fit_w3
testthat::expect_equal(eval1$weights,eval2$weights,tolerance = 1e-6)
testthat::expect_equal(eval1$support,eval2$support,tolerance = 1e-6)
