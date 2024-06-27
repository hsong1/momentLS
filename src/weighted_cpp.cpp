// we only include RcppEigen.h which pulls Rcpp.h in for us
#include <RcppEigen.h>
#include <iostream>

// via the depends attribute we tell Rcpp to create hooks for
// RcppEigen so that the build process will know what to do
//
// [[Rcpp::depends(RcppEigen)]]

using namespace Eigen;
using namespace std;


//' @export
 // [[Rcpp::export]]
 Eigen::ArrayXd poissonKernel_cpp(double rho, Eigen::ArrayXd & wseq){
   return (1-std::pow(rho,2)) * (1-2*cos(wseq)*rho + std::pow(rho,2)).inverse();
 }

//' @export
// [[Rcpp::export]]
Eigen::ArrayXd phi_cpp(Eigen::ArrayXd & wseq, Eigen::ArrayXd support, Eigen::ArrayXd weights){
  
  Eigen::ArrayXd out(wseq.size());
  out.setZero();
  
  for(int i=0; i < support.size(); i++){
    out += poissonKernel_cpp(support[i],wseq) * weights[i];
  }
  return out;
}

//' @export
// [[Rcpp::export]]
Eigen::MatrixXd makeXtX_w_cpp(Eigen::ArrayXd & alphaGrid, Eigen::ArrayXd & wseq, Eigen::ArrayXd & phi_wseq, bool diag){
  
  
  Eigen::ArrayXd wei;
  wei = phi_wseq.pow(2).inverse(); // 1/phi(wseq)^2
  wei = wei / wei.size();
  
  Eigen::MatrixXd Mat(alphaGrid.size(),alphaGrid.size());
  Mat.setZero();
  
  for(int i=0; i < alphaGrid.size(); i++){
    for(int j=(i+1); j < alphaGrid.size(); j++){
      Mat(i,j) = (poissonKernel_cpp(alphaGrid[i],wseq) * poissonKernel_cpp(alphaGrid[j],wseq) * wei).sum();
      Mat(j,i) = Mat(i,j);
    }
  }
  
  if(diag){
    for(int i=0; i<alphaGrid.size(); i++){
      Mat(i,i) = (poissonKernel_cpp(alphaGrid[i], wseq).pow(2)*wei).sum();
    }
  }
  
  return(Mat);
  
}


//' @export
// [[Rcpp::export]]
Eigen::ArrayXd computeXtr_w_cpp(Eigen::ArrayXd & alphaGrid, 
                                Eigen::ArrayXd & FT_r_wseq,
                                Eigen::ArrayXd & wseq, 
                                Eigen::ArrayXd & phi_wseq){
  
  
  Eigen::ArrayXd wei;
  wei = phi_wseq.pow(2).inverse(); // 1/phi(wseq)^2
  wei = wei / wei.size();
  
  Eigen::ArrayXd out;
  out.resize(alphaGrid.size());
  out.setZero();
  
  for(int i=0; i<alphaGrid.size(); i++){
    out[i] = (poissonKernel_cpp(alphaGrid[i],wseq) * FT_r_wseq * wei).sum();
  }
  
  return(out);
  
}

/*** R
library(momentLS)
M0=1e4
wseq = 2*pi*(0:(M0-1))/M0
rho = 0.3

max(abs(poissonKernel_cpp(wseq = wseq,rho = rho)- poissonKernel(rho = rho,w = wseq)))

supp = c(0.1,0.3)
weight = c(0.2,0.5)
max(abs(apply(sapply(1:length(supp), function(i) poissonKernel(w = wseq,rho = supp[i])*weight[i]),1,sum)-phi_cpp(wseq,supp,weight)))


alphaGrid  = makeGrid(upper_threshold = 0.9,nX = 101,cm = F,scale = "log")
phi_wseq = phi_cpp(wseq,supp,weight)
a=makeXtX_w_cpp(alphaGrid,wseq = wseq,phi_wseq = phi_wseq,diag = T)

Mat = matrix(ncol=length(alphaGrid), nrow=length(alphaGrid))
for(i in 1:length(alphaGrid)){
  for(j in i:length(alphaGrid)){
    alpha1 = alphaGrid[i]; alpha2 = alphaGrid[j]
    Mat[i,j] = sum(poissonKernel(w = wseq, rho = alpha1)* poissonKernel(w = wseq,rho = alpha2) / phi_wseq^2) / length(wseq)
    Mat[j,i]= Mat[i,j]
  }
}

max(abs(a-Mat))

library(momentLS)
# chainParams = list(type="MH",  M = 1000,  nStates = 100, g = matrix(rnorm(100*1),ncol=1), discreteMC = simulate_discreteMC(nStates = 100), d = NULL)
chainParams = list(type="AR",  M = 10000,  rho=0.9)
ch = generateChain(chainParams)
# fit momentLS
r = autocov(ch$x);
dhat = tune_delta(x = ch$x,nSplits = 5)$delta*0.8
rho_hat = SR1(r = r,delta = dhat)

Xtr_num1 = c()
raug = c(r, rep(0, length(wseq)-length(r)))
Frm_wseq = Re(2*fft(raug)) - r[1]

alpha1 = alphaGrid[1]
v = sapply(alphaGrid, function(alpha1) sum(poissonKernel(w = wseq,rho = alpha1)*Frm_wseq/phi_wseq^2) / length(wseq))-
  computeXtr_w_cpp(alphaGrid = alphaGrid,FT_r_wseq = Frm_wseq,wseq = wseq,phi_wseq = phi_wseq)
max(abs(v))
*/
