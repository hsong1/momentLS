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
Rcpp::List compute_XtXw_Xtrw_cpp(
    Eigen::ArrayXd & alphaGrid,
    Eigen::ArrayXd & FT_r_wseq,
    Eigen::ArrayXd & wseq,
    Eigen::ArrayXd & phi_wseq){
  
  int n = wseq.size();
  int s = alphaGrid.size();
  
  Eigen::MatrixXd Xw(n,s);
  Eigen::ArrayXd invPhi = phi_wseq.inverse(); // 1/phi
  // X[,i] = K(wseq, alpha_i) / phi(wseq)
  for(int i = 0; i < s; i++){
    Xw.col(i) = (poissonKernel_cpp(alphaGrid[i], wseq) * invPhi).matrix();
  }
  
  // XtX_w = (1/n) * Xw'Xw
  Eigen::MatrixXd XtX_w = (Xw.transpose() * Xw) / n;
  
  // Xtr_w = (1/n) * Xw' yw  where yw = FT_r / phi(wseq)
  Eigen::VectorXd yw = (FT_r_wseq * invPhi).matrix();
  Eigen::VectorXd Xtr_w = (Xw.transpose() * yw) / n;
  
  return Rcpp::List::create(
    Rcpp::Named("XtX_w") = XtX_w,
    Rcpp::Named("Xtr_w") = Xtr_w
  );
}

// // [[Rcpp::export]]
// Eigen::MatrixXd makeXtX_w_cpp(Eigen::ArrayXd & alphaGrid, Eigen::ArrayXd & wseq, Eigen::ArrayXd & phi_wseq, bool diag){
// 
// 
//   Eigen::ArrayXd wei;
//   wei = phi_wseq.pow(2).inverse(); // 1/phi(wseq)^2
//   wei = wei / wei.size();
// 
//   Eigen::MatrixXd Mat(alphaGrid.size(),alphaGrid.size());
//   Mat.setZero();
// 
//   for(int i=0; i < alphaGrid.size(); i++){
//     for(int j=(i+1); j < alphaGrid.size(); j++){
//       Mat(i,j) = (poissonKernel_cpp(alphaGrid[i],wseq) * poissonKernel_cpp(alphaGrid[j],wseq) * wei).sum();
//       Mat(j,i) = Mat(i,j);
//     }
//   }
// 
//   if(diag){
//     for(int i=0; i<alphaGrid.size(); i++){
//       Mat(i,i) = (poissonKernel_cpp(alphaGrid[i], wseq).pow(2)*wei).sum();
//     }
//   }
// 
//   return(Mat);
// 
// }


// // [[Rcpp::export]]
// Eigen::ArrayXd computeXtr_w_cpp(Eigen::ArrayXd & alphaGrid, 
//                                 Eigen::ArrayXd & FT_r_wseq,
//                                 Eigen::ArrayXd & wseq, 
//                                 Eigen::ArrayXd & phi_wseq){
//   
//   
//   Eigen::ArrayXd wei;
//   wei = phi_wseq.pow(2).inverse(); // 1/phi(wseq)^2
//   wei = wei / wei.size();
//   
//   Eigen::ArrayXd out;
//   out.resize(alphaGrid.size());
//   out.setZero();
//   
//   for(int i=0; i<alphaGrid.size(); i++){
//     out[i] = (poissonKernel_cpp(alphaGrid[i],wseq) * FT_r_wseq * wei).sum();
//   }
//   
//   return(out);
//   
// }

/*** R
library(momentLS)
# data
# chainParams = list(type="MH",  M = 1000,  nStates = 100, g = matrix(rnorm(100*1),ncol=1), discreteMC = simulate_discreteMC(nStates = 100), d = NULL)
M0 = 1e4
chainParams = list(type="AR",  M = M0,  rho=0.9)
ch = generateChain(chainParams)
r = autocov(ch$x);

# Frequency domain transform 
wseq = 2*pi*(0:(M0-1))/M0
raug = c(r, rep(0, length(wseq)-length(r)))
Frm_wseq = Re(2*fft(raug)) - r[1]

# some weight
dhat = tune_delta(x = ch$x,nSplits = 5)$delta*0.8
m_hat = SR1(r = r,delta = dhat)
supp = m_hat$support; weight = m_hat$weights
phi_wseq = phi_cpp(wseq,supp,weight)

# grid
alphaGrid  = makeGrid(upper_threshold = 0.9,nX = 501,cm = F,scale = "log")


# test XtX_w and Xtr_w
out = compute_XtXw_Xtrw_cpp(alphaGrid = alphaGrid,FT_r_wseq = Frm_wseq,wseq = wseq, phi_wseq = phi_wseq)

# test XtX_w
Mat = matrix(ncol=length(alphaGrid), nrow=length(alphaGrid))
for(i in 1:length(alphaGrid)){
  for(j in i:length(alphaGrid)){
    alpha1 = alphaGrid[i]; alpha2 = alphaGrid[j]
    Mat[i,j] = sum(poissonKernel(w = wseq, rho = alpha1)* poissonKernel(w = wseq,rho = alpha2) / phi_wseq^2) / length(wseq)
    Mat[j,i]= Mat[i,j]
  }
}
max(abs(Mat-out$XtX_w))

# test Xtr_w
Xtr_num1 = c()
alpha1 = alphaGrid[1]
v = sapply(alphaGrid, function(alpha1) sum(poissonKernel(w = wseq,rho = alpha1)*Frm_wseq/phi_wseq^2) / length(wseq))-out$Xtr_w
max(abs(v))

# # compare with previous functions
# max(abs(out$XtX_w - makeXtX_w_cpp(alphaGrid,wseq,phi_wseq,diag=T)))
# max(abs(out$Xtr_w - computeXtr_w_cpp(alphaGrid,Frm_wseq,wseq,phi_wseq)))
# 
# test = function(alphaGrid,Frm_wseq,wseq,phi_wseq){
#   XtX_w= makeXtX_w_cpp(alphaGrid,wseq,phi_wseq,diag=T)
#   Xtr_w= computeXtr_w_cpp(alphaGrid,Frm_wseq,wseq,phi_wseq)
#   return(list(XtX_w,Xtr_w))
# }
# 
# microbenchmark::microbenchmark(
#   out = compute_XtXw_Xtrw_cpp(alphaGrid = alphaGrid,FT_r_wseq = Frm_wseq,wseq = wseq, phi_wseq = phi_wseq),
#   out2 = test(alphaGrid,Frm_wseq,wseq,phi_wseq),
#   times=10
# )

*/
