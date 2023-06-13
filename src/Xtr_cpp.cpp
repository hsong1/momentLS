// we only include RcppEigen.h which pulls Rcpp.h in for us
#include <RcppEigen.h>

// via the depends attribute we tell Rcpp to create hooks for
// RcppEigen so that the build process will know what to do
//
// [[Rcpp::depends(RcppEigen)]]

//' @export
// [[Rcpp::export]]
Eigen::VectorXd Xtr_cpp(Eigen::VectorXd x,Eigen::VectorXd a, 
                        bool standardization= true, double epsilon=1e-10){
  int n=a.size();
  int J=std::ceil(std::log2(n))+1;
  
  
  int nRho=x.size();
  int m=std::ceil(std::exp(1+1/std::exp(1))*std::log(n/epsilon));
  
  Eigen::VectorXi index;
  index.setZero(nRho);
  
  
  for (int j=0;j<nRho;j++){
    index(j)=(int) std::floor(std::log2(1/(1-std::abs(x(j)))));
    if (index(j)>(J-1)){index(j)=J-1;};
  };
  
  Eigen::VectorXi occupied;
  occupied.setZero(J);
  
  
  for (int j=0; j<nRho;j++){
    if (x(j)>=0){
      int i=index(j);
      occupied(i)=1;
    };
  };
  
  Eigen::MatrixXd reducedVa;
  reducedVa.setZero(m,J);
  
  
  for (int i=0;i<J;i++){
    if (occupied(i)!=0){
      double rho0=1-std::pow(2,-(i+1));
      int kMax=std::min(n-1,(int) std::floor(std::pow(2,i+1)*std::log(n/epsilon)));
      for (int k=0;k<=kMax;k++){
        double currentVal=std::pow(rho0,k);
        for (int r=0; r<=std::min(m-1,k);r++){
          reducedVa(r,i)=reducedVa(r,i)+currentVal*a(k);
          currentVal=currentVal*(k-r)*(1-rho0)/((r+1)*(rho0));
        };
      };
    };
  };
  
  Eigen::VectorXd vals;
  vals.setZero(nRho);
  for (int j=0;j<nRho;j++){
    if (x(j)>=0){
      int i=index(j);
      
      double normalization=1;
      if (standardization){normalization=std::sqrt((1-std::pow(x(j),2))/(1+std::pow(x(j),2)));};
      
      for (int r=0;r<m;r++){
        vals(j)=vals(j)+normalization*std::pow(((x(j)-(1-std::pow(2,-(i+1))))/(std::pow(2,-(i+1)))),r)*reducedVa(r,i);
      };
      vals(j)=vals(j)*2-normalization*a(0);
    };
  };
  
  
  reducedVa.setZero(m,J);
  
  occupied.setZero(J);
  
  for (int j=0;j<nRho;j++){
    if (x(j)<0){
      int i=index(j);
      occupied(i)=1;
    }
  }
  
  
  for (int i=0;i<J;i++){
    if (occupied(i)!=0){
      double rho0=-(1-std::pow(2,-(i+1)));
      int kMax=std::min(n-1,(int) std::floor(std::pow(2,i+1)*std::log(n/epsilon)));
      for (int k=0;k<=kMax;k++){
        double currentVal=std::pow(rho0,k);
        for (int r=0; r<=std::min(m-1,k);r++){
          reducedVa(r,i)=reducedVa(r,i)+currentVal*a(k);
          currentVal=currentVal*(k-r)*(1-(-rho0))/((r+1)*(rho0));
        };
      };
    };
  };
  
  for (int j=0;j<nRho;j++){
    if (x(j)<0){
      int i=index(j);
      
      double normalization=1;
      if (standardization){normalization=std::sqrt((1-std::pow(x(j),2))/(1+std::pow(x(j),2)));};
      
      for (int r=0;r<m;r++){
        vals(j)=vals(j)+normalization*std::pow(((x(j)-(-(1-std::pow(2,-(i+1)))))/(std::pow(2,-(i+1)))),r)*reducedVa(r,i);
      };
      vals(j)=vals(j)*2-normalization*a(0);
    };
  };
  
  
  return(vals);
}