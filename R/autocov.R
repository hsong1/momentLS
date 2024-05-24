#' Compute the empirical auto (cross) covariances of chain x
#' for a vector x, returns vector r with r[k] = lag (k-1) autocovariance (r[1]=lag 0 autocovariance )
#' for n by p matrix x, returns array r with r[k,i,j]=lag (k-1) cross covariance of x[,i] with x[,j]
#' @param x input chain (vector or matrix)
#' @param center a logical value; TRUE for empirical mean centering
#' @import stats
#' @export
autocov<-function(x,center=TRUE){
  
  if (is.null(dim(x))){
    # vector
    if (center){xbar=x-mean(x)}else{xbar = x}
    n=length(xbar) # length of x
    xbar2=c(xbar,rep(0,n)) # zero-pad
    t1=fft(xbar2)
    norms=Re(t1)^2+Im(t1)^2
    unadj=Re(fft(norms,inverse=TRUE)[1:n])
    adj=unadj/(2*n^2)
    return(adj)
  }else{
    if(length(dim(x))!=2){stop("length(dim(x)) should be 2")}
    p=ncol(x)
    n=nrow(x)
    
    if (center){xbar=t(t(x)-apply(x,2,mean))}else{xbar = x}
    xbar2=rbind(xbar,matrix(0,n,p))
    t1=mvfft(xbar2)
    t1conj=Conj(t1)
    
    r=array(0,dim=c(n,p,p))
    f_ij=complex(n)
    
    #i=1; j = i+1
    for (i in 1:p){
      for (j in i:p){
        f_ij=t1conj[,i]*t1[,j]
        r[,i,j]=Re(fft(f_ij,inverse=TRUE)[1:n])/(2*n^2)
        
        if (i!=j){
          r[,j,i]=Re(fft(f_ij,inverse=FALSE)[1:n])/(2*n^2)
        }
      }
    }
    
    return(r)
  }
  stop("Requires vector or matrix input for x")
}
