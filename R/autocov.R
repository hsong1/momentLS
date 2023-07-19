#'@import stats
#'@export
autocov <-
  function(x,center=TRUE){
    xbar=x
    if (center){xbar=x-mean(x)}
    n=length(xbar)
    xbar2=c(xbar,rep(0,n))
    t1=fft(xbar2)
    norms=Re(t1)^2+Im(t1)^2
    unadj=Re(fft(norms,inverse=TRUE)[1:n])
    adj=unadj/(2*n^2)
    return(adj)
  }

