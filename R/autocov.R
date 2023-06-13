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


correlogram_from_x <-
  function(x,center=TRUE){
    M=length(x)
    r=autocov(x,center=center)
    r_full=c(r,r[M:2])
    return(Re(fft(r_full))) #the imaginary components are all 0
  }

autocov_from_correlogram <-
  function(correlogram){
    M_full=length(correlogram)
    M=(M_full+1)/2
    return((Re(fft(correlogram))/M_full)[1:M]) #imaginary components all 0
  }
