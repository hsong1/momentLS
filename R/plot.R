#' @import graphics
#' @export
#' @method plot SRfit1
plot.SRfit1=function(x,...){
  m = x
  plot(m$alphaGrid, m$gradient, type="l")
  inds = which(m$alphaGrid%in%m$support)
  points(m$alphaGrid[inds], m$gradient[inds])
  abline(h = 0,lty=2)
}
