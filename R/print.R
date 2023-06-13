#' @export
#' @method print SRfit1
print.SRfit1=function(x,...,quote=FALSE){
  print(round(data.frame(support= x$support, weights = x$weights),4),quote=quote)
}
