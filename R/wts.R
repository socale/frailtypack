#' Identify weights
#' 
#' This is a special function used in the context of the joint frailty models
#' for data from nested case-control studies. It specifies weights defined by
#' using 'wts' function, and is used of 'frailtyPenal' formula for fitting
#' joint models.
#' 
#' 
#' @usage wts(x)
#' @param x A numeric variable which is supposed to indicate the weights
#' @return \item{x}{A variable identified as weights}
#' @seealso \code{\link{frailtyPenal}}
#' @keywords misc
#' @export
#' @examples
#' 
#' 
#' \dontrun{
#' 
#' data(dataNCC)
#' modJoint.ncc <- frailtyPenal(Surv(t.start,t.stop,event)~cluster(id)+cov1
#' +cov2+terminal(death)+wts(ncc.wts), formula.terminalEvent=~cov1+cov2,
#' data=dataNCC,n.knots=8,kappa=c(1.6e+10, 5.0e+03),recurrentAG=TRUE, RandDist="LogN") 
#' 
#' 
#' 
#' print(modJoint.ncc)
#' 
#' }
#' 
#' 
"wts" <- function(x)
{
  x
}
