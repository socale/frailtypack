#' Plot difference of EPOCE estimators between two joint frailty models.
#' 
#' Plots values of the difference of two Cross-Validated Prognosis Observed
#' Loss (CVPOL) computed with two joint frailty models. Confidence intervals
#' are allowed.
#' 
#' 
#' @usage \method{plot}{Diffepoce}(x, conf.bands=TRUE, Xlab = "Time", Ylab =
#' "EPOCE difference" , ...)
#' @param x An object inheriting from \code{Diffepoce} class.
#' @param conf.bands Logical value. Determines whether confidence intervals
#' will be plotted. The default is FALSE.
#' @param Xlab Label of x-axis. Default is '"Time"'
#' @param Ylab Label of y-axis. Default is '"EPOCE difference"'
#' @param \dots Other unused arguments.
#' @return Print one plot with one curve and its confidence interval.
#' @seealso \code{\link{Diffepoce}}
#' @keywords file
##' @export
"plot.Diffepoce" <- function (x, conf.bands=TRUE, Xlab = "Time", Ylab = "EPOCE difference",  ...)
{

	if (conf.bands){
		matplot(x$pred.times,cbind(x$DEPOCE,x$TIinf,x$TIsup),type=c("b","l","l"),pch="X",lty=c(1,2,2),col="blue",xlab=Xlab,ylab=Ylab)
	}else{
		plot(x$pred.times,x$DEPOCE,type="b",pch="X",col="blue",xlab=Xlab,ylab=Ylab)
	}

	abline(h=0,lty=1,lwd=2)

	return(invisible())
}
