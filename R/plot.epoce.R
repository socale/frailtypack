#' Plot values of estimators of the Expected Prognostic Observed Cross-Entropy
#' (EPOCE).
#' 
#' Plots values of estimators MPOL and CVPOL for evaluating EPOCE. No
#' confidence interval.
#' 
#' 
#' @usage \method{plot}{epoce}(x, type, pos.legend="topright", cex.legend=0.7,
#' Xlab="Time",Ylab="Epoce", ...)
#' @param x An object inheriting from \code{epoce} class
#' @param type Type of estimator to plot. If new dataset was used only mpol can
#' be plotted (\code{"mpol"}), otherwise mpol and cvpol can be plotted
#' (\code{"mpol"} and \code{"cvpol"}, default is \code{"cvpol"}).
#' @param pos.legend The location of the legend can be specified by setting
#' this argument to a single keyword from the list '"bottomright"', '"bottom"',
#' '"bottomleft"', '"left"', '"topleft"', '"top"', '"topright"', '"right"' and
#' '"center"'. The default is '"topright"'.
#' @param cex.legend size of the legend. Default is 0.7.
#' @param Xlab Label of x-axis. Default is '"Time"'
#' @param Ylab Label of y-axis. Default is '"Epoce"'
#' @param \dots Other unused arguments.
#' @return Print a curve of the estimator of EPOCE using time points defined in
#' \code{epoce}.
#' @seealso \code{\link{epoce}}
#' @export
#' @keywords file
"plot.epoce" <- function (x, type, pos.legend="topright", cex.legend=0.7, Xlab="Time",Ylab="Epoce", ...){
	if (missing(type)) stop("'type' argument is missing. Please refer to the documentation.")
	if (x$new.data){
		if(!missing(type) && type!="mpol") stop("For epoce with new dataset only mpol is calculated and can be plotted")
		plot(x$pred.times,x$mpol,type="b",pch="X",col="blue",xlab=Xlab,ylab=Ylab)
		legend(pos.legend,c("mpol"),lty=1,col="blue")
	}else{
		if(!missing(type) && !type%in%c("cvpol","mpol"))stop("Only mpol or cvpol can be plotted for the epoce estimator")
		if (type == "mpol") plot(x$pred.times,x$mpol,type="b",pch="X",col="red",xlab=Xlab,ylab=Ylab)
		else plot(x$pred.times,x$cvpol,type="b",pch="X",col="red",xlab=Xlab,ylab=Ylab)
		# plot(x$pred.times,noquote(paste("x$",type,sep="")),type="b",pch="X",col="red",xlab=Xlab,ylab=Ylab)
		
		legend(pos.legend,type,lty=1,col="red")
	}
	return(invisible())
}
