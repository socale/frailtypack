
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
