
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
