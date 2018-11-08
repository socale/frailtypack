#' Plot predictions using a Cox or a shared frailty model.
#' 
#' Plots predicted probabilities of event. Confidence intervals are allowed.
#' 
#' 
#' @usage \method{plot}{predFrailty}(x, conf.bands=FALSE, pos.legend="topright",
#' cex.legend=0.7, ylim=c(0,1), Xlab = "Time t", Ylab, ...)
#' @param x An object from the 'prediction' function, i.e. a \code{predFrailty}
#' class object.
#' @param conf.bands Logical value. Determines whether confidence intervals
#' will be plotted. The default is FALSE.
#' @param pos.legend The location of the legend can be specified by setting
#' this argument to a single keyword from the list '"bottomright"', '"bottom"',
#' '"bottomleft"', '"left"', '"topleft"', '"top"', '"topright"', '"right"' and
#' '"center"'. The default is '"topright"'.
#' @param cex.legend size of the legend. Default is 0.7.
#' @param ylim range of y-axis. Default is from 0 to 1.
#' @param Xlab Label of x-axis. Default is '"Time t"'
#' @param Ylab Label of y-axis.
#' @param \dots Other unused arguments.
#' @return Print one plot with as many curves as the number of profiles.
#' @keywords file
##' @export
"plot.predFrailty" <- function (x, conf.bands=FALSE, pos.legend="topright", cex.legend=0.7, ylim=c(0,1), Xlab = "Time t", Ylab, ...)
{
	if ((conf.bands) & (!x$icproba)) stop("Confidence intervals were not calculated. Use the MC.sample argument in the 'prediction' function")
	
	if (x$moving.window){ 
		if(missing(Ylab)) legende <- paste("Predicted cumulative probability of event between",x$t,"and time t")
		else legende <- Ylab
		
		if (x$event=="Recurrent"){
			predtimerec <- x$predtimerec[x$predtimerec<=x$t & x$predtimerec!=0]
			xlim <- c(ifelse(length(predtimerec)!=0,min(predtimerec),x$t-(max(x$x.time)-min(x$x.time))*0.1),max(x$x.time))
		}else{
			xlim <- c(x$t-(max(x$x.time)-min(x$x.time))*0.1,max(x$x.time))
		}		
	}else{ 
		if(missing(Ylab)) legende <- paste("Predicted probability of event in the next",x$window)
		else legende <- Ylab
		xlim <- c(min(x$x.time),max(x$x.time))
	}
	
	# par(mfrow=c(2,1))
	if (is.null(x$type)){
		title <- paste("Dynamic prediction for a Cox model")
	}else{
		if (x$type=="marginal") {
			if (x$event=="Recurrent") title <- paste("Marginal dynamic prediction of a recurrent event\nfor a shared frailty model")
			else title <- paste("Marginal dynamic prediction for a shared frailty model")
		}
		else title <- paste("Conditional dynamic prediction for a shared frailty model")
	}
	
	# par(mfrow=c(2,ceiling(x$npred/2)))
	if (conf.bands){
		# for (i in 1:(x$npred)){
			# matplot(x$x.time,t(x$pred)[,i],type="l",lty=1,xlab=Xlab,ylab=legende,main=title,ylim=ylim,xlim=xlim)
			# matlines(x$x.time,t(x$predLow)[,i],type="l",lty=2)
			# matlines(x$x.time,t(x$predHigh)[,i],type="l",lty=2)	
			# legend(pos.legend, rownames(x$pred)[i],lty=1,bty="n",col=i)
		# }
		matplot(x$x.time,t(x$pred),type="l",lty=1,xlab=Xlab,ylab=legende,main=title,ylim=ylim,xlim=xlim)
		matlines(x$x.time,t(x$predLow),type="l",lty=2)
		matlines(x$x.time,t(x$predHigh),type="l",lty=2)		
	}else{
		matplot(x$x.time,t(x$pred),type="l",lty=1,xlab=Xlab,ylab=legende,main=title,ylim=ylim,xlim=xlim)
		# for (i in 1:(x$npred)){
			# matplot(x$x.time,t(x$pred)[,1],type="l",lty=1,xlab=Xlab,ylab=legende,main=title,ylim=ylim,xlim=xlim)			
		# }
	}
	
	if (x$event=="Recurrent"){
	    for (i in 1:(x$npred)){
			if (x$moving.window) predtimereci <- x$predtimerec[i,][x$predtimerec[i,]<=x$t & x$predtimerec[i,]!=0]
			else predtimereci <- x$predtimerec[i,][x$predtimerec[i,]!=0]
			lines(predtimereci,rep(0,length(predtimereci)),type="p",pch="X",col=i)
		}
	}
	legend(pos.legend, rownames(x$pred),lty=1,bty="n",col=(1:x$npred))

	if (x$moving.window) abline(v=x$t,lty=2)

	return(invisible())
}
