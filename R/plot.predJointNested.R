##' @export
"plot.predJointNested" <- function (x, conf.bands=FALSE, relapses=TRUE, pos.legend="topright", cex.legend=0.7, ylim=c(0,1), Xlab = "Time t", Ylab = "Prediction probability of event", ...){	
	
	if ((conf.bands) & (!x$icproba)) stop("Confidence intervals were not calculated. Use the MC.sample argument in the 'prediction' function")
	
	if (x$npred >10) {
			warning("For a better overview only predictions for a maximum of 10 subjects are plotted (the 10 first ones)")
			x$npred <- 10
		}
	if (x$npred <= 5) par(mfrow=c(1,x$npred))
	else par(mfrow=c(2,ceiling(x$npred/2)))
	
	if (x$moving.window){
		if ( missing(Ylab)) legende <- paste("Predicted cumulative probability of terminal event between",x$t,"and time t")
		
		if (relapses){
			predtimerec <- x$predtimerec[x$predtimerec<=x$t & x$predtimerec!=0]
			xlim <- c(ifelse(length(predtimerec)!=0,min(predtimerec),x$t-(max(x$x.time)-min(x$x.time))*0.1),max(x$x.time))
		}else{
			xlim <- c(x$t-(max(x$x.time)-min(x$x.time))*0.1,max(x$x.time))
		}
	}else{ 
		if ( missing(Ylab)) legende <- paste("Predicted probability of terminal event in the next",x$window)
		xlim <- c(min(x$x.time),max(x$x.time))
	}
	for (i in 1:(x$npred)){
		if (conf.bands) {			
			matplot(x$x.time,cbind(x$pred[i,],x$predLow[i,],x$predHigh[i,]),col="blue",type="l",lty=c(1,2,2),xlab=Xlab,ylab=legende,main=paste("id patient :",x$group[i]),ylim=ylim,xlim=xlim)
		}else{
			plot(x$x.time,x$pred[i,],col="blue",type="l",lty=c(1,2,2),xlab=Xlab,ylab=legende,main=paste("id patient :",x$group[i]),ylim=ylim,xlim=xlim)
		}
		
		if (x$moving.window) abline(v=x$t,lty=2)
		
		if (relapses) {
			if (x$moving.window) predtimereci <- x$predtimerec[i,][x$predtimerec[i,]<=x$t & x$predtimerec[i,]!=0]
			else predtimereci <- x$predtimerec[i,][x$predtimerec[i,]!=0]
			lines(predtimereci,rep(0,length(predtimereci)),type="p",pch="X")
		}
		
		if (i==1) {
			if (relapses){					
					legend(pos.legend, c("Prediction of terminal event","recurrent event"), lty=c(1,1,1,0), pch=c("","X"), col=c("blue","red"), cex=cex.legend)
			}else{										
				legend(pos.legend, "Prediction of terminal event", lty=1, col="blue", cex=cex.legend)		
			}
		}
	}
	return(invisible())
}
