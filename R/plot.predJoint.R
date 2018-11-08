#' Plot predictions using a joint frailty model.
#' 
#' Plots predicted probabilities of terminal event. Confidence intervals are
#' allowed.
#' 
#' 
#' @usage \method{plot}{predJoint}(x, conf.bands=FALSE,
#' relapses=TRUE,pos.legend="topright", cex.legend=0.7, ylim=c(0,1), Xlab =
#' "Time t", Ylab = "Prediction probability of event", ...)
#' @param x An object from the 'prediction' function, more generaly a
#' \code{predJoint} class object.
#' @param conf.bands Logical value. Determines whether confidence intervals
#' will be plotted. The default is FALSE.
#' @param relapses Logical value. Determines whether observed recurrent events
#' will be plotted. The default is TRUE.
#' @param pos.legend The location of the legend can be specified by setting
#' this argument to a single keyword from the list '"bottomright"', '"bottom"',
#' '"bottomleft"', '"left"', '"topleft"', '"top"', '"topright"', '"right"' and
#' '"center"'. The default is '"topright"'
#' @param cex.legend size of the legend. Default is 0.7
#' @param ylim range of y-axis. Default is from 0 to 1
#' @param Xlab Label of x-axis. Default is '"Time t"'
#' @param Ylab Label of y-axis. Default is '"Prediction probability of event"'
#' @param \dots Other unused arguments
#' @return Print as many plots as the number of subjects.
#' @keywords file
##' @export
"plot.predJoint" <- function (x, conf.bands=FALSE, relapses=TRUE, pos.legend="topright", cex.legend=0.7, ylim=c(0,1), Xlab = "Time t", Ylab = "Prediction probability of event", ...){	
	event.type <- x$event
	
	if ((conf.bands) & (!x$icproba)) stop("Confidence intervals were not calculated. Use the MC.sample argument in the 'prediction' function")
	
	if (event.type == 1){
		if (x$npred >5) {
			warning("For a better overview only predictions for a maximum of 5 subjects are plotted (the 5 first ones)")
			x$npred <- 5
		}
		par(mfcol=c(2,x$npred))
	}else{
		if (x$npred >10) {
			warning("For a better overview only predictions for a maximum of 10 subjects are plotted (the 10 first ones)")
			x$npred <- 10
		}
		if (x$npred <= 5) par(mfrow=c(1,x$npred))
		else par(mfrow=c(2,ceiling(x$npred/2)))
	}
	
	if (x$moving.window){
		legende1 <- paste("Predicted cumulative probability of terminal event between",x$t,"and time t")
		legende2 <- paste("Predicted cumulative probability of recurrent event between",x$t,"and time t")
		if ((missing(Ylab)) && (event.type==2)) legende <- legende1
		if ((missing(Ylab)) && (event.type==3)) legende <- legende2
		
		if (relapses){
			predtimerec <- x$predtimerec[x$predtimerec<=x$t & x$predtimerec!=0]
			xlim <- c(ifelse(length(predtimerec)!=0,min(predtimerec),x$t-(max(x$x.time)-min(x$x.time))*0.1),max(x$x.time))
		}else{
			xlim <- c(x$t-(max(x$x.time)-min(x$x.time))*0.1,max(x$x.time))
		}
	}else{ 
		legende1 <- paste("Predicted probability of terminal event in the next",x$window)
		legende2 <- paste("Predicted probability of recurrent event in the next",x$window)
		if ((missing(Ylab)) && (event.type==2)) legende <- legende1
		if ((missing(Ylab)) && (event.type==3)) legende <- legende2
		
		xlim <- c(min(x$x.time),max(x$x.time))
	}

	if (!x$intcens){ #(x$joint.clust==1){
		for (i in 1:(x$npred)){
			if (conf.bands) {
				switch(as.character(event.type),				
					"2" = {matplot(x$x.time,cbind(x$pred1[i,],x$predlow1[i,],x$predhigh1[i,]),col="blue",type="l",lty=c(1,2,2),xlab=Xlab,ylab=legende,main=paste("id patient :",x$group[i]),ylim=ylim,xlim=xlim)
					matlines(x$x.time,cbind(x$pred2[i,],x$predlow2[i,],x$predhigh2[i,]),col="red",type="l",lty=c(1,2,2))
					matlines(x$x.time,cbind(x$pred3[i,],x$predlow3[i,],x$predhigh3[i,]),col="green",type="l",lty=c(1,2,2))},
					
					"3" = {matplot(x$x.time,cbind(x$pred1_rec[i,],x$predlow1_rec[i,],x$predhigh1_rec[i,]),col="blue",type="l",lty=c(1,2,2),xlab=Xlab,ylab=legende,main=paste("id patient :",x$group[i]),ylim=ylim,xlim=xlim)},
					
					"1" = {matplot(x$x.time,cbind(x$pred1_rec[i,],x$predlow1_rec[i,],x$predhigh1_rec[i,]),col="blue",type="l",lty=c(1,2,2),xlab=Xlab,ylab=legende2,main=paste("id patient :",x$group[i]),ylim=ylim,xlim=xlim)
					matplot(x$x.time,cbind(x$pred1[i,],x$predlow1[i,],x$predhigh1[i,]),col="blue",type="l",lty=c(1,2,2),xlab=Xlab,ylab=legende1,main=paste("id patient :",x$group[i]),ylim=ylim,xlim=xlim)
					matlines(x$x.time,cbind(x$pred2[i,],x$predlow2[i,],x$predhigh2[i,]),col="red",type="l",lty=c(1,2,2))
					matlines(x$x.time,cbind(x$pred3[i,],x$predlow3[i,],x$predhigh3[i,]),col="green",type="l",lty=c(1,2,2))}
                )			
				#if(event.type==2){ # Terminal
					
				#}else{ # Recurrent					
			}else{
				switch(as.character(event.type),
				#if(event.type==3){ # Terminal
					"2" = {plot(x$x.time,x$pred1[i,],col="blue",type="l",lty=c(1,2,2),xlab=Xlab,ylab=legende,main=paste("id patient :",x$group[i]),ylim=ylim,xlim=xlim)
					lines(x$x.time,x$pred2[i,],col="red",type="l",lty=c(1,2,2))
					lines(x$x.time,x$pred3[i,],col="green",type="l",lty=c(1,2,2))},
					
					"3" = {plot(x$x.time,x$pred1_rec[i,],col="blue",type="l",lty=c(1,2,2),xlab=Xlab,ylab=legende,main=paste("id patient :",x$group[i]),ylim=ylim,xlim=xlim)},
					
					"1" = {plot(x$x.time,x$pred1_rec[i,],col="blue",type="l",lty=c(1,2,2),xlab=Xlab,ylab=legende2,main=paste("id patient :",x$group[i]),ylim=ylim,xlim=xlim)
					plot(x$x.time,x$pred1[i,],col="blue",type="l",lty=c(1,2,2),xlab=Xlab,ylab=legende1,main=paste("id patient :",x$group[i]),ylim=ylim,xlim=xlim)
					lines(x$x.time,x$pred2[i,],col="red",type="l",lty=c(1,2,2))
					lines(x$x.time,x$pred3[i,],col="green",type="l",lty=c(1,2,2))}
				)
					
				#}else{
					#plot(x$x.time,x$pred1_rec[i,],col="blue",type="l",lty=c(1,2,2),xlab=Xlab,ylab=legende,main=paste("id patient :",x$group[i]),ylim=ylim,xlim=xlim)
				#}
			}
			if (x$moving.window) abline(v=x$t,lty=2)
			
			if (relapses) {
				if (x$moving.window) predtimereci <- x$predtimerec[i,][x$predtimerec[i,]<=x$t & x$predtimerec[i,]!=0]
				else predtimereci <- x$predtimerec[i,][x$predtimerec[i,]!=0]
				lines(predtimereci,rep(0,length(predtimereci)),type="p",pch="X")
			}
			
			if (i==1) {
				if (relapses) {
					if(event.type==2){ #Terminal
						legend(pos.legend, c("p1: exactly j recurrences","p2: at least j recurrences","p3: ignoring recurrences","recurrent event"), lty=c(1,1,1,0), pch=c("","","","X"), col=c("blue","red","green","black"), cex=cex.legend)
					}else{
						legend(pos.legend, c("p1: exactly j recurrences", "recurrent event"), lty=c(1,0), pch=c("","X"), col=c("blue", "black"), cex=cex.legend)
					}
				}else{
					if(event.type==3){ #Recurrent
						legend(pos.legend, c("p1: exactly j recurrences","p2: at least j recurrences","p3: ignoring recurrences"), lty=1, col=c("blue","red","green"), cex=cex.legend)
					}else{						
						legend(pos.legend, "p1: exactly j recurrences", lty=1, col="blue", cex=cex.legend)
					}
				}
			}
		}
	}else{ #Pas d'intervalle par censure dans un modele conjoint
		for (i in 1:(x$npred)) {
			if (conf.bands) {
				matplot(x$x.time,cbind(x$pred2[i,],x$predlow2[i,],x$predhigh2[i,]),col="red",type="l",lty=c(1,2,2),xlab=Xlab,ylab=Ylab,main=paste("id patient :",x$group[i]),ylim=ylim) 
			}else{
				plot(x$x.time,x$pred2[i,],col="red",type="l",lty=c(1,2,2),xlab=Xlab,ylab=Ylab,main=paste("id patient :",x$group[i]),ylim=ylim)
			}
			if (i==1) legend(pos.legend, c("probability of death"), lty=1, col=c("red"), cex=cex.legend)
		}
	}
	return(invisible())
}
