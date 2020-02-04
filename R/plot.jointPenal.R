#' Plot Method for a Joint frailty model.
#' 
#' Plots estimated baseline survival and hazard functions of a joint frailty
#' model (output from an object of class 'JointPenal' for joint frailty models
#' ) for each type of event (terminal or recurrent). Confidence bands are
#' allowed.
#' 
#' 
#' @aliases plot.jointPenal lines.jointPenal
#' @aliases plot.jointPenal lines.jointPenal
#' @usage
#' 
#' \method{plot}{jointPenal}(x, event = "Both", type.plot = "Hazard", conf.bands
#' = FALSE, pos.legend="topright", cex.legend = 0.7, ylim, main, color = 2, median=TRUE,
#' Xlab = "Time", Ylab = "Hazard function", ...)
#' @param x A joint model, i.e. an object of class \code{frailtyPenal} for
#' Joint frailty model (output from calling \code{frailtyPenal} function).
#' @param event a character string specifying the type of curve. Possible value
#' are "Terminal", "Recurrent", or "Both". The default is "Both".
#' @param type.plot a character string specifying the type of curve. Possible
#' value are "Hazard", or "Survival". The default is "Hazard". Only the first
#' letters are required, e.g "Haz", "Su"
#' @param conf.bands logical value. Determines whether confidence bands will be
#' plotted. The default is to do so.
#' @param pos.legend The location of the legend can be specified by setting
#' this argument to a single keyword from the list '"bottomright"', '"bottom"',
#' '"bottomleft"', '"left"', '"topleft"', '"top"', '"topright"', '"right"' and
#' '"center"'. The default is '"topright"'
#' @param cex.legend character expansion factor *relative* to current
#' 'par("cex")'. Default is 0.7
#' @param ylim y-axis limits
#' @param main plot title
#' @param color curve color (integer)
#' @param median Logical value. Determines whether survival median will be plotted. Default is TRUE.
#' @param Xlab Label of x-axis. Default is '"Time"'
#' @param Ylab Label of y-axis. Default is '"Hazard function"'
#' @param ... other unused arguments
#' @return Print a plot of the baseline survival or hazard functions for each
#' type of event or both with the confidence bands or not (conf.bands argument)
#' @seealso \code{\link{frailtyPenal}}
#' @keywords methods
##' @export
#' @examples
#' 
#' 
#' \dontrun{
#' 
#' data(readmission)
#' 
#' #-- Gap-time
#' modJoint.gap <- frailtyPenal(Surv(time,event)~cluster(id)+sex+dukes+
#' charlson+terminal(death),formula.terminalEvent=~sex+dukes+charlson,
#' data=readmission,n.knots=14,kappa=c(100,100))
#' 
#' #-- It takes around 1 minute to converge --#
#' 
#' plot(modJoint.gap,type.plot="Haz",event="recurrent",conf.bands=TRUE)
#' plot(modJoint.gap,type.plot="Haz",event="terminal",conf.bands=TRUE)
#' plot(modJoint.gap,type.plot="Haz",event="both",conf.bands=TRUE)
#' 
#' plot(modJoint.gap,type.plot="Su",event="recurrent",conf.bands=TRUE)
#' plot(modJoint.gap,type.plot="Su",event="terminal",conf.bands=TRUE)
#' plot(modJoint.gap,type.plot="Su",event="both",conf.bands=TRUE)
#' 
#' 
#' }
#' 
#' 
"plot.jointPenal" <-
  function (x, event="Both", type.plot="Hazard", conf.bands=FALSE, pos.legend="topright", cex.legend=0.7, ylim, main, color=2, median=TRUE, Xlab = "Time", Ylab = "Hazard function", ...)
  {
    event.type <- charmatch(event, c("Both", "Recurrent", "Terminal"), nomatch = 0)
    if (event.type == 0) {
      stop("event must be 'Both', 'Recurrent' or 'Terminal'")
    }
    
    plot.type <- charmatch(type.plot, c("Hazard", "Survival"), nomatch = 0)
    if (plot.type == 0) {
      stop("estimator must be 'Hazard' or 'Survival'")
    }
    
    if(missing(main)) main<-""
    
    if (event.type==1){ # both
      
      if(plot.type==1){
        if (missing(ylim)){
          yymax<-max(c(x$lamR, x$lamD),na.rm=TRUE)
          yymin<-min(c(x$lamR, x$lamD),na.rm=TRUE)
        }else{
          yymax<-ylim[2]
          yymin<-ylim[1]
        }
        
        if (conf.bands){
          matplot(x$xR[,1], x$lamR[,,1], col=color, type="l", lty=c(1,2,2), xlab=Xlab,ylab=Ylab, ylim=c(yymin,yymax), main=main, ...)
          for (i in (1:x$n.strat)[-1]) matlines(x$xR[,i], x$lamR[,,i], col=color+(i-1), type="l", lty=c(1,2,2), ...)
          matlines(x$xD, x$lamD, col=color+x$n.strat, type="l", lty=c(1,2,2), ...)
        }else{
          plot(x$xR[,1], x$lamR[,1,1], col=color, type="l", lty=1, xlab=Xlab,ylab=Ylab, ylim=c(yymin,yymax), main=main, ...)
          for (i in (1:x$n.strat)[-1]) lines(x$xR[,i], x$lamR[,1,i], col=color+(i-1), type="l", lty=1, ...)
          lines(x$xD, x$lamD[,1], col=color+x$n.strat, type="l", lty=1, ...)
        }
      }else{ 
        if (missing(Ylab)) Ylab <- "Baseline survival function"
        
        if (missing(ylim)){
          yymax<-1
          yymin<-0
        }else{
          yymax<-ylim[2]
          yymin<-ylim[1]
        }
        if (x$typeof == 0){
          if (conf.bands){
            matplot(x$xR[,1], x$survR[,,1], col=color, type="l", lty=c(1,2,2), xlab=Xlab,ylab=Ylab, ylim=c(yymin,yymax), main=main, ...)
            for (i in (1:x$n.strat)[-1]) matlines(x$xR[,i], x$survR[,,i], col=color+(i-1), type="l", lty=c(1,2,2), ...)
            matlines(x$xD, x$survD, col=color+x$n.strat, type="l", lty=c(1,2,2), ...)
            if (median){abline(a=0.5,b=0,cex=.5,col=1,lty=3)}
          }else{
            plot(x$xR[,1], x$survR[,1,1], col=color, type="l", lty=1, xlab=Xlab,ylab=Ylab, ylim=c(yymin,yymax), main=main, ...)
            for (i in (1:x$n.strat)[-1]) lines(x$xR[,i], x$survR[,1,i], col=color+(i-1), type="l", lty=1, ...)
            lines(x$xD, x$survD[,1], col=color+x$n.strat, type="l", lty=1, ...)
            if (median){abline(a=0.5,b=0,cex=.5,col=1,lty=3)}
          }
        }else{
          if (conf.bands){
            matplot(x$xSuR[,1], x$survR[,,1], col=color, type="l", lty=c(1,2,2), xlab=Xlab,ylab=Ylab, ylim=c(yymin,yymax), main=main, ...)
            for (i in (1:x$n.strat)[-1]) matlines(x$xSuR[,i], x$survR[,,i], col=color+(i-1), type="l", lty=c(1,2,2), ...)
            matlines(x$xSuD, x$survD, col=color+x$n.strat, type="l", lty=c(1,2,2), ...)
            if (median){abline(a=0.5,b=0,cex=.5,col=1,lty=3)}
          }else{
            plot(x$xSuR[,1], x$survR[,1,1], col=color, type="l", lty=1, xlab=Xlab,ylab=Ylab, ylim=c(yymin,yymax), main=main, ...)
            for (i in (1:x$n.strat)[-1]) lines(x$xSuR[,i], x$survR[,1,i], col=color+(i-1), type="l", lty=1, ...)
            lines(x$xSuD, x$survD[,1], col=color+x$n.strat, type="l", lty=1, ...)
            if (median){abline(a=0.5,b=0,cex=.5,col=1,lty=3)}
          }
        }
      }
      if (x$n.strat > 1) legend(pos.legend, c(paste("recurrent event strata =",1:x$n.strat),"terminal event"), lty=1, col=color+(0:x$n.strat), xjust=1, cex=cex.legend, ...)
      else legend(pos.legend, c("recurrent event","terminal event"), lty=1, col=c(color,color+x$n.strat), xjust=1, cex=cex.legend, ...)
    }
    
    if (event.type==2){ # recurrent
      
      if(plot.type==1){
        if (missing(ylim)){
          yymax<-max(x$lamR,na.rm=TRUE)
          yymin<-min(x$lamR,na.rm=TRUE)
        }else{
          yymax<-ylim[2]
          yymin<-ylim[1]
        }
        
        if (conf.bands){
          matplot(x$xR[,1], x$lamR[,,1], col=color, type="l", lty=c(1,2,2), xlab=Xlab,ylab=Ylab, ylim=c(yymin,yymax), main=main, ...)
          for (i in (1:x$n.strat)[-1]) matlines(x$xR[,i], x$lamR[,,i], col=color+(i-1), type="l", lty=c(1,2,2), ...)
        }else{
          plot(x$xR[,1], x$lamR[,1,1], col=color, type="l", lty=1, xlab=Xlab,ylab=Ylab, ylim=c(yymin,yymax), main=main, ...)
          for (i in (1:x$n.strat)[-1]) lines(x$xR[,i], x$lamR[,1,i], col=color+(i-1), type="l", lty=1, ...)
        }
      }else{
        if (missing(Ylab)) Ylab <- "Baseline survival function"		
        if (missing(ylim)){
          yymax<-1
          yymin<-0
        }else{
          yymax<-ylim[2]
          yymin<-ylim[1]
        }
        if (x$typeof == 0){
          if (conf.bands){
            matplot(x$xR[,1], x$survR[,,1], col=color, type="l", lty=c(1,2,2), xlab=Xlab,ylab=Ylab, ylim=c(yymin,yymax), main=main, ...)
            for (i in (1:x$n.strat)[-1]) matlines(x$xR[,i], x$survR[,,i], col=color+(i-1), type="l", lty=c(1,2,2), ...)
            if (median){abline(a=0.5,b=0,cex=.5,col=1,lty=3)}
          }else{
            plot(x$xR[,1], x$survR[,1,1], col=color, type="l", lty=1, xlab=Xlab,ylab=Ylab, ylim=c(yymin,yymax), main=main, ...)
            for (i in (1:x$n.strat)[-1]) lines(x$xR[,i], x$survR[,1,i], col=color+(i-1), type="l", lty=1, ...)
            if (median){abline(a=0.5,b=0,cex=.5,col=1,lty=3)}
          }
        }else{
          if (conf.bands){
            matplot(x$xSuR[,1], x$survR[,,1], col=color, type="l", lty=c(1,2,2), xlab=Xlab,ylab=Ylab, ylim=c(yymin,yymax), main=main, ...)
            for (i in (1:x$n.strat)[-1]) matlines(x$xSuR[,i], x$survR[,,i], col=color+(i-1), type="l", lty=c(1,2,2), ...)
            if (median){abline(a=0.5,b=0,cex=.5,col=1,lty=3)}
          }else{
            plot(x$xSuR[,1], x$survR[,1,1], col=color, type="l", lty=1, xlab=Xlab,ylab=Ylab, ylim=c(yymin,yymax), main=main, ...)
            for (i in (1:x$n.strat)[-1]) lines(x$xSuR[,i], x$survR[,1,i], col=color+(i-1), type="l", lty=1, ...)
            if (median){abline(a=0.5,b=0,cex=.5,col=1,lty=3)}
          }
        }
      }
      if (x$n.strat > 1) legend(pos.legend, paste("recurrent event strata =",1:x$n.strat), lty=1, col=color+(1:x$n.strat-1), xjust=1, cex=cex.legend, ...)
      else legend(pos.legend, c("recurrent event"), lty=1, col=color, xjust=1, cex=cex.legend, ...)
    }
    
    if (event.type==3){ # terminal
      
      if(plot.type==1){
        if (missing(ylim)){
          yymax<-max(x$lamD,na.rm=TRUE)
          yymin<-min(x$lamD,na.rm=TRUE)
        }else{
          yymax<-ylim[2]
          yymin<-ylim[1]
        }
        
        if (conf.bands){
          matplot(x$xD, x$lamD, col=color+x$n.strat, type="l", lty=c(1,2,2), xlab=Xlab,ylab=Ylab, ylim=c(yymin,yymax), main=main, ...)
        }else{
          plot(x$xD, x$lamD[,1], col=color+x$n.strat, type="l", lty=1, xlab=Xlab,ylab=Ylab, ylim=c(yymin,yymax), main=main, ...)
        }
      }else{
        
        if (missing(ylim)){
          yymax<-1
          yymin<-0
        }else{
          yymax<-ylim[2]
          yymin<-ylim[1]
        }
        if (x$typeof == 0){
          if (conf.bands){
            matplot(x$xD, x$survD, col=color+x$n.strat, type="l", lty=c(1,2,2), xlab=Xlab,ylab=Ylab, ylim=c(yymin,yymax), main=main, ...)
            if (median){abline(a=0.5,b=0,cex=.5,col=1,lty=3)}
          }else{
            plot(x$xD, x$survD[,1], col=color+x$n.strat, type="l", lty=1, xlab=Xlab,ylab=Ylab, ylim=c(yymin,yymax), main=main, ...)
            if (median){abline(a=0.5,b=0,cex=.5,col=1,lty=3)}
          }
        }else{
          if (conf.bands){
            matplot(x$xSuD, x$survD, col=color+x$n.strat, type="l", lty=c(1,2,2), xlab=Xlab,ylab=Ylab, ylim=c(yymin,yymax), main=main, ...)
            if (median){abline(a=0.5,b=0,cex=.5,col=1,lty=3)}
          }else{
            plot(x$xSuD, x$survD[,1], col=color+x$n.strat, type="l", lty=1, xlab=Xlab,ylab=Ylab, ylim=c(yymin,yymax), main=main, ...)
            if (median){abline(a=0.5,b=0,cex=.5,col=1,lty=3)}
          }
        }
      }
      legend(pos.legend, c("terminal event"), lty=1, col=color+x$n.strat, xjust=1, cex=cex.legend, ...)
    }
    
    return(invisible())
  }
