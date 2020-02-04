#' Plot Method for a trivariate joint model for longitudinal data, recurrent
#' events and a terminal event.
#' 
#' Plots estimated baseline survival and hazard functions of a joint model
#' (output from an object of class 'trivPenal') for each type of event
#' (terminal or recurrent). Confidence bands are allowed.
#' 
#' 
#' @aliases plot.trivPenal lines.trivPenal
#' @usage
#' 
#' \method{plot}{trivPenal}(x, event = "Both", type.plot = "Hazard", conf.bands =
#' FALSE, pos.legend="topright", cex.legend = 0.7, ylim, main, color = 2, median=TRUE, Xlab
#' = "Time", Ylab = "Hazard function", ...)
#' @param x A joint model, an object of class \code{trivPenal}.
#' @param event a character string specifying the type of curve. Possible value
#' are "Terminal", "Recurrent", or "Both". The default is "Both".
#' @param type.plot a character string specifying the type of curve. Possible
#' value are "Hazard", or "Survival". The default is "Hazard". Only the first
#' words are required, e.g "Haz", "Su"
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
#' @seealso \code{\link{trivPenal}}
#' @keywords methods
##' @export
#' @examples
#' 
#' \dontrun{
#' ###--- Trivariate joint model for longitudinal data, ---###
#' ###--- recurrent events and a terminal event ---###
#' 
#' data(colorectal)
#' data(colorectalLongi)
#' 
#' # Weibull baseline hazard function
#' # Random effects as the link function, Gap timescale
#' # (computation takes around 30 minutes)
#' model.weib.RE.gap <-trivPenal(Surv(gap.time, new.lesions) ~ cluster(id)
#' + age + treatment + who.PS + prev.resection + terminal(state),
#' formula.terminalEvent =~ age + treatment + who.PS + prev.resection, 
#' tumor.size ~ year * treatment + age + who.PS, data = colorectal,
#' data.Longi = colorectalLongi, random = c("1", "year"), id = "id", 
#' link = "Random-effects", left.censoring = -3.33, recurrentAG = FALSE,
#' hazard = "Weibull", method.GH="Pseudo-adaptive", n.nodes = 7)
#' 
#' plot(model.weib.RE.gap)
#' plot(model.weib.RE.gap, type = "survival")
#' }
"plot.trivPenal" <-
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
    
    
    if(missing(main))
      main<-""
    
    if (event.type==1){ # both
      
      if(plot.type==1){
        if (missing(ylim)){
          yymax<-max(c(x$lamR[-1], x$lamD[-1]),na.rm=TRUE)
          yymin<-min(c(x$lamR[-1], x$lamD[-1]),na.rm=TRUE)
        }else{
          yymax<-ylim[2]
          yymin<-ylim[1]
        }
        
        if (conf.bands){
          matplot(x$xR[,1], x$lamR[,,1], col=color, type="l", lty=c(1,2,2), xlab=Xlab,ylab=Ylab, ylim=c(yymin,yymax), main=main, ...)
          matlines(x$xD, x$lamD, col=color+x$n.strat, type="l", lty=c(1,2,2), ...)
        }else{
          plot(x$xR[,1], x$lamR[,1,1], col=color, type="l", lty=1, xlab=Xlab,ylab=Ylab, ylim=c(yymin,yymax), main=main, ...)
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
            matlines(x$xD, x$survD, col=color+x$n.strat, type="l", lty=c(1,2,2), ...)
            if (median){abline(a=0.5,b=0,cex=.5,col=1,lty=3)}
          }else{
            plot(x$xR[,1], x$survR[,1,1], col=color, type="l", lty=1, xlab=Xlab,ylab=Ylab, ylim=c(yymin,yymax), main=main, ...)
            lines(x$xD, x$survD[,1], col=color+x$n.strat, type="l", lty=1, ...)
            if (median){abline(a=0.5,b=0,cex=.5,col=1,lty=3)}
          }
        }else{
          if (conf.bands){
            matplot(x$xSuR[,1], x$survR[,,1], col=color, type="l", lty=c(1,2,2), xlab=Xlab,ylab=Ylab, ylim=c(yymin,yymax), main=main, ...)
            matlines(x$xSuD, x$survD, col=color+x$n.strat, type="l", lty=c(1,2,2), ...)
            if (median){abline(a=0.5,b=0,cex=.5,col=1,lty=3)}
          }else{
            plot(x$xSuR[,1], x$survR[,1,1], col=color, type="l", lty=1, xlab=Xlab,ylab=Ylab, ylim=c(yymin,yymax), main=main, ...)
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
          yymax<-max(x$lamR[-1],na.rm=TRUE)
          yymin<-min(x$lamR[-1],na.rm=TRUE)
        }else{
          yymax<-ylim[2]
          yymin<-ylim[1]
        }
        
        if (conf.bands){
          matplot(x$xR[,1], x$lamR[,,1], col=color, type="l", lty=c(1,2,2), xlab=Xlab,ylab=Ylab, ylim=c(yymin,yymax), main=main, ...)
        }else{
          plot(x$xR[,1], x$lamR[,1,1], col=color, type="l", lty=1, xlab=Xlab,ylab=Ylab, ylim=c(yymin,yymax), main=main, ...)
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
            if (median){abline(a=0.5,b=0,cex=.5,col=1,lty=3)}
          }else{
            plot(x$xR[,1], x$survR[,1,1], col=color, type="l", lty=1, xlab=Xlab,ylab=Ylab, ylim=c(yymin,yymax), main=main, ...)
            if (median){abline(a=0.5,b=0,cex=.5,col=1,lty=3)}
          }
        }else{
          if (conf.bands){
            matplot(x$xSuR[,1], x$survR[,,1], col=color, type="l", lty=c(1,2,2), xlab=Xlab,ylab=Ylab, ylim=c(yymin,yymax), main=main, ...)
            if (median){abline(a=0.5,b=0,cex=.5,col=1,lty=3)}
          }else{
            plot(x$xSuR[,1], x$survR[,1,1], col=color, type="l", lty=1, xlab=Xlab,ylab=Ylab, ylim=c(yymin,yymax), main=main, ...)
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
          yymax<-max(x$lamD[-1],na.rm=TRUE)
          yymin<-min(x$lamD[-1],na.rm=TRUE)
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
