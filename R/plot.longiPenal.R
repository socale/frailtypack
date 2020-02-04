#' Plot Method for a joint model for longitudinal data and a terminal event.
#' 
#' Plots estimated baseline survival and hazard functions for a terminal
#' outcome from an object of class 'longiPenal'. Confidence bands are allowed.
#' 
#' 
#' @aliases plot.longiPenal lines.longiPenal
#' @usage
#' 
#' \method{plot}{longiPenal}(x, type.plot = "Hazard", conf.bands=TRUE,
#' pos.legend= "topright", cex.legend=0.7, main, color, median=TRUE, Xlab = "Time", Ylab =
#' "Hazard function", ...)
#' @param x A joint model for longitudinal outcome and a terminal event, i.e. a
#' \code{longiPenal} class object (output from calling \code{longiPenal}
#' function).
#' @param type.plot a character string specifying the type of curve for the
#' terminal event. Possible value are "Hazard", or "Survival". The default is
#' "Hazard". Only the first words are required, e.g "Haz", "Su"
#' @param conf.bands Logical value. Determines whether confidence bands will be
#' plotted.  The default is to do so.
#' @param pos.legend The location of the legend can be specified by setting
#' this argument to a single keyword from the list '"bottomright"', '"bottom"',
#' '"bottomleft"', '"left"', '"topleft"', '"top"', '"topright"', '"right"' and
#' '"center"'. The default is '"topright"'
#' @param cex.legend character expansion factor *relative* to current
#' 'par("cex")'. Default is 0.7
#' @param main title of plot
#' @param color color of the curve (integer)
#' @param median Logical value. Determines whether survival median will be plotted. Default is TRUE.
#' @param Xlab Label of x-axis. Default is '"Time"'
#' @param Ylab Label of y-axis. Default is '"Hazard function"'
#' @param ... other unused arguments
#' @return Print a plot for the terminal event of the joint model for a
#' longitudinal and survival data.
#' @seealso \code{\link{longiPenal}}
#' @keywords file
##' @export
#' @examples
#' 
#' 
#' \dontrun{
#' ###--- Joint model for longitudinal data and a terminal event ---###
#' 
#' data(colorectal)
#' data(colorectalLongi)
#' 
#' # Survival data preparation - only terminal events 
#' colorectalSurv <- subset(colorectal, new.lesions == 0)
#' 
#' # Baseline hazard function approximated with splines
#' # Random effects as the link function
#' 
#' model.spli.RE <- longiPenal(Surv(time1, state) ~ age + treatment + who.PS 
#' + prev.resection, tumor.size ~  year * treatment + age + who.PS ,
#' colorectalSurv,	data.Longi = colorectalLongi, random = c("1", "year"),
#' id = "id", link = "Random-effects", left.censoring = -3.33, 
#' n.knots = 7, kappa = 2)
#' pdf(file = "/home/agareb1/etudiants/al10/newpack/test/plot_longi.pdf")
#' 
#' # Plot the estimated baseline hazard function with the confidence intervals
#' plot(model.spli.RE)	
#' 
#' # Plot the estimated baseline hazard function with the confidence intervals
#' plot(model.spli.RE, type = "Survival")	
#' }
#' 
#' 
"plot.longiPenal" <- function (x, type.plot="Hazard", conf.bands=TRUE, pos.legend="topright", cex.legend=0.7, main, color=2, median=TRUE, Xlab = "Time", Ylab = "Hazard function", ...)
{
  
  plot.type <- charmatch(type.plot, c("Hazard", "Survival"),nomatch = 0)
  if (plot.type == 0) {
    stop("estimator must be Hazard or Survival")
  }
  
  if(missing(main))
    main<-""
  
  if(plot.type==1){ # hazard
    
    if(conf.bands){
      matplot(x$xD[-1,1], x$lamD[-1,,1], col=color, type="l", lty=c(1,2,2), xlab=Xlab,ylab=Ylab, main=main, ...)
    }else{
      plot(x$xD[-1,1], x$lamD[-1,1,1], col=color, type="l", lty=1, xlab=Xlab,ylab=Ylab, main=main, ...)
    }
    
  }else{ # survival
    if (missing(Ylab)) Ylab <- "Baseline survival function"
    if (x$typeof == 0){
      if (conf.bands){
        matplot(x$xD[,1], x$survD[,,1], col=color, type="l", lty=c(1,2,2), xlab=Xlab,ylab=Ylab, main=main, ...)
        if (median){abline(a=0.5,b=0,cex=.5,col=1,lty=3)}
      }else{
        plot(x$xD[,1], x$survD[,1,1], col=color, type="l", lty=1, xlab=Xlab,ylab=Ylab, main=main, ...)
        if (median){abline(a=0.5,b=0,cex=.5,col=1,lty=3)}
      }
    }else{
      if (conf.bands){
        matplot(x$xSuD[,1], x$survD[,,1], col=color, type="l", lty=c(1,2,2), xlab=Xlab,ylab=Ylab, main=main, ...)
        if (median){abline(a=0.5,b=0,cex=.5,col=1,lty=3)}
      }else{
        plot(x$xSuD[,1], x$survD[,1,1], col=color, type="l", lty=1, xlab=Xlab,ylab=Ylab, main=main, ...)
        if (median){abline(a=0.5,b=0,cex=.5,col=1,lty=3)}
      }
    }
    
  }
  
  legend(pos.legend, c("event"), lty=1, col=color, cex=cex.legend, ...)
  
  return(invisible())
}
