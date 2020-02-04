#' Plot Method for a Nested frailty model.
#' 
#' Plots estimated baseline survival and hazard functions (output from an
#' object of class 'NestedPenal' for nested frailty models). 
#' Confidence bands are allowed.
#' 
#' 
#' @aliases plot.nestedPenal lines.nestedPenal
#' @usage
#' \method{plot}{nestedPenal}(x, type.plot="Hazard", conf.bands=TRUE,
#' pos.legend="topright", cex.legend=0.7, main, color=2, median=TRUE, Xlab = "Time", Ylab =
#' "Hazard function", ...)
#' 
#' @param x A nested model, i.e. an object of class \code{frailtyPenal} for
#' Nested frailty models (output from calling \code{frailtyPenal} function).
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
#' @param main plot title
#' @param color curve color (integer)
#' @param median Logical value. Determines whether survival median will be plotted. Default is TRUE.
#' @param Xlab Label of x-axis. Default is '"Time"'
#' @param Ylab Label of y-axis. Default is '"Hazard function"'
#' @param \dots Other graphical parameters like those in
#' \code{\link{plot.frailtyPenal}}
#' @return Print a plot of the baseline survival or hazard functions with the
#' confidence bands or not (conf.bands argument)
#' @seealso \code{\link{frailtyPenal}}
#' @keywords methods
##' @export
#' @examples
#' 
#' 
#' \dontrun{
#' 
#' data(dataNested)
#' modNested <- frailtyPenal(Surv(t1,t2,event)~cluster(group)+
#' subcluster(subgroup)+cov1+cov2,data=dataNested,n.knots=8,
#' kappa=50000,hazard="Splines")
#' 
#' plot(modNested,conf.bands=FALSE)
#' 
#' }
#' 
#' 
"plot.nestedPenal" <- function (x, type.plot="Hazard", conf.bands=TRUE, pos.legend="topright", cex.legend=0.7, main, color=2, median=TRUE, Xlab = "Time", Ylab = "Hazard function", ...)
{
  
  plot.type <- charmatch(type.plot, c("Hazard", "Survival"),nomatch = 0)
  if (plot.type == 0) {
    stop("estimator must be Hazard or Survival")
  }
  
  if(missing(main)) main<-""
  
  if(plot.type==1){
    if(x$n.strat==1){
      if(conf.bands){
        matplot(x$x[,1], x$lam[,,1], col=color, type="l", lty=c(1,2,2), xlab=Xlab,ylab=Ylab, main=main, ...)
      }else{
        plot(x$x[,1], x$lam[,1,1], col=color, type="l", lty=1, xlab=Xlab,ylab=Ylab, main=main, ...)
      }
    }else{
      if(conf.bands){
        matplot(x$x[,1], x$lam[,,1], col=color, type="l", lty=c(1,2,2), xlab=Xlab,ylab=Ylab, main=main, ...)
        matlines(x$x[,2], x$lam[,,2], col=color+1, type="l", lty=c(1,2,2), ...)
      }else{
        plot(x$x[,1], x$lam[,1,1], col=color, type="l", lty=1, xlab=Xlab,ylab=Ylab, main=main, ...)
        lines(x$x[,2], x$lam[,1,2], col=color+1, type="l", lty=1, ...)
      }
      legend(pos.legend, c("strata = 1", "strata = 2"), lty=1, col=c(color,color+1), cex=cex.legend, ...)
    }
    
  }else{
    if (missing(Ylab)) Ylab <- 'Baseline survival function'
    if (x$n.strat==1){
      if (x$typeof == 0){
        if (conf.bands){
          matplot(x$x[,1], x$surv[,,1], col=color, type="l", lty=c(1,2,2), xlab=Xlab,ylab=Ylab, main=main, ...)
          if (median){abline(a=0.5,b=0,cex=.5,col=1,lty=3)}
        }else{
          plot(x$x[,1], x$surv[,1,1], col=color, type="l", lty=1, xlab=Xlab,ylab=Ylab, main=main, ...)
          if (median){abline(a=0.5,b=0,cex=.5,col=1,lty=3)}
        }
      }else{
        if (conf.bands){
          matplot(x$xSu[,1], x$surv[,,1], col=color, type="l", lty=c(1,2,2), xlab=Xlab,ylab=Ylab, main=main, ...)
          if (median){abline(a=0.5,b=0,cex=.5,col=1,lty=3)}
        }else{
          plot(x$xSu[,1], x$surv[,1,1], col=color, type="l", lty=1, xlab=Xlab,ylab=Ylab, main=main, ...)
          if (median){abline(a=0.5,b=0,cex=.5,col=1,lty=3)}
        }
      }
    }else{
      if (x$typeof == 0){
        if (conf.bands){
          matplot(x$x[,1], x$surv[,,1], col=color, type="l", lty=c(1,2,2), xlab=Xlab,ylab=Ylab, main=main, ...)
          matlines(x$x[,2], x$surv[,,2], col=color+1, type="l", lty=c(1,2,2), ...)
          if (median){abline(a=0.5,b=0,cex=.5,col=1,lty=3)}
        }else{
          plot(x$x[,1], x$surv[,1,1], col=color, type="l", lty=1, xlab=Xlab,ylab=Ylab, main=main, ...)
          lines(x$x[,2], x$surv[,1,2], col=color+1, type="l", lty=1, ...)
          if (median){abline(a=0.5,b=0,cex=.5,col=1,lty=3)}
        }
      }else{
        if (conf.bands){
          matplot(x$xSu[,1], x$surv[,,1], col=color, type="l", lty=c(1,2,2), xlab=Xlab,ylab=Ylab, main=main, ...)
          matlines(x$xSu[,2], x$surv[,,2], col=color+1, type="l", lty=c(1,2,2), ...)
          if (median){abline(a=0.5,b=0,cex=.5,col=1,lty=3)}
        }else{
          plot(x$xSu[,1], x$surv[,1,1], col=color, type="l", lty=1, xlab=Xlab,ylab=Ylab, main=main, ...)
          lines(x$xSu[,2], x$surv[,1,2], col=color+1, type="l", lty=1, ...)
          if (median){abline(a=0.5,b=0,cex=.5,col=1,lty=3)}
        }
      }
      legend(pos.legend, c("strata = 1", "strata = 2"), lty=1, col=c(color,color+1), cex=cex.legend, ...)
    }
  }
  return(invisible())
}
