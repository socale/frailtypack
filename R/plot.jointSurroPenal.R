##' Plot Method for the one-step Joint surrogate model for the evaluation of a 
##' canditate surrogate endpoint.
##' 
##' Plots estimated baseline survival and hazard functions for the surrogate 
##' endpoint and the true endpoint from an object of class 'jointSurroPenal'. 
##' Confidence bands are allowed.
##' 
##' 
##' @aliases plot.jointSurroPenal lines.jointSurroPenal
##' @usage
##' 
##' \method{plot}{jointSurroPenal}(x, type.plot = "Hazard", conf.bands=TRUE,
##' pos.legend = "topright", cex.legend=0.7, main, Xlab = "Time", 
##' Ylab = "Baseline hazard function", median = TRUE, xmin = 0, xmax = NULL, 
##' ylim = c(0,1), endpoint = 2, scale = 1, ...)
##' @param x An object inheriting from \code{jointSurroPenal} class
##' (output from calling the function \code{jointSurroPenal} ).
##' @param type.plot A character string specifying the type of curve. Possible
##' value are "Hazard", or "Survival". The default is "Hazard". Only the first
##' letters are required, e.g "Haz", "Su".
##' @param conf.bands Logical value. Determines whether confidence bands will
##' be plotted.  The default is to do so.
##' @param pos.legend The location of the legend can be specified by setting
##' this argument to a single keyword from the list '"bottomright"',
##' '"bottom"', '"bottomleft"', '"left"', '"topleft"', '"top"', '"topright"',
##' '"right"' and '"center"'. The default is '"topright"'.
##' @param cex.legend Character expansion factor *relative* to current
##' 'par("cex")'. Default is 0.7.
##' @param main Title of plot.
##' @param median Logical value. Determines whether survival median will be plotted. Default is TRUE.
##' @param Xlab Label of x-axis. Default is '"Time"'.
##' @param Ylab Label of y-axis. Default is '"Baseline hazard function"'.
##' @param xmin Minimum value for x-axis, the default is \code{0}. 
##' @param xmax Maximum value for x-axis, the default is \code{NULL}.
##' @param ylim Range of y-axis. Default is from 0 to 1.
##' @param endpoint A binary that indicates the endpoint to represent. \code{0} for
##' the surrogate endpoint, \code{1} for the true endpoint, and \code{2} for both
##' surrogate endpoint and true endpoint. The default is \code{2}.
##' @param scale A numeric that allows to rescale (by multiplication) the survival times. If no change is need the
##' argument is set to 1, the default value. eg: 1/365 aims to convert days to years .
##' @param ... other unused arguments.
##' @return Print a plot of the baseline survival or hazard functions for each
##' type of event or both with the confidence bands or not (conf.bands
##' argument)
##' @seealso \code{\link{jointSurroPenal}, \link{jointSurroCopPenal}}
##' @keywords surrogate
##' @export
##' @examples
##' 
##' 
##' \dontrun{
##' 
##' 
##' ###--- Joint surrogate model ---###
##' ###---evaluation of surrogate endpoints---###
##' 
##' data(dataOvarian)
##' joint.surro.ovar <- jointSurroPenal(data = dataOvarian, n.knots = 8, 
##'                 init.kappa = c(2000,1000), indicator.alpha = 0, 
##'                 nb.mc = 200, scale = 1/365)
##' 
##' # Baseline Hazards fonctions for both the surrogate endpoint 
##' # and the true endpoint
##' plot(joint.surro.ovar,endpoint = 2,type.plot = "Haz", conf.bands = T)   
##' 
##' # Baseline survival fonctions for both the surrogate endpoint 
##' # and the true endpoint
##' plot(joint.surro.ovar,endpoint = 2,type.plot = "Su", conf.bands = T)  
##'              
##' 
##' 
##' }
##' 
##' 
"plot.jointSurroPenal" <- function (x, type.plot="Hazard", conf.bands=TRUE, pos.legend = "topright", 
                                    cex.legend = 0.7, main, Xlab = "Time", Ylab = "Baseline hazard function", 
                                    median = TRUE, xmin = 0, xmax = NULL, ylim = c(0,1), endpoint = 2, scale = 1, ...)
{
  color = 2
  # gestion de l'echelle des temps
  x$xT <- x$xT*scale
  x$xS <- x$xS*scale
  xlim = c(min(x$xT), max(x$xT))
  
  if(is.null(xmax)) xmax <- max(x$xT)
  if((xmax > min(x$xT)) & (xmax < max(x$xT))) xlim[2] <- xmax
  if((xmin > 0) & (xmin < max(x$xT))) xlim[1] <- xmin
  
	plot.type <- charmatch(type.plot, c("Hazard", "Survival"),nomatch = 0)
	if (plot.type == 0) {
		stop("estimator must be 'Hazard' or 'Survival'")
	}	

	if(missing(main)) main<-""  
	
	i <- 1
	if(endpoint == 1){ # true endpoint
  	if(plot.type==1){ # hazard
  		if(conf.bands){
  			matplot(x$xT, x$lamT, col=color, type="l", lty=c(1,2,2), xlab=Xlab,ylab=Ylab, main=main, xlim = xlim, ylim = ylim, ...)
  			matlines(x$xT, x$lamT, col=color+(i-1), type="l", lty=c(1,2,2), ...)
  		}else{
  		  if (missing (Ylab)) Ylab <- "Baseline hazard function"
  			plot(x$xT, x$lamT[,1], col=color, type="l", lty=1, xlab=Xlab,ylab=Ylab, main=main, xlim = xlim, ylim = ylim, ...)
  			lines(x$xT, x$lamT[,1], col=color+(i-1), type="l", lty=1, ...)
  		}
  
  	}else{ # survival
  	 if (missing (Ylab)) Ylab <- "Baseline survival function"
  			if (conf.bands){
  				matplot(x$xT, x$survT, col=color, type="l", lty=c(1,2,2), xlab=Xlab,ylab=Ylab, main=main, xlim = xlim, ylim = ylim, ...)
  				matlines(x$xT, x$survT, col=color+(i-1), type="l", lty=c(1,2,2), ...)
  				if (median){abline(a=0.5,b=0,cex=.5,col=1,lty=3)}
  			}else{
  				plot(x$xT, x$survT[,1], col=color, type="l", lty=1, xlab=Xlab,ylab=Ylab, main=main, xlim = xlim, ylim = ylim, ...)
  				lines(x$xT, x$survT[,1], col=color+(i-1), type="l", lty=1, ...)
  				if (median){abline(a=0.5,b=0,cex=.5,col=1,lty=3)}
  			}
  	}
  	
  	legend(pos.legend, c("For the true endpoint"), lty=1, col=color+(0), cex=cex.legend, ...)
	}
	
	if(endpoint == 0){ # surrogate endpoint
	  if(plot.type==1){ # hazard
	    if(conf.bands){
	      matplot(x$xS, x$lamS, col=color + 1, type="l", lty=c(1,2,2), xlab=Xlab,ylab=Ylab, main=main, xlim = xlim, ylim = ylim, ...)
	      matlines(x$xS, x$lamS, col=color + 1, type="l", lty=c(1,2,2), ...)
	    }else{
	      plot(x$xS, x$lamS[,1], col=color + 1, type="l", lty=1, xlab=Xlab,ylab=Ylab, main=main, xlim = xlim, ylim = ylim, ...)
	      lines(x$xS, x$lamS[,1], col=color + 1, type="l", lty=1, ...)
	    }
	    
	  }else{ # survival
	    if (missing (Ylab)) Ylab <- "Baseline survival function"
	    if (conf.bands){
	      matplot(x$xS, x$survS, col=color + 1, type="l", lty=c(1,2,2), xlab=Xlab,ylab=Ylab, main=main, xlim = xlim, ylim = ylim, ...)
	      matlines(x$xS, x$survS, col=color + 1, type="l", lty=c(1,2,2), ...)
	      if (median){abline(a=0.5,b=0,cex=.5,col=1,lty=3)}
	    }else{
	      plot(x$xS, x$survS[,1], col=color + 1, type="l", lty=1, xlab=Xlab,ylab=Ylab, main=main, xlim = xlim, ylim = ylim, ...)
	      lines(x$xS, x$survS[,1], col=color + 1, type="l", lty=1, ...)
	      if (median){abline(a=0.5,b=0,cex=.5,col=1,lty=3)}
	    }
	  }
	  
	  legend(pos.legend, c("For the surrogate endpoint"), lty=1, col=color + 1, cex=cex.legend, ...)
	}
	
	if(endpoint == 2){ # both surrogate and true endpoint
	  mat <- matrix(0, nrow = nrow(x$lamT), ncol = 6)
	  mat[,c(1:3)] <- x$lamT
	  mat[,c(4:6)] <- x$lamS
	  
	  matS <- matrix(0, nrow = nrow(x$lamT), ncol = 6)
	  matS[,c(1:3)] <- x$survT
	  matS[,c(4:6)] <- x$survS
	  if(plot.type==1){ # hazard
	    if (missing (Ylab)) Ylab <- "Baseline hazard functions"
	    if(conf.bands){
	      matplot(x$xT, mat, col=c(color,color,color,color+1,color+1,color+1), type="l", lty=c(1,2,2,1,2,2), xlab=Xlab,ylab=Ylab, main=main, xlim = xlim, ylim = ylim, ...)
	      matlines(x$xT, mat, col=c(color,color,color,color+1,color+1,color+1), type="l", lty=c(1,2,2,1,2,2), ...)
	    }else{
	      mat2 <- mat[,c(1,4)]
	      matplot(x$xT, mat2, col=c(color,color+1), type="l", lty=c(1,1), xlab=Xlab,ylab=Ylab, main=main, xlim = xlim, ylim = ylim, ...)
	      matlines(x$xT, mat2, col=c(color,color+1), type="l", lty=c(1,1), ...)
	    }
	    
	  }else{ # survival
	    if (missing (Ylab)) Ylab <- "Baseline survival functions"
	    if(conf.bands){
	      matplot(x$xT, matS, col=c(color,color,color,color+1,color+1,color+1), type="l", lty=c(1,2,2,1,2,2), xlab=Xlab,ylab=Ylab, main=main, xlim = xlim, ylim = ylim, ...)
	      matlines(x$xT, matS, col=c(color,color,color,color+1,color+1,color+1), type="l", lty=c(1,2,2,1,2,2), ...)
	      if (median){abline(a=0.5,b=0,cex=.5,col=1,lty=3)}
	    }else{
	      matS2 <- matS[,c(1,4)]
	      matplot(x$xT, matS2, col=c(color,color+1), type="l", lty=c(1,1), xlab = Xlab,ylab = Ylab, main = main, xlim = xlim, ylim = ylim, ...)
	      matlines(x$xT, matS2, col=c(color,color+1), type="l", lty=c(1,1), ...)
	      if (median){abline(a=0.5,b=0,cex=.5,col=1,lty=3)}
	    }
	  }
	  
	  legend(pos.legend, c("For the true endpoint", "For the surrogate endpoint"), lty=c(1,1), col=c(color,color+1), cex=cex.legend, ...)
	}
    return(invisible())
}
