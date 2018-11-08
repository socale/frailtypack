#' Plot predictions using a joint model for longitudinal data and a terminal
#' event or a trivariate joint model for longitudinal data, recurrent events
#' and a terminal event.
#' 
#' Plots predicted probabilities of the event. Confidence intervals are
#' allowed.
#' 
#' 
#' @usage \method{plot}{predLongi}(x, conf.bands=FALSE, pos.legend="topright",
#' cex.legend=0.7, ylim=c(0,1), Xlab = "Time t", Ylab, ...)
#' @param x An object inheriting from \code{predLongi}.
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
"plot.predLongi" <- function (x, conf.bands=FALSE, pos.legend="topright", cex.legend=0.7, ylim=c(0,1), Xlab = "Time t", Ylab, ...)
{
        if ((conf.bands) & (!x$icproba)) stop("Confidence intervals were not calculated. Use the MC.sample argument in the 'prediction' function")

        if (x$moving.window){
			if (missing (Ylab)) legende <- paste("Predicted cumulative probability of event between",x$t,"and time t")
			else legende <- Ylab
            xlim <- c(x$t-(max(x$x.time)-min(x$x.time))*0.1,max(x$x.time))
        }else{
            if (missing (Ylab)) legende <- paste("Predicted probability of event in the next",x$window)
            else legende <- Ylab
			xlim <- c(min(x$x.time),max(x$x.time))
        }


#               if(x$trivariate){
                  title <- paste("Dynamic predictions for survival \n using", x$name.fit)

   #   }else{
  #       title <- paste("Dynamic predictions for survival \n longitudinal outcome and terminal event")
  #    }

        if (conf.bands){
                matplot(x$x.time,t(x$pred),type="l",lty=1,xlab=Xlab,ylab=legende,main=title,ylim=ylim,xlim=xlim)
                matlines(x$x.time,t(x$predLow),type="l",lty=2)
                matlines(x$x.time,t(x$predHigh),type="l",lty=2)
        }else{
                matplot(x$x.time,t(x$pred),type="l",lty=1,xlab=Xlab,ylab=legende,main=title,ylim=ylim,xlim=xlim)
        }

        legend(pos.legend, paste("profile",(1:x$npred)),lty=1,bty="n",col=(1:x$npred))

        if (x$moving.window) abline(v=x$t,lty=2)

        return(invisible())
}
