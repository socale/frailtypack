#' summary of parameter estimates of an additive frailty model
#' 
#' This function returns hazard ratios (HR) and its confidence intervals
#' 
#' 
#' @aliases summary.additivePenal print.summary.additivePenal
#' @usage \method{summary}{additivePenal}(object, level = 0.95, len = 6, d = 2,
#' lab="hr", \dots)
#' @param object output from a call to additivePenal.
#' @param level significance level of confidence interval. Default is 95\%.
#' @param d the desired number of digits after the decimal point. Default of 6
#' digits is used.
#' @param len the total field width. Default is 6.
#' @param lab label of printed results.
#' @param \dots other unused arguments.
#' @return Prints HR and its confidence intervals for each covariate.
#' Confidence level is allowed (level argument)
#' @seealso \code{\link{additivePenal}}
#' @keywords methods
##' @export
#' @examples
#' 
#' 
#' \dontrun{
#' 
#' data(dataAdditive)
#' 
#' modAdd <- additivePenal(Surv(t1,t2,event)~cluster(group)+var1+slope(var1),
#' correlation=TRUE,data=dataAdditive,n.knots=8,kappa=862,hazard="Splines")
#' 
#' #- 'var1' is boolean as a treatment variable.
#' 
#' summary(modAdd)
#' 
#' }
#' 
#' 
"summary.additivePenal"<-
 function(object,level=.95, len=6, d=2, lab="hr", ...)
{
	x <- object
	if (!inherits(x, "additivePenal")) 
		stop("Object must be of class 'additivePenal'")
	
	z<-abs(qnorm((1-level)/2))
	co <- x$coef

	if(is.matrix(x$varH)){
		se <- sqrt(diag(x$varH))
	}else{
		se <- sqrt(x$varH)
	}

	or <- exp(co)
	li <- exp(co-z * se)
	ls <- exp(co+z * se)
	r <- cbind(or, li, ls)
	dimnames(r) <- list(names(co), c(lab, paste(level*100,"%",sep=""), "C.I."))
	
	n<-r
	
	dd <- dim(n)
	n[n > 999.99] <- Inf
	a <- formatC(n, d, len,format="f")
	
	dim(a) <- dd
	if(length(dd) == 1){
		dd<-c(1,dd)
		dim(a)<-dd
		lab<-" "
		}
	else
		lab <- dimnames(n)[[1]]
	
	mx <- max(nchar(lab)) + 1
	cat(paste(rep(" ",mx),collapse=""),paste("   ",dimnames(n)[[2]]),"\n")
	for(i in (1):dd[1]){
		lab[i] <- paste(c(rep(" ", mx - nchar(lab[i])), lab[i]),collapse = "")
		cat(lab[i], a[i, 1], "(", a[i, 2], "-", a[i, 3], ") \n")
	}
      
}


