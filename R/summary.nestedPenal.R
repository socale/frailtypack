#' summary of regression coefficient estimates of a nested frailty model
#' 
#' This function returns hazard rations (HR) and its confidence intervals for
#' each regression coefficient.
#' 
#' 
#' @aliases summary.nestedPenal print.summary.nestedPenal
#' @usage \method{summary}{nestedPenal}(object, level = 0.95, len = 6, d = 2,
#' lab="hr", ...)
#' @param object output from a call to nestedPenal.
#' @param level significance level of confidence interval. Default is 95\%.
#' @param d the desired number of digits after the decimal point. Default of 6
#' digits is used.
#' @param len the total field width. Default is 6.
#' @param lab label of printed results.
#' @param \dots other unused arguments.
#' @return Prints HR and its confidence intervals for each regression
#' coefficient. Confidence level is allowed (level argument).
#' @seealso \code{\link{frailtyPenal}}
#' @keywords methods
#' @export
#' @examples
#' 
#' 
#' \dontrun{
#' 
#' data(dataNested)
#' 
#' modNested <- frailtyPenal(Surv(t1,t2,event)~cluster(group)+
#' subcluster(subgroup)+cov1+cov2,data=dataNested,
#' n.knots=8,kappa=c(50000,50000),hazard="Splines")
#' 
#' #- It takes 90 minutes to converge (depends on processor)
#' 
#' summary(modNested)
#' 
#' }
#' 
#' 
"summary.nestedPenal"<-
 function(object,level=.95, len=6, d=2, lab="hr", ...)
{
	x <- object
	if (!inherits(x, "nestedPenal")) 
		stop("Object must be of class 'frailtyPenal'")
	if (is.null(x$coef)){
		cat("     Nested Frailty model: No covariates and no confidence interval\n")
	}else{
		z<-abs(qnorm((1-level)/2))
		co <- x$coef
		se <- sqrt(diag(x$varH))[-c(1:2)]
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
		}else{
			lab <- dimnames(n)[[1]]
		}
	
		mx <- max(nchar(lab)) + 1
		cat(paste(rep(" ",mx),collapse=""),paste("   ",dimnames(n)[[2]]),"\n")
	
		for(i in (1):dd[1]){
			lab[i] <- paste(c(rep(" ", mx - nchar(lab[i])), lab[i]),collapse = "")
			cat(lab[i], a[i, 1], "(", a[i, 2], "-", a[i, 3], ") \n")
		}
	}

}


