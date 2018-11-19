#' summary of parameter estimates of a joint nested frailty model
#' 
#' This function returns hazard rations (HR) and its confidence intervals.
#' 
#' 
#' @aliases summary.jointNestedPenal print.summary.jointNestedPenal
#' @usage \method{summary}{jointNestedPenal}(object, level = 0.95, len = 6, d =
#' 2, lab="hr", ...)
#' @param object output from a call to frailtyPenal for joint nested models
#' @param level significance level of confidence interval. Default is 95\%.
#' @param d the desired number of digits after the decimal point. Default of 6
#' digits is used.
#' @param len the total field width. Default is 6.
#' @param lab label of printed results.
#' @param \dots other unused arguments.
#' @return Prints HR and its confidence intervals for each covariate.
#' Confidence level is allowed (level argument).
#' @seealso \code{\link{frailtyPenal}}
#' @keywords methods
##' @export
#' @examples
#' 
#' 
#' \dontrun{
#' 
#' #-- here is generated cluster (30 clusters)
#' readmissionNested <- transform(readmission,group=id%%30+1)
#' 
#' # Baseline hazard function approximated with splines with calendar-timescale
#' 
#' model.spli.AG <- frailtyPenal(formula = Surv(t.start, t.stop, event)
#'  ~ subcluster(id) + cluster(group) + dukes + terminal(death), 
#'  formula.terminalEvent = ~dukes, data = readmissionNested, 
#'  recurrentAG = TRUE, n.knots = 8, kappa = c(9.55e+9, 1.41e+12),
#'  initialize = TRUE)
#' 
#' summary(model.spli.AG)
#' 
#' }
#' 
 "summary.jointNestedPenal"<-
 function(object,level=.95, len=6, d=2, lab="hr", ...)
    {
        x <- object
        if (!inherits(x, "jointNestedPenal")) 
         stop("Object must be of class 'frailtyPenal'")
		 
		if (is.null(x$coef)) cat("Nested Frailty model: No covariates and no confidence interval\n")
		else {      
			z<-abs(qnorm((1-level)/2))
			co <- x$coef
			
			frail1 <- !is.null(x$theta)
			frail2 <- !is.null(x$eta)
			indic_alpha <- x$indic_alpha
			indic_xi <- x$indic_ksi
			sum_indic <- sum(frail1,frail2,indic_alpha,indic_xi)
			
			se <- sqrt(diag(x$varH))[-c(1:sum_indic)]
			or <- exp(co)
			li <- exp(co-z * se)
			ls <- exp(co+z * se)
			r <- cbind(or, li, ls)
			}

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

        cat("Recurrences:\n")
        cat("------------- \n")
        cat(paste(rep(" ",mx),collapse=""),paste("   ",dimnames(n)[[2]]),"\n")
        for(i in 1:x$nvar[1]) 
         {
          lab[i] <- paste(c(rep(" ", mx - nchar(lab[i])), lab[i]),collapse = "")
          cat(lab[i], a[i, 1], "(", a[i, 2], "-", a[i, 3], ") \n")
        }


        cat("\n") 
        cat("Terminal event:\n")
        cat("--------------- \n")
        cat(paste(rep(" ",mx),collapse=""),paste("   ",dimnames(n)[[2]]),"\n")
        for(i in (x$nvar[1]+1):dd[1]) 
         {
          lab[i] <- paste(c(rep(" ", mx - nchar(lab[i])), lab[i]),collapse = "")
          cat(lab[i], a[i, 1], "(", a[i, 2], "-", a[i, 3], ") \n")
        }
  }
