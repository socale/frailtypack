#' summary of parameter estimates of a shared frailty model
#' 
#' This function returns hazard rations (HR) and its confidence intervals
#' 
#' 
#' @aliases summary.frailtyPenal print.summary.frailtyPenal
#' @usage \method{summary}{frailtyPenal}(object, level = 0.95, len = 6, d = 2,
#' lab="hr", ...)
#' @param object output from a call to frailtyPenal.
#' @param level significance level of confidence interval. Default is 95\%.
#' @param d the desired number of digits after the decimal point. Default of 6
#' digits is used.
#' @param len the total field width. Default is 6.
#' @param lab label of printed results.
#' @param \dots other unused arguments.
#' @return Prints HR and its confidence intervals. Confidence level is allowed
#' (level argument).
#' @seealso \code{\link{frailtyPenal}}
#' @keywords methods
##' @export
#' @examples
#' 
#' 
#' \dontrun{
#' 
#' data(kidney)
#' 
#' ##-- Shared frailty model --##
#' 
#' modSha <- frailtyPenal(Surv(time,status)~age+sex+cluster(id),
#' n.knots=8,kappa=10000,data=kidney,hazard="Splines")
#' 
#' ##-- Cox proportional hazard model --##
#' 
#' modCox <- frailtyPenal(Surv(time,status)~age+sex,
#' n.knots=8,kappa=10000,data=kidney,hazard="Splines")
#' 
#' #-- confidence interval at 95% level (default)
#' 
#' summary(modSha)
#' summary(modCox)
#' 
#' #-- confidence interval at 99% level
#' 
#' summary(modSha,level=0.99)
#' summary(modCox,level=0.99)
#' 
#' }
#' 
#' 
"summary.frailtyPenal" <- function(object,level=.95, len=6, d=2, lab="hr", ...)
{

  
# frailtyPENAL original (NON generalise)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if(is.null(object$family)){
	x <- object
	if (!inherits(x, "frailtyPenal")) 
		stop("Object must be of class 'frailtyPenal'")
	
	nvar <- length(x$coef)
	
	if (is.null(x$coef)){
		cat("     Shared Gamma Frailty model: No covariates and no confidence interval\n")
	}else{
		if (nvar == 0){
			cat("No constant coefficients, only time-varying effects of the covariates \n")
		}else{
			z<-abs(qnorm((1-level)/2))
			co <- x$coef
			if(is.matrix(x$varH)){
				se <- sqrt(diag(x$varH))#[-1]
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
			}else{
				lab <- dimnames(n)[[1]]
			}
	
			mx <- max(nchar(lab)) + 1
			cat(paste(rep(" ",mx),collapse=""),paste("   ",dimnames(n)[[2]]),"\n")
			
			for(i in (1):dd[1]) {
				lab[i] <- paste(c(rep(" ", mx - nchar(lab[i])), lab[i]),collapse = "")
				cat(lab[i], a[i, 1], "(", a[i, 2], "-", a[i, 3], ") \n")
			}
		}
	}
}
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  
# GENfrailtyPENAL (cas du PHM, ie object$family==0)  
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
else if(object$family == 0){
    x <- object
    if (!inherits(x, "frailtyPenal")) 
      stop("Object must be of class 'frailtyPenal'")
    
    nvar <- length(x$coef)
    
    if (is.null(x$coef)){
      cat("     Shared Gamma Frailty model: No covariates and no confidence interval\n")
    }else{
      if (nvar == 0){
        cat("No constant coefficients, only time-varying effects of the covariates \n")
      }else{
        z<-abs(qnorm((1-level)/2))
        co <- x$coef
        if(is.matrix(x$varH)){
          se <- sqrt(diag(x$varH))#[-1]
        }else{
          se <- sqrt(x$varH)
        }
        
        
        # COEFS ORIGINAUX
        #####################################################################
        or <- co
        li <- co - z*se
        ls <- co + z*se
        r <- cbind(or, li, ls)
        dimnames(r) <- list(names(co), c("coefs", paste(level*100,"%",sep=""), "C.I."))
        
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
        cat(  
          paste(rep(" ",mx-2),collapse=""),
          paste("  ",dimnames(n)[[2]][1]),
          paste("",dimnames(n)[[2]][2]),
          paste("",dimnames(n)[[2]][3]),
          "\n"  
        ) 
        
        for(i in (1):dd[1]){
          lab[i] <- paste(c(rep(" ", mx - nchar(lab[i])), lab[i]),collapse = "")
          cat(  lab[i], " " ,a[i,1], "  (", a[i,2], ",", a[i,3], " )", sep="")
          cat("\n")
        }
        cat("\n")
        #####################################################################
        
        
        # HAZARDS RATIOS
        #####################################################################
        or <- exp(co)
        li <- exp(co-z * se)
        ls <- exp(co+z * se)
        r <- cbind(or, li, ls)
        dimnames(r) <- list(names(co), c("hr", paste(level*100,"%",sep=""), "C.I."))
        
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
        cat(  
          paste(rep(" ",mx),collapse=""),
          paste("   ",dimnames(n)[[2]][1]),
          paste("",dimnames(n)[[2]][2]),
          paste("",dimnames(n)[[2]][3]),
          "\n"  
          ) 
        
        for(i in (1):dd[1]) {
          lab[i] <- paste(c(rep(" ", mx - nchar(lab[i])), lab[i]),collapse = "")
          cat(  lab[i], " " ,a[i,1], "  (", a[i,2], ",", a[i,3], " )", sep="")
          cat("\n")
        }
        #####################################################################
      }
    }
}
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  
  
# GENfrailtyPENAL (cas des POM et probit, ie object$family==1 ou object$family==2)  
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
else if(object$family %in% c(1,2,3,4)){
    x <- object
    if (!inherits(x, "frailtyPenal")) 
      stop("Object must be of class 'frailtyPenal'")
    
    nvar <- length(x$coef)
    
    if (is.null(x$coef)){
      cat("     Shared Gamma Frailty model: No covariates and no confidence interval\n")
    }else{
      if (nvar == 0){
        cat("No constant coefficients, only time-varying effects of the covariates \n")
      }else{
        z<-abs(qnorm((1-level)/2))
        co <- x$coef
        if(is.matrix(x$varH)){
          se <- sqrt(diag(x$varH))#[-1]
        }else{
          se <- sqrt(x$varH)
        }
        
        
        # COEFS ORIGINAUX
        #####################################################################
        or <- co
        li <- co - z*se
        ls <- co + z*se
        r <- cbind(or, li, ls)
        dimnames(r) <- list(names(co), c("coefs", paste(level*100,"%",sep=""), "C.I."))
        
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
        cat(  
          paste(rep(" ",mx-2),collapse=""),
          paste("  ",dimnames(n)[[2]][1]),
          paste("",dimnames(n)[[2]][2]),
          paste("",dimnames(n)[[2]][3]),
          "\n"  
        ) 
        
        for(i in (1):dd[1]){
          lab[i] <- paste(c(rep(" ", mx - nchar(lab[i])), lab[i]),collapse = "")
          cat(  lab[i], " " ,a[i,1], "  (", a[i,2], ",", a[i,3], " )", sep="")
          cat("\n")
        }
        #####################################################################
      }
    }
  }
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  
  
	
	
	
	
	}


