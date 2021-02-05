#' summary of parameter estimates of a joint frailty model
#' 
#' This function returns hazard rations (HR) and its confidence intervals.
#' 
#' 
#' @aliases summary.jointPenal print.summary.jointPenal
#' @usage \method{summary}{jointPenal}(object, level = 0.95, len = 6, d = 2,
#' lab="hr", ...)
#' @param object output from a call to frailtyPenal for joint models
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
#' data(readmission)
#' 
#' #-- gap-time
#' modJoint.gap <- frailtyPenal(Surv(time,event)~cluster(id)+sex+dukes+
#' charlson+terminal(death),formula.terminalEvent=~sex+dukes+charlson,
#' data=readmission,n.knots=14,kappa=c(9.55e+9,1.41e+12))
#' 
#' #-- calendar time
#' modJoint.calendar <- frailtyPenal(Surv(t.start,t.stop,event)~cluster(id)+
#' sex+dukes+charlson+terminal(death),formula.terminalEvent=~sex+dukes+charlson,
#' data=readmission,n.knots=10,kappa=c(9.55e+9,1.41e+12),recurrentAG=TRUE)
#' 
#' #-- It takes around 1 minute to converge
#' 
#' summary(modJoint.gap)
#' summary(modJoint.calendar)
#' 
#' }
#' 
#' 
"summary.jointPenal"<-
 function(object,level=.95, len=6, d=2, lab="hr", ...)
    {
   
   
   
  if(is.null(object$family)){
        x <- object
        if (!inherits(x, "jointPenal")) 
         stop("Object must be of class 'frailtyPenal'")
      
        z<-abs(qnorm((1-level)/2))
        co <- x$coef
        se <- sqrt(diag(x$varH))[-c(1:(1+x$indic_alpha))]
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

        cat("Recurrences:\n")
        cat("------------- \n")
        cat(paste(rep(" ",mx),collapse=""),paste("   ",dimnames(n)[[2]]),"\n")
        for(i in 1:x$nvarnotdep[1]) 
         {
          lab[i] <- paste(c(rep(" ", mx - nchar(lab[i])), lab[i]),collapse = "")
          cat(lab[i], a[i, 1], "(", a[i, 2], "-", a[i, 3], ") \n")
        }


        cat("\n") 
        cat("Terminal event:\n")
        cat("--------------- \n")
        cat(paste(rep(" ",mx),collapse=""),paste("   ",dimnames(n)[[2]]),"\n")
        for(i in (x$nvarnotdep[1]+1):dd[1]) 
         {
          lab[i] <- paste(c(rep(" ", mx - nchar(lab[i])), lab[i]),collapse = "")
          cat(lab[i], a[i, 1], "(", a[i, 2], "-", a[i, 3], ") \n")
        }
        
  }
  else{
    x <- object
    if (!inherits(x, "jointPenal")) 
      stop("Object must be of class 'frailtyPenal'")
    
    z<-abs(qnorm((1-level)/2))
    co <- x$coef
    se <- sqrt(diag(x$varH))[-c(1:(1+x$indic_alpha))]
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
    }
    else
      lab <- dimnames(n)[[1]]
    
    mx <- max(nchar(lab)) + 1
    
    cat("Recurrences:\n")
    cat("------------- \n")
    cat(  
      paste(rep(" ",mx-2),collapse=""),
      paste("  ",dimnames(n)[[2]][1]),
      paste("",dimnames(n)[[2]][2]),
      paste("",dimnames(n)[[2]][3]),
      "\n"  
    ) 
    #cat(paste(rep(" ",mx),collapse=""),paste("   ",dimnames(n)[[2]]),"\n")
    for(i in 1:x$nvarnotdep[1]) 
    {
      lab[i] <- paste(c(rep(" ", mx - nchar(lab[i])), lab[i]),collapse = "")
      cat(  lab[i], " " ,a[i,1], "  (", a[i,2], ",", a[i,3], " )", sep="")
      cat("\n")
      #cat(lab[i], a[i, 1], "(", a[i, 2], "-", a[i, 3], ") \n")
    }
    
    
    cat("\n") 
    cat("Terminal event:\n")
    cat("--------------- \n")
    cat(  
      paste(rep(" ",mx-2),collapse=""),
      paste("  ",dimnames(n)[[2]][1]),
      paste("",dimnames(n)[[2]][2]),
      paste("",dimnames(n)[[2]][3]),
      "\n"  
    ) 
    #cat(paste(rep(" ",mx),collapse=""),paste("   ",dimnames(n)[[2]]),"\n")
    for(i in (x$nvarnotdep[1]+1):dd[1]) 
    {
      lab[i] <- paste(c(rep(" ", mx - nchar(lab[i])), lab[i]),collapse = "")
      cat(  lab[i], " " ,a[i,1], "  (", a[i,2], ",", a[i,3], " )", sep="")
      cat("\n")
      #cat(lab[i], a[i, 1], "(", a[i, 2], "-", a[i, 3], ") \n")
    }
  }
   
   
   
}
