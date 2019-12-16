#' Short summary of fixed covariates estimates of a joint model for
#' longitudinal data and a terminal event
#' 
#' This function returns coefficients estimates and their standard error with
#' p-values of the Wald test for the longitudinal outcome and hazard ratios
#' (HR) and their confidence intervals for the terminal event.
#' 
#' 
#' @aliases summary.longiPenal print.summary.longiPenal
#' @usage \method{summary}{longiPenal}(object, level = 0.95, len = 6, d = 2,
#' lab=c("coef","hr"), ...)
#' @param object an object inheriting from \code{longiPenal} class
#' @param level significance level of confidence interval. Default is 95\%.
#' @param d the desired number of digits after the decimal point. Default of 6
#' digits is used.
#' @param len the total field width for the terminal part. Default is 6.
#' @param lab labels of printed results for the longitudinal outcome and the
#' terminal event respectively.
#' @param \dots other unused arguments.
#' @return For the longitudinal outcome it prints the estimates of coefficients
#' of the fixed covariates with their standard error and p-values of the Wald
#' test. For the terminal event it prints HR and its confidence intervals for
#' each covariate. Confidence level is allowed (level argument).
#' @seealso \code{\link{longiPenal}}
#' @keywords methods
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
#' 
#' # Weibull baseline hazard function
#' # Current level of the biomarker as the link function
#' 
#' model.weib.CL <- longiPenal(Surv(time1, state) ~ age + treatment + who.PS
#' + prev.resection, tumor.size ~  year * treatment + age + who.PS , 
#' colorectalSurv, data.Longi = colorectalLongi, random = c("1", "year"),
#' id = "id", link = "Current-level", left.censoring = -3.33, hazard = "Weibull")
#' 	
#' summary(model.spli.RE)
#' summary(model.weib.CL)
#' }
#' 
#' 
"summary.longiPenal"<-
  function(object,level=.95, len=6, d=2, lab=c("coef","hr"), ...)
  {
    x <- object
    if (!inherits(x, "longiPenal"))
      stop("Object must be of class 'longiPenal'")

    z<-abs(qnorm((1-level)/2))
    co <- x$coef
    se <- sqrt(diag(x$varH))
    or <- exp(co)
    li <- exp(co-z * se)
    ls <- exp(co+z * se)
    p <-  signif(1 - pchisq((co/se)^2, 1), 5)

    rl <- cbind(co, se, p)

    dimnames(rl) <- list(names(co), c(lab[1], "SE","p"))

    ddl <- dim(rl)

    al <- formatC(rl, d, len,format="f")

    dim(al) <- ddl
    if(length(ddl) == 1){
      ddl<-c(1,ddl)
      dim(al)<-ddl
      labl<-" "
    }
    else
      labl <- dimnames(rl)[[1]]

    mxl <- max(nchar(labl)) + 1

    al[which(al[,3]==formatC(0, d, len,format="f")),3]<-formatC("<1e-16", d, len,format="f")

    if(x$TwoPart==0){
    cat("Longitudinal outcome:\n")
    cat("------------- \n")
    cat(paste(rep(" ",mxl),collapse=""),paste("  ",dimnames(rl)[[2]]),"\n")
    for(i in (x$nvar[1]+1):ddl[1])
    {
      labl[i] <- paste(c(rep(" ", mxl - nchar(labl[i])), labl[i]),collapse = "")
      cat(labl[i], al[i, 1:3]," \n")
    }



    r <- cbind(or, li, ls)

    dimnames(r) <- list(names(co), c(lab[2], paste(level*100,"%",sep=""), "C.I."))

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
    }
    else{
      cat("Binary outcome:\n")
      cat("------------- \n")
      cat(paste(rep(" ",mxl),collapse=""),paste("  ",dimnames(rl)[[2]]),"\n")
      for(i in (x$nvar[1]+x$nvarY+1):(ddl[1]))
      {
        labl[i] <- paste(c(rep(" ", mxl - nchar(labl[i])), labl[i]),collapse = "")
        cat(labl[i], al[i, 1:3]," \n")
      }
      cat("\n")
      cat("Continuous outcome:\n")
      cat("------------- \n")
      cat(paste(rep(" ",mxl),collapse=""),paste("  ",dimnames(rl)[[2]]),"\n")
      for(i in (x$nvar[1]+1):(ddl[1]-x$nvarB))
      {
        labl[i] <- paste(c(rep(" ", mxl - nchar(labl[i])), labl[i]),collapse = "")
        cat(labl[i], al[i, 1:3]," \n")
      }
      r <- cbind(or, li, ls)
      
      dimnames(r) <- list(names(co), c(lab[2], paste(level*100,"%",sep=""), "C.I."))
      
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
      
    }
    cat("\n")
    cat("Terminal event:\n")
    cat("--------------- \n")
    cat(paste(rep(" ",mx),collapse=""),paste("   ",dimnames(n)[[2]]),"\n")
    for(i in 1:x$nvar[1])
    {
      lab[i] <- paste(c(rep(" ", mx - nchar(lab[i])), lab[i]),collapse = "")
      cat(lab[i], a[i, 1], "(", a[i, 2], "-", a[i, 3], ") \n")
    }


  }
