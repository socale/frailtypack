#' Print a Summary of parameter estimates of a joint model for longitudinal
#' data and a terminal event
#' 
#' Prints a short summary of parameter estimates of a joint model for
#' longitudinal data and a terminal event, an object inheriting from class
#' 'longiPenal'.
#' 
#' 
#' @usage
#' 
#' \method{print}{longiPenal}(x, digits = max(options()$digits - 4, 6), ...)
#' @param x an object inheriting from \code{longiPenal} class
#' @param digits number of digits to print
#' @param \dots other unused arguments
#' @return
#' 
#' Print, separately for each part of the model (longitudinal and terminal) the
#' parameter estimates and details on the estimation.
#' @seealso \code{\link{longiPenal}}
##' @export
#' @keywords methods
"print.longiPenal" <- function (x, digits = max(options()$digits - 4, 6), ...)
{
  
  #for further developments
  #if (x$istop == 1){
  # plot des coefficient dependant du temps
  # if (any(x$nvartimedep != 0)) par(mfrow=c(sum(as.logical(x$nvartimedep)),max(x$nvartimedep)))
  #  if ((x$nvartimedep[1] != 0) & (x$istop == 1)){
  #   for (i in 0:(x$nvartimedep[1]-1)){
  #    matplot(x$BetaTpsMat[,1],x$BetaTpsMat[,(2:4)+4*i],col="blue",type="l",lty=c(1,2,2),xlab="t",ylab="beta(t)",main=paste("Recurrent : ",x$Names.vardep[i+1]),ylim=c(min(x$BetaTpsMat[,-1]),max(x$BetaTpsMat[,-1])))
  #    }
  #  }
  #  if ((x$nvartimedep[2] != 0) & (x$istop == 1)){
  #    for (i in 0:(x$nvartimedep[2]-1)){
  #      matplot(x$BetaTpsMatDc[,1],x$BetaTpsMatDc[,(2:4)+4*i],col="blue",type="l",lty=c(1,2,2),xlab="t",ylab="beta(t)",main=paste("Death : ",x$Names.vardepdc[i+1]),ylim=c(min(x$BetaTpsMatDc[,-1]),max(x$BetaTpsMatDc[,-1])))
  #    }
  #  }
  #}
  
  if (!is.null(cl <- x$call)){
    cat("Call:\n")
    dput(cl)
    
    cat("\n")
  }
  
  
  if (!is.null(x$fail)) {
    cat(" longiPenal failed.", x$fail, "\n")
    return()
  }
  savedig <- options(digits = digits)
  
  on.exit(options(savedig))
  
  coef <- x$coef
  nvar <- length(x$coef) #+x$nvartimedep
  
  #     if (is.null(coef))
  #     {
  #             x$varH<-matrix(x$varH)
  #             x$varHIH<-matrix(x$varHIH)
  #     }
  #AD:
  if (x$typeof == 0){
    if (x$n.knots.temp < 4){
      cat("\n")
      cat("  The minimum number of knots is 4","\n")
      cat("\n")
    }
    if (x$n.knots.temp > 20){
      cat("\n")
      cat("  The maximum number of knots is 20","\n")
    }
  }
  #AD
  
  
  
  if (x$istop == 1){
    if (!is.null(coef)){
      if(nvar != 1){
        seH <- sqrt(diag(x$varH))
        seHIH <- sqrt(diag(x$varHIH))
      }else{
        seH <- sqrt(x$varH)
        seHIH <- sqrt(x$varHIH)
      }
      if (x$typeof == 0){
        tmp <- cbind(coef, exp(coef), seH, seHIH, coef/seH, signif(1 - pchisq((coef/seH)^2, 1), digits - 1))
        if(x$global_chisq.test==1) tmpwald <- cbind(x$global_chisq, x$dof_chisq, ifelse(x$p.global_chisq == 0, "< 1e-16", x$p.global_chisq))
        if(x$global_chisq.test_d==1) tmpwalddc <- cbind(x$global_chisq_d, x$dof_chisq_d, ifelse(x$p.global_chisq_d == 0, "< 1e-16", x$p.global_chisq_d))
       if(x$TwoPart==1)if(x$global_chisq.test_B==1) tmpwaldB <- cbind(x$global_chisq_B,x$dof_chisq_B,x$p.global_chisq_B)#add TwoPart             
      }else{
        tmp <- cbind(coef, exp(coef), seH, coef/seH, signif(1 - pchisq((coef/seH)^2, 1), digits - 1))
        if(x$global_chisq.test==1) tmpwald <- cbind(x$global_chisq, x$dof_chisq, ifelse(x$p.global_chisq == 0, "< 1e-16", x$p.global_chisq))
        if(x$global_chisq.test_d==1) tmpwalddc <- cbind(x$global_chisq_d, x$dof_chisq_d, ifelse(x$p.global_chisq_d == 0, "< 1e-16", x$p.global_chisq_d))
        if(x$TwoPart==1)if(x$global_chisq.test_B==1) tmpwaldB <- cbind(x$global_chisq_B,x$dof_chisq_B,x$p.global_chisq_B)#add TwoPart                   
      }
      
      cat("\n")

      if(x$TwoPart==0){
          cat("  Joint Model for Longitudinal Data and a Terminal Event","\n")

        if (x$typeof == 0){
          cat("  Parameter estimates using a Penalized Likelihood on the hazard function","\n")
        }else{
          cat("  Parameter estimates using a Parametrical approach for the hazard function","\n")
        }
        
      if(x$leftCensoring==TRUE)cat("  and assuming left-censored longitudinal outcome \n")
      }else if(x$TwoPart==1){
      if(x$TwoPart==1)cat("Joint Model for Semi-continuous Data and a Terminal Event\n")
        if (x$typeof == 0){
          cat("  Parameter estimates using a Penalized Likelihood on the hazard function","\n")
        }else{
          cat("  Parameter estimates using a Parametrical approach for the hazard function","\n")
        }
        cat("Conditional Two-Part model with correlated random-effects for Semi-continuous outcome\n")
      }

      if(x$link=='Random-effects') cat("  Association function: random effects","\n")
      if(x$link=='current-level') cat("  Association function: current level of the longitudinal outcome","\n")
      
      
      cat("\n")
      #AL:
      if(x$typeof==0){ind <- 6}
      else {ind <- 5}
      
      
      if(sum(tmp[,ind]<1e-16)>=1){
        d1 <- dim(tmp)[1]
        d2 <- dim(tmp)[2]
        which <- which(tmp[,ind]<1e-16)
        
        sprint<-paste("%.",digits-2,"e",sep="")
        tmp2 <- matrix(c(sprintf(sprint,tmp[,ind])),d1,1)
        
        tmp2[which,1]<-"<1e-16"
        
        sprint<-paste("%.",digits,"f",sep="")
        tmp <- matrix(c(sprintf(sprint,tmp[,-ind])),d1,d2-1)
        tmp <- cbind(tmp,tmp2)
        
      }
      
      if (x$typeof == 0){
        if(x$global_chisq.test==1){
          dimnames(tmpwald) <- list(x$names.factor,c("chisq", "df", "global p"))
          
        }
        if(x$global_chisq.test_d==1){
          dimnames(tmpwalddc) <- list(x$names.factordc,c("chisq", "df", "global p"))
          
        }
        if(x$TwoPart==1){ #add TwoPart
            if(x$global_chisq.test_B==1){
              dimnames(tmpwaldB) <- list(x$names.factorB,c("chisq", "df", "global p"))

            }
        }      
        dimnames(tmp) <- list(names(coef), c("coef", "exp(coef)",
                                             "SE coef (H)", "SE coef (HIH)", "z", "p"))
      }else{
        if(x$global_chisq.test==1){
          dimnames(tmpwald) <- list(x$names.factor,c("chisq", "df", "global p"))
          
        }
        if(x$global_chisq.test_d==1){
          dimnames(tmpwalddc) <- list(x$names.factordc,c("chisq", "df", "global p"))
          

        }
        if(x$TwoPart==1){ #add TwoPart
            if(x$global_chisq.test_B==1){
          dimnames(tmpwaldB) <- list(x$names.factorB,c("chisq", "df", "global p"))

            }  
        }
        dimnames(tmp) <- list(names(coef), c("coef", "exp(coef)",
                                             "SE coef (H)", "z", "p"))
      }
    if(x$TwoPart==1){
        if (x$noVarB == 0){
          cat("Binary outcome:\n")
          cat("---------------- \n")
          prmatrix(tmp[-c(1:(x$nvarY+x$nvarEnd)),-2 ,drop=FALSE],quote=FALSE,right=TRUE)
          if(x$global_chisq.test_B==1){
            cat("\n")
            prmatrix(tmpwaldB)
          }
        }
      cat("\n")
        if (x$noVarY == 0){
          cat("Semi-continuous outcome:\n")
          cat("---------------- \n")
          if(!x$GLMlog){
            prmatrix(tmp[c((x$nvarEnd+1):(x$nvarY+x$nvarEnd)),-2 ,drop=FALSE],quote=FALSE,right=TRUE)
          }else if(x$GLMlog){
            prmatrix(tmp[c((x$nvarEnd+1):(x$nvarY+x$nvarEnd)),,drop=FALSE],quote=FALSE,right=TRUE)
          }
          if(x$global_chisq.test==1){
            cat("\n")
            prmatrix(tmpwald)
          }
        }
    }else if(x$TwoPart==0){
      
        if (x$noVarY == 0){
          cat("Longitudinal outcome:\n")
          cat("---------------- \n")
          if(!x$GLMlog){
          prmatrix(tmp[-c(1:(x$nvarEnd)),-2 ,drop=FALSE],quote=FALSE,right=TRUE)
          }else if(x$GLMlog){
            prmatrix(tmp[-c(1:(x$nvarEnd)),,drop=FALSE],quote=FALSE,right=TRUE)
          }
          if(x$global_chisq.test==1){
            cat("\n")
            prmatrix(tmpwald)
          }
        }
    }   
      cat("\n")
      #for further developments
      #      if (x$nvarnotdep[1] == 0){
      #       cat("Terminal event:\n")
      #        cat("------------- \n")
      #        cat("No constant coefficients, only time-varying effects of the covariates \n")
      #      }else{
      if (x$noVarEnd == 0){
        cat("Terminal event:\n")
        cat("------------- \n")
        prmatrix(tmp[1:x$nvarEnd, ,drop=FALSE],quote=FALSE,right=TRUE)
        if(x$global_chisq.test_d==1){
          cat("\n")
          prmatrix(tmpwalddc)
           
        }
      }
      #  }
      cat("\n")
    }
    
    
    if (x$noVarY == 1){
      cat("\n")
      cat("    Longitudinal outcome: No fixed covariates \n")
      cat("    ----------- \n")
    }
    
    
    
    if (x$noVarEnd == 1){
      cat("\n")
      cat("    Terminal event: No covariates \n")
      cat("    -------------- \n")
      cat("\n")
    }
    
    
    cat(" \n")
    
    cat("Components of Random-effects covariance matrix B1: \n")
    tab.B1 <- round(x$B1,6)
    dimnames(tab.B1) <- list(x$names.re,rep("",dim(x$B1)[1]))
    prmatrix(tab.B1)
    
    
    cat("\n")
    
    
    cat("Association parameters: \n")
    if(x$link=='Random-effects'){
      tab.Asso <- cbind(x$eta, x$se.eta, x$eta/x$se.eta, signif(1 - pchisq((x$eta/x$se.eta)^2, 1), digits - 1))
      
      if(sum(tab.Asso[,4]<1e-16)>1){
        d1 <- dim(tab.Asso)[1]
        d2 <- dim(tab.Asso)[2]
        which <- which(tab.Asso[,4]<1e-16)
        
        sprint<-paste("%.",digits-2,"e",sep="")
        tmp2 <- matrix(c(sprintf(sprint,tab.Asso[,4])),d1,1)
        
        tmp2[which,1]<-"<1e-16"
        sprint<-paste("%.",digits,"f",sep="")
        tab.Asso <- matrix(c(sprintf(sprint,tab.Asso[,-4])),d1,d2-1)
        tab.Asso <- cbind(tab.Asso,tmp2)
        
      }
      
      #  which <- which(tab.Asso[,4]<1e-16)
      #   sprint<-paste("%.",digits,"f",sep="")
      #   tab.Asso <- matrix(c(sprintf(sprint,tab.Asso)),d1,d2)
      
      #  tab.Asso[which,4]<-"<1e-16"
      # tab.Asso <- round(tab.Asso,digits)
      # tab.Asso[which(tab.Asso[,4]==0),4]<-noquote("<1e-16")
     if(x$TwoPart==0) dimnames(tab.Asso) <- list(x$names.re,c("coef",  "SE", "z", "p"))
      if(x$TwoPart==1) dimnames(tab.Asso) <- list(x$names.re,c("coef",  "SE", "z", "p"))
      prmatrix(tab.Asso,quote=FALSE,right=TRUE)
    }else{
      tab.Asso <- cbind(x$eta, x$se.eta, x$eta/x$se.eta, signif(1 - pchisq((x$eta/x$se.eta)^2, 1), digits - 1))
      # d1 <- dim(tab.Asso)[1]
      #  d2 <- dim(tab.Asso)[2]
      # which <- which(tab.Asso[,4]<1e-16)
      # sprint<-paste("%.",digits,"f",sep="")
      # tab.Asso <- matrix(c(sprintf(sprint,tab.Asso)),d1,d2)
      if(sum(tab.Asso[,4]<1e-16)>1){
        d1 <- dim(tab.Asso)[1]
        d2 <- dim(tab.Asso)[2]
        which <- which(tab.Asso[,4]<1e-16)
        
        sprint<-paste("%.",digits-2,"e",sep="")
        tmp2 <- matrix(c(sprintf(sprint,tab.Asso[,4])),d1,1)
        
        tmp2[which,1]<-"<1e-16"
        sprint<-paste("%.",digits,"f",sep="")
        tab.Asso <- matrix(c(sprintf(sprint,tab.Asso[,-4])),d1,d2-1)
        print(tab.Asso)
        print(tmp2)
        tab.Asso <- cbind(tab.Asso,tmp2)
        
      }
      # tab.Asso[which,4]<-"<1e-16"
      # tab.Asso <- round(tab.Asso,digits)
      # tab.Asso[which(tab.Asso[,4]==0),4]<-noquote("<1e-16")
if(x$link=='Current-level'){
     dimnames(tab.Asso) <- list("Current level",c("coef",  "SE", "z", "p"))
     }else if (x$link=="Two-part"){
     dimnames(tab.Asso) <- list(c("Binary part", "Continuous part"),c("coef",  "SE", "z", "p"))
     }
      prmatrix(tab.Asso,quote=FALSE,right=TRUE)
    }
    
    cat("\n")
    
    cat("Residual standard error: ",round(x$ResidualSE,6), " (SE (H): ", round(x$se.ResidualSE,6), ") \n \n")
    
    if (x$typeof == 0){
      cat(paste("      penalized marginal log-likelihood =", round(x$logLikPenal,2)))
      cat("\n")
      cat("      Convergence criteria: \n")
      cat("      parameters =",signif(x$EPS[1],3),"likelihood =",signif(x$EPS[2],3),"gradient =",signif(x$EPS[3],3),"\n")
      cat("\n")
      cat("      LCV = the approximate likelihood cross-validation criterion\n")
      cat("            in the semi parametrical case     =",x$LCV,"\n")
    }else{
      cat(paste("      marginal log-likelihood =", round(x$logLik,2)))
      cat("\n")
      cat("      Convergence criteria: \n")
      cat("      parameters =",signif(x$EPS[1],3),"likelihood =",signif(x$EPS[2],3),"gradient =",signif(x$EPS[3],3),"\n")
      cat("\n")
      cat("      AIC = Aikaike information Criterion     =",x$AIC,"\n")
      cat("\n")
      cat("The expression of the Aikaike Criterion is:","\n")
      cat("        'AIC = (1/n)[np - l(.)]'","\n")
      if (x$typeof == 2){
        cat("\n")
        
        cat("      Scale for the weibull hazard function is :",round(x$scale.weib,2),"\n")
        cat("      Shape for the weibull hazard function is :",round(x$shape.weib,2),"\n")
        
        cat("\n")
        cat("The expression of the Weibull hazard function is:","\n")
        cat("        'lambda(t) = (shape.(t^(shape-1)))/(scale^shape)'","\n")
        cat("The expression of the Weibull survival function is:","\n")
        cat("        'S(t) = exp[- (t/scale)^shape]'")
        cat("\n")
      }
    }
    #AD:
    cat("\n")
    cat("      n=", x$groups)
    cat("\n")
    if(x$TwoPart==0) cat("      n repeated measurements=", x$n.measurements)
    if(x$TwoPart==1) cat("\n      n repeated measurements (Semi-continuous)=", x$n.measurements)
    if(x$TwoPart==1) cat("\n      n repeated measurements (Binary)=", x$n.measurementsB)                   
    if(x$leftCensoring)cat("\n      Percentage of left-censored measurements=",round(x$prop.censored,4)*100,"%\n      Censoring threshold s=",x$leftCensoring.threshold)
    if (length(x$na.action)){
      cat("  (", length(x$na.action), " observation deleted due to missing) \n")
    }else{
      cat("\n")
    }
    
    cat("      n events=", x$n.deaths)
    
    cat( "\n")
    cat("      number of iterations: ", x$n.iter,"\n")
    
    if (x$typeof == 0){
      cat("\n")
      cat("      Exact number of knots used: ", x$n.knots, "\n")
      
      
      cat("      Value of the smoothing parameter: ", x$kappa, sep=" ")
      
      
      #     cat(", DoF: ", formatC(-x$DoF, format="f",digits=2))
    }
    cat("\n")
    cat("      Gaussian quadrature method: ",x$methodGH,"with",x$n.nodes, "nodes",sep=" ", "\n")
    
  }else{ #did not converge
    
    
    cat("   Longitudinal Data and Terminal Event Joint Model parameter estimates ","\n")
    
    if (x$typeof == 0){
      cat("  using a Penalized Likelihood on the hazard function","\n")
    }else{
      cat("  using a Parametrical approach for the hazard function","\n")
    }
    
    if (x$noVarY == 1){
      cat("\n")
      cat("    Longitudinal Outcome: No fixed covariates \n")
      cat("    ----------- \n")
    }
    
    if (x$noVarEnd == 1){
      cat("\n")
      cat("    Terminal event: No covariates \n")
      cat("    -------------- \n")
      cat("\n")
    }
    
    
    
    
    cat("      Convergence criteria: \n")
    cat("      parameters =",signif(x$EPS[1],3),"likelihood =",signif(x$EPS[2],3),"gradient =",signif(x$EPS[3],3),"\n")
    
    cat("\n")
    cat("      n=", x$groups)
    
    if (length(x$na.action)){
      cat("  (", length(x$na.action), " observation deleted due to missing) \n")
    }else{
      cat("\n")
    }
    
    cat("      n events=", x$n.deaths)
    
    cat( "\n")
    cat("      number of iterations: ", x$n.iter)
  }
  cat("\n")
  invisible()
}
