#' Print a Short Summary of parameter estimates of a joint frailty model
#' 
#' Prints a short summary of parameter estimates of a joint frailty model, or
#' more generally an object of class 'frailtyPenal' for joint frailty models.
#' 
#' @importFrom utils type.convert
#' 
#' @usage
#' 
#' \method{print}{jointPenal}(x, digits = max(options()$digits - 4, 6), ...)
#' @param x the result of a call to the jointPenal function
#' @param digits number of digits to print
#' @param \dots other unused arguments
#' @return
#' 
#' Print, separately for each type of event (recurrent and terminal), the
#' parameter estimates of the survival or hazard functions.
#' @seealso \code{\link{frailtyPenal}}
#' @keywords methods
##' @export
"print.jointPenal" <- function (x, digits = max(options()$digits - 4, 6), ...) 
{
  
  

if(is.null(x$family)){ 
  
  if (x$istop == 1){
    # plot des coefficient dependant du temps
    if (any(x$nvartimedep != 0)) par(mfrow=c(1,2))#sum(as.logical(x$nvartimedep)),max(x$nvartimedep)))
    if ((x$nvartimedep[1] != 0) & (x$istop == 1)){
      for (i in 0:(x$nvartimedep[1]-1)){
        matplot(x$BetaTpsMat[,1],x$BetaTpsMat[,(2:4)+4*i],col="blue",type="l",lty=c(1,2,2),xlab="t",ylab="beta(t)",main=paste("Recurrent : ",x$Names.vardep[i+1]),ylim=c(min(x$BetaTpsMat[,-1]),max(x$BetaTpsMat[,-1])))
      }
    }
    if ((x$nvartimedep[2] != 0) & (x$istop == 1)){
      for (i in 0:(x$nvartimedep[2]-1)){
        matplot(x$BetaTpsMatDc[,1],x$BetaTpsMatDc[,(2:4)+4*i],col="blue",type="l",lty=c(1,2,2),xlab="t",ylab="beta(t)",main=paste("Death : ",x$Names.vardepdc[i+1]),ylim=c(min(x$BetaTpsMatDc[,-1]),max(x$BetaTpsMatDc[,-1])))
      }
    }
  }  
  if (!is.null(cl <- x$call)){
    cat("Call:\n")
    dput(cl)
    #if (x$type == "counting" & x$AG == FALSE){ # pour l'instant joint n'accepte pas la vraie troncature a gauche
    #	cat("\n      left truncated structure used")
    #}
    if (x$AG == TRUE){
      cat("\n      Calendar timescale")
    }
    if (x$intcens == TRUE){
      cat("\n      interval censored data used")
    }
    cat("\n")
  }
  if (!is.null(x$fail)) {
    cat(" frailtyPenal failed.", x$fail, "\n")
    return()
  }
  savedig <- options(digits = digits)
  on.exit(options(savedig))
  coef <- x$coef
  nvar <- length(x$coef)
  
  if (is.null(coef)){
    x$varH<-matrix(x$varH)
    x$varHIH<-matrix(x$varHIH)
  }
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
  }else{
    if ((x$typeof == 1) & (x$indic.nb.intR == 1)) cat("  The maximum number of time intervals is 20","\n")
    if ((x$typeof == 1) & (x$indic.nb.intD == 1)) cat("  The maximum number of time intervals is 20","\n")
  }
  #AD
  
  if (x$logNormal == 0) frail <- x$theta
  else frail <- x$sigma2
  indic_alpha <- x$indic_alpha
  
  if (x$istop == 1){
    if (!is.null(coef)){
      if (indic_alpha == 1 || x$joint.clust==2) {
        seH <- sqrt(diag(x$varH))[-c(1,2)]
        seHIH <- sqrt(diag(x$varHIH))[-c(1,2)]
      }else{
        seH <- sqrt(diag(x$varH))[-1]
        seHIH <- sqrt(diag(x$varHIH))[-1]
      }
      if (x$typeof == 0){
        tmp <- cbind(coef, exp(coef), seH, seHIH, coef/seH, ifelse(signif(1 - pchisq((coef/seH)^2, 1), digits - 1) == 0, "< 1e-16", signif(1 - pchisq((coef/seH)^2, 1), digits - 1)))
        if(x$global_chisq.test==1) tmpwald <- cbind(x$global_chisq, x$dof_chisq, ifelse(x$p.global_chisq == 0, "< 1e-16", x$p.global_chisq))
        if(x$global_chisq.test_d==1) tmpwalddc <- cbind(x$global_chisq_d, x$dof_chisq_d, ifelse(x$p.global_chisq_d == 0, "< 1e-16", x$p.global_chisq_d))
      }else{
        tmp <- cbind(coef, exp(coef), seH, coef/seH, ifelse(signif(1 - pchisq((coef/seH)^2, 1), digits - 1) == 0, "< 1e-16", signif(1 - pchisq((coef/seH)^2, 1), digits - 1)))
        if(x$global_chisq.test==1) tmpwald <- cbind(x$global_chisq, x$dof_chisq, ifelse(x$p.global_chisq == 0, "< 1e-16", x$p.global_chisq))
        if(x$global_chisq.test_d==1) tmpwalddc <- cbind(x$global_chisq_d, x$dof_chisq_d, ifelse(x$p.global_chisq_d == 0, "< 1e-16", x$p.global_chisq_d))
      }
      cat("\n")
      if (x$joint.clust == 0) cat("  For clustered data","\n")
      if (x$joint.clust == 0){
        if (x$logNormal == 0){
          cat("  Joint gamma frailty model for a survival and a terminal event processes","\n")
        }else{
          cat("  Joint Log-Normal frailty model for a survival and a terminal event processes","\n")
        }
      }else{
        if ((x$logNormal == 0)&(x$joint.clust==1)){
          cat("  Joint gamma frailty model for recurrent and a terminal event processes","\n")
        }
        else if ((x$logNormal == 0)&(x$joint.clust==2)){
          cat("  General Joint gamma frailty model for recurrent and a terminal event processes","\n")
        }
        else{
          cat("  Joint Log-Normal frailty model for recurrent and a terminal event processes","\n")
        }
      }
      if (x$typeof == 0){
        cat("  using a Penalized Likelihood on the hazard function","\n")
      }else{
        cat("  using a Parametrical approach for the hazard function","\n")
      }
      if (any(x$nvartimedep != 0)) cat("  and some time-dependant covariates","\n")
      if (x$n.strat>1) cat("  (Stratification structure used for recurrences) :",x$n.strat,"strata \n")
      if(x$ncc==TRUE)cat("  and considering weights for the nested case-control design \n")
      if (x$typeof == 0){
        if(x$global_chisq.test==1){
          dimnames(tmpwald) <- list(x$names.factor,c("chisq", "df", "global p"))
          
        }
        if(x$global_chisq.test_d==1){
          dimnames(tmpwalddc) <- list(x$names.factordc,c("chisq", "df", "global p"))
          
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
        dimnames(tmp) <- list(names(coef), c("coef", "exp(coef)",
                                             "SE coef (H)", "z", "p"))
      }
      cat("\n")
      
      if (x$nvarnotdep[1] == 0){
        if (x$joint.clust == 0) cat("Survival event:\n")
        if ((x$joint.clust == 1) | (x$joint.clust == 2)) cat("Recurrences:\n")
        cat("------------- \n")
        cat("No constant coefficients, only time-varying effects of the covariates \n")
      }else{
        if (x$noVar1 == 0){
          if (x$joint.clust == 0) cat("Survival event:\n")
          if (x$joint.clust >= 1) cat("Recurrences:\n")
          cat("------------- \n")
          prmatrix(tmp[1:x$nvarnotdep[1], ,drop=FALSE])
          if(x$global_chisq.test==1){
            cat("\n")
            prmatrix(tmpwald)
          }
        }
      }
      cat("\n")
      
      if (x$nvarnotdep[2] == 0){
        cat("Terminal event:\n")
        cat("---------------- \n")
        cat("No constant coefficients, only time-varying effects of the covariates \n")
      }else{
        if (x$noVar2 == 0){
          cat("Terminal event:\n")
          cat("---------------- \n")
          prmatrix(tmp[-c(1:x$nvarnotdep[1]), ,drop=FALSE])
          if(x$global_chisq.test_d==1){
            cat("\n")
            prmatrix(tmpwalddc)
          }
        }
      }
      cat("\n")
    }
    #theta <- x$theta
    temp <- diag(x$varH)[1]
    seH.frail <- sqrt(((2 * (frail^0.5))^2) * temp) # delta methode
    temp <- diag(x$varHIH)[1]
    seHIH.frail <- sqrt(((2 * (frail^0.5))^2) * temp) # delta methode
    #AD:
    if (x$noVar1 == 1){
      cat("\n")
      if (x$joint.clust == 0) cat("    Survival event: No covariates \n")
      if (x$joint.clust >= 1) cat("    Recurrences: No covariates \n")
      cat("    ----------- \n")
    }
    
    if (x$noVar2 == 1){
      cat("\n")
      cat("    Terminal event: No covariates \n")
      cat("    -------------- \n")
      cat("\n")
    }
    #AD:  
    cat(" Frailty parameters: \n")
    # 		if (x$typeof == 0){
    # 			if (x$logNormal == 0){
    # 				cat("   theta (variance of Frailties, w):", frail, "(SE (H):",seH.frail, ")", "(SE (HIH):", seHIH.frail, ")", "\n")
    # 				if (indic_alpha == 1) cat("   alpha (w^alpha for terminal event):", x$alpha, "(SE (H):",sqrt(diag(x$varH))[2], ")","(SE (HIH):",sqrt(diag(x$varHIH))[2],")","\n")
    # 				else cat("   alpha is fixed (=1) \n")
    # 			}else{
    # 				cat("   sigma square (variance of Frailties, eta):", frail, "(SE (H):",seH.frail, ")", "(SE (HIH):", seHIH.frail, ")", "\n")
    # 				if (indic_alpha == 1) cat("   alpha (exp(alpha.eta) for terminal event):", x$alpha, "(SE (H):",sqrt(diag(x$varH))[2], ")","(SE (HIH):",sqrt(diag(x$varHIH))[2],")","\n")
    # 				else cat("   alpha is fixed (=1) \n")
    # 			}
    # 			cat(" \n")
    # 		}else{
    # 			if (x$logNormal == 0){
    # 				cat("   theta (variance of Frailties, Z):", frail, "(SE (H):",seH.frail, ")", "\n")
    # 				if (indic_alpha == 1) cat("   alpha (w^alpha for terminal event):", x$alpha, "(SE (H):",sqrt(diag(x$varH))[2], ")", "\n")
    # 				else cat("   alpha is fixed (=1) \n")
    # 			}else{
    # 				cat("   sigma square (variance of Frailties, Z):", frail, "(SE (H):",seH.frail, ")", "\n")
    # 				if (indic_alpha == 1) cat("   alpha (exp(alpha.eta) for terminal event):", x$alpha, "(SE (H):",sqrt(diag(x$varH))[2], ")", "\n")
    # 				else cat("   alpha is fixed (=1) \n")
    # 			}
    # 			cat(" \n")
    # 		}
    if (x$logNormal == 0){
      if (indic_alpha == 1 & x$joint.clust<=1){
        cat("   theta (variance of Frailties, w):", frail, "(SE (H):",seH.frail, ")", "p =", ifelse(signif(1 - pnorm(frail/seH.frail), digits - 1) == 0, "< 1e-16", signif(1 - pnorm(frail/seH.frail), digits - 1)), "\n")
        cat("   alpha (w^alpha for terminal event):", x$alpha, "(SE (H):",sqrt(diag(x$varH))[2], ")", "p =", ifelse(signif(1 - pchisq((x$alpha/sqrt(diag(x$varH))[2])^2,1), digits - 1) == 0, "< 1e-16", signif(1 - pchisq((x$alpha/sqrt(diag(x$varH))[2])^2,1), digits - 1)), "\n")
      }else if (x$joint.clust ==2) {
        cat("   theta (variance of u, association between recurrences and terminal event):", frail, "(SE (H):",seH.frail, ")", "p =", ifelse(signif(1 - pnorm(frail/seH.frail), digits - 1) == 0, "< 1e-16", signif(1 - pnorm(frail/seH.frail), digits - 1)), "\n")
        cat("   eta (variance of v, intra-subject correlation):", x$eta, "(SE (H):",sqrt(((2 * (x$eta^0.5))^2) * diag(x$varH)[2]), ")", "p =", ifelse(signif(1 - pnorm (x$eta/sqrt(((2 * (x$eta^0.5))^2) * diag(x$varH)[2]),1), digits - 1) == 0, "< 1e-16", signif(1 - pnorm (x$eta/sqrt(((2 * (x$eta^0.5))^2) * diag(x$varH)[2]),1), digits - 1)), "\n")
      } else {
        cat("   theta (variance of Frailties, w):", frail, "(SE (H):",seH.frail, ")", "p =", ifelse(signif(1 - pnorm(frail/seH.frail), digits - 1) == 0, "< 1e-16", signif(1 - pnorm(frail/seH.frail), digits - 1)), "\n")
        cat("   alpha is fixed (=1) \n")
      }
    }else{
      cat("   sigma square (variance of Frailties, eta):", frail, "(SE (H):",seH.frail, ")", "p =", ifelse(signif(1 - pnorm(frail/seH.frail), digits - 1) == 0, "< 1e-16", signif(1 - pnorm(frail/seH.frail), digits - 1)), "\n")
      if (indic_alpha == 1) cat("   alpha (exp(alpha.eta) for terminal event):", x$alpha, "(SE (H):",sqrt(diag(x$varH))[2], ")", "p =", ifelse(signif(1 - pchisq((x$alpha/sqrt(diag(x$varH))[2])^2,1), digits - 1) == 0, "< 1e-16", signif(1 - pchisq((x$alpha/sqrt(diag(x$varH))[2])^2,1), digits - 1)), "\n")
      else cat("   alpha is fixed (=1) \n")
    }
    cat(" \n")
    
    if (x$typeof == 0){
      cat(paste("Penalized marginal log-likelihood =", round(x$logLikPenal,2)))
      cat("\n")
      cat("   Convergence criteria: \n")
      cat("   parameters =",signif(x$EPS[1],3),"likelihood =",signif(x$EPS[2],3),"gradient =",signif(x$EPS[3],3),"\n")
      cat("\n")
      cat("Likelihood Cross-Validation (LCV) criterion in the semi parametrical case:\n")
      cat("   approximate LCV =",x$LCV,"\n")
    }else{
      cat(paste("   marginal log-likelihood =", round(x$logLik,2)))
      cat("\n")
      cat("   Convergence criteria: \n")
      cat("   parameters =",signif(x$EPS[1],3),"likelihood =",signif(x$EPS[2],3),"gradient =",signif(x$EPS[3],3),"\n")
      cat("\n")
      #			cat("   LCV = the approximate likelihood cross-validation criterion\n")
      #			cat("         in the parametric case     =",x$LCV,"\n")
      cat("   AIC = Aikaike information Criterion     =",x$AIC,"\n")
      cat("\n")
      cat("The expression of the Aikaike Criterion is:","\n")
      cat("        'AIC = (1/n)[np - l(.)]'","\n")
      if (x$typeof == 2){
        cat("\n")
        cat("      Scale for the weibull hazard function is :",round(x$scale.weib[1],2),round(x$scale.weib[2],2),"\n")	
        cat("      Shape for the weibull hazard function is :",round(x$shape.weib[1],2),round(x$shape.weib[2],2),"\n")
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
    if (x$joint.clust == 0){
      cat("   n observations=", x$n, " n subjects=", x$ind, " n groups=", x$groups)
    }else{
      cat("   n observations=", x$n, " n subjects=", x$groups)
    }
    if (length(x$na.action)){
      cat("      (", length(x$na.action), " observation deleted due to missing) \n")
    }else{ 
      cat("\n")
    }
    if (x$joint.clust == 0){
      cat("   n events=", x$n.events)
    }else{
      cat("   n recurrent events=", x$n.events)
    }
    cat("\n")
    cat("   n terminal events=", x$n.deaths)
    cat("\n")
    cat("   n censored events=" ,x$n.censored)
    cat("\n")
    cat("   number of iterations: ", x$n.iter,"\n")
    if (x$logNormal == 0) {
      cat("   Number of nodes for the Gauss-Laguerre quadrature: ", x$nb.gl,"\n")
    }
    else {cat("   Number of nodes for the Gauss-Hermite quadrature: ", x$nb.gh,"\n")}
    
    if ((x$typeof == 1) & (x$indic.nb.intR == 1)){
      cat("   Exact number of time intervals used: 20","\n")
    }else{
      if (x$typeof == 1) cat("   Exact number of time intervals used: ",x$nbintervR,"\n")
    }
    if ((x$typeof == 1) & (x$indic.nb.intD == 1)){
      cat("   Exact number of time intervals used: 20","\n")
    }else{
      if (x$typeof == 1) cat("   Exact number of time intervals used: ",x$nbintervDC,"\n")
    }
    
    if (x$typeof == 0){ 
      cat("\n")
      cat("   Exact number of knots used: ", x$n.knots, "\n")
      cat("   Value of the smoothing parameters: ", x$kappa, sep=" ")
      cat("\n")
    }
  }else{
    if (!is.null(coef)){ 
      cat("\n")
      if (x$joint.clust == 0) cat("  For clustered data","\n")
      if (x$joint.clust == 0){
        if (x$logNormal == 0){
          cat("  Joint gamma frailty model for a survival and a terminal event processes","\n")
        }else{
          cat("  Joint Log-Normal frailty model for a survival and a terminal event processes","\n")
        }
      }else{
        if ((x$logNormal == 0)&(x$joint.clust==1)){
          cat("  Joint gamma frailty model for recurrent and a terminal event processes","\n")
        }
        else if ((x$logNormal == 0)&(x$joint.clust==2)){
          cat("  General Joint gamma frailty model for recurrent and a terminal event processes","\n")
        }
        else{
          cat("  Joint Log-Normal frailty model for recurrent and a terminal event processes","\n")
        }
      }
      if (x$typeof == 0){
        cat("  using a Penalized Likelihood on the hazard function","\n")
      }else{
        cat("  using a Parametrical approach for the hazard function","\n")
      }
      if (any(x$nvartimedep != 0)) cat("  and some time-dependant covariates","\n")
      if (x$noVar1 == 1){
        cat("\n")
        if (x$joint.clust == 0) cat("    Survival event: No covariates \n")
        if (x$joint.clust >= 1) cat("    Recurrences: No covariates \n")
        cat("    ----------- \n")
      }
      
      if (x$noVar2 == 1){
        cat("\n")
        cat("    Terminal event: No covariates \n")
        cat("    -------------- \n")
        cat("\n")
      }
      
      cat("\n")
      
      cat("   Convergence criteria: \n")
      cat("   parameters =",signif(x$EPS[1],3),"likelihood =",signif(x$EPS[2],3),"gradient =",signif(x$EPS[3],3),"\n")
      
      cat("\n")
      cat("   n=", x$n)
      if (length(x$na.action)){
        cat("      (", length(x$na.action), " observation deleted due to missing) \n")
      }else{
        cat("\n")
      }
      if (x$joint.clust == 0){
        cat("   n events=", x$n.events)
      }else{
        cat("   n recurrent events=", x$n.events)
      }
      cat("\n")
      cat("   n terminal events=", x$n.deaths)
      cat("\n")
      if (x$logNormal == 0) {
        cat("   Number of nodes for the Gauss-Laguerre quadrature: ", x$nb.gl,"\n")
      }
      else {cat("   Number of nodes for the Gauss-Hermite quadrature: ", x$nb.gh,"\n")}
    }
  }
  invisible()
  
}
  
  
  
else{
  
  if (x$istop == 1){
    # plot des coefficient dependant du temps
    # if (any(x$nvartimedep != 0)) par(mfrow=c(1,2))#sum(as.logical(x$nvartimedep)),max(x$nvartimedep)))
    if(x$family[2] != 3){
      if ((x$nvartimedep[1] != 0) & (x$istop == 1)){
        for (i in 0:(x$nvartimedep[1]-1)){
          matplot(x$BetaTpsMat[,1],x$BetaTpsMat[,(2:4)+4*i],col="blue",type="l",lty=c(1,2,2),xlab="t",ylab="beta(t)",main=paste("Recurrent : ",x$Names.vardep[i+1]),ylim=c(min(x$BetaTpsMat[,-1]),max(x$BetaTpsMat[,-1])))
        }
      }
    }else{
      if ((x$nvartimedep[1] != 0) & (x$istop == 1)){
        par(mfrow=c(1,2))
        trapz <- function(x,y){
          idx = 2:length(x)
          return (as.double( (x[idx] - x[idx-1]) %*% (y[idx] + y[idx-1])) / 2)
        }
        for (i in 0:(x$nvartimedep[1]-1)){
          matplot(x$BetaTpsMat[,1],x$BetaTpsMat[,(2:4)+4*i],col="blue",type="l",lty=c(1,2,2),xlab="t",ylab="beta(t)",main=paste("Recurrent : ",x$Names.vardep[i+1]),ylim=c(min(x$BetaTpsMat[,-1]),max(x$BetaTpsMat[,-1])))
          nblignes = nrow(x$BetaTpsMat)-1
          matcumul = matrix(NA, nrow = nblignes, ncol = 3)
          abs = x$BetaTpsMat[, 1]
          ord = x$BetaTpsMat[,(2:4)+4*i]
          for(j in 1:nblignes){
            matcumul[j, ] = c(
              trapz(abs[1:(j+1)], ord[1:(j+1), 1]),
              trapz(abs[1:(j+1)], ord[1:(j+1), 2]),
              trapz(abs[1:(j+1)], ord[1:(j+1), 3])
            )
          }
          matplot(x=abs[-1],
                  y=matcumul[,(1:3)],
                  col="blue", type="l", lty=c(1,2,2),
                  xlab="t", 
                  ylab="Cumulative effect", 
                  main=paste("Recurrent : ",x$Names.vardep[i+1])
          )
        }
      }
    }
    if(x$family[1] != 3){
      if ((x$nvartimedep[2] != 0) & (x$istop == 1)){
        for (i in 0:(x$nvartimedep[2]-1)){
          matplot(x$BetaTpsMatDc[,1],x$BetaTpsMatDc[,(2:4)+4*i],col="blue",type="l",lty=c(1,2,2),xlab="t",ylab="beta(t)",main=paste("Death : ",x$Names.vardepdc[i+1]),ylim=c(min(x$BetaTpsMatDc[,-1]),max(x$BetaTpsMatDc[,-1])))
        }
      }
    }else{
      if ((x$nvartimedep[2] != 0) & (x$istop == 1)){
        par(mfrow=c(1,2))
        trapz <- function(x,y){
          idx = 2:length(x)
          return (as.double( (x[idx] - x[idx-1]) %*% (y[idx] + y[idx-1])) / 2)
        }
        for (i in 0:(x$nvartimedep[2]-1)){
          matplot(x$BetaTpsMatDc[,1],x$BetaTpsMatDc[,(2:4)+4*i],col="blue",type="l",lty=c(1,2,2),xlab="t",ylab="beta(t)",main=paste("Death : ",x$Names.vardepdc[i+1]),ylim=c(min(x$BetaTpsMatDc[,-1]),max(x$BetaTpsMatDc[,-1])))
          nblignes = nrow(x$BetaTpsMatDc)-1
          matcumul = matrix(NA, nrow = nblignes, ncol = 3)
          abs = x$BetaTpsMatDc[, 1]
          ord = x$BetaTpsMatDc[,(2:4)+4*i]
          for(j in 1:nblignes){
            matcumul[j, ] = c(
              trapz(abs[1:(j+1)], ord[1:(j+1), 1]),
              trapz(abs[1:(j+1)], ord[1:(j+1), 2]),
              trapz(abs[1:(j+1)], ord[1:(j+1), 3])
            )
          }
          matplot(x=abs[-1],
                  y=matcumul[,(1:3)],
                  col="blue", type="l", lty=c(1,2,2),
                  main=paste("Death : ",x$Names.vardepdc[i+1]),
                  xlab="t", 
                  ylab="Cumulative effect"
          )
        }
      }
    }
  }  
  if (!is.null(cl <- x$call)){
    cat("Call:\n")
    dput(cl)
    #if (x$type == "counting" & x$AG == FALSE){ # pour l'instant joint n'accepte pas la vraie troncature a gauche
    #	cat("\n      left truncated structure used")
    #}
    if (x$AG == TRUE){
      cat("\n      Calendar timescale")
    }
    if (x$intcens == TRUE){
      cat("\n      interval censored data used")
    }
    cat("\n")
  }
  if (!is.null(x$fail)) {
    cat(" frailtyPenal failed.", x$fail, "\n")
    return()
  }
  savedig <- options(digits = digits)
  on.exit(options(savedig))
  coef <- x$coef
  nvar <- length(x$coef)
  
  if (is.null(coef)){
    x$varH<-matrix(x$varH)
    x$varHIH<-matrix(x$varHIH)
  }
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
  }else{
    if ((x$typeof == 1) & (x$indic.nb.intR == 1)) cat("  The maximum number of time intervals is 20","\n")
    if ((x$typeof == 1) & (x$indic.nb.intD == 1)) cat("  The maximum number of time intervals is 20","\n")
  }
  #AD
  
  if (x$logNormal == 0) frail <- x$theta
  else frail <- x$sigma2
  indic_alpha <- x$indic_alpha
  
  if (x$istop == 1){
    if (!is.null(coef)){
      if (indic_alpha == 1 || x$joint.clust==2) {
        seH <- sqrt(diag(x$varH))[-c(1,2)]
        seHIH <- sqrt(diag(x$varHIH))[-c(1,2)]
      }else{
        seH <- sqrt(diag(x$varH))[-1]
        seHIH <- sqrt(diag(x$varHIH))[-1]
      }
      if (x$typeof == 0){
        tmp <- cbind(coef, exp(coef), seH, seHIH, coef/seH, ifelse(signif(1 - pchisq((coef/seH)^2, 1), digits - 1) == 0, "< 1e-16", signif(1 - pchisq((coef/seH)^2, 1), digits - 1)))
        if(x$global_chisq.test==1) tmpwald <- cbind(x$global_chisq, x$dof_chisq, ifelse(x$p.global_chisq == 0, "< 1e-16", x$p.global_chisq))
        if(x$global_chisq.test_d==1) tmpwalddc <- cbind(x$global_chisq_d, x$dof_chisq_d, ifelse(x$p.global_chisq_d == 0, "< 1e-16", x$p.global_chisq_d))
      }else{
        tmp <- cbind(coef, exp(coef), seH, coef/seH, ifelse(signif(1 - pchisq((coef/seH)^2, 1), digits - 1) == 0, "< 1e-16", signif(1 - pchisq((coef/seH)^2, 1), digits - 1)))
        if(x$global_chisq.test==1) tmpwald <- cbind(x$global_chisq, x$dof_chisq, ifelse(x$p.global_chisq == 0, "< 1e-16", x$p.global_chisq))
        if(x$global_chisq.test_d==1) tmpwalddc <- cbind(x$global_chisq_d, x$dof_chisq_d, ifelse(x$p.global_chisq_d == 0, "< 1e-16", x$p.global_chisq_d))
      }
      
      #cat("\n")
      if (x$joint.clust == 0) cat("  For clustered data","\n")
      if (x$joint.clust == 0){
        if (x$logNormal == 0){
          #cat("  Joint gamma frailty model for a survival and a terminal event processes","\n")
          cat("  Generalized Joint Survival Model with Shared Gamma Frailty","\n")
          cat("  for a survival and a terminal event processes")
        }else{
          #cat("  Joint Log-Normal frailty model for a survival and a terminal event processes","\n")
          cat("  Generalized Joint Survival Model with Shared Log-Normal Frailty","\n")
          cat("  a survival and a terminal event processes")
        }
      }else{
        if ((x$logNormal == 0)&(x$joint.clust==1)){
          #cat("  Joint gamma frailty model for recurrent and a terminal event processes","\n")
          cat("  Generalized Joint Survival Model with Shared Gamma Frailty","\n")
          cat("  for recurrent events and a terminal event")
        }
        else if ((x$logNormal == 0)&(x$joint.clust==2)){
          cat("  General Joint gamma frailty model for recurrent and a terminal event processes","\n")
        }
        else{
          #cat("  Joint Log-Normal frailty model for recurrent and a terminal event processes","\n")
          cat("  Generalized Joint Survival Model with Shared Log-Normal Frailty","\n")
          cat("  for recurrent events and a terminal event")
        }
      }
      if (x$typeof == 0){
        cat("  using a Penalized Likelihood on the hazard function","\n")
      }else{
        #cat("  using a Parametrical approach for the hazard function","\n")
        cat("  using parametrical approaches","\n")
      }
      if (any(x$nvartimedep != 0)) cat("  and some time-dependent covariates","\n")
      if (x$n.strat>1) cat("  (Stratification structure used for recurrences) :",x$n.strat,"strata \n")
      if(x$ncc==TRUE)cat("  and considering weights for the nested case-control design \n")
      if (x$typeof == 0){
        if(x$global_chisq.test==1){
          dimnames(tmpwald) <- list(x$names.factor,c("chisq", "df", "global p"))
          
        }
        if(x$global_chisq.test_d==1){
          dimnames(tmpwalddc) <- list(x$names.factordc,c("chisq", "df", "global p"))
          
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
        dimnames(tmp) <- list(names(coef), c("coef", "exp(coef)",
                                             "SE coef (H)", "z", "p"))
      }
      cat("\n")
      
      if (x$nvarnotdep[1] == 0){
        if (x$joint.clust == 0) cat("Survival event:\n")
        if ((x$joint.clust == 1) | (x$joint.clust == 2)) cat("Recurrences:\n")
        cat("------------- \n")
        cat("No constant coefficients, only time-varying effects of the covariates \n")
      }else{
        if (x$noVar1 == 0){
          if (x$joint.clust == 0) cat("Survival event:\n")
          if (x$joint.clust >= 1) cat("\n Recurrences:\n")
          cat("------------- \n")
          if (x$family[2]==0){
            if(!(x$typeof==0)){
              cat("  Parametrical approach with link   g() = log(-log()) ", "\n")
              cat("  S(t) = [ g^-1(eta) ]^frailty ", "\n")
              cat("  eta = shape.log(t) - shape.log(scale) + beta'X", "\n")
              cat("  (Proportional Hazards Frailty Model with a Weibull distribution) ", "\n")
              cat("  Expression of the Weibull hazard function:     'lambda(t) = (shape.(t^(shape-1)))/(scale^shape)'", "\n")
              cat("  Expression of the Weibull survival function:   'S(t) = exp[- (t/scale)^shape]'")
              cat("\n\n")
            }else{
              cat("  Semi-Parametrical approach ", "\n")
              # cat("  S(t) = [ g^-1(eta) ]^frailty ", "\n")
              # cat("  eta = shape.log(t) - shape.log(scale) + beta'X", "\n")
              # cat("  (Proportional Hazards Frailty Model with a Weibull distribution) ", "\n")
              cat("  Expression of the hazard function:     'lambda(t) = lambda_0(t) * exp(beta(t)'X)'", "\n")
              cat("  (Baseline hazard function lambda_0(.) estimated using M-splines)")
              cat("\n\n")
            }
          }
          else if (x$family[2]==1){
            cat("  Parametrical approach with link   g() = -logit() ", "\n")
            cat("  S(t) = [ g^-1(eta) ]^frailty ", "\n")
            cat("  eta = shape.log(t) - shape.log(scale) + beta'X ", "\n")
            cat("  (Proportional Odds Frailty Model with a log-logistic distribution) ", "\n")
            cat("  Expression of the log-logistic hazard function:     'lambda(t) = 1 / [ 1+exp(-eta) ] * d.eta/d.t'", "\n")
            cat("  Expression of the log-logistic survival function:   'S(t) = 1 / [ 1 + (t/scale)^shape ]'")
            cat("\n\n")
          }
          else if (x$family[2]==2){
            cat("  Parametrical approach with link   g() = -PHI^-1() ", "\n")
            cat("  S(t) = [ g^-1(eta) ]^frailty ", "\n")
            cat("  eta = shape.log(t) - shape.log(scale) + beta'X ", "\n")
            cat("  (Probit Frailty Model with a log-normal distribution) ", "\n")
            cat("  Expression of the log-normal hazard function:     'lambda(t) = phi(-eta)/PHI(-eta) * d.eta/d.t'", "\n")
            cat("  Expression of the log-normal survival function:   'S(t) = PHI(-eta)'")
            cat("\n\n")
          }
          else if (x$family[2]==3){
            if(!(x$typeof==0)){
              cat("  Parametrical approach with link   g() = -log() ", "\n")
              cat("  S(t) = [ g^-1(eta) ]^frailty ", "\n")
              cat("  eta = (t/scale)^shape + t*beta'X", "\n")
              cat("  (Additive Hazards Frailty Model with a Weibull distribution) ", "\n")
              cat("  Expression of the Weibull hazard function:     'lambda(t) = (shape.(t^(shape-1)))/(scale^shape)'", "\n")
              cat("  Expression of the Weibull survival function:   'S(t) = exp[- (t/scale)^shape]'")
              cat("\n\n")
            }else{
              cat("  Semi-Parametrical approach ", "\n")
              # cat("  S(t) = [ g^-1(eta) ]^frailty ", "\n")
              # cat("  eta = shape.log(t) - shape.log(scale) + beta'X", "\n")
              # cat("  (Proportional Hazards Frailty Model with a Weibull distribution) ", "\n")
              cat("  Expression of the hazard function:     'lambda(t) = lambda_0(t) + beta(t)'X'", "\n")
              cat("  (Baseline hazard function lambda_0(.) estimated using M-splines)")
              cat("\n\n")
            }
          }
          #browser()
          tmp.mieux = data.frame(tmp)
          tmp.mieux[] = lapply(tmp.mieux, type.convert)
          tmp.mieux.rec = tmp.mieux[1:x$nvarnotdep[1], -2]
          if(x$typeof == 0){
            names(tmp.mieux.rec) = c("coef", "SE coef (H)", "SE coef (HIH)", "z", "p")
          }
          if(x$typeof == 2){
            names(tmp.mieux.rec) = c("coef", "SE coef (H)", "z", "p")
          }
          prmatrix(tmp.mieux.rec, quote=FALSE)
          
          
          if(x$global_chisq.test==1){
            cat("\n")
            tmpwald.mieux = data.frame(tmpwald)
            tmpwald.mieux[] = lapply(tmpwald.mieux, type.convert)
            prmatrix(tmpwald.mieux, quote=FALSE)
          }
          cat("\n")
          if ( (x$family[2]==0)&(!(x$typeof==0)) ){
          cat("Scale for the Weibull hazard function:", round(x$scale.param[1],2), "\n")	
          cat("Shape for the Weibull hazard function:", round(x$shape.param[1],2), "\n")
          }
          else if ( (x$family[2]==1)&(!(x$typeof==0)) ){
            cat("Scale for the log-logistic hazard function:", round(x$scale.param[1],2), "\n")	
            cat("Shape for the log-logistic hazard function:", round(x$shape.param[1],2), "\n")
          }
          else if ( (x$family[2]==2)&(!(x$typeof==0)) ){
            cat("Scale for the log-normal hazard function:", round(x$scale.param[1],2), "\n")	
            cat("Shape for the log-normal hazard function:", round(x$shape.param[1],2), "\n")
          }
          else if ( (x$family[2]==3)&(!(x$typeof==0)) ){
            cat("Scale for the Weibull hazard function:", round(x$scale.param[1],2), "\n")	
            cat("Shape for the Weibull hazard function:", round(x$shape.param[1],2), "\n")
          }
        }
      }
      cat("\n")
      
      if (x$nvarnotdep[2] == 0){
        cat("\n Terminal event:\n")
        cat("---------------- \n")
        cat("No constant coefficients, only time-varying effects of the covariates \n")
      }else{
        if (x$noVar2 == 0){
          cat("\n Terminal event:\n")
          cat("---------------- \n")
          if (x$family[1]==0){
            if(!(x$typeof==0)){
              cat("  Parametrical approach with link   g() = log(-log()) ", "\n")
              cat("  S(t) = [ g^-1(eta) ] ^ (frailty^gamma) ", "\n")
              cat("  eta = shape.log(t) - shape.log(scale) + beta'X", "\n")
              cat("  (Proportional Hazards Frailty Model with a Weibull distribution) ", "\n")
              cat("  Expression of the Weibull hazard function:     'lambda(t) = (shape.(t^(shape-1)))/(scale^shape)'", "\n")
              cat("  Expression of the Weibull survival function:   'S(t) = exp[- (t/scale)^shape]'")
              cat("\n\n")
            }else{
              cat("  Semi-Parametrical approach ", "\n")
              # cat("  S(t) = [ g^-1(eta) ]^frailty ", "\n")
              # cat("  eta = shape.log(t) - shape.log(scale) + beta'X", "\n")
              # cat("  (Proportional Hazards Frailty Model with a Weibull distribution) ", "\n")
              cat("  Expression of the hazard function:     'lambda(t) = lambda_0(t) * exp(beta(t)'X)'", "\n")
              cat("  (Baseline hazard function lambda_0(.) estimated using M-splines)")
              cat("\n\n")
            }
          }
          else if (x$family[1]==1){
            cat("  Parametrical approach with link   g() = -logit() ", "\n")
            cat("  S(t) = [ g^-1(eta) ] ^ (frailty^gamma) ", "\n")
            cat("  eta = shape.log(t) - shape.log(scale) + beta'X ", "\n")
            cat("  (Proportional Odds Frailty Model with a log-logistic distribution) ", "\n")
            cat("  Expression of the log-logistic hazard function:     'lambda(t) = 1 / [ 1+exp(-eta) ] * d.eta/d.t'", "\n")
            cat("  Expression of the log-logistic survival function:   'S(t) = 1 / [ 1 + (t/scale)^shape ]'")
            cat("\n\n")
          }
          else if (x$family[1]==2){
            cat("  Parametrical approach with link   g() = -PHI^-1() ", "\n")
            cat("  S(t) = [ g^-1(eta) ] ^ (frailty^gamma) ", "\n")
            cat("  eta = shape.log(t) - shape.log(scale) + beta'X ", "\n")
            cat("  (Probit Frailty Model with a log-normal distribution) ", "\n")
            cat("  Expression of the log-normal hazard function:     'lambda(t) = phi(-eta)/PHI(-eta) * d.eta/d.t'", "\n")
            cat("  Expression of the log-normal survival function:   'S(t) = PHI(-eta)'")
            cat("\n\n")
          }
          else if (x$family[1]==3){
            if(!(x$typeof==0)){
              cat("  Parametrical approach with link   g() = -log() ", "\n")
              cat("  S(t) = [ g^-1(eta) ]^frailty ", "\n")
              cat("  eta = (t/scale)^shape + t*beta'X", "\n")
              cat("  (Additive Hazards Frailty Model with a Weibull distribution) ", "\n")
              cat("  Expression of the Weibull hazard function:     'lambda(t) = (shape.(t^(shape-1)))/(scale^shape)'", "\n")
              cat("  Expression of the Weibull survival function:   'S(t) = exp[- (t/scale)^shape]'")
              cat("\n\n")
            }else{
              cat("  Semi-Parametrical approach ", "\n")
              # cat("  S(t) = [ g^-1(eta) ]^frailty ", "\n")
              # cat("  eta = shape.log(t) - shape.log(scale) + beta'X", "\n")
              # cat("  (Proportional Hazards Frailty Model with a Weibull distribution) ", "\n")
              cat("  Expression of the hazard function:     'lambda(t) = lambda_0(t) + beta(t)'X'", "\n")
              cat("  (Baseline hazard function lambda_0(.) estimated using M-splines)")
              cat("\n\n")
            }
          }
          tmp.mieux.dc = tmp.mieux[-c(1:x$nvarnotdep[1]), -2]
          if(x$typeof==0){
            names(tmp.mieux.dc) = c("coef", "SE coef (H)", "SE coef (HIH)", "z", "p")
          }
          if(x$typeof==2){
            names(tmp.mieux.dc) = c("coef", "SE coef (H)", "z", "p")
          }
          
          row.names(tmp.mieux.dc) = names(coef)[-c(1:x$nvarnotdep[1])]
          prmatrix(tmp.mieux.dc, quote=FALSE)
          
          
          if(x$global_chisq.test_d==1){
            cat("\n")
            tmpwalddc.mieux = data.frame(tmpwalddc)
            tmpwalddc.mieux[] = lapply(tmpwalddc.mieux, type.convert)
            prmatrix(tmpwalddc.mieux, quote=FALSE)
          }
          cat("\n")
          if ( (x$family[1]==0)&(!(x$typeof==0)) ){
            cat("Scale for the Weibull hazard function:", round(x$scale.param[2],2), "\n")	
            cat("Shape for the Weibull hazard function:", round(x$shape.param[2],2), "\n")
          }
          else if ( (x$family[1]==1)&(!(x$typeof==0)) ){
            cat("Scale for the log-logistic hazard function:", round(x$scale.param[2],2), "\n")	
            cat("Shape for the log-logistic hazard function:", round(x$shape.param[2],2), "\n")
          }
          else if ( (x$family[1]==2)&(!(x$typeof==0)) ){
            cat("Scale for the log-normal hazard function:", round(x$scale.param[2],2), "\n")	
            cat("Shape for the log-normal hazard function:", round(x$shape.param[2],2), "\n")
          }
          else if ( (x$family[1]==3)&(!(x$typeof==0)) ){
            cat("Scale for the Weibull hazard function:", round(x$scale.param[2],2), "\n")	
            cat("Shape for the Weibull hazard function:", round(x$shape.param[2],2), "\n")
          }
        }
      }
      cat("\n")
    }
    #theta <- x$theta
    temp <- diag(x$varH)[1]
    seH.frail <- sqrt(((2 * (frail^0.5))^2) * temp) # delta methode
    temp <- diag(x$varHIH)[1]
    seHIH.frail <- sqrt(((2 * (frail^0.5))^2) * temp) # delta methode
    #AD:
    if (x$noVar1 == 1){
      cat("\n")
      if (x$joint.clust == 0) cat("    Survival event: No covariates \n")
      if (x$joint.clust >= 1) cat("    Recurrences: No covariates \n")
      cat("    ----------- \n")
    }
    
    if (x$noVar2 == 1){
      cat("\n")
      cat("    Terminal event: No covariates \n")
      cat("    -------------- \n")
      cat("\n")
    }
    #AD:  
    cat("\nFrailty parameters: \n")
    # 		if (x$typeof == 0){
    # 			if (x$logNormal == 0){
    # 				cat("   theta (variance of Frailties, w):", frail, "(SE (H):",seH.frail, ")", "(SE (HIH):", seHIH.frail, ")", "\n")
    # 				if (indic_alpha == 1) cat("   alpha (w^alpha for terminal event):", x$alpha, "(SE (H):",sqrt(diag(x$varH))[2], ")","(SE (HIH):",sqrt(diag(x$varHIH))[2],")","\n")
    # 				else cat("   alpha is fixed (=1) \n")
    # 			}else{
    # 				cat("   sigma square (variance of Frailties, eta):", frail, "(SE (H):",seH.frail, ")", "(SE (HIH):", seHIH.frail, ")", "\n")
    # 				if (indic_alpha == 1) cat("   alpha (exp(alpha.eta) for terminal event):", x$alpha, "(SE (H):",sqrt(diag(x$varH))[2], ")","(SE (HIH):",sqrt(diag(x$varHIH))[2],")","\n")
    # 				else cat("   alpha is fixed (=1) \n")
    # 			}
    # 			cat(" \n")
    # 		}else{
    # 			if (x$logNormal == 0){
    # 				cat("   theta (variance of Frailties, Z):", frail, "(SE (H):",seH.frail, ")", "\n")
    # 				if (indic_alpha == 1) cat("   alpha (w^alpha for terminal event):", x$alpha, "(SE (H):",sqrt(diag(x$varH))[2], ")", "\n")
    # 				else cat("   alpha is fixed (=1) \n")
    # 			}else{
    # 				cat("   sigma square (variance of Frailties, Z):", frail, "(SE (H):",seH.frail, ")", "\n")
    # 				if (indic_alpha == 1) cat("   alpha (exp(alpha.eta) for terminal event):", x$alpha, "(SE (H):",sqrt(diag(x$varH))[2], ")", "\n")
    # 				else cat("   alpha is fixed (=1) \n")
    # 			}
    # 			cat(" \n")
    # 		}
    if (x$logNormal == 0){
      if (indic_alpha == 1 & x$joint.clust<=1){
        cat("   theta (variance of frailties, u_i):     ", frail, " (SE(H): ",seH.frail, "), ", 
            "p = ", ifelse(signif(1 - pnorm(frail/seH.frail), digits - 1) == 0, "< 1e-16", signif(1 - pnorm(frail/seH.frail), digits - 1)), sep="", "\n")
        cat("   alpha ((u_i)^alpha for terminal event): ", x$alpha, " (SE(H): ",sqrt(diag(x$varH))[2], "), ", 
            "p = ", ifelse(signif(1 - pchisq((x$alpha/sqrt(diag(x$varH))[2])^2,1), digits - 1) == 0, "< 1e-16", signif(1 - pchisq((x$alpha/sqrt(diag(x$varH))[2])^2,1), digits - 1)),  sep="", "\n")
      }else if (x$joint.clust ==2) {
        cat("   theta (variance of u, association between recurrences and terminal event):", frail, "(SE (H):",seH.frail, ")", "p =", ifelse(signif(1 - pnorm(frail/seH.frail), digits - 1) == 0, "< 1e-16", signif(1 - pnorm(frail/seH.frail), digits - 1)), "\n")
        cat("   eta (variance of v, intra-subject correlation):", x$eta, "(SE (H):",sqrt(((2 * (x$eta^0.5))^2) * diag(x$varH)[2]), ")", "p =", ifelse(signif(1 - pnorm (x$eta/sqrt(((2 * (x$eta^0.5))^2) * diag(x$varH)[2]),1), digits - 1) == 0, "< 1e-16", signif(1 - pnorm (x$eta/sqrt(((2 * (x$eta^0.5))^2) * diag(x$varH)[2]),1), digits - 1)), "\n")
      } else {
        cat("   theta (variance of Frailties, w):", frail, "(SE (H):",seH.frail, ")", "p =", ifelse(signif(1 - pnorm(frail/seH.frail), digits - 1) == 0, "< 1e-16", signif(1 - pnorm(frail/seH.frail), digits - 1)), "\n")
        cat("   gamma is fixed (=1) \n")
      }
    }else{
      cat("   sigma square (variance of Frailties, eta):", frail, "(SE (H):",seH.frail, ")", "p =", ifelse(signif(1 - pnorm(frail/seH.frail), digits - 1) == 0, "< 1e-16", signif(1 - pnorm(frail/seH.frail), digits - 1)), "\n")
      if (indic_alpha == 1) cat("   alpha (exp(alpha.eta) for terminal event):", x$alpha, "(SE (H):",sqrt(diag(x$varH))[2], ")", "p =", ifelse(signif(1 - pchisq((x$alpha/sqrt(diag(x$varH))[2])^2,1), digits - 1) == 0, "< 1e-16", signif(1 - pchisq((x$alpha/sqrt(diag(x$varH))[2])^2,1), digits - 1)), "\n")
      else cat("   alpha is fixed (=1) \n")
    }
    cat(" \n")
    
    if (x$typeof == 0){
      cat(paste("Penalized marginal log-likelihood =", round(x$logLikPenal,2)))
      cat("\n")
      cat("   Convergence criteria: \n")
      cat("   parameters =",signif(x$EPS[1],3),"; likelihood =",signif(x$EPS[2],3),"; gradient =",signif(x$EPS[3],3),"\n")
      cat("\n")
      cat("Likelihood Cross-Validation (LCV) criterion in the semi parametrical case:\n")
      cat("   approximate LCV =",x$LCV,"\n")
    }else{
      cat(paste("Marginal log-likelihood =", round(x$logLik,2)))
      cat("\n")
      cat("   Convergence criteria: \n")
      cat("   parameters =",signif(x$EPS[1],3),"; likelihood =",signif(x$EPS[2],3),"; gradient =",signif(x$EPS[3],3),"\n")
      cat("\n")
      #			cat("   LCV = the approximate likelihood cross-validation criterion\n")
      #			cat("         in the parametric case     =",x$LCV,"\n")
      cat("AIC (Aikaike Information Criterion) =",x$AIC,"\n")
      cat("The expression of the Aikaike Information Criterion is:","\n")
      cat("        'AIC = (1/n)[np - l(.)]'","\n")
      #if (x$typeof == 2){
        #cat("\n")
        #cat("      Scale for the weibull hazard function is :",round(x$scale.weib[1],2),round(x$scale.weib[2],2),"\n")	
        #cat("      Shape for the weibull hazard function is :",round(x$shape.weib[1],2),round(x$shape.weib[2],2),"\n")
        #cat("\n")
        #cat("The expression of the Weibull hazard function is:","\n")
        #cat("        'lambda(t) = (shape.(t^(shape-1)))/(scale^shape)'","\n")
        #cat("The expression of the Weibull survival function is:","\n")
        #cat("        'S(t) = exp[- (t/scale)^shape]'")
        #cat("\n")
      #}
    }
    #AD:
    cat("\n")
    if (x$joint.clust == 0){
      cat("   n.observations =", x$n, " n.subjects =", x$ind, " n.groups =", x$groups)
    }else{
      cat("   n.observations = ", x$n, ", n.subjects = ", x$groups, sep="")
    }
    if (length(x$na.action)){
      cat("      (", length(x$na.action), " observation deleted due to missing) \n")
    }else{ 
      cat("\n")
    }
    if (x$joint.clust == 0){
      cat("   n.events =", x$n.events)
    }else{
      cat("   n.recurrent events =", x$n.events)
    }
    cat("\n")
    cat("   n.terminal events =", x$n.deaths)
    cat("\n")
    cat("   n.censored events =" ,x$n.censored)
    cat("\n")
    cat("   number of iterations:", x$n.iter,"\n")
    if (x$logNormal == 0) {
      cat("   Number of nodes for the Gauss-Laguerre quadrature:", x$nb.gl,"\n")
    }
    else {cat("   Number of nodes for the Gauss-Hermite quadrature:", x$nb.gh,"\n")}
    
    if ((x$typeof == 1) & (x$indic.nb.intR == 1)){
      cat("   Exact number of time intervals used: 20","\n")
    }else{
      if (x$typeof == 1) cat("   Exact number of time intervals used: ",x$nbintervR,"\n")
    }
    if ((x$typeof == 1) & (x$indic.nb.intD == 1)){
      cat("   Exact number of time intervals used: 20","\n")
    }else{
      if (x$typeof == 1) cat("   Exact number of time intervals used:",x$nbintervDC,"\n")
    }
    
    if (x$typeof == 0){ 
      cat("\n")
      cat("   Exact number of knots used:", x$n.knots, "\n")
      cat("   Value of the smoothing parameters:", x$kappa, sep=" ")
      cat("\n")
    }
  }else{
    if (!is.null(coef)){ 
      cat("\n")
      if (x$joint.clust == 0) cat("  For clustered data","\n")
      if (x$joint.clust == 0){
        if (x$logNormal == 0){
          cat("  Joint gamma frailty model for a survival and a terminal event processes","\n")
        }else{
          cat("  Joint Log-Normal frailty model for a survival and a terminal event processes","\n")
        }
      }else{
        if ((x$logNormal == 0)&(x$joint.clust==1)){
          cat("  Joint gamma frailty model for recurrent and a terminal event processes","\n")
        }
        else if ((x$logNormal == 0)&(x$joint.clust==2)){
          cat("  General Joint gamma frailty model for recurrent and a terminal event processes","\n")
        }
        else{
          cat("  Joint Log-Normal frailty model for recurrent and a terminal event processes","\n")
        }
      }
      if (x$typeof == 0){
        cat("  using a Penalized Likelihood on the hazard function","\n")
      }else{
        cat("  using a Parametrical approach for the hazard function","\n")
      }
      if (any(x$nvartimedep != 0)) cat("  and some time-dependant covariates","\n")
      if (x$noVar1 == 1){
        cat("\n")
        if (x$joint.clust == 0) cat("    Survival event: No covariates \n")
        if (x$joint.clust >= 1) cat("    Recurrences: No covariates \n")
        cat("    ----------- \n")
      }
      
      if (x$noVar2 == 1){
        cat("\n")
        cat("    Terminal event: No covariates \n")
        cat("    -------------- \n")
        cat("\n")
      }
      
      cat("\n")
      
      cat("   Convergence criteria: \n")
      cat("   parameters =",signif(x$EPS[1],3),"likelihood =",signif(x$EPS[2],3),"gradient =",signif(x$EPS[3],3),"\n")
      
      cat("\n")
      cat("   n=", x$n)
      if (length(x$na.action)){
        cat("      (", length(x$na.action), " observation deleted due to missing) \n")
      }else{
        cat("\n")
      }
      if (x$joint.clust == 0){
        cat("   n events=", x$n.events)
      }else{
        cat("   n recurrent events=", x$n.events)
      }
      cat("\n")
      cat("   n terminal events=", x$n.deaths)
      cat("\n")
      if (x$logNormal == 0) {
        cat("   Number of nodes for the Gauss-Laguerre quadrature: ", x$nb.gl,"\n")
      }
      else {cat("   Number of nodes for the Gauss-Hermite quadrature: ", x$nb.gh,"\n")}
    }
  }
  invisible()
  
}
  
  
  
}
