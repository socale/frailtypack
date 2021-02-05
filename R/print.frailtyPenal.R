#' Print a Short Summary of parameter estimates of a shared frailty model
#' 
#' Prints a short summary of parameter estimates of a 'frailtyPenal' object
#' 
#' 
#' @usage \method{print}{frailtyPenal}(x, digits = max(options()$digits - 4, 6),
#' ...)
#' @param x the result of a call to the frailtyPenal function.
#' @param digits number of digits to print.
#' @param \dots other unused arguments.
#' @return
#' 
#' Print the parameter estimates of the survival or hazard functions.
#' @seealso \code{\link{frailtyPenal}}
#' @keywords methods
##' @export
"print.frailtyPenal" <- function (x, digits = max(options()$digits - 4, 6), ...) 
{

  
  
# frailtyPENAL original (NON generalise)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if(is.null(x$family)){
  if (x$istop == 1){
# plot des coefficient dependant du temps
	if ((x$nvartimedep != 0) & (x$istop == 1)){
		par(mfrow=c(1,x$nvartimedep))
		for (i in 0:(x$nvartimedep-1)){
			matplot(x$BetaTpsMat[,1],x$BetaTpsMat[,(2:4)+4*i],col="blue",type="l",lty=c(1,2,2),xlab="t",ylab="beta(t)",main=x$Names.vardep[i+1],ylim=c(min(x$BetaTpsMat[,-1]),max(x$BetaTpsMat[,-1])))
		}
	}
  }

	if (!is.null(cl <- x$call)){
		cat("Call:\n")
		dput(cl)
		if ((x$type == "counting" | x$type == "intervaltronc") & (x$AG == FALSE)){
			cat("\n      left truncated structure used")
		}
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
	nvar <- length(x$coef) #+x$nvartimedep
 
#	if (is.null(coef))
#	{
#		x$varH<-matrix(x$varH) 
#		x$varHIH<-matrix(x$varHIH)
#	}
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
		if ((x$typeof == 1) & (x$indic.nb.int == 1)) cat("  The maximum number of time intervals is 20","\n")
	}
#AD

	if (x$logNormal == 0) frail <- x$theta
	else frail <- x$sigma2

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
		    tmp <- cbind(coef, exp(coef), seH, seHIH, coef/seH, ifelse(signif(1 - pchisq((coef/seH)^2, 1), digits - 1) == 0, "< 1e-16", signif(1 - pchisq((coef/seH)^2, 1), digits - 1))) # ifelse pour eviter que R n'affiche pval = 0 lorsque celle-ci est tres petite
		    if(x$global_chisq.test==1) tmpwald <- cbind(x$global_chisq,x$dof_chisq,ifelse(x$p.global_chisq == 0, "< 1e-16", x$p.global_chisq))
		  }else{
		    tmp <- cbind(coef, exp(coef), seH, coef/seH, ifelse(signif(1 - pchisq((coef/seH)^2, 1), digits - 1) == 0, "< 1e-16", signif(1 - pchisq((coef/seH)^2, 1), digits - 1)))
		    if(x$global_chisq.test==1) tmpwald <- cbind(x$global_chisq,x$dof_chisq,ifelse(x$p.global_chisq == 0, "< 1e-16", x$p.global_chisq))
		  }
		
			cat("\n")
			if (!is.null(frail)){
				if (x$logNormal == 0){
					cat("  Shared Gamma Frailty model parameter estimates ","\n")
				}else{
					cat("  Shared Log-Normal Frailty model parameter estimates ","\n")
				}
				if (x$typeof == 0){
					cat("  using a Penalized Likelihood on the hazard function","\n")
				}else{
					cat("  using a Parametrical approach for the hazard function","\n")
				}
				if (x$nvartimedep != 0) cat("  and some time-dependant covariates","\n")
				if (x$n.strat>1) cat("  (Stratification structure used) :",x$n.strat,"strata \n")
			}else{
				if (x$typeof == 0){
					cat("  Cox proportional hazards model parameter estimates ","\n")
					cat("  using a Penalized Likelihood on the hazard function","\n")
				}else{
					cat("  Cox proportional hazards model parameter estimates ","\n")
					cat("  using a Parametrical approach for the hazard function","\n")
				}
				if (x$nvartimedep != 0) cat("  and some time-dependant covariates","\n")
				if (x$n.strat>1) cat("  (Stratification structure used) :",x$n.strat,"strata \n")
			}
			
			if (x$typeof == 0){
				if(x$global_chisq.test==1){
					dimnames(tmpwald) <- list(x$names.factor,c("chisq", "df", "global p"))
					
				}
				dimnames(tmp) <- list(names(coef), c("coef", "exp(coef)",
				"SE coef (H)", "SE coef (HIH)", "z", "p"))

			}else{
				if(x$global_chisq.test==1){
					dimnames(tmpwald) <- list(x$names.factor,c("chisq", "df", "global p"))
					
				}
				dimnames(tmp) <- list(names(coef), c("coef", "exp(coef)",
				"SE coef (H)", "z", "p"))

			}
			cat("\n")
#AL:
			if (nvar == 0){
				cat("No constant coefficients, only time-varying effects of the covariates \n")
			}else{
				prmatrix(tmp)
				if(x$global_chisq.test==1){
					cat("\n")
					prmatrix(tmpwald)
				}
			}
#AL
			cat("\n")
		}
	
		if (!is.null(frail)) {
			#tetha <- x$theta
			#temp <- x$varTheta[1]
			seH <- sqrt(x$varTheta[1]) #sqrt(((2 * (frail^0.5))^2) * temp)
			#temp <- x$varTheta[2]
			seHIH <- sqrt(x$varTheta[2]) #sqrt(((2 * (frail^0.5))^2) * temp)
			
			
#AD:
			if (x$noVar1 == 1){
				cat("\n")
				if (x$logNormal == 0){ cat("    Shared Gamma Frailty model: No covariates \n") }
				else { cat("    Shared Log-Normal Frailty model: No covariates \n") }
				cat("    -------------------------- \n")
				cat("\n")
			}
#AD:
# 			if (x$typeof == 0){
# 				if (x$logNormal == 0){
# 					cat("    Frailty parameter, Theta:", frail, "(SE (H):",
# 					seH, ")", "(SE (HIH):", seHIH, ")", "\n")
# 				}else{
# 					cat("    Frailty parameter, Sigma Square:", frail, "(SE (H):",
# 					seH, ")", "(SE (HIH):", seHIH, ")", "\n")
# 				}
# 			}else{
# 				if (x$logNormal == 0){
# 					cat("    Frailty parameter, Theta:", frail, "(SE (H):",
# 					seH, ")", "\n")
# 				}else{
# 					cat("    Frailty parameter, Sigma Square:", frail, "(SE (H):",
# 					seH, ")", "\n")
# 				}
# 			}

			if (x$logNormal == 0){
			  cat("    Frailty parameter, Theta:", frail, "(SE (H):",
			  seH, ")", "p =", ifelse(signif(1 - pnorm(frail/seH), digits - 1) == 0, "< 1e-16", signif(1 - pnorm(frail/seH), digits - 1)), "\n")
			}else{
			  cat("    Frailty parameter, Sigma Square:", frail, "(SE (H):",
			  seH, ")", "p =", ifelse(signif(1 - pnorm(frail/seH), digits - 1) == 0, "< 1e-16", signif(1 - pnorm(frail/seH), digits - 1)), "\n")
			}
			
        }

		cat(" \n")
#AD:	
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
				if (x$n.strat == 1){
					cat("      Scale for the weibull hazard function is :",round(x$scale.weib[1],2),"\n")
					cat("      Shape for the weibull hazard function is :",round(x$shape.weib[1],2),"\n")
				}else{
					cat("      Scale for the weibull hazard function is :",round(x$scale.weib[1],2),round(x$scale.weib[2],2),"\n")
					cat("      Shape for the weibull hazard function is :",round(x$shape.weib[1],2),round(x$shape.weib[2],2),"\n")
				}
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
		cat("      n=", x$n)
		
		if (length(x$na.action)){
			cat("  (", length(x$na.action), " observation deleted due to missing) \n")
		}else{
			cat("\n")
		}
		
		if (!is.null(frail)) cat("      n events=", x$n.events, " n groups=", x$groups)
		else cat("      n events=", x$n.events)
		
		cat( "\n")
		cat("      number of iterations: ", x$n.iter,"\n")
		if ((x$typeof == 1) & (x$indic.nb.int == 1)){
			cat("      Exact number of time intervals used: 20","\n")
		 }else{
		 	if (x$typeof == 1) cat("      Exact number of time intervals used: ",x$nbintervR,"\n")
		 }	
		if (x$typeof == 0){ 
			cat("\n")
			cat("      Exact number of knots used: ", x$n.knots, "\n")
	
			if (!x$cross.Val){
				cat("      Value of the smoothing parameter: ", x$kappa, sep=" ")
			}
		
			if (x$cross.Val){
				if (is.null(frail)){
					cat("      Smoothing parameter estimated by Cross validation: ", x$kappa, sep=" ")
				}else{ 
					cat("      Best smoothing parameter estimated by")
					cat("\n")
					cat("      an approximated Cross validation: ", x$kappa, sep=" ")
				}
			}
			cat(", DoF: ", formatC(-x$DoF, format="f",digits=2))
		}
	}else{
		if (!is.null(frail)){
			if (x$logNormal == 0){
				cat("  Shared Gamma Frailty model parameter estimates ","\n")
			}else{
				cat("  Shared Log-Normal Frailty model parameter estimates ","\n")
			}
			if (x$typeof == 0){
				cat("  using a Penalized Likelihood on the hazard function","\n")
			}else{
				cat("  using a Parametrical approach for the hazard function","\n")
			}
			if (x$nvartimedep != 0) cat("  and some time-dependant covariates","\n")
			if (x$n.strat>1) cat("  (Stratification structure used)", "\n")
			if (x$noVar1 == 1){
				cat("\n")
				if (x$logNormal == 0){ cat("    Shared Gamma Frailty model: No covariates \n") }
				else { cat("    Shared Log-Normal Frailty model: No covariates \n") }
				cat("    -------------------------- \n")
				cat("\n")
			}

		}else{
			if (x$typeof == 0){
				cat("  Cox proportional hazards model parameter estimates ","\n")
				cat("  using a Penalized Likelihood on the hazard function","\n")
			}else{
				cat("  Cox proportional hazards model parameter estimates ","\n")
				cat("  using a Parametrical approach for the hazard function","\n")
			}
			if (x$nvartimedep != 0) cat("  and some time-dependant covariates","\n")
			if (x$n.strat>1) cat("  (Stratification structure used)", "\n")
		}
		cat("\n")
		
		cat("      Convergence criteria: \n")
		cat("      parameters =",signif(x$EPS[1],3),"likelihood =",signif(x$EPS[2],3),"gradient =",signif(x$EPS[3],3),"\n")
		
		cat("\n")
		cat("      n=", x$n)
		
		if (length(x$na.action)){
			cat("  (", length(x$na.action), " observation deleted due to missing) \n")
		}else{
			cat("\n")
		}

		if (!is.null(frail)) cat("      n events=", x$n.events, " n groups=", x$groups)
		else cat("      n events=", x$n.events)
		
		cat( "\n")
		cat("      number of iterations: ", x$n.iter)
	}
	cat("\n")
	invisible()
}	
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~	
	
	
	
	
	
# GENfrailtyPENAL (cas du PHM, ie object$family==0)  
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
else if(x$family == 0){	
  if (x$istop == 1){
    # plot des coefficient dependant du temps
    if ((x$nvartimedep != 0) & (x$istop == 1)){
      par(mfrow=c(1,x$nvartimedep))
      for (i in 0:(x$nvartimedep-1)){
        matplot(x$BetaTpsMat[,1],x$BetaTpsMat[,(2:4)+4*i],col="blue",type="l",lty=c(1,2,2),xlab="t",ylab="beta(t)",main=x$Names.vardep[i+1],ylim=c(min(x$BetaTpsMat[,-1]),max(x$BetaTpsMat[,-1])))
      }
    }
  }
  
  if (!is.null(cl <- x$call)){
    cat("Call:\n")
    dput(cl)
    if ((x$type == "counting" | x$type == "intervaltronc") & (x$AG == FALSE)){
      cat("\n      left truncated structure used")
    }
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
  nvar <- length(x$coef) #+x$nvartimedep
  
  #	if (is.null(coef))
  #	{
  #		x$varH<-matrix(x$varH) 
  #		x$varHIH<-matrix(x$varHIH)
  #	}
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
    if ((x$typeof == 1) & (x$indic.nb.int == 1)) cat("  The maximum number of time intervals is 20","\n")
  }
  #AD
  
  if (x$logNormal == 0) frail <- x$theta
  else frail <- x$sigma2
  
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
        tmp <- cbind(coef, exp(coef), seH, seHIH, coef/seH, ifelse(signif(1 - pchisq((coef/seH)^2, 1), digits - 1) == 0, "< 1e-16", signif(1 - pchisq((coef/seH)^2, 1), digits - 1))) # ifelse pour eviter que R n'affiche pval = 0 lorsque celle-ci est tres petite
        if(x$global_chisq.test==1) tmpwald <- cbind(x$global_chisq,x$dof_chisq,ifelse(x$p.global_chisq == 0, "< 1e-16", x$p.global_chisq))
      }else{
        tmp <- cbind(coef, exp(coef), seH, coef/seH, ifelse(signif(1 - pchisq((coef/seH)^2, 1), digits - 1) == 0, "< 1e-16", signif(1 - pchisq((coef/seH)^2, 1), digits - 1)))
        if(x$global_chisq.test==1) tmpwald <- cbind(x$global_chisq,x$dof_chisq,ifelse(x$p.global_chisq == 0, "< 1e-16", x$p.global_chisq))
      }
      
      cat("\n")
      if (!is.null(frail)){
        if (x$logNormal == 0){
          #cat("  Shared Gamma Frailty model parameter estimates ","\n")
          cat("  Generalized Survival Model with Shared Gamma Frailty ", "\n")
        }else{
          #cat("  Shared Log-Normal Frailty model parameter estimates ","\n")
          cat("  Generalized Survival Model with Shared Log-Normal Frailty ", "\n")
        }
        if (x$typeof == 0){
          cat("  using a Penalized Likelihood on the hazard function","\n")
        }else{
          #cat("  using a Parametrical approach for the hazard function","\n")
          cat("  Parametrical approach with link   g() = log(-log()) ", "\n")
          cat("  S(t) = [ g^-1(eta) ]^frailty ", "\n")
          cat("  eta = shape.log(t) - shape.log(scale) + beta'X", "\n")
          cat("  (Proportional Hazards Frailty Model with a Weibull distribution) ", "\n\n")
        }
        if (x$nvartimedep != 0) cat("  and some time-dependant covariates","\n")
        if (x$n.strat>1) cat("  (Stratification structure used) :",x$n.strat,"strata \n")
      }else{
        if (x$typeof == 0){
          #cat("  Cox proportional hazards model parameter estimates ","\n")
          cat("  Generalized Survival Model ","\n")
          cat("  using a Penalized Likelihood on the hazard function ","\n")
        }else{
          #cat("  Cox proportional hazards model parameter estimates ","\n")
          cat("  Generalized Survival Model ","\n")
          cat("  Parametrical approach with link   g() = log(-log()) ", "\n")
          cat("  g(S(t)) = eta = shape.log(t) - shape.log(scale) + beta'X ", "\n")
          cat("  (Proportional Hazards Model with a Weibull distribution) ", "\n\n")
          #cat("  using a Parametrical approach for the hazard function","\n")
          #cat("  using a Parametrical approach ","\n")
        }
        if (x$nvartimedep != 0) cat("  and some time-dependant covariates","\n")
        if (x$n.strat>1) cat("  (Stratification structure used) :",x$n.strat,"strata \n")
      }
      
      if (x$typeof == 0){
        if(x$global_chisq.test==1){
          dimnames(tmpwald) <- list(x$names.factor,c("chisq", "df", "global p"))
          
        }
        dimnames(tmp) <- list(names(coef), c("coef", "exp(coef)",
                                             "SE coef (H)", "SE coef (HIH)", "z", "p"))
        
      }else{
        if(x$global_chisq.test==1){
          dimnames(tmpwald) <- list(x$names.factor,c("chisq", "df", "global p"))
          
        }
        dimnames(tmp) <- list(names(coef), c("coef", "exp(coef)",
                                             "SE coef (H)", "z", "p"))
        
      }
      cat("\n")
      #AL:
      if (nvar == 0){
        cat("No constant coefficients, only time-varying effects of the covariates \n")
      }else{
        prmatrix(tmp)
        if(x$global_chisq.test==1){
          cat("\n")
          prmatrix(tmpwald)
        }
      }
      #AL
      cat("\n")
    }
    
    if (!is.null(frail)) {
      #tetha <- x$theta
      #temp <- x$varTheta[1]
      seH <- sqrt(x$varTheta[1]) #sqrt(((2 * (frail^0.5))^2) * temp)
      #temp <- x$varTheta[2]
      seHIH <- sqrt(x$varTheta[2]) #sqrt(((2 * (frail^0.5))^2) * temp)
      
      
      #AD:
      if (x$noVar1 == 1){
        cat("\n")
        if (x$logNormal == 0){ cat("    Shared Gamma Frailty model: No covariates \n") }
        else { cat("    Shared Log-Normal Frailty model: No covariates \n") }
        cat("    -------------------------- \n")
        cat("\n")
      }
      #AD:
      # 			if (x$typeof == 0){
      # 				if (x$logNormal == 0){
      # 					cat("    Frailty parameter, Theta:", frail, "(SE (H):",
      # 					seH, ")", "(SE (HIH):", seHIH, ")", "\n")
      # 				}else{
      # 					cat("    Frailty parameter, Sigma Square:", frail, "(SE (H):",
      # 					seH, ")", "(SE (HIH):", seHIH, ")", "\n")
      # 				}
      # 			}else{
      # 				if (x$logNormal == 0){
      # 					cat("    Frailty parameter, Theta:", frail, "(SE (H):",
      # 					seH, ")", "\n")
      # 				}else{
      # 					cat("    Frailty parameter, Sigma Square:", frail, "(SE (H):",
      # 					seH, ")", "\n")
      # 				}
      # 			}
      
      if (x$logNormal == 0){
        cat("    Frailty parameter, Theta:", frail, "(SE (H):",
            seH, ")", "p =", ifelse(signif(1 - pnorm(frail/seH), digits - 1) == 0, "< 1e-16", signif(1 - pnorm(frail/seH), digits - 1)), "\n")
      }else{
        cat("    Frailty parameter, Sigma Square:", frail, "(SE (H):",
            seH, ")", "p =", ifelse(signif(1 - pnorm(frail/seH), digits - 1) == 0, "< 1e-16", signif(1 - pnorm(frail/seH), digits - 1)), "\n")
      }
      
    }
    
    cat(" \n")
    #AD:	
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
        if (x$n.strat == 1){
          cat("      Scale for the Weibull hazard function is :",round(x$scale.weib[1],2),"\n")
          cat("      Shape for the Weibull hazard function is :",round(x$shape.weib[1],2),"\n")
        }else{
          cat("      Scale for the Weibull hazard function is :",round(x$scale.weib[1],2),round(x$scale.weib[2],2),"\n")
          cat("      Shape for the Weibull hazard function is :",round(x$shape.weib[1],2),round(x$shape.weib[2],2),"\n")
        }
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
    cat("      n=", x$n)
    
    if (length(x$na.action)){
      cat("  (", length(x$na.action), " observation deleted due to missing) \n")
    }else{
      cat("\n")
    }
    
    if (!is.null(frail)) cat("      n events=", x$n.events, " n groups=", x$groups)
    else cat("      n events=", x$n.events)
    
    cat( "\n")
    cat("      number of iterations: ", x$n.iter,"\n")
    if ((x$typeof == 1) & (x$indic.nb.int == 1)){
      cat("      Exact number of time intervals used: 20","\n")
    }else{
      if (x$typeof == 1) cat("      Exact number of time intervals used: ",x$nbintervR,"\n")
    }	
    if (x$typeof == 0){ 
      cat("\n")
      cat("      Exact number of knots used: ", x$n.knots, "\n")
      
      if (!x$cross.Val){
        cat("      Value of the smoothing parameter: ", x$kappa, sep=" ")
      }
      
      if (x$cross.Val){
        if (is.null(frail)){
          cat("      Smoothing parameter estimated by Cross validation: ", x$kappa, sep=" ")
        }else{ 
          cat("      Best smoothing parameter estimated by")
          cat("\n")
          cat("      an approximated Cross validation: ", x$kappa, sep=" ")
        }
      }
      cat(", DoF: ", formatC(-x$DoF, format="f",digits=2))
    }
  }else{
    if (!is.null(frail)){
      if (x$logNormal == 0){
        #cat("  Shared Gamma Frailty model parameter estimates ","\n")
        cat("  Generalized Survival Model with Shared Gamma Frailty ", "\n")
      }else{
        #cat("  Shared Log-Normal Frailty model parameter estimates ","\n")
        cat("  Generalized Survival Model with Shared Log-Normal Frailty ", "\n")
      }
      if (x$typeof == 0){
        cat("  using a Penalized Likelihood on the hazard function","\n")
      }else{
        #cat("  using a Parametrical approach for the hazard function","\n")
        cat("  using a Parametrical approach ","\n")
      }
      if (x$nvartimedep != 0) cat("  and some time-dependant covariates","\n")
      if (x$n.strat>1) cat("  (Stratification structure used)", "\n")
      if (x$noVar1 == 1){
        cat("\n")
        if (x$logNormal == 0){ cat("    Shared Gamma Frailty model: No covariates \n") }
        else { cat("    Shared Log-Normal Frailty model: No covariates \n") }
        cat("    -------------------------- \n")
        cat("\n")
      }
      
    }else{
      if (x$typeof == 0){
        #cat("  Cox proportional hazards model parameter estimates ","\n")
        cat("  Generalized Survival Model ","\n")
        cat("  using a Penalized Likelihood on the hazard function","\n")
      }else{
        #cat("  Cox proportional hazards model parameter estimates ","\n")
        cat("  Generalized Survival Model ","\n")
        #cat("  using a Parametrical approach for the hazard function","\n")
        cat("  using a Parametrical ","\n")
      }
      if (x$nvartimedep != 0) cat("  and some time-dependant covariates","\n")
      if (x$n.strat>1) cat("  (Stratification structure used)", "\n")
    }
    cat("\n")
    
    cat("      Convergence criteria: \n")
    cat("      parameters =",signif(x$EPS[1],3),"likelihood =",signif(x$EPS[2],3),"gradient =",signif(x$EPS[3],3),"\n")
    
    cat("\n")
    cat("      n=", x$n)
    
    if (length(x$na.action)){
      cat("  (", length(x$na.action), " observation deleted due to missing) \n")
    }else{
      cat("\n")
    }
    
    if (!is.null(frail)) cat("      n events=", x$n.events, " n groups=", x$groups)
    else cat("      n events=", x$n.events)
    
    cat( "\n")
    cat("      number of iterations: ", x$n.iter)
  }
  cat("\n")
  invisible()	
}
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~    
  
  
  
 
   
# GENfrailtyPENAL (cas du POM, ie object$family==1)  
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
else if(x$family == 1){	
    if (x$istop == 1){
      # plot des coefficient dependant du temps
      if ((x$nvartimedep != 0) & (x$istop == 1)){
        par(mfrow=c(1,x$nvartimedep))
        for (i in 0:(x$nvartimedep-1)){
          matplot(x$BetaTpsMat[,1],x$BetaTpsMat[,(2:4)+4*i],col="blue",type="l",lty=c(1,2,2),xlab="t",ylab="beta(t)",main=x$Names.vardep[i+1],ylim=c(min(x$BetaTpsMat[,-1]),max(x$BetaTpsMat[,-1])))
        }
      }
    }
    
    if (!is.null(cl <- x$call)){
      cat("Call:\n")
      dput(cl)
      if ((x$type == "counting" | x$type == "intervaltronc") & (x$AG == FALSE)){
        cat("\n      left truncated structure used")
      }
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
    nvar <- length(x$coef) #+x$nvartimedep
    
    #	if (is.null(coef))
    #	{
    #		x$varH<-matrix(x$varH) 
    #		x$varHIH<-matrix(x$varHIH)
    #	}
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
      if ((x$typeof == 1) & (x$indic.nb.int == 1)) cat("  The maximum number of time intervals is 20","\n")
    }
    #AD
    
    if (x$logNormal == 0) frail <- x$theta
    else frail <- x$sigma2
    
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
          tmp <- cbind(coef, exp(coef), seH, seHIH, coef/seH, ifelse(signif(1 - pchisq((coef/seH)^2, 1), digits - 1) == 0, "< 1e-16", signif(1 - pchisq((coef/seH)^2, 1), digits - 1))) # ifelse pour eviter que R n'affiche pval = 0 lorsque celle-ci est tres petite
          if(x$global_chisq.test==1) tmpwald <- cbind(x$global_chisq,x$dof_chisq,ifelse(x$p.global_chisq == 0, "< 1e-16", x$p.global_chisq))
        }else{
          tmp <- cbind(coef, exp(coef), seH, coef/seH, ifelse(signif(1 - pchisq((coef/seH)^2, 1), digits - 1) == 0, "< 1e-16", signif(1 - pchisq((coef/seH)^2, 1), digits - 1)))
          if(x$global_chisq.test==1) tmpwald <- cbind(x$global_chisq,x$dof_chisq,ifelse(x$p.global_chisq == 0, "< 1e-16", x$p.global_chisq))
        }
        
        cat("\n")
        if (!is.null(frail)){
          if (x$logNormal == 0){
            cat("  Generalized Survival Model with Shared Gamma Frailty ", "\n")
            #cat("  Shared Gamma Frailty model parameter estimates ","\n")
          }else{
            cat("  Generalized Survival Model with Shared Log-Normal Frailty ", "\n")
            #cat("  Shared Log-Normal Frailty model parameter estimates ","\n")
          }
          if (x$typeof == 0){
            cat("  using a Penalized Likelihood on the hazard function","\n")
          }else{
            cat("  Parametrical approach with link   g() = -logit() ", "\n")
            cat("  S(t) = [ g^-1(eta) ]^frailty ", "\n")
            cat("  eta = shape.log(t) - shape.log(scale) + beta'X ", "\n")
            cat("  (Proportional Odds Frailty Model with a log-logistic distribution) ", "\n\n")
            #cat("  using a Parametrical approach for the hazard function","\n")
          }
          if (x$nvartimedep != 0) cat("  and some time-dependant covariates","\n")
          if (x$n.strat>1) cat("  (Stratification structure used) :",x$n.strat,"strata \n")
        }else{
          if (x$typeof == 0){
            cat("  Cox proportional hazards model parameter estimates ","\n")
            cat("  using a Penalized Likelihood on the hazard function","\n")
          }else{
            cat("  Generalized Survival Model ","\n")
            cat("  Parametrical approach with link   g() = -logit() ", "\n")
            cat("  g(S(t)) = eta = shape.log(t) - shape.log(scale) + beta'X ", "\n")
            cat("  (Proportional Odds Model with a log-logistic distribution) ", "\n\n")
            #cat("  Cox proportional hazards model parameter estimates ","\n")
            #cat("  using a Parametrical approach for the hazard function","\n")
          }
          if (x$nvartimedep != 0) cat("  and some time-dependant covariates","\n")
          if (x$n.strat>1) cat("  (Stratification structure used) :",x$n.strat,"strata \n")
        }
        
        if (x$typeof == 0){
          if(x$global_chisq.test==1){
            dimnames(tmpwald) <- list(x$names.factor,c("chisq", "df", "global p"))
            
          }
          dimnames(tmp) <- list(names(coef), c("coef", "exp(coef)",
                                               "SE coef (H)", "SE coef (HIH)", "z", "p"))
          
        }else{
          if(x$global_chisq.test==1){
            dimnames(tmpwald) <- list(x$names.factor,c("chisq", "df", "global p"))
            
          }
          dimnames(tmp) <- list(names(coef), c("coef", "exp(coef)",
                                               "SE coef (H)", "z", "p"))
          
        }
        cat("\n")
        #AL:
        if (nvar == 0){
          cat("No constant coefficients, only time-varying effects of the covariates \n")
        }else{
          prmatrix(tmp[, -2])
          if(x$global_chisq.test==1){
            cat("\n")
            prmatrix(tmpwald)
          }
        }
        #AL
        cat("\n")
      }
      
      if (!is.null(frail)) {
        #tetha <- x$theta
        #temp <- x$varTheta[1]
        seH <- sqrt(x$varTheta[1]) #sqrt(((2 * (frail^0.5))^2) * temp)
        #temp <- x$varTheta[2]
        seHIH <- sqrt(x$varTheta[2]) #sqrt(((2 * (frail^0.5))^2) * temp)
        
        
        #AD:
        if (x$noVar1 == 1){
          cat("\n")
          if (x$logNormal == 0){ cat("    Shared Gamma Frailty model: No covariates \n") }
          else { cat("    Shared Log-Normal Frailty model: No covariates \n") }
          cat("    -------------------------- \n")
          cat("\n")
        }
        #AD:
        # 			if (x$typeof == 0){
        # 				if (x$logNormal == 0){
        # 					cat("    Frailty parameter, Theta:", frail, "(SE (H):",
        # 					seH, ")", "(SE (HIH):", seHIH, ")", "\n")
        # 				}else{
        # 					cat("    Frailty parameter, Sigma Square:", frail, "(SE (H):",
        # 					seH, ")", "(SE (HIH):", seHIH, ")", "\n")
        # 				}
        # 			}else{
        # 				if (x$logNormal == 0){
        # 					cat("    Frailty parameter, Theta:", frail, "(SE (H):",
        # 					seH, ")", "\n")
        # 				}else{
        # 					cat("    Frailty parameter, Sigma Square:", frail, "(SE (H):",
        # 					seH, ")", "\n")
        # 				}
        # 			}
        
        if (x$logNormal == 0){
          cat("    Frailty parameter, Theta:", frail, "(SE (H):",
              seH, ")", "p =", ifelse(signif(1 - pnorm(frail/seH), digits - 1) == 0, "< 1e-16", signif(1 - pnorm(frail/seH), digits - 1)), "\n")
        }else{
          cat("    Frailty parameter, Sigma Square:", frail, "(SE (H):",
              seH, ")", "p =", ifelse(signif(1 - pnorm(frail/seH), digits - 1) == 0, "< 1e-16", signif(1 - pnorm(frail/seH), digits - 1)), "\n")
        }
        
      }
      
      cat(" \n")
      #AD:	
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
          if (x$n.strat == 1){
            cat("      Scale for the log-logistic survival function is :",round(x$scale.weib[1],2),"\n")
            cat("      Shape for the log-logistic survival function is :",round(x$shape.weib[1],2),"\n")
          }else{
            cat("      Scale for the log-logistic survival function is :",round(x$scale.weib[1],2),round(x$scale.weib[2],2),"\n")
            cat("      Shape for the log-logistic survival function is :",round(x$shape.weib[1],2),round(x$shape.weib[2],2),"\n")
          }
          cat("\n")
          cat("The expression of the log-logistic survival function is:","\n")
          cat("        'S(t) = 1 / [ 1 + (t/scale)^shape ]'", "\n")
          cat("The expression of the log-logistic hazard function is:","\n")
          cat("        'lambda(t) = 1 / [ 1+exp(-eta) ] * d.eta/d.t' ","\n")
          cat("\n")
        }
      }
      #AD:
      cat("\n")
      cat("      n=", x$n)
      
      if (length(x$na.action)){
        cat("  (", length(x$na.action), " observation deleted due to missing) \n")
      }else{
        cat("\n")
      }
      
      if (!is.null(frail)) cat("      n events=", x$n.events, " n groups=", x$groups)
      else cat("      n events=", x$n.events)
      
      cat( "\n")
      cat("      number of iterations: ", x$n.iter,"\n")
      if ((x$typeof == 1) & (x$indic.nb.int == 1)){
        cat("      Exact number of time intervals used: 20","\n")
      }else{
        if (x$typeof == 1) cat("      Exact number of time intervals used: ",x$nbintervR,"\n")
      }	
      if (x$typeof == 0){ 
        cat("\n")
        cat("      Exact number of knots used: ", x$n.knots, "\n")
        
        if (!x$cross.Val){
          cat("      Value of the smoothing parameter: ", x$kappa, sep=" ")
        }
        
        if (x$cross.Val){
          if (is.null(frail)){
            cat("      Smoothing parameter estimated by Cross validation: ", x$kappa, sep=" ")
          }else{ 
            cat("      Best smoothing parameter estimated by")
            cat("\n")
            cat("      an approximated Cross validation: ", x$kappa, sep=" ")
          }
        }
        cat(", DoF: ", formatC(-x$DoF, format="f",digits=2))
      }
    }else{
      if (!is.null(frail)){
        if (x$logNormal == 0){
          cat("  Shared Gamma Frailty model parameter estimates ","\n")
        }else{
          cat("  Shared Log-Normal Frailty model parameter estimates ","\n")
        }
        if (x$typeof == 0){
          cat("  using a Penalized Likelihood on the hazard function","\n")
        }else{
          cat("  using a Parametrical approach for the hazard function","\n")
        }
        if (x$nvartimedep != 0) cat("  and some time-dependant covariates","\n")
        if (x$n.strat>1) cat("  (Stratification structure used)", "\n")
        if (x$noVar1 == 1){
          cat("\n")
          if (x$logNormal == 0){ cat("    Shared Gamma Frailty model: No covariates \n") }
          else { cat("    Shared Log-Normal Frailty model: No covariates \n") }
          cat("    -------------------------- \n")
          cat("\n")
        }
        
      }else{
        if (x$typeof == 0){
          cat("  Cox proportional hazards model parameter estimates ","\n")
          cat("  using a Penalized Likelihood on the hazard function","\n")
        }else{
          cat("  Cox proportional hazards model parameter estimates ","\n")
          cat("  using a Parametrical approach for the hazard function","\n")
        }
        if (x$nvartimedep != 0) cat("  and some time-dependant covariates","\n")
        if (x$n.strat>1) cat("  (Stratification structure used)", "\n")
      }
      cat("\n")
      
      cat("      Convergence criteria: \n")
      cat("      parameters =",signif(x$EPS[1],3),"likelihood =",signif(x$EPS[2],3),"gradient =",signif(x$EPS[3],3),"\n")
      
      cat("\n")
      cat("      n=", x$n)
      
      if (length(x$na.action)){
        cat("  (", length(x$na.action), " observation deleted due to missing) \n")
      }else{
        cat("\n")
      }
      
      if (!is.null(frail)) cat("      n events=", x$n.events, " n groups=", x$groups)
      else cat("      n events=", x$n.events)
      
      cat( "\n")
      cat("      number of iterations: ", x$n.iter)
    }
    cat("\n")
    invisible()	
}
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~     
  
  
  
    
    
# GENfrailtyPENAL (cas du probit, ie object$family==2)  
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
else if(x$family == 2){	
  if (x$istop == 1){
    # plot des coefficient dependant du temps
    if ((x$nvartimedep != 0) & (x$istop == 1)){
      par(mfrow=c(1,x$nvartimedep))
      for (i in 0:(x$nvartimedep-1)){
        matplot(x$BetaTpsMat[,1],x$BetaTpsMat[,(2:4)+4*i],col="blue",type="l",lty=c(1,2,2),xlab="t",ylab="beta(t)",main=x$Names.vardep[i+1],ylim=c(min(x$BetaTpsMat[,-1]),max(x$BetaTpsMat[,-1])))
      }
    }
  }
  
  if (!is.null(cl <- x$call)){
    cat("Call:\n")
    dput(cl)
    if ((x$type == "counting" | x$type == "intervaltronc") & (x$AG == FALSE)){
      cat("\n      left truncated structure used")
    }
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
  nvar <- length(x$coef) #+x$nvartimedep
  
  #	if (is.null(coef))
  #	{
  #		x$varH<-matrix(x$varH) 
  #		x$varHIH<-matrix(x$varHIH)
  #	}
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
    if ((x$typeof == 1) & (x$indic.nb.int == 1)) cat("  The maximum number of time intervals is 20","\n")
  }
  #AD
  
  if (x$logNormal == 0) frail <- x$theta
  else frail <- x$sigma2
  
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
        tmp <- cbind(coef, exp(coef), seH, seHIH, coef/seH, ifelse(signif(1 - pchisq((coef/seH)^2, 1), digits - 1) == 0, "< 1e-16", signif(1 - pchisq((coef/seH)^2, 1), digits - 1))) # ifelse pour eviter que R n'affiche pval = 0 lorsque celle-ci est tres petite
        if(x$global_chisq.test==1) tmpwald <- cbind(x$global_chisq,x$dof_chisq,ifelse(x$p.global_chisq == 0, "< 1e-16", x$p.global_chisq))
      }else{
        tmp <- cbind(coef, exp(coef), seH, coef/seH, ifelse(signif(1 - pchisq((coef/seH)^2, 1), digits - 1) == 0, "< 1e-16", signif(1 - pchisq((coef/seH)^2, 1), digits - 1)))
        if(x$global_chisq.test==1) tmpwald <- cbind(x$global_chisq,x$dof_chisq,ifelse(x$p.global_chisq == 0, "< 1e-16", x$p.global_chisq))
      }
      
      cat("\n")
      if (!is.null(frail)){
        if (x$logNormal == 0){
          cat("  Generalized Survival Model with Shared Gamma Frailty ", "\n")
          #cat("  Shared Gamma Frailty model parameter estimates ","\n")
        }else{
          cat("  Generalized Survival Model with Shared Log-Normal Frailty ", "\n")
          #cat("  Shared Log-Normal Frailty model parameter estimates ","\n")
        }
        if (x$typeof == 0){
          cat("  using a Penalized Likelihood on the hazard function","\n")
        }else{
          cat("  Parametrical approach with link   g() = -PHI^-1() ", "\n")
          cat("  S(t) = [ g^-1(eta) ]^frailty ", "\n")
          cat("  eta = shape.log(t) - shape.log(scale) + beta'X ", "\n")
          cat("  (Probit Frailty Model with a log-normal distribution) ", "\n\n")
          #cat("  using a Parametrical approach for the hazard function","\n")
        }
        if (x$nvartimedep != 0) cat("  and some time-dependant covariates","\n")
        if (x$n.strat>1) cat("  (Stratification structure used) :",x$n.strat,"strata \n")
      }else{
        if (x$typeof == 0){
          cat("  Cox proportional hazards model parameter estimates ","\n")
          cat("  using a Penalized Likelihood on the hazard function","\n")
        }else{
          cat("  Generalized Survival Model ","\n")
          cat("  Parametrical approach with link   g() = -PHI^-1() ", "\n")
          cat("  g(S(t)) = eta = shape.log(t) - shape.log(scale) + beta'X ", "\n")
          cat("  (Probit Model with a log-normal distribution) ", "\n\n")
          #cat("  Cox proportional hazards model parameter estimates ","\n")
          #cat("  using a Parametrical approach for the hazard function","\n")
        }
        if (x$nvartimedep != 0) cat("  and some time-dependant covariates","\n")
        if (x$n.strat>1) cat("  (Stratification structure used) :",x$n.strat,"strata \n")
      }
      
      if (x$typeof == 0){
        if(x$global_chisq.test==1){
          dimnames(tmpwald) <- list(x$names.factor,c("chisq", "df", "global p"))
          
        }
        dimnames(tmp) <- list(names(coef), c("coef", "exp(coef)",
                                             "SE coef (H)", "SE coef (HIH)", "z", "p"))
        
      }else{
        if(x$global_chisq.test==1){
          dimnames(tmpwald) <- list(x$names.factor,c("chisq", "df", "global p"))
          
        }
        dimnames(tmp) <- list(names(coef), c("coef", "exp(coef)",
                                             "SE coef (H)", "z", "p"))
        
      }
      cat("\n")
      #AL:
      if (nvar == 0){
        cat("No constant coefficients, only time-varying effects of the covariates \n")
      }else{
        prmatrix(tmp[, -2])
        if(x$global_chisq.test==1){
          cat("\n")
          prmatrix(tmpwald)
        }
      }
      #AL
      cat("\n")
    }
    
    if (!is.null(frail)) {
      #tetha <- x$theta
      #temp <- x$varTheta[1]
      seH <- sqrt(x$varTheta[1]) #sqrt(((2 * (frail^0.5))^2) * temp)
      #temp <- x$varTheta[2]
      seHIH <- sqrt(x$varTheta[2]) #sqrt(((2 * (frail^0.5))^2) * temp)
      
      
      #AD:
      if (x$noVar1 == 1){
        cat("\n")
        if (x$logNormal == 0){ cat("    Shared Gamma Frailty model: No covariates \n") }
        else { cat("    Shared Log-Normal Frailty model: No covariates \n") }
        cat("    -------------------------- \n")
        cat("\n")
      }
      #AD:
      # 			if (x$typeof == 0){
      # 				if (x$logNormal == 0){
      # 					cat("    Frailty parameter, Theta:", frail, "(SE (H):",
      # 					seH, ")", "(SE (HIH):", seHIH, ")", "\n")
      # 				}else{
      # 					cat("    Frailty parameter, Sigma Square:", frail, "(SE (H):",
      # 					seH, ")", "(SE (HIH):", seHIH, ")", "\n")
      # 				}
      # 			}else{
      # 				if (x$logNormal == 0){
      # 					cat("    Frailty parameter, Theta:", frail, "(SE (H):",
      # 					seH, ")", "\n")
      # 				}else{
      # 					cat("    Frailty parameter, Sigma Square:", frail, "(SE (H):",
      # 					seH, ")", "\n")
      # 				}
      # 			}
      
      if (x$logNormal == 0){
        cat("    Frailty parameter, Theta:", frail, "(SE (H):",
            seH, ")", "p =", ifelse(signif(1 - pnorm(frail/seH), digits - 1) == 0, "< 1e-16", signif(1 - pnorm(frail/seH), digits - 1)), "\n")
      }else{
        cat("    Frailty parameter, Sigma Square:", frail, "(SE (H):",
            seH, ")", "p =", ifelse(signif(1 - pnorm(frail/seH), digits - 1) == 0, "< 1e-16", signif(1 - pnorm(frail/seH), digits - 1)), "\n")
      }
      
    }
    
    cat(" \n")
    #AD:	
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
        if (x$n.strat == 1){
          cat("      Scale for the log-normal survival function is :",round(x$scale.weib[1],2),"\n")
          cat("      Shape for the log-normal survival function is :",round(x$shape.weib[1],2),"\n")
        }else{
          cat("      Scale for the log-normal survival function is :",round(x$scale.weib[1],2),round(x$scale.weib[2],2),"\n")
          cat("      Shape for the log-normal survival function is :",round(x$shape.weib[1],2),round(x$shape.weib[2],2),"\n")
        }
        cat("\n")
        cat("The expression of the log-normal survival function is:","\n")
        cat("        'S(t) = PHI(-eta)'", "\n")
        cat("The expression of the log-normal hazard function is:","\n")
        cat("        'lambda(t) = phi(-eta)/PHI(-eta) * d.eta/d.t' ","\n")
        cat("\n")
      }
    }
    #AD:
    cat("\n")
    cat("      n=", x$n)
    
    if (length(x$na.action)){
      cat("  (", length(x$na.action), " observation deleted due to missing) \n")
    }else{
      cat("\n")
    }
    
    if (!is.null(frail)) cat("      n events=", x$n.events, " n groups=", x$groups)
    else cat("      n events=", x$n.events)
    
    cat( "\n")
    cat("      number of iterations: ", x$n.iter,"\n")
    if ((x$typeof == 1) & (x$indic.nb.int == 1)){
      cat("      Exact number of time intervals used: 20","\n")
    }else{
      if (x$typeof == 1) cat("      Exact number of time intervals used: ",x$nbintervR,"\n")
    }	
    if (x$typeof == 0){ 
      cat("\n")
      cat("      Exact number of knots used: ", x$n.knots, "\n")
      
      if (!x$cross.Val){
        cat("      Value of the smoothing parameter: ", x$kappa, sep=" ")
      }
      
      if (x$cross.Val){
        if (is.null(frail)){
          cat("      Smoothing parameter estimated by Cross validation: ", x$kappa, sep=" ")
        }else{ 
          cat("      Best smoothing parameter estimated by")
          cat("\n")
          cat("      an approximated Cross validation: ", x$kappa, sep=" ")
        }
      }
      cat(", DoF: ", formatC(-x$DoF, format="f",digits=2))
    }
  }else{
    if (!is.null(frail)){
      if (x$logNormal == 0){
        cat("  Shared Gamma Frailty model parameter estimates ","\n")
      }else{
        cat("  Shared Log-Normal Frailty model parameter estimates ","\n")
      }
      if (x$typeof == 0){
        cat("  using a Penalized Likelihood on the hazard function","\n")
      }else{
        cat("  using a Parametrical approach for the hazard function","\n")
      }
      if (x$nvartimedep != 0) cat("  and some time-dependant covariates","\n")
      if (x$n.strat>1) cat("  (Stratification structure used)", "\n")
      if (x$noVar1 == 1){
        cat("\n")
        if (x$logNormal == 0){ cat("    Shared Gamma Frailty model: No covariates \n") }
        else { cat("    Shared Log-Normal Frailty model: No covariates \n") }
        cat("    -------------------------- \n")
        cat("\n")
      }
      
    }else{
      if (x$typeof == 0){
        cat("  Cox proportional hazards model parameter estimates ","\n")
        cat("  using a Penalized Likelihood on the hazard function","\n")
      }else{
        cat("  Cox proportional hazards model parameter estimates ","\n")
        cat("  using a Parametrical approach for the hazard function","\n")
      }
      if (x$nvartimedep != 0) cat("  and some time-dependant covariates","\n")
      if (x$n.strat>1) cat("  (Stratification structure used)", "\n")
    }
    cat("\n")
    
    cat("      Convergence criteria: \n")
    cat("      parameters =",signif(x$EPS[1],3),"likelihood =",signif(x$EPS[2],3),"gradient =",signif(x$EPS[3],3),"\n")
    
    cat("\n")
    cat("      n=", x$n)
    
    if (length(x$na.action)){
      cat("  (", length(x$na.action), " observation deleted due to missing) \n")
    }else{
      cat("\n")
    }
    
    if (!is.null(frail)) cat("      n events=", x$n.events, " n groups=", x$groups)
    else cat("      n events=", x$n.events)
    
    cat( "\n")
    cat("      number of iterations: ", x$n.iter)
  }
  cat("\n")
  invisible()	
}
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~    

  
  
  
# GENfrailtyPENAL (cas du modele additif, ie object$family==3)  
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
else if(x$family %in% c(3,4)){	
    if (x$istop == 1){
      # plot des coefficient dependant du temps
      if ((x$nvartimedep != 0) & (x$istop == 1)){
        
        trapz <- function(x,y){
          idx = 2:length(x)
          return (as.double( (x[idx] - x[idx-1]) %*% (y[idx] + y[idx-1])) / 2)
        }
        par(mfrow=c(x$nvartimedep, 2))
        
        for (i in 0:(x$nvartimedep-1)){
          matplot(x$BetaTpsMat[,1],
                  x$BetaTpsMat[,(2:4)+4*i],
                  col="blue",type="l",lty=c(1,2,2),
                  xlab="t",ylab="beta(t)",
                  main=x$Names.vardep[i+1],
                  ...
                  )
          
          nblignes = nrow(x$BetaTpsMat)-1
          matcumul = matrix(NA, nrow = nblignes, ncol = 3)
          abs = x$BetaTpsMat[, 1]
          ord = x$BetaTpsMat[, (2:4)+4*i] #exp(x$BetaTpsMat[, (2:4)+4*i]) 
          #ordinf = x$BetaTpsMat[,3]
          #ordsup = x$BetaTpsMat[,4]
          for(j in 1:nblignes){
            matcumul[j, ] = c(
              trapz(abs[1:(j+1)], ord[1:(j+1), 1]),
              trapz(abs[1:(j+1)], ord[1:(j+1), 2]),
              trapz(abs[1:(j+1)], ord[1:(j+1), 3])
            )
          }
          
          matplot(abs[-1],
                  matcumul[,(1:3)],
                  col="blue",type="l",lty=c(1,2,2),
                  xlab="t",ylab="Cumulative effect",
                  main=x$Names.vardep,
                  ...
                  )
          
        }
      }
    }
    
    if (!is.null(cl <- x$call)){
      cat("Call:\n")
      dput(cl)
      if ((x$type == "counting" | x$type == "intervaltronc") & (x$AG == FALSE)){
        cat("\n      left truncated structure used")
      }
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
    nvar <- length(x$coef) #+x$nvartimedep
    
    #	if (is.null(coef))
    #	{
    #		x$varH<-matrix(x$varH) 
    #		x$varHIH<-matrix(x$varHIH)
    #	}
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
      if ((x$typeof == 1) & (x$indic.nb.int == 1)) cat("  The maximum number of time intervals is 20","\n")
    }
    #AD
    
    if (x$logNormal == 0) frail <- x$theta
    else frail <- x$sigma2
    
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
          tmp <- cbind(coef, exp(coef), seH, seHIH, coef/seH, ifelse(signif(1 - pchisq((coef/seH)^2, 1), digits - 1) == 0, "< 1e-16", signif(1 - pchisq((coef/seH)^2, 1), digits - 1))) # ifelse pour eviter que R n'affiche pval = 0 lorsque celle-ci est tres petite
          if(x$global_chisq.test==1) tmpwald <- cbind(x$global_chisq,x$dof_chisq,ifelse(x$p.global_chisq == 0, "< 1e-16", x$p.global_chisq))
        }else{
          tmp <- cbind(coef, exp(coef), seH, coef/seH, ifelse(signif(1 - pchisq((coef/seH)^2, 1), digits - 1) == 0, "< 1e-16", signif(1 - pchisq((coef/seH)^2, 1), digits - 1)))
          if(x$global_chisq.test==1) tmpwald <- cbind(x$global_chisq,x$dof_chisq,ifelse(x$p.global_chisq == 0, "< 1e-16", x$p.global_chisq))
        }
        
        cat("\n")
        if (!is.null(frail)){
          if (x$logNormal == 0){
            cat("  Generalized Survival Model with Shared Gamma Frailty ", "\n")
            #cat("  Shared Gamma Frailty model parameter estimates ","\n")
          }else{
            cat("  Generalized Survival Model with Shared Log-Normal Frailty ", "\n")
            #cat("  Shared Log-Normal Frailty model parameter estimates ","\n")
          }
          if (x$typeof == 0){
            cat("  using a Penalized Likelihood on the hazard function","\n")
          }else{
            cat("  Parametrical approach with link   g() = -log() ", "\n")
            cat("  S(t) = [ g^-1(eta) ]^frailty ", "\n")
            cat("  eta = (t/scale)^shape + t*beta'X", "\n")
            cat("  (Additive Hazards Frailty Model with a Weibull distribution) ", "\n\n")
            #cat("  using a Parametrical approach for the hazard function","\n")
          }
          if (x$nvartimedep != 0) cat("  and some time-dependant covariates","\n")
          if (x$n.strat>1) cat("  (Stratification structure used) :",x$n.strat,"strata \n")
        }else{
          if (x$typeof == 0){
            cat("  Cox proportional hazards model parameter estimates ","\n")
            cat("  using a Penalized Likelihood on the hazard function","\n")
          }else{
            cat("  Generalized Survival Model ","\n")
            cat("  Parametrical approach with link   g() = -log() ", "\n")
            cat("  g(S(t)) = eta =(t/scale)^shape + t*beta'X ", "\n")
            cat("  (Additive Hazards Model with a Weibull distribution) ", "\n\n")
            #cat("  Cox proportional hazards model parameter estimates ","\n")
            #cat("  using a Parametrical approach for the hazard function","\n")
          }
          if (x$nvartimedep != 0) cat("  and some time-dependant covariates","\n")
          if (x$n.strat>1) cat("  (Stratification structure used) :",x$n.strat,"strata \n")
        }
        
        if (x$typeof == 0){
          if(x$global_chisq.test==1){
            dimnames(tmpwald) <- list(x$names.factor,c("chisq", "df", "global p"))
            
          }
          dimnames(tmp) <- list(names(coef), c("coef", "exp(coef)",
                                               "SE coef (H)", "SE coef (HIH)", "z", "p"))
          
        }else{
          if(x$global_chisq.test==1){
            dimnames(tmpwald) <- list(x$names.factor,c("chisq", "df", "global p"))
            
          }
          dimnames(tmp) <- list(names(coef), c("coef", "exp(coef)",
                                               "SE coef (H)", "z", "p"))
          
        }
        cat("\n")
        #AL:
        if (nvar == 0){
          cat("No constant coefficients, only time-varying effects of the covariates \n")
        }else{
          prmatrix(tmp[, -2])
          if(x$global_chisq.test==1){
            cat("\n")
            prmatrix(tmpwald)
          }
        }
        #AL
        cat("\n")
      }
      
      if (!is.null(frail)) {
        #tetha <- x$theta
        #temp <- x$varTheta[1]
        seH <- sqrt(x$varTheta[1]) #sqrt(((2 * (frail^0.5))^2) * temp)
        #temp <- x$varTheta[2]
        seHIH <- sqrt(x$varTheta[2]) #sqrt(((2 * (frail^0.5))^2) * temp)
        
        
        #AD:
        if (x$noVar1 == 1){
          cat("\n")
          if (x$logNormal == 0){ cat("    Shared Gamma Frailty model: No covariates \n") }
          else { cat("    Shared Log-Normal Frailty model: No covariates \n") }
          cat("    -------------------------- \n")
          cat("\n")
        }
        #AD:
        # 			if (x$typeof == 0){
        # 				if (x$logNormal == 0){
        # 					cat("    Frailty parameter, Theta:", frail, "(SE (H):",
        # 					seH, ")", "(SE (HIH):", seHIH, ")", "\n")
        # 				}else{
        # 					cat("    Frailty parameter, Sigma Square:", frail, "(SE (H):",
        # 					seH, ")", "(SE (HIH):", seHIH, ")", "\n")
        # 				}
        # 			}else{
        # 				if (x$logNormal == 0){
        # 					cat("    Frailty parameter, Theta:", frail, "(SE (H):",
        # 					seH, ")", "\n")
        # 				}else{
        # 					cat("    Frailty parameter, Sigma Square:", frail, "(SE (H):",
        # 					seH, ")", "\n")
        # 				}
        # 			}
        
        if (x$logNormal == 0){
          cat("    Frailty parameter, Theta:", frail, "(SE (H):",
              seH, ")", "p =", ifelse(signif(1 - pnorm(frail/seH), digits - 1) == 0, "< 1e-16", signif(1 - pnorm(frail/seH), digits - 1)), "\n")
        }else{
          cat("    Frailty parameter, Sigma Square:", frail, "(SE (H):",
              seH, ")", "p =", ifelse(signif(1 - pnorm(frail/seH), digits - 1) == 0, "< 1e-16", signif(1 - pnorm(frail/seH), digits - 1)), "\n")
        }
        
      }
      
      cat(" \n")
      #AD:	
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
          if (x$n.strat == 1){
            cat("      Scale for the Weibull survival function is :",round(x$scale.weib[1],2),"\n")
            cat("      Shape for the Weibull survival function is :",round(x$shape.weib[1],2),"\n")
          }else{
            cat("      Scale for the Weibull survival function is :",round(x$scale.weib[1],2),round(x$scale.weib[2],2),"\n")
            cat("      Shape for the Weibull survival function is :",round(x$shape.weib[1],2),round(x$shape.weib[2],2),"\n")
          }
          cat("\n")
          cat("The expression of the Weibull survival function is:","\n")
          cat("        'S(t) = exp[- (t/scale)^shape]'", "\n")
          cat("The expression of the log-normal hazard function is:","\n")
          cat("        'lambda(t) = lambda(t) = (shape.(t^(shape-1)))/(scale^shape)' ","\n")
          cat("\n")
        }
      }
      #AD:
      cat("\n")
      cat("      n=", x$n)
      
      if (length(x$na.action)){
        cat("  (", length(x$na.action), " observation deleted due to missing) \n")
      }else{
        cat("\n")
      }
      
      if (!is.null(frail)) cat("      n events=", x$n.events, " n groups=", x$groups)
      else cat("      n events=", x$n.events)
      
      cat( "\n")
      cat("      number of iterations: ", x$n.iter,"\n")
      if ((x$typeof == 1) & (x$indic.nb.int == 1)){
        cat("      Exact number of time intervals used: 20","\n")
      }else{
        if (x$typeof == 1) cat("      Exact number of time intervals used: ",x$nbintervR,"\n")
      }	
      if (x$typeof == 0){ 
        cat("\n")
        cat("      Exact number of knots used: ", x$n.knots, "\n")
        
        if (!x$cross.Val){
          cat("      Value of the smoothing parameter: ", x$kappa, sep=" ")
        }
        
        if (x$cross.Val){
          if (is.null(frail)){
            cat("      Smoothing parameter estimated by Cross validation: ", x$kappa, sep=" ")
          }else{ 
            cat("      Best smoothing parameter estimated by")
            cat("\n")
            cat("      an approximated Cross validation: ", x$kappa, sep=" ")
          }
        }
        cat(", DoF: ", formatC(-x$DoF, format="f",digits=2))
      }
    }else{
      if (!is.null(frail)){
        if (x$logNormal == 0){
          cat("  Shared Gamma Frailty model parameter estimates ","\n")
        }else{
          cat("  Shared Log-Normal Frailty model parameter estimates ","\n")
        }
        if (x$typeof == 0){
          cat("  using a Penalized Likelihood on the hazard function","\n")
        }else{
          cat("  using a Parametrical approach for the hazard function","\n")
        }
        if (x$nvartimedep != 0) cat("  and some time-dependant covariates","\n")
        if (x$n.strat>1) cat("  (Stratification structure used)", "\n")
        if (x$noVar1 == 1){
          cat("\n")
          if (x$logNormal == 0){ cat("    Shared Gamma Frailty model: No covariates \n") }
          else { cat("    Shared Log-Normal Frailty model: No covariates \n") }
          cat("    -------------------------- \n")
          cat("\n")
        }
        
      }else{
        if (x$typeof == 0){
          cat("  Cox proportional hazards model parameter estimates ","\n")
          cat("  using a Penalized Likelihood on the hazard function","\n")
        }else{
          cat("  Cox proportional hazards model parameter estimates ","\n")
          cat("  using a Parametrical approach for the hazard function","\n")
        }
        if (x$nvartimedep != 0) cat("  and some time-dependant covariates","\n")
        if (x$n.strat>1) cat("  (Stratification structure used)", "\n")
      }
      cat("\n")
      
      cat("      Convergence criteria: \n")
      cat("      parameters =",signif(x$EPS[1],3),"likelihood =",signif(x$EPS[2],3),"gradient =",signif(x$EPS[3],3),"\n")
      
      cat("\n")
      cat("      n=", x$n)
      
      if (length(x$na.action)){
        cat("  (", length(x$na.action), " observation deleted due to missing) \n")
      }else{
        cat("\n")
      }
      
      if (!is.null(frail)) cat("      n events=", x$n.events, " n groups=", x$groups)
      else cat("      n events=", x$n.events)
      
      cat( "\n")
      cat("      number of iterations: ", x$n.iter)
    }
    cat("\n")
    invisible()	
  }
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~      
  
  
}

