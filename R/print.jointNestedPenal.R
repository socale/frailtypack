#' Print a Short Summary of parameter estimates of a joint nested frailty model
#' 
#' Prints a short summary of parameter estimates of a joint nested frailty
#' model, or more generally an object of class 'jointNestedPenal' for joint
#' nested frailty models.
#' 
#' 
#' @usage \method{print}{jointNestedPenal}(x, digits = max(options()$digits - 4,
#' 6), ...)
#' @param x the result of a call to the jointNestedPenal function
#' @param digits number of digits to print
#' @param \dots other unused arguments
#' @return
#' 
#' Print, separately for each type of event (recurrent and terminal), the
#' parameter estimates of the survival or hazard functions.
#' @seealso \code{\link{frailtyPenal}}
#' @keywords methods
##' @export
"print.jointNestedPenal" <- function (x, digits = max(options()$digits - 4, 6), ...){  
	if (!is.null(cl <- x$call)){
		cat("Call:\n")
		dput(cl)
		if (x$AG == TRUE) cat("\n      Calendar timescale")
		#if (x$intcens == TRUE) cat("\n      interval censored data used")
	}
    cat("\n")
	savedig <- options(digits = digits)
	on.exit(options(savedig))
	coef <- x$coef
	nvar <- sum(x$nvar)  
	if (is.null(coef)){ 
		x$varH<-matrix(x$varH)
		x$varHIH<-matrix(x$varHIH)
	}
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
  
	frail1 <- x$theta
	frail2 <- x$eta
	indic_alpha <- x$indic_alpha
	indic_xi <- x$indic_ksi
  
	if (x$istop == 1){
		if (!is.null(coef)){
			if (indic_alpha == 1 & indic_xi == 1){
				seH <- sqrt(diag(x$varH))[-c(1:4)]
				seHIH <- sqrt(diag(x$varHIH))[-c(1:4)]
			}
			if (indic_alpha == 0 | indic_xi == 0) {
				seH <- sqrt(diag(x$varH))[-c(1:3)]
				seHIH <- sqrt(diag(x$varHIH))[-c(1:3)]
			}
			if(indic_alpha == 0 & indic_xi == 0) {
				seH <- sqrt(diag(x$varH))[-c(1:2)]
				seHIH <- sqrt(diag(x$varHIH))[-c(1:2)]
			}
		  if (x$typeof == 0){
		    tmp <- cbind(coef, exp(coef), seH, seHIH, coef/seH, ifelse(signif(1 - pchisq((coef/seH)^2, 1), digits - 1) == 0, "< 1e-16", signif(1 - pchisq((coef/seH)^2, 1), digits - 1)))
		    if(x$global_chisq.test==1) tmpwald <- cbind(x$global_chisq, x$dof_chisq, ifelse(x$p.global_chisq == 0, "< 1e-16", x$p.global_chisq))
		    if(x$global_chisq.test_d==1) tmpwalddc <- cbind(x$global_chisq_d, x$dof_chisq_d, ifelse(x$p.global_chisq_d == 0, "< 1e-16", x$p.global_chisq_d))
		  }
		  else{
		    tmp <- cbind(coef, exp(coef), seH, coef/seH, ifelse(signif(1 - pchisq((coef/seH)^2, 1), digits - 1) == 0, "< 1e-16", signif(1 - pchisq((coef/seH)^2, 1), digits - 1)))
		    if(x$global_chisq.test==1) tmpwald <- cbind(x$global_chisq, x$dof_chisq, ifelse(x$p.global_chisq == 0, "< 1e-16", x$p.global_chisq))
		    if(x$global_chisq.test_d==1) tmpwalddc <- cbind(x$global_chisq_d, x$dof_chisq_d, ifelse(x$p.global_chisq_d == 0, "< 1e-16", x$p.global_chisq_d))
		  }
			cat("\n")
			
			cat("  Joint nested gamma frailty model for recurrent and a terminal event processes","\n")
						
			if (x$typeof == 0) cat("  using a Penalized Likelihood on the hazard function","\n")
			else cat("  using a Parametrical approach for the hazard function","\n")
			if (x$n.strat>1) cat("  (Stratification structure used for recurrences) :",x$n.strat,"strata \n")
			
			if (x$typeof == 0){
				if(x$global_chisq.test==1) dimnames(tmpwald) <- list(x$names.factor,c("chisq", "df", "global p"))
				if(x$global_chisq.test_d==1) dimnames(tmpwalddc) <- list(x$names.factordc,c("chisq", "df", "global p"))
                dimnames(tmp) <- list(names(coef), c("coef", "exp(coef)","SE coef (H)", "SE coef (HIH)", "z", "p"))
			}else{
				if(x$global_chisq.test==1) dimnames(tmpwald) <- list(x$names.factor,c("chisq", "df", "global p"))          
				if(x$global_chisq.test_d==1) dimnames(tmpwalddc) <- list(x$names.factordc,c("chisq", "df", "global p"))          
				dimnames(tmp) <- list(names(coef), c("coef", "exp(coef)","SE coef (H)", "z", "p")) 
			}
			cat("\n")      
			
			if (x$noVar1 == 0){
				cat("Recurrences:\n")
				cat("------------- \n")
				prmatrix(tmp[1:x$nvar[1], ,drop=FALSE])
				if(x$global_chisq.test==1){
					cat("\n")
					prmatrix(tmpwald)
				}
			}
			
			cat("\n")      
			if (x$noVar2 == 0){
				cat("Terminal event:\n")
				cat("---------------- \n")
				prmatrix(tmp[-c(1:x$nvar[1]), ,drop=FALSE]) 
				if(x$global_chisq.test_d==1){
					cat("\n")
					prmatrix(tmpwalddc)
				}
			}
			cat("\n")
		}	
		temp1 <- diag(x$varH)[1]
		temp2 <- diag(x$varH)[2]
		seH.frail1 <- sqrt(((2 * (frail1^0.5))^2) * temp1) # delta methode 
		seH.frail2 <- sqrt(((2 * (frail2^0.5))^2) * temp2)
		temp1 <- diag(x$varHIH)[1]
		temp2 <- diag(x$varHIH)[2]
		seHIH.frail1 <- sqrt(((2 * (frail1^0.5))^2) * temp1) # delta methode
		seHIH.frail1 <- sqrt(((2 * (frail2^0.5))^2) * temp2)
		if (x$noVar1 == 1){
			cat("\n")
			cat("    Recurrences: No covariates \n")
			cat("    ----------- \n")
		}		
		if (x$noVar2 == 1){
			cat("\n")
			cat("    Terminal event: No covariates \n")
			cat("    -------------- \n")
			cat("\n")
		}
	
		cat(" Frailty parameters: \n")
		cat("   theta (variance of Frailties, u):", frail1, "(SE (H):",seH.frail1, ")", "p =", ifelse(signif(1 - pnorm(frail1/seH.frail1), digits - 1) == 0, "< 1e-16", signif(1 - pnorm(frail1/seH.frail1), digits - 1)), "\n")
		cat("	eta (variance of Frailties, w) :", frail2, "(SE (H) :", seH.frail2,")", "p =", ifelse(signif(1 - pnorm(frail2/seH.frail2), digits - 1) == 0, "< 1e-16", signif(1 - pnorm(frail2/seH.frail2), digits - 1)), "\n")
		if (indic_alpha == 1) cat("   alpha (u^alpha for terminal event):", x$alpha, "(SE (H):",sqrt(diag(x$varH))[3], ")", "p =", ifelse(signif(1 - pchisq((x$alpha/sqrt(diag(x$varH))[3])^2,1), digits - 1) == 0,"< 1e-16", signif(1 - pchisq((x$alpha/sqrt(diag(x$varH))[3])^2,1), digits - 1)), "\n")
		if (indic_xi == 1 & indic_alpha == 1) cat("   ksi (w^ksi for recurrent event):", x$ksi, "(SE (H):",sqrt(diag(x$varH))[4], ")", "p =", ifelse(signif(1 - pchisq((x$ksi/sqrt(diag(x$varH))[4])^2,1), digits - 1) == 0,"< 1e-16",signif(1 - pchisq((x$ksi/sqrt(diag(x$varH))[4])^2,1), digits - 1)), "\n")
		if (indic_xi == 1 & indic_alpha == 0) cat("   ksi (w^ksi for recurrent event):", x$ksi, "(SE (H):",sqrt(diag(x$varH))[3], ")", "p =", ifelse(signif(1 - pchisq((x$ksi/sqrt(diag(x$varH))[3])^2,1), digits - 1) == 0,"< 1e-16", signif(1 - pchisq((x$ksi/sqrt(diag(x$varH))[3])^2,1), digits - 1)), "\n")
		cat(" \n")
    
		if (x$typeof == 0){
			cat(paste("   penalized marginal log-likelihood =", round(x$logLikPenal,2)))
			cat("\n")
			cat("   Convergence criteria: \n")
			cat("   parameters =",signif(x$EPS[1],3),"likelihood =",signif(x$EPS[2],3),"gradient =",signif(x$EPS[3],3),"\n")
			cat("\n")
			cat("   LCV = the approximate likelihood cross-validation criterion\n")
			cat("         in the semi parametric case     =",x$LCV,"\n")
		}
		else{ 
			cat(paste("   marginal log-likelihood =", round(x$logLik,2)))
			cat("\n")
			cat("   Convergence criteria: \n")
			cat("   parameters =",signif(x$EPS[1],3),"likelihood =",signif(x$EPS[2],3),"gradient =",signif(x$EPS[3],3),"\n")
			cat("\n")
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
		cat("\n")
		cat("   n observations=", x$n, " n subjects=", x$subgroups, " n groups=", x$groups)   
		
		if (length(x$na.action)) cat("      (", length(x$na.action), " observation deleted due to missing) \n")
		else cat("\n")
		
		cat("   n recurrent events=", x$n.events)
		cat("\n")
		cat("   n terminal events=", x$n.death)
		cat("\n")
		cat("   n censored events=" ,x$n.censored)
		cat("\n")
		cat("   number of iterations: ", x$n.iter,"\n")
		if (x$logNormal == 0) {
		  cat("   Number of nodes for the Gauss-Laguerre quadrature: ", x$nb.gl,"\n")
		}
    else {cat("   Number of nodes for the Gauss-Hermite quadrature: ", x$nb.gh,"\n")}
		
		if (x$typeof == 0){ # splines
			cat("\n")
			cat("   Exact number of knots used: ", x$n.knots, "\n")
			cat("   Value of the smoothing parameters: ", x$kappa, sep=" ")
			cat("\n")
		}
	}else{ # s'il y a eu un pb de convergence et istop != 1
		if (!is.null(coef)){ 
			cat("\n")				
			cat("  Joint nested gamma frailty model for recurrent and a terminal event processes","\n")
			
			if (x$typeof == 0) cat("  using a Penalized Likelihood on the hazard function","\n")
			else cat("  using a Parametrical approach for the hazard function","\n")
			
			if (x$noVar1 == 1){
				cat("\n")
				cat("    Recurrences: No covariates \n")
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
			cat("   n observations=", x$n, " n subjects=", x$subgroups, " n families=", x$groups)
			if (length(x$na.action)) cat("      (", length(x$na.action), " observation deleted due to missing) \n")
			else cat("\n")
			
			cat("   n recurrent events=", x$n.events)
			cat("\n")
			cat("   n terminal events=", x$n.death)
			cat("\n")			
			cat("   n censored events=" ,x$n.censored)
			cat("\n")
			cat("   number of iterations: ", x$n.iter,"\n")
			if (x$logNormal == 0) {
			  cat("   Number of nodes for the Gauss-Laguerre quadrature: ", x$nb.gl,"\n")
			}
			else {cat("   Number of nodes for the Gauss-Hermite quadrature: ", x$nb.gh,"\n")}
		}
	}
	invisible()
}
