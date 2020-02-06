##' Short summary of the random effects parameters, the fixed treatment 
##' effects, and the surrogacy evaluation criteria estimated from a joint surrogate model
##' 
##' This function returns the estimate of the coefficients and their standard error with p-values 
##' of the Wald test for the joint surrogate model, also hazard ratios (HR) and their 
##' confidence intervals for the fixed treatment effects, and finaly an estimate of the 
##' surrogacy evaluation criterian (Kendall's \eqn{\tau} and \if{latex}{\eqn{R^2_{trial}}}
#'    \if{html}{\code{R}\out{<sup>2</sup><sub>trial</sub>}})
##' 
##' 
##' @aliases summary.jointSurroPenal print.summary.jointSurroPenal
##' @usage \method{summary}{jointSurroPenal}(object, d = 4, len = 3, int.method.kt = 0, 
##' nb.gh = 32, ...)
##' 
##' @param object An object inheriting from \code{jointSurroPenal} class.
##' @param d The desired number of digits after the decimal point for parameters. 
##' The maximum of 4 digits is required for the estimates. Default of 3 digits is used.
##' @param len The desired number of digits after the decimal point for p-value and convergence 
##' criteria. Default of 4 digits is used.
##' @param int.method.kt A binary, indicates the integration method for Kendall's \eqn{\tau}
##' estimation : \code{0} for Monte carlo, and \code{1} for Gaussian Hermite quadrature. 
##' the default is \code{0}. 
##' @param  nb.gh Number of nodes for the Gaussian-Hermite quadrature.  The default is \code{32}
##' \code{1} for Gaussian-Hermite quadrature.
##' @param \dots other unused arguments.
##' 
##' @return For the variances parameters of the random effects, it prints the estimate of
##' the coefficients with their standard error, Z-statistics and p-values
##' of the Wald test. For the fixed treatment effects, it also prints HR and its confidence
##' intervals for each covariate. For the surrogacy evaluation criteria, its prints the estimated 
##' Kendall's \eqn{\tau} with its 95\% Confidence interval obtained by the parametric bootstrap
##'  or Delta-method, 
##' the estimated \if{latex}{\eqn{R^2_{trial}}}
#'    \if{html}{\code{R}\out{<sup>2</sup><sub>trial</sub>}}(R2trial) with standard error and the 95\% Confidence interval 
##' obtained by Delta-method (Dowd \emph{et al.}, 2014), \if{latex}{\eqn{R^2_{trial}}}
#'    \if{html}{\code{R}\out{<sup>2</sup><sub>trial</sub>}}(R2.boot) and its 95\% 
##' Confidence interval obtained by the parametric bootstrap. 
##' We notice that, using bootstrap, 
##' the standard error of the point estimate is not available. We propose a classification of \if{latex}{\eqn{R^2_{trial}}}
#'    \if{html}{\code{R}\out{<sup>2</sup><sub>trial</sub>}} according to 
##' the suggestion of the Institute of Quality and Efficiency in Health Care 
##' (Prasad \emph{et al.}, 2015). 
##' We also display the surrogate threshold effect (\code{\link[=ste]{ste}}) with the associated hazard risk.
##' The rest of parameters concerns the convergence characteristics and 
##' included: the penalized marginal log-likelihood, the number of iterations, the LCV and the Convergence criteria.
##' @seealso \code{\link{jointSurroPenal}, \link{jointSurroCopPenal}, \link{jointSurroTKendall}}
##' 
##' @author Casimir Ledoux Sofeu \email{casimir.sofeu@u-bordeaux.fr}, \email{scl.ledoux@gmail.com} and 
##' Virginie Rondeau \email{virginie.rondeau@inserm.fr}
##' 
##' @keywords methods
##' 
##' @references 
##' Dowd BE, Greene WH, Norton EC (2014). "Computation of Standard Errors." Health Services
##' Research, 49(2), 731-750.
##' 
##' Prasad V, Kim C, Burotto M, Vandross A (2015). "The strength of association between
##' surrogate end points and survival in oncology: A systematic review of trial-level meta-
##' alyses." JAMA Internal Medicine, 175(8), 1389-1398.
##' @export
##' @importFrom stats sd
##' @examples
##' 
##' 
##' 
##' ###---Data generation---###
##' data.sim <-jointSurrSimul(n.obs=400, n.trial = 20,cens.adm=549, 
##'           alpha = 1.5, theta = 3.5, gamma = 2.5, zeta = 1, 
##'           sigma.s = 0.7, sigma.t = 0.7, cor = 0.8, betas = -1.25, 
##'           betat = -1.25, full.data = 0, random.generator = 1, 
##'           seed = 0, nb.reject.data = 0)
##' \dontrun{
##' ###---Estimation---###
##' joint.surrogate <- jointSurroPenal(data = data.sim, nb.mc = 300, 
##'                    nb.gh = 20, indicator.alpha = 1, n.knots = 6)
##'                             
##' summary(joint.surrogate)
##' summary(joint.surrogate, d = 4, len = 3, int.method.kt = 1, nb.gh = 25)
##' }
##' 
##' 
"summary.jointSurroPenal"<-
  function(object, d = 4, len = 3, int.method.kt = 0, nb.gh = 32, ...){
    x <- object
    if (!inherits(x, "jointSurroPenal"))
      stop("Object must be of class 'jointSurroPenal'")
    
    coef <- data.frame(x$Coefficients)
    
    beta <- coef[(1 : (nrow(coef)-2)),1:2]
    names(beta)[2] <- "SE"
    beta$"z" <- round(beta$Estimate/beta$SE,len)
    beta$P <- signif(1 - pchisq((beta$Estimate/beta$SE)^2, 1), 5)
    beta$" " <- ifelse(beta$P < 0.001,"***",ifelse(beta$P < 0.01,"**",
                      ifelse(beta$P < 0.05,"*",ifelse(beta$P < 0.1,"."," "))))
    beta$P <- formatC(beta$P, d, format = "g")
    beta$Estimate <- round(beta$Estimate,min(4,len))
    names(beta)[2] <-"Std Error"
    
    # travail des P
    p <- NULL
    p <- ifelse(as.numeric(beta$P) < 10^-10, "< e-10", beta$P)
    beta$P <- p
    
    cat("Estimates for variances parameters of the random effects", "\n")
    rownames(beta)[(nrow(beta) - 4)] <- "sigma2_S"
    rownames(beta)[(nrow(beta) - 3)] <- "sigma2_T"
    rownames(beta)[(nrow(beta) - 2)] <- "sigma_ST"
    print(beta[1:(nrow(beta) - 2),])
    
    beta2 <- beta[((nrow(beta) - 1) : nrow(beta)),]
    beta2$Estimate <- round(beta2$Estimate,min(4,len))
   
    cat(" ", "\n")
    cat("Estimates for the fixed treatment effects", "\n")
    print(beta2)
    
    cat("---","\n")
    cat("Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1")
    
    cat(" ", "\n")
    cat(" ","\n")
    
    cat("hazard ratios (HR) and confidence intervals for the fixed treatment effects", "\n")
    HR <- round(exp(coef[((nrow(coef) - 3) : (nrow(coef)-2)),-2]), len)
    names(HR)[1] <- c("exp(coef)")
    print(HR)
    
    # 95%CI of R2 using Delta-method
    coef <- rbind(coef,coef[nrow(coef)-1,])
    coef[nrow(coef),c(3,4)] <- object$R2.boot[-1]
    coef[nrow(coef),1] <- object$R2.boot[1]
    coef[nrow(coef),2] <- NA
    rownames(coef)[nrow(coef)] <- "R2.boot"
    
    
    validation <- coef[c(nrow(coef) - 1, nrow(coef) - 2,nrow(coef)),]
    validation[,2] <- as.character(validation[,2])
    if(x$type.joint == 1) {
      validation[1,2] <- "--"
    }else{
      validation[1,2] <- as.character(round(as.numeric(validation[1,2]), len))
    }
    validation[3,2] <- "--"
    validation[,1] <- round(validation[,1], len)
    validation[,3] <- round(validation[,3], len)
    validation[,4] <- round(validation[,4], len)
    validation[2,2] <- as.character(round(as.numeric(validation[2,2]), len))
    if(int.method.kt == 1){
      validation[1,1] <- jointSurroTKendall(theta = object$Coefficients["Theta",1],
                                            gamma = ifelse(is.na(object$Coefficients["gamma",1]), 0, object$Coefficients["gamma",1]),
                                            alpha = ifelse(is.na(object$Coefficients["alpha",1]), 1, object$Coefficients["alpha",1]),
                                            zeta = ifelse(is.na(object$Coefficients["zeta",1]), 1, object$Coefficients["zeta",1]),
                                            nb.gh, ui = ifelse(is.na(object$Coefficients["gamma",1]), 0, 1))
    }
    
    names(validation)[2] <- "Std Error"
    
    # More precision on the surrogacy evaluation
    validation2 <- data.frame(matrix(rep(NA,18), nrow = 3, ncol = 6))
    names(validation2) <- c("Level", names(validation), "Strength")
    rownames(validation2) <- rownames(validation)
    validation2[,1] <- c("Individual", "Trial", "Trial")
    validation2[,2:5] <- validation
    validation2[2,6] <- ifelse(validation2[2,4] <= 0.49,"Low",ifelse(validation2[2,4]<0.72,
                               "Medium","High"))
    validation2[3,6] <- ifelse(validation2[3,4] <= 0.49,"Low",ifelse(validation2[2,4]<0.72,
                                "Medium","High"))
    validation2[1,6] <- " "
      
    cat(" ", "\n")
    cat("Surrogacy evaluation criterion", "\n")
    print(validation2)
    cat("---","\n")
    cat("Correlation strength: <= 0.49 'Low'; ]0.49 - 0.72[ 'Medium'; >= 0.72 'High' ","\n")
    cat("---","\n")
    
    cat(c("Surrogate threshold effect (STE) :",round(ste(object),len),"(HR =",round(exp(ste(object)),len),")"),"\n")
    
    cat(" ", "\n")
    cat("Convergence parameters", "\n")
    cat(c("Penalized marginal log-likelihood = ", round(object$loglikPenal, len)), "\n")
    cat(c("Number of iterations = ", object$n.iter),"\n")
    cat("LCV = the approximate likelihood cross-validation criterion", "\n")
    cat(c("      in the semi parametrical case     = ", round(object$LCV, len)),"\n")
    cat("Convergence criteria:", "\n")
    EPS <- formatC(object$EPS, d, format = "g")
    cat(c("  parameters = ",EPS[1], "likelihood = ", EPS[2], "gradient = ", EPS[3]), "\n")
  }

