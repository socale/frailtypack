##' Surrogate threshold effect for the one-step Joint surrogate model for the evaluation of a 
##' canditate surrogate endpoint.
##' 
##' @description{
##' This function compute the surrogate threshold effect (STE) from the one-step joint 
##' surrogate \code{\link[=jointSurroPenal]{model}}. The STE is defined as the minimum treament effect on the surrogate 
##' necessary to predict a non-zero effect on the true endpoint (Burzykowski \emph{et al.}, 2006).
##' }
##' 
##' @details{  
##' The STE is obtained by solving the equation \eqn{l(\alpha_0) = 0} (resp. \eqn{u(\alpha_0) = 0}), where \eqn{\alpha_0} represents
##' the corresponding STE, and \eqn{l(\alpha_0)} (resp. \eqn{u(\alpha_0}) is the lower (resp. upper) bound of the prediction interval 
##' of the treatment effect on the true endpoint (\eqn{\beta + b_0}). Thereby,
##' 
##' \eqn{l(\alpha_0) \equiv 
##' E(\beta + b_0|\alpha_0, \vartheta) - Z_{1-(\gamma/2)} \sqrt(Var(\beta + b_0|\alpha_0, \vartheta))}
##' 
##' and 
##' 
##' \eqn{u(\alpha_0) \equiv 
##' E(\beta + b_0|\alpha_0, \vartheta) + Z_{1-(\gamma/2)} \sqrt(Var(\beta + b_0|\alpha_0, \vartheta))}
##' 
##' where \eqn{\vartheta} represents the set of estimates for the fixed-effects and the 
##' variance-covariance parameters of the random effects obtained from the joint surrogate 
##' \code{\link[=jointSurroPenal]{model}} 
##' (Sofeu \emph{et al.}, 2018). 
##' 
##' Given that negative values of treatment effect indicate a reduction of the risk 
##' of failure and are considered beneficial, STE is recommended to be computed from 
##' the upper prediction
##' limit \eqn{u(\alpha_0)}.
##' 
##' The details on the computation of STE is describes in 
##' Burzykowski \emph{et al.} (2006).
##' }
##' 
##' @aliases ste 
##' @usage
##' 
##' ste(object, var.used = "error.estim", alpha. = 0.05, 
##'     pred.int.use = "up")
##' @param object An object inheriting from \code{jointSurroPenal} class
##' (output from calling the function \code{jointSurroPenal}).
##' @param var.used This argument takes two values. The first one is \code{"error.estim"}
##' and indicates if the prediction error take into account
##' the estimation error of the estimates of the parameters. If the estimates 
##' are supposed to be known or if the dataset includes a high number of trials with 
##' a high number of subject per trial, value \code{No.error} can be used. 
##' The default is \code{error.estim}.
##' @param alpha. The confidence level for the prediction interval. The default is \code{0.05}
##' @param pred.int.use A character string that indicates the bound of the prediction interval 
##' to use to compute the STE. Possible values are \code{up} for the upper bound (the default)
##' or \code{lw} for the lower bound.
##' 
##' @return Returns and displays the STE.
##' @seealso \code{\link{jointSurroPenal}} \code{\link[=predict.jointSurroPenal]{predict}}
##' 
##' @author Casimir Ledoux Sofeu \email{casimir.sofeu@u-bordeaux.fr}, \email{scl.ledoux@gmail.com} and 
##' Virginie Rondeau \email{virginie.rondeau@inserm.fr}
##' 
##' @references 
##' Burzykowski T, Buyse M (2006). "Surrogate threshold effect: an alternative 
##' measure for meta-analytic surrogate endpoint validation." Pharmaceutical 
##' Statistics, 5(3), 173-186.ISSN 1539-1612.
##' 
##' Sofeu C.L., Emura T. and Rondeau V. (2018). One-step validation method for surrogate 
##' endpoints in multiple randomized cancer clinical trials with failure-time endpoints. 
##' 
##' @keywords surrogate prediction ste Surrogate threshold effect
##' @export
##' @importFrom stats optimize
##' @examples
##' 
##' 
##' \dontrun{
##' 
##' 
##' ###--- Joint surrogate model ---###
##' ###---evaluation of surrogate endpoints---###
##' 
##' data(dataOvarian)
##' joint.surro.ovar <- jointSurroPenal(data = dataOvarian, n.knots = 8, 
##'                 init.kappa = c(2000,1000), indicator.alpha = 0, 
##'                 nb.mc = 200, scale = 1/365)
##' 
##' # ======STE=====
##' # ste(joint.surro.ovar, var.used = "error.estim")
##' # Assuming no errors on the estimates
##' # ste(joint.surro.ovar, var.used = "No.error", pred.int.use = "up")
##' 
##' }
##' 
##' 
ste <- function (object, var.used = "error.estim", alpha. = 0.05, pred.int.use = "up")
{
  if (!inherits(object, "jointSurroPenal"))
    stop("object must be of class 'jointSurroPenal'")
  
  if(! var.used %in% c("error.estim","No.error"))
    stop("Argument 'var.used' must be specified to 'error.estim' or 'No.error' ")
  
  if(! pred.int.use %in% c("up","lw"))
    stop("Argument 'pred.int.use' must be specified to 'up' or 'lw' ")
  
  # beta  <- object$beta.t
  # alpha <- object$beta.s
  # dab   <- object$Coefficients$Estimate[nrow(object$Coefficients)-4]
  # daa   <- object$Coefficients$Estimate[nrow(object$Coefficients)-6]
  # dbb   <- object$Coefficients$Estimate[nrow(object$Coefficients)-5]
  # 
  # #alpha0 <- matrixPred$beta.S[i]
  # x     <- t(matrix(c(1, -dab/daa),1,2))
  # 
  # # Vmu (sigma_ST, sigma_SS). on utilise la matrice obtenu par delta methode a partir de la hessienne
  # Vmu   <- matrix(c(object$varcov.Sigma[3,3], object$varcov.Sigma[3,1],
  #                   object$varcov.Sigma[3,1], object$varcov.Sigma[1,1]),2,2)
  # nparam <- nrow(object$varH)
  # 
  # # VD (bete_T, beta_S). on utilise la hesienne directement car pas de changement de variable
  # VD     <- matrix(c(object$varH[nparam,nparam], object$varH[nparam,nparam-1],
  #                    object$varH[nparam-1,nparam], object$varH[nparam -1,nparam - 1]),2,2)
  # R2trial <- object$Coefficients$Estimate[nrow(object$Coefficients)-1]

  # moyenne conditionnelle de beta + b_0
  # beta.T.i <- beta + (dab/daa) * (alpha0 - alpha)
  # variance.inf <- dbb * (1 - R2trial) 
  # variance.N <- t(x) %*% (Vmu + (((alpha0 - alpha)/daa)**2) * VD) %*% x
  # + variance.inf
  
  # if(var.used == "error.estim") 
  #   variance <- variance.N
  # else 
  #   variance <- variance.inf
  # 
  # l.alpha0 <- beta.T.i - qnorm(1-alpha./2) * sqrt(variance) # borne inferieure
  # u.alpha0 <- beta.T.i + qnorm(1-alpha./2) * sqrt(variance) # borne superieure
  
  # resolution de l'equation l.alpha0 =0 a l'aide de la fonction optimize. cette fonction retourne
  # la solution minimale de l'equation a resoudre
  
  f <- function(x, object, var.used, alpha., pred.int.use){
    beta  <- object$beta.t
    alpha <- object$beta.s
    dab   <- object$Coefficients$Estimate[nrow(object$Coefficients)-4]
    daa   <- object$Coefficients$Estimate[nrow(object$Coefficients)-6]
    dbb   <- object$Coefficients$Estimate[nrow(object$Coefficients)-5]
    
    #x <- matrixPred$beta.S[i]
    x.     <- t(matrix(c(1, -dab/daa),1,2))
    
    # Vmu (sigma_ST, sigma_SS). on utilise la matrice obtenu par delta methode a partir de la hessienne
    Vmu   <- matrix(c(object$varcov.Sigma[3,3], object$varcov.Sigma[3,1],
                      object$varcov.Sigma[3,1], object$varcov.Sigma[1,1]),2,2)
    nparam <- nrow(object$varH)
    
    # VD (bete_T, beta_S). on utilise la hesienne directement car pas de changement de variable
    VD     <- matrix(c(object$varH[nparam,nparam], object$varH[nparam,nparam-1],
                       object$varH[nparam-1,nparam], object$varH[nparam -1,nparam - 1]),2,2)
    R2trial <- object$Coefficients$Estimate[nrow(object$Coefficients)-1]
    
    if(var.used == "error.estim") {
      if(pred.int.use == "lw"){
        return((beta + (dab/daa) * (x - alpha) - qnorm(1-alpha./2) * 
          sqrt(t(x.) %*% (Vmu + (((x - alpha)/daa)**2) * VD) %*% x.
               + dbb * (1 - R2trial)))^2)
       }
      else{
         return((beta + (dab/daa) * (x - alpha) + qnorm(1-alpha./2) * 
                   sqrt(t(x.) %*% (Vmu + (((x - alpha)/daa)**2) * VD) %*% x.
                        + dbb * (1 - R2trial)))^2)
      }
    }
    else{
      if(pred.int.use == "lw"){
        return((beta + (dab/daa) * (x - alpha) - qnorm(1-alpha./2) * 
          sqrt(dbb * (1 - R2trial)))^2)
      }
      else{
        return((beta + (dab/daa) * (x - alpha) + qnorm(1-alpha./2) *
                  sqrt(dbb * (1 - R2trial)))^2)
      }
    }
    
  }
  ste <- optimize(f, c(-1e8,1e8), object = object, var.used = var.used, alpha. = alpha.,
                  pred.int.use = pred.int.use)
  
  #cat(c("STE = ", ste$minimum, "objective = ",ste$objective))
  
  return(ste$minimum)
}


