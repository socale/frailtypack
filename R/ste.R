##' Surrogate threshold effect for the one-step Joint surrogate model for the evaluation of a 
##' canditate surrogate endpoint.
##' 
##' @description{
##' This function compute the surrogate threshold effect (STE) from the one-step joint frailty 
##' \code{\link[=jointSurroPenal]{model}} or joint frailty-copula 
##' \code{\link[=jointSurroCopPenal]{model}}. The STE is defined as the minimum treament effect 
##' on surrogate endpoint, necessary to predict a non-zero effect on 
##' true endpoint (Burzykowski \emph{et al.}, 2006).
##' }
##' 
##' @details{  
##' The STE is obtained by solving the equation \if{latex}{\eqn{l(\alpha_0)}} 
#' \if{html}{\code{l}(\eqn{\alpha}\out{<sub>0</sub>})} \code{= 0}
#'   (resp. \if{latex}{\eqn{u(\alpha_0)}} 
#' \if{html}{\code{u}(\eqn{\alpha}\out{<sub>0</sub>})} \code{= 0}), where 
#' \if{latex}{\eqn{\alpha_0}} \if{html}{\eqn{\alpha}\out{<sub>0</sub>}} represents
##' the corresponding STE, and \if{latex}{\eqn{l(\alpha_0)}} 
#' \if{html}{\code{l}(\eqn{\alpha}\out{<sub>0</sub>})} (resp. \if{latex}{\eqn{u(\alpha_0)}} 
#' \if{html}{\code{u}(\eqn{\alpha}\out{<sub>0</sub>})}) is the lower (resp. upper) bound of the prediction interval 
##' of the treatment effect on the true endpoint \if{html}{(\eqn{\beta} + b\out{<sub>0</sub>})} \if{latex}{(\eqn{\beta} + \eqn{b_0})}. Thereby,
##' 
##' \if{latex}{\eqn{l(\alpha_0) \equiv 
##' E(\beta + b_0|\alpha_0, \vartheta) - Z_{1-(\gamma/2)} \sqrt(Var(\beta + b_0|\alpha_0, \vartheta))}
##' 
##' and 
##' 
##' \eqn{u(\alpha_0) \equiv 
##' E(\beta + b_0|\alpha_0, \vartheta) + Z_{1-(\gamma/2)} \sqrt(Var(\beta + b_0|\alpha_0, \vartheta))}
##' }
##' \if{html}{\figure{ste.png}{options: width="60\%"}}
##' 
##' where \if{latex}{\eqn{\vartheta}}\if{html}{\figure{vartheta.png}{options: width="3\%"}}
##' represents the set of estimates for the fixed-effects and the 
##' variance-covariance parameters of the random effects obtained from the joint surrogate 
##' \code{\link[=jointSurroPenal]{model}} 
##' (Sofeu \emph{et al.}, 2018). 
##' 
##' If the previous equations gives two solutions, STE can be the 
##' minimum (resp. the maximum) value or both of them, according to the shape of the function. 
##' If the concavity of the function is turned upwards, STE is the first value and
##' the second value represents the maximum (res. the minimum) treament value observable 
##' on the surrogate that can predict a non zero treatment effect on true endpoint. 
##' If the concavity of the function is turned down, both the two solutions
##' represent the STE and the interpretation is such that accepted values of the 
##' treatment effects on \code{S} predict a no zero treatment effects on \code{T}
##' 
##' Given that negative values of treatment effect indicate a reduction of the risk 
##' of failure and are considered beneficial, STE is recommended to be computed from 
##' the upper prediction
##' limit  \if{latex}{\eqn{u(\alpha_0)}} 
#' \if{html}{\code{u}(\eqn{\alpha}\out{<sub>0</sub>})}.
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
##' (output from calling the \code{jointSurroPenal} or \code{jointSurroCopPenal} function ).
##' @param var.used This argument takes two values. The first one is \code{"error.estim"}
##' and indicates if the prediction error take into account
##' the estimation error of the estimates of the parameters. If the estimates 
##' are supposed to be known or if the dataset includes a high number of trials with 
##' a high number of subject per trial, value \code{No.error} can be used. 
##' The default is \code{error.estim}.
##' @param alpha. The confidence level for the prediction interval. The default is \code{0.05}
##' @param pred.int.use A character string that indicates the bound of the prediction interval 
##' to use to compute the STE. Possible values are \code{up} for the upper bound (the default)
##' or \code{lw} for the lower bound. \code{up} induces protective treatment effects and \code{lw}
#' induces risk factors.
##' 
##' @return Returns and displays the STE.
##' @seealso \code{\link{jointSurroPenal}, \link{jointSurroCopPenal}}, \code{\link[=predict.jointSurroPenal]{predict}}
##' 
##' @author Casimir Ledoux Sofeu \email{casimir.sofeu@u-bordeaux.fr}, \email{scl.ledoux@gmail.com} and 
##' Virginie Rondeau \email{virginie.rondeau@inserm.fr}
##' 
##' @references 
##' Burzykowski T, Buyse M (2006). "Surrogate threshold effect: an alternative 
##' measure for meta-analytic surrogate endpoint validation." Pharmaceutical 
##' Statistics, 5(3), 173-186.ISSN 1539-1612.
##' 
##' Sofeu, C. L., Emura, T., and Rondeau, V. (2019). One-step validation method for surrogate 
##' endpoints using data from multiple randomized cancer clinical trials with failure-time endpoints. 
##' Statistics in Medicine 38, 2928-2942. 
##' 
##' Sofeu, C. L. and Rondeau, V. (2020). How to use frailtypack for validating failure-time surrogate 
##' endpoints using individual patient data from meta-analyses of randomized controlled trials. 
##' PLOS ONE; 15, 1-25.
##' 
##' @keywords surrogate prediction ste Surrogate threshold effect
##' @export
##' @importFrom stats optimize
##' @importFrom rootSolve uniroot.all
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
    
    #===proble dans l'indexation de la matrice des coefficient. scl: 11/04/2020
    # dab   <- object$Coefficients$Estimate[nrow(object$Coefficients)-4]
    # daa   <- object$Coefficients$Estimate[nrow(object$Coefficients)-6]
    # dbb   <- object$Coefficients$Estimate[nrow(object$Coefficients)-5]
    dab   <- object$Coefficients[rownames(object$Coefficients)=="sigma_sT",1]
    daa   <- object$Coefficients[rownames(object$Coefficients)=="sigma_s",1]
    dbb   <- object$Coefficients[rownames(object$Coefficients)=="sigma_t",1]
    #====

    
    #x <- matrixPred$beta.S[i]
    x.     <- t(matrix(c(1, -dab/daa),1,2))
    
    # Vmu (sigma_ST, sigma_SS). on utilise la matrice obtenu par delta methode a partir de la hessienne
    Vmu   <- matrix(c(object$varcov.Sigma[3,3], object$varcov.Sigma[3,1],
                      object$varcov.Sigma[3,1], object$varcov.Sigma[1,1]),2,2)
    nparam <- nrow(object$varH)
    
    # VD (bete_T, beta_S). on utilise la hesienne directement car pas de changement de variable
    VD     <- matrix(c(object$varH[nparam,nparam], object$varH[nparam,nparam-1],
                       object$varH[nparam-1,nparam], object$varH[nparam -1,nparam - 1]),2,2)
    R2trial <- object$Coefficients[rownames(object$Coefficients)=="R2trial",1] # scl: 11/04/2020
    
    if(var.used == "error.estim") {
      if(pred.int.use == "lw"){
        return(
          sapply(x, function(x, beta, dab, daa, dbb, alpha, alpha., x., Vmu, VD, R2trial) 
            beta + (dab/daa) * (x - alpha) - qnorm(1-alpha./2) * 
              sqrt(t(x.) %*% (Vmu + (((x - alpha)/daa)**2) * VD) %*% x.
                   + dbb * (1 - R2trial)), 
            beta = beta, dab = dab, daa = daa, dbb = dbb, alpha = alpha, alpha. = alpha., 
            x. = x., Vmu = Vmu, VD = VD, R2trial = R2trial
          )
        )
       }
      else{
        return(
          sapply(x, function(x, beta, dab, daa, dbb, alpha, alpha., x., Vmu, VD, R2trial) 
            beta + (dab/daa) * (x - alpha) + qnorm(1-alpha./2) * sqrt(t(x.) %*% 
            (Vmu + (((x - alpha)/daa)**2) * VD) %*% x. + dbb * (1 - R2trial)), 
            beta = beta, dab = dab, daa = daa, dbb = dbb, alpha = alpha, alpha. = alpha., 
            x. = x., Vmu = Vmu, VD = VD, R2trial = R2trial
            )
        )
      }
    }
    else{
      if(pred.int.use == "lw"){
        return(
          sapply(x, function(x, beta, dab, daa, dbb, alpha, alpha., R2trial) 
            beta + (dab/daa) * (x - alpha) - qnorm(1-alpha./2) * 
              sqrt(dbb * (1 - R2trial)),
            beta = beta, dab = dab, daa = daa, dbb = dbb, alpha = alpha, alpha. = alpha., 
            R2trial = R2trial
          )
        )
      }
      else{
        return(
          
          sapply(x, function(x, beta, dab, daa, dbb, alpha, alpha., R2trial) 
            beta + (dab/daa) * (x - alpha) + qnorm(1-alpha./2) *
              sqrt(dbb * (1 - R2trial)),
            beta = beta, dab = dab, daa = daa, dbb = dbb, alpha = alpha, alpha. = alpha., 
            R2trial = R2trial
          )
        )
      }
    }
    
  }
  # ste <- optimize(f, c(-1e8,1e8), object = object, var.used = var.used, alpha. = alpha.,
  #                 pred.int.use = pred.int.use)
  ste <- uniroot.all(f, lower = -10, upper = 10, object = object, var.used = var.used, alpha. = alpha.,
                  pred.int.use = pred.int.use)
  # suivant le nombre de solutions de l'equation, le "STE" peut etre un point ou alors un ensemble de deux points
  # definissants l'ensemble des valeurs prises par "bateS" pour une prediction significativement differente de 0
  # de "betaT"
  
  if(length(ste) == 0){# on est dans le cas Delta = 0, pas de solution entire pour cette equation
    message("Warning : STE does not exist for this intermediate endpoint. Therefore, 
            regarding the values of R2trial and Kendall tau, the observed treatment effect on the candidate 
            surrogate endpoint can not permitted to predict a non zero treatment effect on true endpoint
            using the considered joint surrogate model and the meta-analysis")
  }else{
    if(length(ste) == 2){
      # recherche du sens de la concavite (bref, signe de "a" dans l'equation "ax^2 + bx + c")
      # je prends un pont au hazard dans l'intervalle [x1,x2] et je regarde le signe de son image
      if(f(sum(ste)/2, object = object, var.used = var.used, alpha. = alpha.,
           pred.int.use = pred.int.use) < 0){ # concavite tournee vers le haut
        message("The treatement effects on the surrogate (beta_S) that can predict a non zero 
                treatment effect on true endpoint (beta_T) belong to the intervall: ]",
                round(ste[1], 3), " ; ", round(ste[2], 3), "[ : HR= ]", round(exp(ste[1]), 3), " ; ", round(exp(ste[2]), 3), "[")
      }
      else { # concavite tournee vers le bas
        message("The treatement effects on the surrogate (beta_S) that can predict a non zero 
                treatment effect on true endpoint (beta_T) belong to the intervall: ]-Inf ; ",
                round(ste[1], 3), "[ U ]", round(ste[2], 3), " ; +Inf[ : HR= ]-Inf ; ",
                round(exp(ste[1]), 3), "[ U ]", round(exp(ste[2]), 3), " ; +Inf[")
      }
    } else{ # une seule solution
      message("The treatement effects on the surrogate (beta_S) that can predict a non zero 
                treatment effect on true endpoint (beta_T) belong to the intervall: ]-Inf ; ",
              round(ste, 3), "[ : HR= ]-Inf ; ",round(exp(ste, 3)), "[")
    }
  }
  return(ste)
}


