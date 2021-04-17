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
##' (Sofeu \emph{et al.}, 2019). 
##' 
##' If the previous equations gives two solutions, STE can be the 
##' minimum (resp. the maximum) value or both of them, according to the shape of the function. 
##' If the concavity of the function is turned upwards, STE is the first value and
##' the second value represents the maximum (res. the minimum) treament value observable 
##' on the surrogate that can predict a nonzero treatment effect on true endpoint. 
##' If the concavity of the function is turned down, both of the solutions
##' represent the STE and the interpretation is such that accepted values of the 
##' treatment effects on \code{S} predict a nonzero treatment effects on \code{T}
##' 
##' Given that negative values of treatment effect indicate a reduction of the risk 
##' of failure and are considered beneficial, STE is recommended to be computed from 
##' the upper prediction
##' limit  \if{latex}{\eqn{u(\alpha_0)}} 
#' \if{html}{\code{u}(\eqn{\alpha}\out{<sub>0</sub>})}.
##' 
##' The details on the computation of STE are described in 
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
##' and indicates if the prediction error takes into account
##' the estimation error of the estimates of the parameters. If the estimates 
##' are supposed to be known or if the dataset includes a high number of trials with 
##' a high number of subject per trial, value \code{No.error} can be used. 
##' The default is \code{error.estim}, which is highly recommended in practice.
##' @param alpha. The confidence level for the prediction interval. The default is \code{0.05}
##' @param pred.int.use A character string that indicates the bound of the prediction interval 
##' to use to compute the STE. Possible values are \code{up} for the upper bound (the default)
##' or \code{lw} for the lower bound. \code{up} when we have a protective treatment effect and \code{lw} 
##' when we have a deleterious treatment effect (see details).
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
##' # Assuming errors on the estimates
##' ste(joint.surro.ovar, var.used = "error.estim")
##' # Assuming no errors on the estimates
##' ste(joint.surro.ovar, var.used = "No.error", pred.int.use = "up")
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
  
  #++La fonction f() est definie dans le fichier autresFonctions.R
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
            surrogate endpoint is not able to predict a nonzero treatment effect on the true endpoint
            using this model and this dataset")
  }else{
    if(length(ste) == 2){
      # recherche du sens de la concavite (bref, signe de "a" dans l'equation "ax^2 + bx + c")
      # je prends un pont au hazard dans l'interval [x1,x2] et je regarde le signe de son image
      if(f(sum(ste)/2, object = object, var.used = var.used, alpha. = alpha.,
           pred.int.use = pred.int.use) < 0){ # concavite tournee vers le haut
        message("The treatment effects on the surrogate endpoint (beta_S) that can predict a nonzero 
treatment effect on the true endpoint (beta_T) belongs to the interval: ]",
                round(ste[1], 3), " ; ", round(ste[2], 3), "[ : HR= ]", round(exp(ste[1]), 3), " ; ", round(exp(ste[2]), 3), "[")
      }
      else { # concavite tournee vers le bas
        message("The treatment effects on the surrogate endpoint (beta_S) that can predict a nonzero 
treatment effect on the true endpoint (beta_T) belongs to the interval: ]-Inf ; ",
                round(ste[1], 3), "[ U ]", round(ste[2], 3), " ; +Inf[ : HR= ]0 ; ",
                round(exp(ste[1]), 3), "[ U ]", round(exp(ste[2]), 3), " ; +Inf[")
      }
    } else{ # une seule solution
      message("The treatment effects on the surrogate endpoint (beta_S) that can predict a nonzero 
treatment effect on the true endpoint (beta_T) belongs to the interval: ]-Inf ; ",
              round(ste, 3), "[ : HR= ]0 ; ",round(exp(ste), 3), "[")
    }
  }
  return(ste)
}


