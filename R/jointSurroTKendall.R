
#' Kendall's \eqn{\tau} estimation using numerical integration methods
#' 
#' This function estimate the Kendall's \eqn{\tau} based on the joint surrogate model 
#' described in \link{jointSurroPenal} (Sofeu \emph{et al.}, 2018), for the evaluation of 
#' a candidate surrogate endpoints, at the individual-level . We used the Monte-carlo and the gaussian Hermite 
#' quadrature methods for numerical integration. in the event of Gaussian Hermite quadrature, 
#' it is better to choose at least \code{20} quadature nodes for better results. 
#' The actual value of nodes used is the maximum between \code{20} and \code{nb.gh}
#'
#' @aliases jointSurroTKendall
#' @usage 
#' jointSurroTKendall(object = NULL, theta, gamma, alpha = 1, zeta = 1, 
#'                    sigma.v = matrix(rep(0,4),2,2), int.method = 0, 
#'                    nb.MC.kendall = 10000, nb.gh = 32, 
#'                    random.generator = 1, random = 0, 
#'                    random.nb.sim = 0, seed = 0, ui = 1)
#'
#' @param object An object inheriting from \code{jointSurroPenal} class. The default is \code{NULL}
#' @param theta Variance of the individual-level random effect, \if{latex}{\eqn{\omega_{ij}}} 
#' \if{html}{\eqn{\omega}\out{<sub>ij</sub>}}. 
#' Required if \code{object} is set to \code{NULL}
#' @param gamma Variance of the trial-level random effect associated with the baseline risk, 
#' \if{latex}{\eqn{u_i}} \if{html}{\code{u}\out{<sub>i</sub>}}. 
#' Required if \code{object} is set to \code{NULL}. The default is \code{3.5}.
#' @param alpha Power parameter associated with \if{latex}{\eqn{u_i}} 
#' \if{html}{\code{u}\out{<sub>i</sub>}}. Required if \code{object} is set to \code{NULL}.
#'  The default is \code{1}.
#' @param zeta Power parameter associated with \if{latex}{\eqn{\omega_{ij}}} 
#' \if{html}{\eqn{\omega}\out{<sub>ij</sub>}}. Required if \code{object} is set to \code{NULL} 
#' The default is \code{1}.
#' @param sigma.v Covariance matrix  of the random effects treatment-by-trial interaction 
#' (\if{latex}{\eqn{v_{S_i}}, \eqn{v_{T_i}}}\if{html}{v\out{<sub>S<sub>i</sub></sub>},
#' v\out{<sub>T<sub>i</sub></sub>}})
#' @param int.method A numeric, indicates the integration method: \code{0} for Monte carlo and 
#' \code{1} for Gaussian-Hermite quadrature. The default is \code{0}
#' @param nb.MC.kendall Number of generated points used with the Monte-Carlo to estimate
#' integrals in the Kendall's \eqn{\tau} formulation. Beter to use at least 4000 points for
#' stable results. The default is \code{10000}.
# @param method.int.kendall A numeric, indicates in the event of the Monte-carlo integration, if only one 
# kendall's \eqn{\tau} should be considered (\code{1}), or four kendall's \eqn{\tau}, according to the 
# randomization group of considered two patiens used for kendall's \eqn{\tau} estimation (\code{0}).
# The default is \code{1}
#' @param nb.gh Number of nodes for the Gaussian-Hermite quadrature.  The default is \code{32}.
#' @param random.generator Random number generator to use by the Fortran compiler, 
#' \code{1} for the intrinsec subroutine \code{Random_number} and \code{2} for the 
#' subroutine \code{uniran()}. The default is \code{1}.
#' @param random A binary that says if we reset the random number generation with a different environment 
#' at each call \code{(1)} or not \code{(0)}. If it is set to \code{1}, we use the computer clock 
#' as a seed. In the last case, it is not possible to reproduce the generated datesets". 
#' The default is \code{0}.
#' @param random.nb.sim If \code{random} is set to \code{1}, a binary that indicates the number 
#' of generations that will be made.
#' @param seed The seed to use for data (or samples) generation. required if \code{random} is set to \code{0}. 
#' The default is \code{0}.
#' @param ui A binary, indicates whether one considered trial random effect associated with 
#' the baseline risk (\code{1}) or not (\code{0}). The default is \code{1}.
#'
#' @return This function return the estimated Kendall's \eqn{\tau} 
#' 
#' @seealso \code{\link{jointSurrSimul}}, \code{\link{summary.jointSurroPenal}}
#' @author Casimir Ledoux Sofeu \email{casimir.sofeu@u-bordeaux.fr}, \email{scl.ledoux@gmail.com} and 
#' Virginie Rondeau \email{virginie.rondeau@inserm.fr}
#' 
#' @references
#'
#' Sofeu C.L., Emura T. and Rondeau V. (2018). One-step validation method for surrogate 
#' endpoints in multiple randomized cancer clinical trials with failure-time endpoints. 
#' \code{Under review}
#' 
#' @export
#'
#' @examples
#' Ktau1 <- jointSurroTKendall(theta = 3.5, gamma = 2.5, nb.gh = 32)
#' Ktau2 <- jointSurroTKendall(theta = 1, gamma = 0.8, alpha = 1, zeta = 1, 
#'          nb.gh = 32)
#' 
##' ###---Kendall's \eqn{\tau} from a joint surrogate model ---###
##' 
##' \dontrun{
##' data.sim <-jointSurrSimul(n.obs=400, n.trial = 20,cens.adm=549, 
##'           alpha = 1.5, theta = 3.5, gamma = 2.5, zeta = 1, 
##'           sigma.s = 0.7, sigma.t = 0.7,cor = 0.8, betas = -1.25, 
##'           betat = -1.25, full.data = 0, random.generator = 1, 
##'           seed = 0, nb.reject.data = 0)
##'           
##' 
##' ###---Estimation---###
##' joint.surrogate <- jointSurroPenal(data = data.sim, nb.mc = 300, 
##'                    nb.gh = 20, indicator.alpha = 1, n.knots = 6)
##'                    
#'  Ktau3 <- jointSurroTKendall(joint.surrogate)
#'  Ktau4 <- jointSurroTKendall(joint.surrogate,nb.MC.kendall = 4000,
#'           seed = 1)
#' }
#' 
jointSurroTKendall <- function(object = NULL, theta, gamma, alpha = 1, zeta = 1,
                               sigma.v = matrix(rep(0,4),2,2), int.method = 0, 
                               nb.MC.kendall = 10000, nb.gh = 32, random.generator = 1, random = 0, 
                               random.nb.sim = 0, seed = 0, ui = 1){
  
  
  # estimation par quadrature de gauss hermite non adaptative
  adaptative = 0
  ui_chap_Essai=rep(0,4)
  invBi_chol_Essai_k = matrix(rep(0,16), nrow = 4, ncol = 4)
  method.int.kendall <- 1
  
  # si on a un seul taux de kendall, sans vouloir prendre en compte l'heterogeneite sur les risques de base
  if(int.method == 0 & ui == 0){
    if(method.int.kendall==1) 
      method.int.kendall <- 5
  }
  
  # si on a un seul taux de kendall, et l'on voudrait prendre en compte l'heterogeneite sur les risques de base
  if(int.method == 0 & ui == 1){
    if(method.int.kendall==1) 
      method.int.kendall <- 4
  }
  
  if(!is.null(object)){
    theta <- object$theta
    gamma <- object$gamma
    alpha <- object$alpha
    zeta <- object$zeta
    ui <- object$ui
    sigma.v[1,] <- c(object$Coefficients["sigma_s",1],object$Coefficients["sigma_sT",1])
    sigma.v[2,] <- c(object$Coefficients["sigma_sT",1],object$Coefficients["sigma_t",1])
  }
  
  # recherche des points et poids de quadrature
  npg <- max(nb.gh,20)
  xx1 <- statmod::gauss.quad(npg,kind="hermite")$nodes
  ww1 <- statmod::gauss.quad(npg,kind="hermite")$weights * exp(xx1**2) # w = w*exp(x**2)
  
  if(ui==1){
    ndim <- 4
  }
  else{
    ndim <- 2
    ui_chap_Essai <- ui_chap_Essai[1:2]
    invBi_chol_Essai_k <- invBi_chol_Essai_k[1:2,1:2]
  }
  
  ss <- as.double(0)
  tau.kendal.00 <- as.double(0)
  tau.kendal.11 <- as.double(0)
  tau.kendal.10 <- as.double(0)
  tau.kendal.01 <- as.double(0)
  ans <- .Fortran(C_jointsurrokendall,
                  as.double(theta),
                  as.double(gamma),
                  as.double(alpha),
                  as.double(zeta),
                  as.integer(adaptative),
                  as.integer(npg),
                  as.integer(ndim),
                  ui_chap_Essai,
                  as.matrix(invBi_chol_Essai_k),
                  xx1,
                  ww1,
                  as.double(sigma.v), 
                  as.integer(int.method), 
                  as.integer(nb.MC.kendall),
                  as.integer(method.int.kendall), 
                  as.integer(random.generator),
                  as.integer(random),
                  as.integer(random.nb.sim),
                  as.integer(seed),
                  tau.kendal.00 = tau.kendal.00,
                  tau.kendal.01 = tau.kendal.01, 
                  tau.kendal.10 = tau.kendal.10, 
                  tau.kendal.11 = tau.kendal.11,
                  ss = ss,
                  PACKAGE="frailtypack"
  )
  # if((int.method == 1) | ((int.method == 0) & (method.int.kendall == 1))){
    return(ans$ss)
  # }else{
  #   result <- NULL
  #   result$tau.kendal.00 <- ans$tau.kendal.00
  #   result$tau.kendal.11 <- ans$tau.kendal.11
  #   result$tau.kendal.10 <- ans$tau.kendal.10
  #   result$tau.kendal.01 <- ans$tau.kendal.01
  #   return(result)
  # }
}