#' Plot of the prediction of the treatment effect on the true endpoint
#' 
#' Plot the prediction of the treatment effect on the true endpoint based on the observed treatment effect
#' on the surrogate endpoint, with the prediction interval: results from the one-step Joint surrogate model  
#' for evaluating a canditate surrogate endpoint. The graphic also includes a vertical line that cut 
#' the x axis to the value of \link{ste}.
#'
#'
#' @aliases plotTreatPredJointSurro
#' @usage
#' plotTreatPredJointSurro(object, from = -2, to = 2, type = "Coef", 
#'    var.used = "error.estim", alpha. = 0.05, n = 1000, lty = 2, d = 3, 
#'    colCI = "blue", xlab = "beta.S", ylab = "beta.T.predict", 
#'    pred.int.use = "up")
#' 
#' @param object An object inheriting from \code{jointSurroPenal} class
##' (output from calling the function \code{jointSurroPenal} ).
#' @param from The range (with \code{to}) over which the function will be plotted. The default is 
#' \code{from -2 to 2}
#' @param to The range (with \code{from}) over which the function will be plotted. The default is 
#' \code{from -2 to 2}
#' @param type The type of graphic, \code{"Coef"} for the \code{log HR} or \code{"HR"} for hazard ratio.
#' If set to \code{HR}, the arguments \code{from} and \code{to} must take positive values.
#' The default is \code{"Coef"}.
#' @param var.used This argument can take two values. The first one is \code{"error.estim"}
##' and indicates if the prediction error take into account
##' the estimation error of the estimates of the parameters. If the estimates 
##' are supposed to be known or if the dataset includes a high number of trials with 
##' a high number of subject per trial, value \code{No.error} can be used. 
##' The default is \code{error.estim}.
#' @param alpha. The confidence level for the prediction interval. The default is \code{0.05}
#' @param n Integer; the number of \if{latex}{\eqn{\beta_S}} \if{html}{\eqn{\beta}\out{<sub>S</sub>}} 
#' values at which to evaluate. The default is \code{1000}.
#' @param lty The line type. Line types can either be specified as an integer 
#' (0=blank, 1=solid (default), 2=dashed, 3=dotted, 4=dotdash, 5=longdash, 6=twodash) or as one 
#' of the character strings \code{"blank", "solid", "dashed", "dotted", "dotdash"}, \code{"longdash",
#'  or "twodash"}, where \code{"blank"} uses "invisible lines" (i.e., does not draw them). 
#'  The default is \code{2}.
#' @param d The desired number of digits after the decimal point for parameters
##' and confidence intervals. Default of 3 digits is used.
#' @param colCI The color used to display the confidence interval.
#' @param xlab A title for the x axis.
#' @param ylab A title for the y axis.
#' @param pred.int.use A character string that indicates the bound of the prediction interval 
#' to use to compute the STE. Possible values are \code{up} for the upper bound (the default)
#' or \code{lw} for the lower bound. \code{up} induces protective treatment effects and \code{lw}
#' induces risk factors.
#'
#' @return For a considered treatment effects on the surrogate enpoint, plot the
#' associated treatment effects on the true endpoint predicted from the joint surrogate model
#' with the prediction interval.
#' 
#' @seealso \code{\link{jointSurroPenal}, \link{jointSurroCopPenal}, \link{predict.jointSurroPenal}}
#' 
#' @author Casimir Ledoux Sofeu \email{casimir.sofeu@u-bordeaux.fr}, \email{scl.ledoux@gmail.com} and 
#' Virginie Rondeau \email{virginie.rondeau@inserm.fr}
#' 
#' @references 
#' Burzykowski T, Buyse M (2006). "Surrogate threshold effect: an alternative 
#' measure for meta-analytic surrogate endpoint validation." Pharmaceutical 
#' Statistics, 5(3), 173-186.ISSN 1539-1612.
#' 
#' Sofeu, C. L. and Rondeau, V. (2020). How to use frailtypack for validating failure-time surrogate 
#' endpoints using individual patient data from meta-analyses of randomized controlled trials. 
#' PLOS ONE; 15, 1-25.
#' @importFrom graphics points curve
#' @keywords surrogate prediction
#' @export
#'
#' @examples
#' 
#' \dontrun{
#' 
#' 
#' ###--- Joint surrogate model ---###
#' ###---evaluation of surrogate endpoints---###
#' 
#' data(dataOvarian)
#' joint.surro.ovar <- jointSurroPenal(data = dataOvarian, n.knots = 8, 
#'                 init.kappa = c(2000,1000), indicator.alpha = 0, 
#'                 nb.mc = 200, scale = 1/365)
#' 
#' ## "HR"
#' plotTreatPredJointSurro(joint.surro.ovar, from = 0, to = 4, 
#'                 type = "HR", var.used = "error.estim", lty = 2)
#'              
#' ## "log HR"
#' plotTreatPredJointSurro(joint.surro.ovar, from = -2, to = 2, 
#'                 type = "Coef", var.used = "error.estim", lty = 2)
#'                 
#' ### For a value of ste greater than 0 (HR > 1), which induces deleterious
#' ### treatment effet, argument "pred.int.use" can be set to "lw"  
#' 
#' plotTreatPredJointSurro(joint.surro.ovar, from = 0, to = 2, 
#'                 type = "HR", var.used = "error.estim", lty = 2,
#'                 pred.int.use = "lw")
#' 
#' }
#' 
plotTreatPredJointSurro <- function(object, from = -2, to = 2, type = "Coef", var.used = "error.estim", 
                                       alpha. = 0.05, n = 1000, lty = 2, d = 3, colCI = "blue", xlab = "beta.S", 
                                       ylab = "beta.T.predict", pred.int.use = "up"){
  # type  = "coef" or "HR"
  # n = number of points for the curve
  # colCI = color Confidence interval
  
  # data control
  if((type == "HR") & (to < 0 | from < 0) ){
    if(to > 0) {
      warning("With argument 'type' = 'HR', arguments 'to' and 'from' must take positive values")
      from <- 0
    }
    else{
      stop("With argument 'type' = 'HR', arguments 'to' and 'from' must take positive values")
    }
  }
    
  beta  <- object$beta.t
  dab   <- object$Coefficients$Estimate[nrow(object$Coefficients)-4]
  daa   <- object$Coefficients$Estimate[nrow(object$Coefficients)-6]
  dbb   <- object$Coefficients$Estimate[nrow(object$Coefficients)-5]
  alpha <- object$beta.s
  # alpha0 <- matrixPred$beta.S[i]
  x1     <- t(matrix(c(1, -dab/daa),1,2))
  
  Vmu   <- matrix(c(object$varcov.Sigma[3,3], object$varcov.Sigma[3,1], 
                    object$varcov.Sigma[3,1], object$varcov.Sigma[1,1]),2,2)
  nparam <- nrow(object$varH)
  
  VD     <- matrix(c(object$varH[nparam,nparam], object$varH[nparam,nparam-1],
                     object$varH[nparam-1,nparam], object$varH[nparam -1,nparam - 1]),2,2)
  R2trial <- object$Coefficients$Estimate[nrow(object$Coefficients)-1] 
  
  if(var.used == "error.estim"){
    if(type == "Coef"){ # log HR
      expressx <- function(x){
        beta + (dab/daa) * (x - alpha)
      }
      
      expressxVect <- Vectorize(expressx)
      curve (expr = expressxVect, from = from, to = to, n = n, xlab = xlab, ylab = ylab, 
             main = paste("STE = ", round(ste(object, var.used = var.used, pred.int.use = pred.int.use), d), 
                          "(HR = ", round(exp(ste(object, var.used = var.used, pred.int.use = pred.int.use)), d), ")"))
      #inf
      expressInf <- function(x){
        beta + (dab/daa) * (x - alpha) - qnorm(1-alpha./2) * 
                    sqrt(t(x1) %*% (Vmu + (((x - alpha)/daa)**2) * VD) %*% x1
                         + dbb * (1 - R2trial))
      }
      
      expressInfVect <- Vectorize(expressInf)
      curve (expr = expressInfVect, from = from, to = to, add = T, col = colCI, n = n, lty = lty)
      #sup
      expressSup <- function(x){
        beta + (dab/daa) * (x - alpha) + qnorm(1-alpha./2) * 
          sqrt(t(x1) %*% (Vmu + (((x - alpha)/daa)**2) * VD) %*% x1
               + dbb * (1 - R2trial))
      }
      
      expressSupVect <- Vectorize(expressSup)
      curve (expr = expressSupVect, from = from, to = to, add = T, col = colCI, n = n, lty = lty)
    }
    else{ # HR
      
      expressx <- function(x){
        x <- log(x) # on suppose que les entrees sont des HR et donc on les converti en log HR
        exp(beta + (dab/daa) * (x - alpha))
      }
      
      expressxVect <- Vectorize(expressx)
      
      curve (expr = expressxVect, from = from, to = to, n = n, xlab = xlab, ylab = ylab, 
             main = paste("STE = ", round(ste(object, var.used = var.used, pred.int.use = pred.int.use), d), 
                          "(HR = ", round(exp(ste(object, var.used = var.used, pred.int.use = pred.int.use)), d), ")"))
      #inf
      expressInf <- function(x){
        x <- log(x) # on suppose que les entrees sont des HR et donc on les converti en log HR
        exp(beta + (dab/daa) * (x - alpha) - qnorm(1-alpha./2) * 
          sqrt(t(x1) %*% (Vmu + (((x - alpha)/daa)**2) * VD) %*% x1
               + dbb * (1 - R2trial)))
      }
      
      expressInfVect <- Vectorize(expressInf)
      curve (expr = expressInfVect, from = from, to = to, add = T, col = colCI, n = n, lty = lty)
      #sup
      expressSup <- function(x){
        x <- log(x) # on suppose que les entrees sont des HR et donc on les converti en log HR
        exp(beta + (dab/daa) * (x - alpha) + qnorm(1-alpha./2) * 
          sqrt(t(x1) %*% (Vmu + (((x - alpha)/daa)**2) * VD) %*% x1
               + dbb * (1 - R2trial)))
      }
      
      expressSupVect <- Vectorize(expressSup)
      curve (expr = expressSupVect, from = from, to = to, add = T, col = colCI, n = n, lty = lty)
    }
  }
  else{
    if(type == "Coef"){ # log HR
      expressx <- function(x){
        beta + (dab/daa) * (x - alpha)
      }
      
      expressxVect <- Vectorize(expressx)
      curve (expr = expressxVect, from = from, to = to, n = n, xlab = xlab, ylab = ylab, 
             main = paste("STE = ", round(ste(object, var.used = var.used, pred.int.use = pred.int.use), d), 
                          "(HR = ", round(exp(ste(object, var.used = var.used, pred.int.use = pred.int.use)), d), ")"))
      #inf
      expressInf <- function(x){
        beta + (dab/daa) * (x - alpha) - qnorm(1-alpha./2) * sqrt(
          dbb * (1 - R2trial))
      }
      
      expressInfVect <- Vectorize(expressInf)
      curve (expr = expressInfVect,
        from = from, to = to, add = T, col = colCI, n = n, lty = lty)
      
      #sup
      expressSup <- function(x){
        beta + (dab/daa) * (x - alpha) + qnorm(1-alpha./2) * sqrt(
          dbb * (1 - R2trial))
      }
      
      expressSupVect <- Vectorize(expressSup)
      curve (expr = expressSupVect,
        from = from, to = to, add = T, col = colCI, n = n, lty = lty)
    }
    else{ # HR
      expressx <- function(x){
        x <- log(x) # on suppose que les entrees sont des HR et donc on les converti en log HR
        exp(beta + (dab/daa) * (x - alpha))
      }
      
      expressxVect <- Vectorize(expressx)
      curve (expr = expressxVect, from = from, to = to, n = n, xlab = xlab, ylab = ylab, 
             main = paste("STE = ", round(ste(object, var.used = var.used, pred.int.use = pred.int.use), d), 
                          "(HR = ", round(exp(ste(object, var.used = var.used, pred.int.use = pred.int.use)), d), ")"))
      #inf
      expressInf <- function(x){
        x <- log(x) # on suppose que les entrees sont des HR et donc on les converti en log HR
        exp(beta + (dab/daa) * (x - alpha) - qnorm(1-alpha./2) * sqrt(
          dbb * (1 - R2trial)))
      }
      
      expressInfVect <- Vectorize(expressInf)
      curve (expr = expressInfVect,
        from = from, to = to, add = T, col = colCI, n = n, lty = lty)
      #sup
      expressSup <- function(x){
        x <- log(x) # on suppose que les entrees sont des HR et donc on les converti en log HR
        exp(beta + (dab/daa) * (x - alpha) + qnorm(1-alpha./2) * sqrt(
          dbb * (1 - R2trial)))
      }
      
      expressSupVect <- Vectorize(expressSup)
      curve (expr = expressSupVect,
        from = from, to = to, add = T, col = colCI, n = n, lty = lty)
    }
  }
  
  #ste
  if(type == "HR"){ # log HR
    abline(h = 1, col = "cyan", lty = 4)
    points(exp(ste(object, var.used = var.used, pred.int.use = pred.int.use)),1)
    abline(v = exp(ste(object, var.used = var.used, pred.int.use = pred.int.use)), col = "cyan", lty = 4)
  }
  else{
    abline(h = 0, col = "cyan", lty = 4)
    points(ste(object, var.used = var.used, pred.int.use = pred.int.use),0)
    abline(v = ste(object, var.used = var.used, pred.int.use = pred.int.use), col = "cyan", lty = 4)
  }
  
}