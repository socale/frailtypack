#' Plot of the prediction of the treatment effect on the true endpoint and the STE
#' 
#' Plot the prediction of the treatment effect on the true endpoint based on the observed treatment effect
#' on the surrogate endpoint, with the prediction interval: results from the one-step Joint surrogate model  
#' for evaluating a canditate surrogate endpoint. The graphic also includes vertical lines that cut 
#' the x axis to the values of \link{ste}. A hatched rectagle/zone indicates the values of 
#' \if{latex}{\eqn{\beta_S}} \if{html}{\eqn{\beta}\out{<sub>S</sub>}} that predict a non zeto 
#' \if{latex}{\eqn{\beta_T}} \if{html}{\eqn{\beta}\out{<sub>T</sub>}}, according to the number of value 
#' for \code{STE} and the shape of the upper confidence limit for the prediction model.
#'
#'
#' @aliases plotTreatPredJointSurro
# 
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
##' The default is \code{error.estim} (highly recommended).
#' @param alpha. The confidence level for the prediction interval. The default is \code{0.05}
#' @param n An integer that indicates the number of values for \if{latex}{\eqn{\beta_S}} \if{html}{\eqn{\beta}\out{<sub>S</sub>}}. The default is \code{1000}.
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
#' or \code{lw} for the lower bound. \code{up} when we have a protective treatment effect and \code{lw} 
##' when we have a deleterious treatment effect.
#' @param main Title of the graphics
#' @param add.accept.area.betaS A boolean that indicates if the plot should add acceptance area for 
#' \if{latex}{\eqn{\beta_S}} \if{html}{\eqn{\beta}\out{<sub>S</sub>}}
#' that predict a nonzero \if{latex}{\eqn{\beta_T}} \if{html}{\eqn{\beta}\out{<sub>T</sub>}}. The default is TRUE
#' @param ybottom A scalar for the left y bottom position of the rectangle on the x-axis associated with acceptable 
#' value for \if{latex}{\eqn{\beta_S}} \if{html}{\eqn{\beta}\out{<sub>S</sub>}} to predict a 
#' non zero \if{latex}{\eqn{\beta_T}} \if{html}{\eqn{\beta}\out{<sub>T</sub>}}. The default is \code{-0.05}.
#' @param ytop A scalar for the top right y position of the rectangle on the x-axis associated with acceptable 
#' value for \if{latex}{\eqn{\beta_S}} \if{html}{\eqn{\beta}\out{<sub>S</sub>}} to predict a 
#' non zero \if{latex}{\eqn{\beta_T}} \if{html}{\eqn{\beta}\out{<sub>T</sub>}}. The default is \code{0.05}.
#' @param density The density of shading lines, in lines per inch.  The default 
#' value of 'NULL' means that no shading lines are drawn.  A 
#' zero value of 'density' means no shading lines whereas
#' negative values (and 'NA') suppress shading (and so allow color filling). The default is \code{20}
#' @param angle Angle (in degrees) of the shading lines. The default is \code{45}
#' @param legend.show A boolean that indicates if the legend should be displayed
#' @param leg.x The x co-ordinate to be used to position the legend.
#' @param leg.y The y co-ordinate to be used to position the legend. The default is \code{4}
#' @param legend A character or expression vector of length >= 1 to appear in the legend
#' @param leg.text.col The color used for the legend text. The default is \code{black}.
#' @param leg.lty The line type, width and color for the legend box (if bty = "o").
#' @param leg.pch = The plotting symbols appearing in the legend, as numeric vector or a 
#' vector of 1-character strings (see \link{points}). Unlike \code{points}, this can all be 
#' specified as a single multi-character string. Must be specified for symbol drawing.
#' @param leg.bg The background color for the legend box. (Note that this is only used if bty \code{!= "n"}.)
#' @param leg.bty The type of box to be drawn around the legend. The allowed values are \code{"o"} 
#' (the default) and \code{"n"}.
#' @param leg.cex Character expansion factor relative to current par(\code{"cex"}). Used 
#' for text as defined in \link{legend}.
#' @param \dots other unused arguments
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
#' @importFrom graphics points curve rect segments
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
#'                 type = "HR", lty = 2, leg.y = 13)
#'                 
#' ## or without acceptance area for betaS:
#' plotTreatPredJointSurro(joint.surro.ovar, from = 0, to = 4, 
#'                 type = "HR", lty = 2, leg.y = 13, 
#'                 add.accept.area.betaS = FALSE)
#'              
#' ## "log HR"
#' plotTreatPredJointSurro(joint.surro.ovar, from = -2, to = 2, 
#'                 type = "Coef", lty = 2, leg.y = 3.5)
#'                 
#' ### For a value of ste greater than 0 (HR > 1), which induces deleterious
#' ### treatment effet, argument "pred.int.use" can be set to "lw"  
#' 
#' plotTreatPredJointSurro(joint.surro.ovar, from = 0, to = 2, 
#'                 type = "HR", lty = 2, leg.y = 4,
#'                 pred.int.use = "lw")
#' 
#' }
#' 
plotTreatPredJointSurro <- function(object, from = -3, to = 2, type = "Coef", var.used = "error.estim", 
                                    alpha. = 0.05, n = 1000, lty = 2, d = 3, colCI = "blue", xlab = "beta.S", 
                                    ylab = "beta.T.predict", pred.int.use = "up", main = NULL, 
                                    add.accept.area.betaS = TRUE, ybottom = -0.05, ytop = 0.05, 
                                    density = 20, angle = 45, legend.show = TRUE, leg.x = NULL, leg.y = 2, 
                                    legend = c("Prediction model", "95% prediction Interval", "Beta.S for nonzero beta.T", "STE"), 
                                    leg.text.col = "black", leg.lty = c(1, 2, 4, NA), 
                                    leg.pch = c(NA, NA, 7, 1), leg.bg = "white", leg.bty = "n", 
                                    leg.cex = 0.85, ...){
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
  
  #ste : il est obtenu a partir de la resolution d'une equation de scond degre 
  # de la forme "ax^2 + bx + c = 0"
  
  if((xlab == "beta.S") & (type == "HR")) xlab = "HR.S (Cox model)"
  if((ylab == "beta.T.predict") & (type == "HR")) ylab = "Predicted HR.T (Joint model)"
  
  STE <- ste(object, var.used = var.used, pred.int.use = pred.int.use)
  
  if(var.used == "error.estim"){
    if(type == "Coef"){ # log HR
      expressx <- function(x){
        beta + (dab/daa) * (x - alpha)
      }
      
      expressxVect <- Vectorize(expressx)
      if(length(STE) == 1){
        curve(expr = expressxVect, from = from, to = to, n = n, xlab = xlab, ylab = ylab, 
               main = if(is.null(main)) paste("STE = ", round(STE, d), "(HR = ", round(exp(STE), d), ")") else main, ...)
      }else{
        if(length(STE) == 2){
          curve(expr = expressxVect, from = from, to = to, n = n, xlab = xlab, ylab = ylab, 
               main = if(is.null(main)) paste("STE = ", round(max(STE), d), "(HR = ", round(exp(max(STE)), d), 
               ") and min(beta_S) = ", round(min(STE), d), "(HR = ", round(exp(min(STE)), d), ")") else main, ...)
        }else{ 
          curve(expr = expressxVect, from = from, to = to, n = n, xlab = xlab, ylab = ylab, 
                 main = if(is.null(main)) paste("No possible value for STE")
                else main, ...)
        }
      }
      
      #inf
      expressInf <- function(x){
        beta + (dab/daa) * (x - alpha) - qnorm(1-alpha./2) * 
                    sqrt(t(x1) %*% (Vmu + (((x - alpha)/daa)**2) * VD) %*% x1
                         + dbb * (1 - R2trial))
      }
      
      expressInfVect <- Vectorize(expressInf)
      curve (expr = expressInfVect, from = from, to = to, add = T, col = colCI, n = n, lty = lty, ...)
      #sup
      expressSup <- function(x){
        beta + (dab/daa) * (x - alpha) + qnorm(1-alpha./2) * 
          sqrt(t(x1) %*% (Vmu + (((x - alpha)/daa)**2) * VD) %*% x1
               + dbb * (1 - R2trial))
      }
      
      expressSupVect <- Vectorize(expressSup)
      curve (expr = expressSupVect, from = from, to = to, add = T, col = colCI, n = n, lty = lty, ...)
    }
    else{ # HR
      
      expressx <- function(x){
        x <- log(x) # on suppose que les entrees sont des HR et donc on les converti en log HR
        exp(beta + (dab/daa) * (x - alpha))
      }
      
      expressxVect <- Vectorize(expressx)
      
      
      if(length(STE) == 1){
        curve(expr = expressxVect, from = from, to = to, n = n, xlab = xlab, ylab = ylab, 
              main = if(is.null(main)) paste("STE = ", round(STE, d), "(HR = ", round(exp(STE), d), ")") else main, ...)
      }else{
        if(length(STE) == 2){
          curve(expr = expressxVect, from = from, to = to, n = n, xlab = xlab, ylab = ylab, 
                main = if(is.null(main)) paste("STE = ", round(max(STE), d), "(HR = ", round(exp(max(STE)), d), 
                                               ") and min(beta_S) = ", round(min(STE), d), "(HR = ", round(exp(min(STE)), d), ")") else main, ...)
        }else{ 
          curve(expr = expressxVect, from = from, to = to, n = n, xlab = xlab, ylab = ylab, 
                main = if(is.null(main)) paste("No possible value for STE")
                else main, ...)
        }
      }
      
      #inf
      expressInf <- function(x){
        x <- log(x) # on suppose que les entrees sont des HR et donc on les converti en log HR
        exp(beta + (dab/daa) * (x - alpha) - qnorm(1-alpha./2) * 
          sqrt(t(x1) %*% (Vmu + (((x - alpha)/daa)**2) * VD) %*% x1
               + dbb * (1 - R2trial)))
      }
      
      expressInfVect <- Vectorize(expressInf)
      curve (expr = expressInfVect, from = from, to = to, add = T, col = colCI, n = n, lty = lty, ...)
      #sup
      expressSup <- function(x){
        x <- log(x) # on suppose que les entrees sont des HR et donc on les converti en log HR
        exp(beta + (dab/daa) * (x - alpha) + qnorm(1-alpha./2) * 
          sqrt(t(x1) %*% (Vmu + (((x - alpha)/daa)**2) * VD) %*% x1
               + dbb * (1 - R2trial)))
      }
      
      expressSupVect <- Vectorize(expressSup)
      curve (expr = expressSupVect, from = from, to = to, add = T, col = colCI, n = n, lty = lty, ...)
    }
  }
  else{
    if(type == "Coef"){ # log HR
      expressx <- function(x){
        beta + (dab/daa) * (x - alpha)
      }
      
      expressxVect <- Vectorize(expressx)
      
      if(length(STE) == 1){
        curve(expr = expressxVect, from = from, to = to, n = n, xlab = xlab, ylab = ylab, 
              main = if(is.null(main)) paste("STE = ", round(STE, d), "(HR = ", round(exp(STE), d), ")") else main, ...)
      }else{
        if(length(STE) == 2){
          curve(expr = expressxVect, from = from, to = to, n = n, xlab = xlab, ylab = ylab, 
                main = if(is.null(main)) paste("STE = ", round(max(STE), d), "(HR = ", round(exp(max(STE)), d), 
                                               ") and min(beta_S) = ", round(min(STE), d), "(HR = ", round(exp(min(STE)), d), ")") else main, ...)
        }else{ 
          curve(expr = expressxVect, from = from, to = to, n = n, xlab = xlab, ylab = ylab, 
                main = if(is.null(main)) paste("No possible value for STE")
                else main, ...)
        }
      }
      
      #inf
      expressInf <- function(x){
        beta + (dab/daa) * (x - alpha) - qnorm(1-alpha./2) * sqrt(
          dbb * (1 - R2trial))
      }
      
      expressInfVect <- Vectorize(expressInf)
      curve (expr = expressInfVect,
        from = from, to = to, add = T, col = colCI, n = n, lty = lty, ...)
      
      #sup
      expressSup <- function(x){
        beta + (dab/daa) * (x - alpha) + qnorm(1-alpha./2) * sqrt(
          dbb * (1 - R2trial))
      }
      
      expressSupVect <- Vectorize(expressSup)
      curve (expr = expressSupVect,
        from = from, to = to, add = T, col = colCI, n = n, lty = lty, ...)
    }
    else{ # HR
      expressx <- function(x){
        x <- log(x) # on suppose que les entrees sont des HR et donc on les converti en log HR
        exp(beta + (dab/daa) * (x - alpha))
      }
      
      expressxVect <- Vectorize(expressx)
      if(length(STE) == 1){
        curve(expr = expressxVect, from = from, to = to, n = n, xlab = xlab, ylab = ylab, 
              main = if(is.null(main)) paste("STE = ", round(STE, d), "(HR = ", round(exp(STE), d), ")") else main, ...)
      }else{
        if(length(STE) == 2){
          curve(expr = expressxVect, from = from, to = to, n = n, xlab = xlab, ylab = ylab, 
                main = if(is.null(main)) paste("STE = ", round(max(STE), d), "(HR = ", round(exp(max(STE)), d), 
                                               ") and min(beta_S) = ", round(min(STE), d), "(HR = ", round(exp(min(STE)), d), ")") else main, ...)
        }else{ 
          curve(expr = expressxVect, from = from, to = to, n = n, xlab = xlab, ylab = ylab, 
                main = if(is.null(main)) paste("No possible value for STE")
                else main, ...)
        }
      }
      
      #inf
      expressInf <- function(x){
        x <- log(x) # on suppose que les entrees sont des HR et donc on les converti en log HR
        exp(beta + (dab/daa) * (x - alpha) - qnorm(1-alpha./2) * sqrt(
          dbb * (1 - R2trial)))
      }
      
      expressInfVect <- Vectorize(expressInf)
      curve (expr = expressInfVect,
        from = from, to = to, add = T, col = colCI, n = n, lty = lty, ...)
      #sup
      expressSup <- function(x){
        x <- log(x) # on suppose que les entrees sont des HR et donc on les converti en log HR
        exp(beta + (dab/daa) * (x - alpha) + qnorm(1-alpha./2) * sqrt(
          dbb * (1 - R2trial)))
      }
      
      expressSupVect <- Vectorize(expressSup)
      curve (expr = expressSupVect,
        from = from, to = to, add = T, col = colCI, n = n, lty = lty, ...)
    }
  }
  
  #ste 
  
  leg.col = c("black", colCI, "red", "black")
  if(is.null(leg.x)) leg.x = from
  
  if(length(STE) == 0){ # on est dans le cas Delta = 0, pas de solution entire pour cette equation
    # message("Warning : STE does not exist for this intermediate endpoint. Therefore, 
    #         regarding the values of R2trial and Kendall tau, the observed treatment effect on the candidate 
    #         surrogate endpoint can not permitted to predict a non zero treatment effect on true endpoint
    #         using the considered joint surrogate model and the meta-analysis")
      if(type == "HR"){ # log HR
        abline(h = 1, col = "cyan", lty = 3)
        
      }
      else{
        abline(h = 0, col = "cyan", lty = 3)
      }
    if(legend.show == TRUE){
      legend(x = leg.x, y = leg.y, legend = legend[c(1,2)], col = leg.col[c(1,2)], text.col = leg.text.col, 
             lty = leg.lty[c(1,2)], pch = leg.pch[c(1,2)], bg = leg.bg, cex = leg.cex, 
             bty = leg.bty)
      
    }
  }else{
    if(length(STE) == 1){ # une seule solution de l'equation 
      if(type == "HR"){ # log HR
        abline(h = 1, col = "cyan", lty = 3)
        points(exp(STE),1)
        #abline(v = exp(STE), col = "cyan", lty = 4) 
        segments(x0 = exp(STE), y0 = 0, x1 = exp(STE), y1 = 1, col = "cyan", lty = 4)
        if(add.accept.area.betaS == TRUE)
          rect(from, ybottom + 1, exp(STE), ytop + 1, col = "red", density = density, angle = angle)
      }
      else{
        abline(h = 0, col = "cyan", lty = 3)
        points(STE,0)
        #abline(v = STE, col = "cyan", lty = 4)
        segments(x0 = STE, y0 = -6, x1 = STE, y1 = 0, col = "cyan", lty = 4)
        if(add.accept.area.betaS == TRUE)
          rect(from, ybottom, STE, ytop, col = "red", density = density, angle = angle)
      }
    } else{ # on a deux valeurs du STE
      # recherche du sens de la concavite (bref, signe de "a" dans l'equation "ax^2 + bx + c")
      # je prends un pont au hazard dans l'intervalle [x1,x2] et je regarde le signe de son image
      if(f(sum(STE)/2, object = object, var.used = var.used, alpha. = alpha.,
           pred.int.use = pred.int.use) < 0){ # concavite tournee vers le haut
        if(type == "HR"){ # log HR
          abline(h = 1, col = "cyan", lty = 3)
          points(exp(STE[1]),1)
          points(exp(STE[2]),1)
          # abline(v = exp(STE), col = "cyan", lty = 4)
          segments(x0 = exp(STE[1]), y0 = 0, x1 = exp(STE[1]), y1 = 1, col = "cyan", lty = 4)
          segments(x0 = exp(STE[2]), y0 = 0, x1 = exp(STE[2]), y1 = 1, col = "cyan", lty = 4)
          if(add.accept.area.betaS == TRUE)
            rect(exp(STE[1]), ybottom + 1, exp(STE[2]), ytop + 1, col = "red", density = density, angle = angle)
        }
        else{
          abline(h = 0, col = "cyan", lty = 3)
          points(STE[1],0)
          points(STE[2],0)
          # abline(v = STE, col = "cyan", lty = 4)
          segments(x0 = STE[1], y0 = -6, x1 = STE[1], y1 = 0, col = "cyan", lty = 4)
          segments(x0 = STE[2], y0 = -6, x1 = STE[2], y1 = 0, col = "cyan", lty = 4)
          if(add.accept.area.betaS == TRUE)
            rect(STE[1], ybottom, STE[2], ytop, col = "red", density = density, angle = angle)
        }
      }else{ # concavite tournee vers le bas
        if(type == "HR"){ # log HR
          abline(h = 1, col = "cyan", lty = 3)
          points(exp(STE[1]),1)
          points(exp(STE[2]),1)
          # abline(v = exp(STE), col = "cyan", lty = 4)
          segments(x0 = exp(STE[1]), y0 = 0, x1 = exp(STE[1]), y1 = 1, col = "cyan", lty = 4)
          segments(x0 = exp(STE[2]), y0 = 0, x1 = exp(STE[2]), y1 = 1, col = "cyan", lty = 4)
          if(add.accept.area.betaS == TRUE){
            rect(from, ybottom + 1, exp(STE[1]), ytop + 1, col = "cyan", density = density, angle = angle)
            rect(exp(STE[2]), ybottom + 1, to, ytop + 1, col = "red", density = density, angle = angle)
          }
        }
        else{
          abline(h = 0, col = "cyan", lty = 3)
          points(STE[1],0)
          points(STE[2],0)
          # abline(v = STE, col = "cyan", lty = 4)
          segments(x0 = STE[1], y0 = -6, x1 = STE[1], y1 = 0, col = "cyan", lty = 4)
          segments(x0 = STE[2], y0 = -6, x1 = STE[2], y1 = 0, col = "cyan", lty = 4)
          if(add.accept.area.betaS == TRUE){
            rect(from, ybottom, STE[1], ytop, col = "cyan", density = density, angle = angle)
            rect(STE[2], ybottom, to, ytop, col = "red", density = density, angle = angle)
          }
        }
      }
    }
    if(legend.show == TRUE){
      if(add.accept.area.betaS == TRUE)
        legend(x = leg.x, y = leg.y, legend = legend, col = leg.col, text.col = leg.text.col, 
               lty = leg.lty, pch = leg.pch, bg = leg.bg, cex = leg.cex, 
               bty = leg.bty)
      else
        legend(x = leg.x, y = leg.y, legend = legend[c(1, 2, 4)], col = leg.col[c(1, 2, 4)], text.col = leg.text.col, 
               lty = leg.lty[c(1, 2, 4)], pch = leg.pch[c(1, 2, 4)], bg = leg.bg, cex = leg.cex, 
               bty = leg.bty)
      
    }
  }
  
}