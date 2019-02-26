##' Predict Method for the one-step Joint surrogate model for the evaluation of a 
##' canditate surrogate endpoint.
##' 
##' @description{
##' Predict the treatment effect on the true endpoint (\eqn{beta_T}), based on the 
##' treatment effect observed on the surrogate endpoint (\eqn{beta_S}).
##' }
##' 
##' @details{
##' Prediction is based on the formulas described in (Burzikwosky \emph{et al.}, 2006).
##' We do not consider the case of prediction which suppose estimation error on 
##' the estimate of the treatment effect on the surrogate endpoint in the new trial.
##' }
##' @aliases predict.jointSurroPenal 
##' @usage
##' 
##' \method{predict}{jointSurroPenal}(object, datapred = NULL, var.used = "error.estim", 
##' alpha. = 0.05, dec = 3, ...)
##' @param object An object inheriting from \code{jointSurroPenal} class
##' (output from calling the function \code{jointSurroPenal}).
##' @param datapred Dataset to use for the prediction. If this argument is specified,
##' the data structure must be the same as the parameter \code{data} in the 
##' function \link{jointSurroPenal}. However, if observation on te true endpoint are
##' not available, columns timeT and \code{statusT} can be absent.
##' @param var.used This argument can take two values. The first one is \code{"error.estim"}
##' and indicates if the prediction error take into account
##' the estimation error of the estimates of the parameters. If the estimates 
##' are supposed to be known or if the dataset includes a high number of trials with 
##' a high number of subject per trial, value \code{No.error} can be used. 
##' The default is \code{error.estim}.
##' @param alpha. The confidence level for the prediction interval. The default is \code{0.05}
##' @param dec The desired number of digits after the decimal point for parameters
##' and confidence intervals. Default of 3 digits is used.
##' @param ... other unused arguments.
##' 
##' @return Return and display a dataframe including for each trial the number of included subjects, 
##' the observed 
##' treatment effect on the surrogate endpoint, the observed treatment effect on
##' the true endpoint (if available) and the predicted treatment effect on the 
##' true enpoint with the associated prediction intervals. If the observed treatment effect on the true 
##' endpoint (if available) is included into the prediction interval, the last columns contains "*".
##' @seealso \code{\link{jointSurroPenal}}
##' 
##' @author Casimir Ledoux Sofeu \email{casimir.sofeu@u-bordeaux.fr}, \email{scl.ledoux@gmail.com} and 
##' Virginie Rondeau \email{virginie.rondeau@inserm.fr}
##' 
##' @references 
##' Burzykowski T, Buyse M (2006). "Surrogate threshold effect: an alternative 
##' measure for meta-analytic surrogate endpoint validation." Pharmaceutical 
##' Statistics, 5(3), 173-186.ISSN 1539-1612.
##' 
##' @keywords surrogate prediction
##' @export
##' 
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
##' # prediction of the treatment effects on the true endpoint in each trial of 
##' the dataOvarian dataset
##' # predict(joint.surro.ovar)
##' 
##' }
##' 
##' 
"predict.jointSurroPenal" <- function (object, datapred = NULL, var.used = "error.estim", alpha. = 0.05
                                       , dec = 3, ...)
{
  if (!inherits(object, "jointSurroPenal"))
    stop("object must be of class 'jointSurroPenal'")
  
  if(! var.used %in% c("error.estim","No.error"))
    stop("Argument 'var.used' must be specified to 'error.estim' or 'No.error' ")
  
  if(is.null(datapred)){ # we used the dataset from the model
    data <- object$data
  }
  else{
    # ================ data checking=======================
      data <- datapred
      # The initial followup time. The default value is 0
      data$initTime <- 0 
      # dataset's names control
      varStatus=(c("initTime","timeS","statusS","trialID","patienID","trt") %in% names(data))
      if(F %in% varStatus){
        stop("Control the names of your variables. They must contain at leat 5 variables named: timeS, statusS, trialID, patienID and trt. seed the help on this function")
      }
      
      # traitement des donnees
      if(max(table(data$patienID)) > 1){
        stop("Control your dataset. You probably have a duplicate on individual (patienID variable)")
      }
      
      if(!is.numeric(data$timeS)|!is.numeric(data$trialID)){
        stop("The variables timeS, and trialID must be numeric") 
      }
      
      if(F %in% c(levels(as.factor(as.character(data$statusS))) %in% c(0,1),levels(as.factor(as.character(data$trt))) %in% c(0,1))){
        stop("The variables statusS, and trt must be coded 0 or 1")
      }
      if(T %in% c((data$timeS - data$initTime) <= 0)){
        stop("Controll the follow up times of your sujects. the is at leat one subjects with intTime > timeS ")
      }
    # ====================== end data checking ==============================
  }
  
  trial <- unique(data$trialID)
  trialtable <- table(data$trialID)
  
  if(F %in% (c("timeT","statusT") %in% names(data))){
    matrixPred <- data.frame(matrix(0, nrow = length(trial), ncol = 6))
    names(matrixPred) <- c("trialID","ntrial","bata.S", "beta.T.i", "Inf.95.CI", "Sup.95.CI" )
    matrixPred$trialID <- trial
    for(i in 1:length(trial)){
      subdata <- data[data$trialID == trial[i],]
      matrixPred$beta.S[i] <- coxph(Surv(timeS, statusS) ~ trt, subdata)$coefficients
    }
  }
  else{
    matrixPred <- data.frame(matrix(0, nrow = length(trial), ncol = 8))
    names(matrixPred) <- c("trialID","ntrial","beta.S", "beta.T", "beta.T.i", "Inf.95.CI", "Sup.95.CI","" )
    matrixPred$trialID <- trial
    for(i in 1:length(trial)){
      subdata <- data[data$trialID == trial[i],]
      matrixPred$beta.S[i] <- coxph(Surv(timeS, statusS) ~ trt, subdata)$coefficients
      matrixPred$beta.T[i] <- coxph(Surv(timeT, statusT) ~ trt, subdata)$coefficients
    }
  }
  
  for(i in 1:length(trial)){
    beta  <- object$beta.t
    dab   <- object$Coefficients$Estimate[nrow(object$Coefficients)-4]
    daa   <- object$Coefficients$Estimate[nrow(object$Coefficients)-6]
    dbb   <- object$Coefficients$Estimate[nrow(object$Coefficients)-5]
    alpha <- object$beta.s
    alpha0 <- matrixPred$beta.S[i]
    x     <- t(matrix(c(1, -dab/daa),1,2))
    
    # Vmu (sigma_ST, sigma_SS). on utilise la matrice obtenu par delta methode a partir de la hessienne
    Vmu   <- matrix(c(object$varcov.Sigma[3,3], object$varcov.Sigma[3,1], 
                      object$varcov.Sigma[3,1], object$varcov.Sigma[1,1]),2,2)
    nparam <- nrow(object$varH)
    
    # VD (bete_T, beta_S). on utilise la hesienne directement car pas de changement de variable
    VD     <- matrix(c(object$varH[nparam,nparam], object$varH[nparam,nparam-1],
                       object$varH[nparam-1,nparam], object$varH[nparam -1,nparam - 1]),2,2)
    R2trial <- object$Coefficients$Estimate[nrow(object$Coefficients)-1] 
    matrixPred$beta.T.i[i] <- beta + (dab/daa) * (alpha0 - alpha)
    variance.inf <- dbb * (1 - R2trial) 
    variance.N <- t(x) %*% (Vmu + (((alpha0 - alpha)/daa)**2) * VD) %*% x
    + variance.inf
    
    if(var.used == "error.estim") 
      variance <- variance.N
    else 
      variance <- variance.inf
    
    matrixPred$Inf.95.CI[i] <- matrixPred$beta.T.i[i] - qnorm(1-alpha./2) * sqrt(variance)
    matrixPred$Sup.95.CI[i] <- matrixPred$beta.T.i[i] + qnorm(1-alpha./2) * sqrt(variance)
    
    # matrixPred[i,-c(1,2)] <- round(matrixPred[i,-c(1,2)],dec)
    # ajout du nombre d'essais
    # matrixPred$ntrial[i] <- trialtable[matrixPred[i,1]]
    matrixPred$ntrial[i] <- trialtable[i]
    
    # je mets une "*" si la valeur observee est incluse dans l'intervalle de prediction
    if(!(F %in% (c("timeT","statusT") %in% names(data)))){
      if((matrixPred$beta.T[i] >= matrixPred$Inf.95.CI[i]) & (matrixPred$beta.T[i] <= matrixPred$Sup.95.CI[i]))
        matrixPred[i,ncol(matrixPred)] <- "*"
      else
        matrixPred[i,ncol(matrixPred)] <- " "
    }
  }
  matrixPred[,-c(1,2,8)] <- round(matrixPred[,-c(1,2,8)],dec)
  
  
 # print(matrixPred)
  return(matrixPred)

}
