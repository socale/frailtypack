##' Plot of leave-one-out crossvalidation Outputs from the one-step Joint surrogate model for evaluating a 
##' canditate surrogate endpoint.
##' 
##' @description{
##' Plot of leave-one-out crossvalidation Outputs for evaluating the joint surrogate model 
##' }
##' 
##' @aliases plot.jointSurroPenalloocv
##' @usage 
##' 
##' \method{plot}{jointSurroPenalloocv}(x, unusedtrial = NULL, xleg = "bottomleft", 
##' yleg = NULL, main = NULL, ...)
##' @param x An object inherent from the \code{jointSurroPenalloocv} Class
##' @param unusedtrial Vector of unconsidered trials, may be due to the fact that the 
##' predicted treatment effects on true endpoint have an outlier. In this case, 
##' one can drop from the data the trials with very hight absolute predicted value 
##' @param xleg X-coordinate for the location of the legend
##' @param yleg Y-coordinate for the location of the legend, the default is \code{NULL}
##' @param main The desired main
##' @param ... other unused arguments.
##' 
##' @return This function displays the boxplots corresponding to the number of trials in the 
##' dataset. Each boxplot includes 3 elements correnponding to the predicted treatment effect on true endpoint
##' with the prediction interval. The circle inside or outside the boxplot represents the observed
##' treatment effect on true endpoint. For all trials with convergence issues or outliers, the boxplot is just represented
##' by a dash. In the last case, we display in the main of the figure a vector of these trials, if argument \code{main} 
##' is set to \code{NULL}. The function returns the list of unused trials.
##' @seealso \code{\link{loocv}}
##' 
##' @author Casimir Ledoux Sofeu \email{casimir.sofeu@u-bordeaux.fr}, \email{scl.ledoux@gmail.com} and 
##' Virginie Rondeau \email{virginie.rondeau@inserm.fr}
##' 
##' @references 
##' Burzykowski T, Buyse M (2006). "Surrogate threshold effect: an alternative 
##' measure for meta-analytic surrogate endpoint validation." Pharmaceutical 
##' Statistics, 5(3), 173-186.ISSN 1539-1612.
##' 
##' @importFrom graphics boxplot points
##' @keywords surrogate prediction loocv plot
##' @export
##' @examples
##' 
##' 
##' \dontrun{
##' library(frailtypack)
##' data(dataOvarian)
##' 
##' joint.surro.Gumbel <- jointSurroCopPenal(data = dataOvarian, int.method = 0, 
##'                       n.knots = 8, maxit=50, kappa.use = 4, nb.mc = 1000, 
##'                       typecopula = 2, print.iter = T, scale = 1/365)
##' summary(joint.surro.Gumbel)
##' 
##' loocv.result <- loocv(joint.surro.Gumbel)
##' loocv.result
##' 
##' plot(x = loocv.result, unusedtrial = c(22, 30, 33, 38, 42, 47, 49), 
##'                           xleg = "bottomleft", y = NULL)
##' }
##' 
"plot.jointSurroPenalloocv" <- function(x, unusedtrial = NULL, xleg = "bottomleft", 
                                        yleg = NULL, main = NULL, ...){
    object <- x
    data.gumbel <- object$result
    data.gumbel <- data.gumbel[!(data.gumbel$trialID %in% unusedtrial),]
    data.plot <- data.frame(matrix(NA, nrow = 3 * nrow(data.gumbel), ncol = 4))
    names(data.plot) <- c("trialID", "beta.T", "ordre", "beta")
    data.plot$ordre <- rep(c(0,1,2),nrow(data.gumbel))
    
    # remove from unusedtrial all trials allready in the vector object$notconvtrial  
    unusedtrial <- unusedtrial[!c(unusedtrial %in% intersect(unusedtrial, object$notconvtrial))]
    
    j <- 1
    for(i in 1:nrow(data.gumbel)){
      if(data.plot$trialID[i] %in% unusedtrial){
        data.plot$trialID[j:(j+2)] <- data.gumbel$trial[i]
        data.plot$beta.T[j:(j+2)] <- 0
        data.plot$beta[j:(j+2)] <- c(0,0,0)
      }else{
        data.plot$trialID[j:(j+2)] <- data.gumbel$trial[i]
        data.plot$beta.T[j:(j+2)] <- data.gumbel$beta.T[i]
        data.plot$beta[j:(j+2)] <- as.numeric(data.gumbel[i,c("beta.T.i", "Inf.95.CI", "Sup.95.CI")])
      }
      j = j+3
    }
    
    # I allow 0 to all trials belonging in the vector unusedtrial
    data.plot2 <- data.frame(matrix(NA, nrow = 3 * object$ntrial, ncol = 4))
    names(data.plot2) <- c("trialID", "beta.T", "ordre", "beta")
    data.plot2$ordre <- rep(c(0,1,2),object$ntrial)
    for(i in 1:nrow(data.plot))
      data.plot2[i,] <- data.plot[i,] 
    k<-1
    while (k <= length(unusedtrial)){
      data.plot2$trialID[j:(j+2)] <- unusedtrial[k]
      data.plot2$beta.T[j:(j+2)] <- 0
      data.plot2$beta[j:(j+2)] <- c(0,0,0)
      j <- j+3
      k <- k+1
    }
    
    k<-1
    while (k <= length(object$notconvtrial)){
      data.plot2$trialID[j:(j+2)] <- object$notconvtrial[k]
      data.plot2$beta.T[j:(j+2)] <- 0
      data.plot2$beta[j:(j+2)] <- c(0, 0,0)
      j <- j+3
      k <- k+1
    }
    
    doBy::orderBy(~ trialID+ordre,data.plot2)
    trialnotused <- sort(c(unusedtrial, object$notconvtrial))
    if(length(trialnotused)>=1) mainlabel <- trialnotused[1]
    if(length(trialnotused)>1){
      for(i in 2:length(trialnotused))
        mainlabel <- paste(mainlabel,trialnotused[i], sep = ", ")
    }
    
    if(is.null(main)){
      if(length(mainlabel) > 0){ # to avoid to display the main in cases where all trials have been used
        boxplot(data.plot2$beta ~ data.plot2$trialID, xlab = "Trials", 
                ylab = "Log Hazard ration of the true endpoint",
                main = paste("Unused trials = ", mainlabel, sep = ""))
      }
      else{
        boxplot(data.plot2$beta ~ data.plot2$trialID, xlab = "Trials",
                ylab = "Log Hazard ration of the true endpoint")
      }
    }else{
      boxplot(data.plot2$beta ~ data.plot2$trialID, xlab = "Trials", 
              ylab = "Log Hazard ration of the true endpoint",
              main = main)
    }
    points(data.plot2$trialID[!(data.plot2$trialID %in% unusedtrial)],data.plot2$beta.T[!(data.plot2$trialID %in% unusedtrial)])
    legend(x = xleg, y = yleg, c("Beta observed", "Beta predict"), cex = 1, pch= c(1,15))
    abline(h = 0)
    return(paste("Unused trials = ", mainlabel, sep = ""))
  }