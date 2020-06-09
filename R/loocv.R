##' The trials leave-one-out crossvalidation for the one-step Joint surrogate model for evaluating a 
##' canditate surrogate endpoint.
##' 
##' @description{
##' The trials leave-one-out crossvalidation for evaluating the joint surrogate model 
##' }
##' 
##' @aliases loocv 
##' @usage
##' 
##' loocv(object, unusedtrial, var.used = "error.estim", alpha. = 0.05, 
##' dec = 3, print.times = TRUE)
##' 
##' @param object An object inheriting from \code{jointSurroPenal} class
##' (output from calling the function \code{jointSurroPenal} or \code{jointSurroCopPenal}).
##' @param unusedtrial A list of trial not to be taken into account in the cross-validation.
##' This parameter is useful when after excluding some trials, the model is facing 
##' convergence problem.
##' @param var.used This argument takes two values. The first one is \code{"error.estim"}
##' and indicates if the prediction variance takes into account
##' the estimation errors from the estimates of the parameters. If estimates 
##' are supposed to be known or if the dataset includes a high number of trials with 
##' a high number of subject per trial, value \code{"No.error"} can be used. 
##' The default is \code{error.estim}.
##' @param alpha. The confidence level for the prediction interval. The default is \code{0.05}
##' @param dec The desired number of digits after the decimal point for parameters
##' and confidence intervals. Default of 3 digits is used.
##' @param print.times a logical parameter to print estimation time. Default is TRUE.
##' 
##' @return This function returns an object of class \code{jointSurroPenalloocv} containing:
##' \item{result}{A dataframe 
##' including for each trial the number of included subjects, the observed 
##' treatment effect on the surrogate endpoint, the observed treatment effect on
##' the true endpoint and the predicted treatment effect on the 
##' true enpoint with the associated prediction intervals. If the observed treatment effect on the true 
##' endpoint is included into the prediction interval, the last columns contains "*".} 
##' \item{ntrial}{The number of trials in the meta-analysis}
##' \item{notconvtrial}{The vector of trials that have not converged}
##' \item{pred.error}{The prediction error, corresponding to the number of cases where the prediction interval does not included the observed treatment effect on T}
##' \item{different.models}{The list of the \code{G} models obtained after excuded for the \code{i-th} trial}
##' \item{loocv.summary}{A dataframe of the estimates for the \code{G} models; each raw including the results
##'  without the subjects of the given trial}
##' @seealso \code{\link{jointSurroPenal}, \link{jointSurroCopPenal}}
##' 
##' @author Casimir Ledoux Sofeu \email{casimir.sofeu@u-bordeaux.fr}, \email{scl.ledoux@gmail.com} and 
##' Virginie Rondeau \email{virginie.rondeau@inserm.fr}
##' 
##' @references 
##' Burzykowski T, Buyse M (2006). "Surrogate threshold effect: an alternative 
##' measure for meta-analytic surrogate endpoint validation." Pharmaceutical 
##' Statistics, 5(3), 173-186.ISSN 1539-1612.
##' 
##' @keywords surrogate prediction loocv
##' @export
##' @examples
##' 
##' 
##' \dontrun{
##' # Generation of data to use 
##'  data.sim <- jointSurrSimul(n.obs=300, n.trial = 10,cens.adm=549.24,
##'              alpha = 1.5, theta = 3.5, gamma = 2.5, zeta = 1, sigma.s = 0.7,
##'              sigma.t = 0.7, cor = 0.8, betas = -1.25, betat = -1.25,
##'              full.data = 0, random.generator = 1, seed = 0,
##'              nb.reject.data = 0)
##' 
##' ###--- Joint surrogate model ---###
##'  
##' joint.surro.sim.MCGH <- jointSurroPenal(data = data.sim, int.method = 2,
##'                         nb.mc = 300, nb.gh = 20, print.iter = F)
##' 
##' # Example of loocv taking into accountn ony trial 2 trials (1 and 3)
##' dloocv <- loocv(joint.surro.sim.MCGH, unusedtrial = c(2,4:10))
##' dloocv$result
##' dloocv$loocv.summary
##' 
##' # In order to summarize all the estimated models during the loocv proccess:
##' dloocv$different.models
##' 
##' }
##' 
##' 
loocv <- function (object, unusedtrial = NULL, var.used = "error.estim", alpha. = 0.05,
                   dec = 3, print.times = TRUE)
{
  if (!inherits(object, "jointSurroPenal"))
    stop("object must be of class 'jointSurroPenal'")
  
  if(! var.used %in% c("error.estim","No.error"))
    stop("Argument 'var.used' must be specified to 'error.estim' or 'No.error' ")
  
  # gestion de l'affichage a l'ecran
  flush.console()
  if (print.times){
    ptm<-proc.time()
    cat("\n")
    cat("Be patient. The program is computing ... \n")
  }
  
    
  dataUse <- object$data
  
  trial <- unique(dataUse$trialID)
  # init of the result
  d <- data.frame(matrix(rep(NA,8), nrow = 1, ncol = 8))[-1,]
  names(d) <- c("trialID","ntrial","beta.S", "beta.T", "beta.T.i", "Inf.95.CI", "Sup.95.CI","" )
  notconvtrial = unusedtrial
  trialused = NULL
  lloocv = list()
  for(i in 1:length(trial)){
    if(!(i %in% unusedtrial)){ # one can identifie trials that pose problem when they are removed, and then ignore them
      message("Trial: ", i)
      dataUseloo <- dataUse[!(dataUse$trialID %in% trial[i]),]
      if(object$type.joint == 1){ # joint surrogate model
        # Estimation
        if(!is.na(object$parameter["init.kappa1"])){
          joint.surro <- jointSurroPenal(dataUseloo, maxit = object$parameter["maxit"][[1]],indicator.zeta = object$parameter["indicator.zeta"][[1]], 
                        indicator.alpha = object$parameter["indicator.alpha"][[1]], frail.base = object$parameter["frail.base"][[1]], 
                        n.knots = object$parameter["n.knots"][[1]], LIMparam = object$parameter["LIMparam"][[1]], LIMlogl = object$parameter["LIMlogl"][[1]], 
                        LIMderiv = object$parameter["LIMderiv"][[1]], nb.mc = object$parameter["nb.mc"][[1]], nb.gh = object$parameter["nb.gh"][[1]], 
                        nb.gh2 = object$parameter["nb.gh2"][[1]], adaptatif = object$parameter["adaptatif"][[1]], 
                        int.method = object$parameter["int.method"][[1]], nb.iterPGH = object$parameter["nb.iterPGH"][[1]], 
                        nb.MC.kendall = object$parameter["nb.MC.kendall"][[1]], nboot.kendall = object$parameter["nboot.kendall"][[1]], 
                        true.init.val = object$parameter["true.init.val"][[1]], theta.init = object$parameter["theta.init"][[1]], 
                        sigma.ss.init = object$parameter["sigma.ss.init"][[1]], sigma.tt.init = object$parameter["sigma.tt.init"][[1]], 
                        sigma.st.init = object$parameter["sigma.st.init"][[1]], gamma.init = object$parameter["gamma.init"][[1]], 
                        alpha.init = object$parameter["alpha.init"][[1]], zeta.init = object$parameter["zeta.init"][[1]], 
                        betas.init = object$parameter["betas.init"][[1]], betat.init = object$parameter["betat.init"][[1]], 
                        scale = object$parameter["scale"][[1]], random.generator = object$parameter["random.generator"][[1]], 
                        kappa.use = object$parameter["kappa.use"][[1]], random = object$parameter["random"][[1]], 
                        random.nb.sim = object$parameter["random.nb.sim"][[1]], seed = object$parameter["seed"][[1]], 
                        init.kappa = c(object$parameter["init.kappa1"][[1]],object$parameter["init.kappa2"][[1]]), 
                        nb.decimal = object$parameter["nb.decimal"][[1]], print.times = object$parameter["print.times"][[1]], 
                        print.iter = object$parameter["print.iter"][[1]])
        }
        
        if(is.na(object$parameter["init.kappa1"])){
          joint.surro <- jointSurroPenal(dataUseloo, maxit = object$parameter["maxit"][[1]],indicator.zeta = object$parameter["indicator.zeta"][[1]], 
                         indicator.alpha = object$parameter["indicator.alpha"][[1]], frail.base = object$parameter["frail.base"][[1]], 
                         n.knots = object$parameter["n.knots"][[1]], LIMparam = object$parameter["LIMparam"][[1]], LIMlogl = object$parameter["LIMlogl"][[1]], 
                         LIMderiv = object$parameter["LIMderiv"][[1]], nb.mc = object$parameter["nb.mc"][[1]], nb.gh = object$parameter["nb.gh"][[1]], 
                         nb.gh2 = object$parameter["nb.gh2"][[1]], adaptatif = object$parameter["adaptatif"][[1]], 
                         int.method = object$parameter["int.method"][[1]], nb.iterPGH = object$parameter["nb.iterPGH"][[1]], 
                         nb.MC.kendall = object$parameter["nb.MC.kendall"][[1]], nboot.kendall = object$parameter["nboot.kendall"][[1]], 
                         true.init.val = object$parameter["true.init.val"][[1]], theta.init = object$parameter["theta.init"][[1]], 
                         sigma.ss.init = object$parameter["sigma.ss.init"][[1]], sigma.tt.init = object$parameter["sigma.tt.init"][[1]], 
                         sigma.st.init = object$parameter["sigma.st.init"], gamma.init = object$parameter["gamma.init"][[1]], 
                         alpha.init = object$parameter["alpha.init"][[1]], zeta.init = object$parameter["zeta.init"][[1]], 
                         betas.init = object$parameter["betas.init"][[1]], betat.init = object$parameter["betat.init"][[1]], 
                         scale = object$parameter["scale"][[1]], random.generator = object$parameter["random.generator"][[1]], 
                         kappa.use = object$parameter["kappa.use"][[1]], random = object$parameter["random"][[1]], 
                         random.nb.sim = object$parameter["random.nb.sim"][[1]], seed = object$parameter["seed"][[1]], 
                         init.kappa = NULL, 
                         nb.decimal = object$parameter["nb.decimal"][[1]], print.times = object$parameter["print.times"][[1]], 
                         print.iter = object$parameter["print.iter"][[1]])
        }
      } else{ # joint frailty copula model
          if(!is.na(object$parameter["init.kappa1"])){
            joint.surro <- jointSurroCopPenal(dataUseloo, maxit = object$parameter["maxit"][[1]], 
                             indicator.alpha = object$parameter["indicator.alpha"][[1]], frail.base = object$parameter["frail.base"][[1]], 
                             n.knots = object$parameter["n.knots"][[1]], LIMparam = object$parameter["LIMparam"][[1]], LIMlogl = object$parameter["LIMlogl"][[1]], 
                             LIMderiv = object$parameter["LIMderiv"][[1]], nb.mc = object$parameter["nb.mc"][[1]], nb.gh = object$parameter["nb.gh"][[1]], 
                             nb.gh2 = object$parameter["nb.gh2"][[1]], adaptatif = object$parameter["adaptatif"][[1]], 
                             int.method = object$parameter["int.method"][[1]], nb.iterPGH = object$parameter["nb.iterPGH"][[1]], 
                             #nboot.kendall = object$parameter["nboot.kendall"][[1]], 
                             true.init.val = object$parameter["true.init.val"][[1]], thetacopula.init = object$parameter["theta.init"][[1]], 
                             sigma.ss.init = object$parameter["sigma.ss.init"][[1]], sigma.tt.init = object$parameter["sigma.tt.init"][[1]], 
                             sigma.st.init = object$parameter["sigma.st.init"][[1]], gamma.init = object$parameter["gamma.init"][[1]], 
                             alpha.init = object$parameter["alpha.init"][[1]], 
                             betas.init = object$parameter["betas.init"][[1]], betat.init = object$parameter["betat.init"][[1]], 
                             scale = object$parameter["scale"], random.generator = object$parameter["random.generator"][[1]], 
                             kappa.use = object$parameter["kappa.use"][[1]], random = object$parameter["random"][[1]], 
                             random.nb.sim = object$parameter["random.nb.sim"][[1]], seed = object$parameter["seed"][[1]], 
                             init.kappa = c(object$parameter["init.kappa1"][[1]],object$parameter["init.kappa2"][[1]]),
                             typecopula = object$parameter["typecopula"][[1]],
                             nb.decimal = object$parameter["nb.decimal"][[1]], print.times = object$parameter["print.times"][[1]], 
                             print.iter = object$parameter["print.iter"][[1]])                 
          }
          
          if(is.na(object$parameter["init.kappa1"])){
            joint.surro <- jointSurroCopPenal(dataUseloo, maxit = object$parameter["maxit"][[1]], 
                             indicator.alpha = object$parameter["indicator.alpha"][[1]], frail.base = object$parameter["frail.base"][[1]], 
                             n.knots = object$parameter["n.knots"][[1]], LIMparam = object$parameter["LIMparam"][[1]], LIMlogl = object$parameter["LIMlogl"][[1]], 
                             LIMderiv = object$parameter["LIMderiv"][[1]], nb.mc = object$parameter["nb.mc"][[1]], nb.gh = object$parameter["nb.gh"][[1]], 
                             nb.gh2 = object$parameter["nb.gh2"][[1]], adaptatif = object$parameter["adaptatif"][[1]], 
                             int.method = object$parameter["int.method"][[1]], nb.iterPGH = object$parameter["nb.iterPGH"][[1]], 
                             #nboot.kendall = object$parameter["nboot.kendall"][[1]], 
                             true.init.val = object$parameter["true.init.val"][[1]], thetacopula.init = object$parameter["theta.init"][[1]], 
                             sigma.ss.init = object$parameter["sigma.ss.init"][[1]], sigma.tt.init = object$parameter["sigma.tt.init"][[1]], 
                             sigma.st.init = object$parameter["sigma.st.init"][[1]], gamma.init = object$parameter["gamma.init"][[1]], 
                             alpha.init = object$parameter["alpha.init"][[1]], 
                             betas.init = object$parameter["betas.init"][[1]], betat.init = object$parameter["betat.init"][[1]], 
                             scale = object$parameter["scale"][[1]], random.generator = object$parameter["random.generator"][[1]], 
                             kappa.use = object$parameter["kappa.use"][[1]], random = object$parameter["random"][[1]], 
                             random.nb.sim = object$parameter["random.nb.sim"][[1]], seed = object$parameter["seed"][[1]], 
                             init.kappa = NULL, typecopula = object$parameter["typecopula"][[1]], 
                             nb.decimal = object$parameter["nb.decimal"][[1]], print.times = object$parameter["print.times"][[1]], 
                             print.iter = object$parameter["print.iter"][[1]])              
          }
        }
    }else{
      joint.surro =NULL
    }
      
    # Prediction
    if(is.null(joint.surro)){ 
      if(!(i %in% unusedtrial)) 
        {
        cat(c("===Model without trial", i, "did not converge===: \n"))
        notconvtrial[length(notconvtrial)+1] <- i
      }
    }else{
      trialused[length(trialused)+1] <- i
      d1 <- predict.jointSurroPenal(joint.surro,datapred = dataUse[dataUse$trialID %in% trial[i],], dec = dec)
      # Merger of the results
      d <- rbind(d,d1)
      # impression du temps de calcul
      if (print.times){
        cost<-(proc.time()-ptm)/60
        cat("The program took", round(cost[3],2), "minutes \n")
      }
      lloocv[[length(lloocv)+1]] <- joint.surro
    }
  }
  
  if(!is.null(d)){
    result <- NULL
    result$result <- d
    result$ntrial <- length(trial)
    result$notconvtrial <- notconvtrial
    result$pred.error <- round(prop.table(table(result$result[,ncol(result$result)]))[1],dec)
    result$different.models <- lloocv
    result$loocv.summary <- loocv.summary(loocv.object = result, trialused = trialused,
    #result$loocv.summary <- loocv.summary(loocv.object = result,
                                          nb.parameters = nrow(object$Coefficients),
                                          names.parameters = rownames(object$Coefficients))
    class(result) <- "jointSurroPenalloocv"
  }
  
  return(result)
}
