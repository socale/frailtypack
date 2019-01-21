##' Leave-one-out crossvalidation for the one-step Joint surrogate model for the evaluation of a 
##' canditate surrogate endpoint.
##' 
##' @description{
##' Leave-one-out crossvalidation for the evaluation of the joint surrogate model 
##' }
##' 
##' @details{
##' 
##' }
##' @aliases loocv 
##' @usage
##' 
##' loocv(object, var.used = "error.meta", alpha. = 0.05, print.times = T)
##' 
##' @param object An object inheriting from \code{jointSurroPenal} class
##' (output from calling \code{jointSurroPenal} function).
##' @param var.used This argument takes two values. The first one is \code{"error.meta"}
##' and indicates if the prediction error take into account
##' the estimation error of the estimates of the parameters. If the estimates 
##' are suppose knew or if the dataset includes a high number of trials with 
##' a high number of subject per trial, value \code{No.error} can be used. 
##' The default is \code{error.meta}.
##' @param alpha. The confidence level for the prediction interval. The default is \code{0.05}
##' @param print.times a logical parameter to print estimation time. Default is TRUE.
##' 
##' @return Returns an object of class \code{jointSurroPenalloocv} containing a dataframe 
##' including for each trial the observed 
##' treatment effect on the surrogate endpoint, the observed treatment effect on
##' the true endpoint and the predicted treatment effect on the 
##' true enpoint with the associated prediction intervalls. If the observed treatment effect on the true 
##' endpoint is included into the prediction interval, the last columns contains "*".
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
##' @keywords surrogate prediction loocv
##' @export
##' @examples
##' 
##' 
##' \dontrun{
##' # Generation of data to use 
##'  data.sim <- jointSurrSimul(n.obs=600, n.trial = 30,cens.adm=549.24, 
##'          alpha = 1.5, theta = 3.5, gamma = 2.5, zeta = 1, sigma.s = 0.7, 
##'          sigma.t = 0.7, rsqrt = 0.8, betas = -1.25, betat = -1.25, 
##'          full.data = 0, random.generator = 1, seed = 0, nb.reject.data = 0)
##' 
##' ###--- Joint surrogate model ---###
##'  (Computation takes around 5 minutes)
##' 
##' joint.surro.sim.MCGH <- jointSurroPenal(data = data.sim, int.method = 2, 
##'                    nb.mc = 300, nb.gh = 20)
##'                 
##' dloocv <- loocv(joint.surro.sim.MCGH)
##' 
##' }
##' 
##' 
loocv <- function (object, var.used = "error.meta", alpha. = 0.05, print.times = T)
{
  if (!inherits(object, "jointSurroPenal"))
    stop("object must be of class 'jointSurroPenal'")
  
  if(! var.used %in% c("error.meta","No.error"))
    stop("Argument 'var.used' must be specified to 'error.meta' or 'No.error' ")
  
  # gestion de l'affichage a l'ecran
  flush.console()
  if (print.times){
    ptm<-proc.time()
    cat("\n")
    cat("Be patient. The program is computing ... \n")
  }
  
    
  dataUse <- object$data
  
  trial <- unique(dataUse$trialID)
  d <- NULL
  for(i in 1:length(trial)){
    dataUseloo <- dataUse[!(dataUse$trialID %in% trial[i]),]
    # Estimation
    joint.surro <- jointSurroPenal(dataUseloo, maxit = object$parameter["maxit"],indicator.zeta = object$parameter["indicator.zeta"], 
                    indicator.alpha = object$parameter["indicator.alpha"], frail.base = object$parameter["frail.base"], 
                    n.knots = object$parameter["n.knots"], LIMparam = object$parameter["LIMparam"], LIMlogl = object$parameter["LIMlogl"], 
                    LIMderiv = object$parameter["LIMderiv"], nb.mc = object$parameter["nb.mc"], nb.gh = object$parameter["nb.gh"], 
                    nb.gh2 = object$parameter["nb.gh2"], adaptatif = object$parameter["adaptatif"], 
                    int.method = object$parameter["int.method"], nb.iterPGH = object$parameter["nb.iterPGH"], 
                    nb.MC.kendall = object$parameter["nb.MC.kendall"], nboot.kendall = object$parameter["nboot.kendall"], 
                    true.init.val = object$parameter["true.init.val"], theta.init = object$parameter["theta.init"], 
                    sigma.ss.init = object$parameter["sigma.ss.init"], sigma.tt.init = object$parameter["sigma.tt.init"], 
                    sigma.st.init = object$parameter["sigma.st.init"], gamma.init = object$parameter["gamma.init"], 
                    alpha.init = object$parameter["alpha.init"], zeta.init = object$parameter["zeta.init"], 
                    betas.init = object$parameter["betas.init"], betat.init = object$parameter["betat.init"], 
                    scale = object$parameter["scale"], random.generator = object$parameter["random.generator"], 
                    kappa.use = object$parameter["kappa.use"], random = object$parameter["random"], 
                    random.nb.sim = object$parameter["random.nb.sim"], seed = object$parameter["seed"], 
                    init.kappa = object$parameter["init.kappa"], nb.decimal = object$parameter["nb.decimal"], 
                    print.times = object$parameter["print.times"], print.iter = object$parameter["print.iter"])
    # Prediction
    if(is.null(joint.surro)) 
      cat(c("===Model without trial", i, "did not converged!!! please try to modified initial values or others parameters===: \n"))
    else{
      d1 <- predict(joint.surro,datapred = dataUse[dataUse$trialID %in% trial[i],])
      # Merger of the results
      d <- rbind(d,d1)
    }
  }
  
  # impression du temps de calcul
  if (print.times){
    cost<-(proc.time()-ptm)/60
    cat("The program took", round(cost[3],2), "minutes \n")
  }
  
  if(!is.null(d)) class(d) <- "jointSurroPenalloocv"
  
  return(d)
}
