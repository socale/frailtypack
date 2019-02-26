##' Short summary of the simulation studies based on a joint surrogate model
##' 
##' This function returns the true value, the mean of the estimates, 
##' the empirical standard error, the mean of the estimated standard errors 
##' (Mean SE), and the coverage probability for model parameters
##' 
##' 
##' @aliases summary.jointSurroPenalSimul print.summary.jointSurroPenalSimul
##' @usage \method{summary}{jointSurroPenalSimul}(object, d = 3, R2boot = 0,  ...)
##' 
##' @param object an object inheriting from \code{jointSurroPenalSimul} class.
##' @param d The desired number of digits after the decimal point f. Default of 3 
##' @param R2boot A binary that specifies whether the confidence interval of \eqn{R^2_{trial}} 
##' should be computed using parametric bootstrap (\code{1}) or Delta-method (\code{0}). 
##' The default is \code{0}
##' @param \dots other unused arguments.
##'  
##' @return For each parameter of the joint surrogate model , we print the true simulation value,  
##' the empirical standard error (empirical SE), the mean of the estimated standard errors 
##' (Mean SE), and the coverate probability (CP) for each model parameters. 
##' For the Kendall's \eqn{\tau}, the 95\% Confidence interval is obtained by the 
##' parametric bootstrap. For \eqn{R^2_{trial}}(R2trial), the standard error is obtained
##' by Delta-method and the 95\% Confidence interval could be obtained directly or by 
##' parametric bootstrap. We also display the total number of non convergence case with 
##' the associated percentage (R : n(\%)), the mean number of iterations to reach convergence,
##' and other estimation and simulation parameters
##' .
##' @seealso \code{\link{jointSurroPenalSimul}}
##' 
##' @author Casimir Ledoux Sofeu \email{casimir.sofeu@u-bordeaux.fr}, \email{scl.ledoux@gmail.com} and 
##' Virginie Rondeau \email{virginie.rondeau@inserm.fr}
##' 
##' @keywords methods
##' @export
##' @examples
##' 
##' # Studies simulation
#' \dontrun{
#' # (Computation takes around 45 minutes using a processor including 40
#' # cores and a read only memory of 378 Go)
#' joint.simul <- jointSurroPenalSimul(nb.dataset = 10, nbSubSimul=600, 
#'                    ntrialSimul=30, LIMparam = 0.001, LIMlogl = 0.001, 
#'                    LIMderiv = 0.001, nb.mc = 200, nb.gh = 20, 
#'                    nb.gh2 = 32, true.init.val = 1, print.iter=F)
#'
#' # results
#' summary(joint.simul, d = 3, R2boot = 1) # bootstrap
#' summary(joint.simul, d = 3, R2boot = 0) # Delta-method
#' 
#' }
##' 
##' 
"summary.jointSurroPenalSimul"<-
  function(object, d =3, R2boot = 0,  ...){
    
    ick=1 # on tien compte (1) ou non (0) du calcul de l'IC du tau de kendall
    n_bootstrap <- 1000 # mais pas utilise car tout le calcul sed fait dans la subroutine jointSurogate
    nb.paquet <- 1
    nb.decimal <- d
    R2parboot <- R2boot
    
    x <- object
    if (!inherits(x, "jointSurroPenalSimul"))
      stop("Object must be of class 'jointSurroPenalSimul'")
    
    cat("Simulation and estimation pamareters", "\n")
    
    cat(c("nb.subject = ", object$nb.subject), "\n")
    cat(c("nb.trials = ", object$nb.trials), "\n")
    cat(c("nb.simul = ", object$nb.simul), "\n")
    cat(c("int.method = ", object$int.method), "\n")
    if(object$int.method != 0) {
      cat(c("nb.gh = ", object$nb.gh), "\n")
      cat(c("nb.gh2 = ", object$nb.gh2), "\n")
    }
    if(object$int.method %in% c(0, 2, 4)) cat(c("nb.mc = ", object$nb.mc), "\n")
    cat(c("kappa.use = ", object$kappa.use), "\n")
    cat(c("n.knots = ", object$n.knots), "\n")
    cat(c("n.iter = ", object$n.iter), "\n")
    cat(" ", "\n")
    
    tau <- jointSurroTKendall(theta = object$theta2, gamma = object$gamma.ui, alpha = object$alpha.ui, zeta = object$zeta)
    
    resultSimul <- synthese_result_modele_reduit(object$dataParamEstim, object$dataTkendall, 
                                                object$dataR2boot, nb.paquet, nb.decimal, object$nb.simul,
                                                object$theta2, object$zeta, object$gamma.ui, object$alpha.ui, 
                                                object$sigma.s, object$sigma.t, object$sigma.st, object$betas,
                                                object$betat, object$R2, tau, n_bootstrap, ick, R2parboot)
    
    resultSimul[-nrow(resultSimul),1] <- substr(resultSimul[-nrow(resultSimul),1],1,nchar(resultSimul[-nrow(resultSimul),1])-6)
    if(is.na(resultSimul[nrow(resultSimul)-2,ncol(resultSimul)-1]))resultSimul[nrow(resultSimul)-2,ncol(resultSimul)-1] <- "-"
    cat("Simulation results", "\n")
    print(resultSimul[-nrow(resultSimul),])
    cat(c("Rejected datasets : n(%) = ",resultSimul[nrow(resultSimul),3]), "\n")
    
  }