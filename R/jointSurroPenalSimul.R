
#' Simulation studies based on the one-step Joint surrogate model for the evaluation of a canditate 
#' surrogate endpoint
#'
#'@description{
#' This function aims to allow simulation studies, based on the joint frailty surrogate model, 
#' described in \link{jointSurroPenal}
#' }
#' 
#' @details{
#' The estimated parameter are obtained using the robust Marquardt algorithm
#' (Marquardt, 1963) which is a combination between a Newton-Raphson algorithm
#' and a steepest descent algorithm. The iterations are stopped when the
#' difference between two consecutive log-likelihoods was small
#' \eqn{(<10^{-3})}, the estimated coefficients were stable (consecutive
#' values \eqn{(<10^{-3})}, and the gradient small enough \eqn{(<10^{-3})}, by default.
#' Cubic M-splines of order 4 are used for the hazard function, and I-splines (integrated M-splines) are
#' used for the cumulative hazard function.
#' 
#' The inverse of the Hessian matrix is the variance estimator and to deal
#' with the positivity constraint of the variance component and the spline
#' coefficients, a squared transformation is used and the standard errors are
#' computed by the \eqn{\Delta}-method (Knight & Xekalaki, 2000). The smooth
#' parameter can be chosen by maximizing a likelihood cross validation
#' criterion (Joly and other, 1998). 
#' 
#' We proposed based on the joint surrogate model a new definition 
#' of the Kendall's \eqn{\tau}. Moreover, distinct numerical integration methods are available to approximate the 
#' integrals in the marginal log-likelihood.
#' 
#' \bold{Non-convergence case management procedure}
#' 
#' Special attention must be given to initializing model parameters, the choice of the number of 
#' spline knots, the smoothing parameters and the number of quadrature points to solve convergence 
#' issues. We first initialized parameters using the user's desired strategy, as specified 
#' by the option \code{true.init.val}. When numerical or convergence problems are encountered, 
#' with \code{kappa.use} set to \code{4}, the model is fitted again using a combination of the following strategies: 
#' vary the number of quadrature point (\code{nb.gh} to \code{nb.gh2} or \code{nb.gh2} to \code{nb.gh})
#' in case of the use of the Gaussian Hermite quadrature integration (see \code{int.method}); 
#' divided or multiplied the smoothing parameters (\code{k_1}, \code{k_2}) by 10 or 100 according to 
#' their preceding values, or used parameter vectors obtained during the last iteration (with a 
#' modification of the number of quadrature points and smoothing parameters). Using this strategy, 
#' we usually obtained during simulation the rejection rate less than 3\%. A sensitivity analysis 
#' was conducted without this strategy, and similar results were obtained on the converged samples, 
#' with about a 23\% rejection rate. 
#' }
#' 
#' @aliases jointSurroPenalSimul
#' @usage 
#' jointSurroPenalSimul(maxit=40, indicator.zeta = 1, 
#'    indicator.alpha = 1, frail.base = 1, n.knots = 6, nb.dataset = 1, 
#'    nbSubSimul=1000, ntrialSimul=30, LIMparam = 0.001, 
#'    LIMlogl = 0.001, LIMderiv = 0.001, nb.mc = 300, nb.gh = 32, 
#'    nb.gh2 = 20, adaptatif = 0, int.method = 2, nb.iterPGH = 5, 
#'    nb.MC.kendall = 10000, nboot.kendall = 1000, true.init.val = 0, 
#'    theta.init = 1, sigma.ss.init = 0.5, sigma.tt.init = 0.5, 
#'    sigma.st.init = 0.48, gamma.init = 0.5, alpha.init = 1, 
#'    zeta.init = 1, betas.init = 0.5, betat.init = 0.5, 
#'    random.generator = 1, equi.subj.trial = 1, prop.subj.trial = NULL, 
#'    equi.subj.trt = 1, prop.subj.trt = NULL, 
#'    theta2 = 3.5, zeta = 1, gamma.ui = 2.5, alpha.ui = 1, 
#'    betas = -1.25, betat = -1.25, lambdas = 1.8, nus = 0.0045, 
#'    lambdat = 3, nut = 0.0025, time.cens = 549, R2 = 0.81,
#'    sigma.s = 0.7, sigma.t = 0.7, kappa.use = 4, random = 0, 
#'    random.nb.sim = 0, seed = 0, init.kappa = NULL, 
#'    nb.decimal = 4, print.times = TRUE, print.iter=FALSE)
#'
#' @param maxit maximum number of iterations for the Marquardt algorithm.
#' Default is \code{40}. 
#' @param indicator.zeta A binary, indicates whether the power's parameter \eqn{\zeta} should 
#' be estimated (1) or not (0). If \code{0}, \eqn{\zeta} will be set to \code{1} during estimation. 
#' The default is \code{1}. This parameter can be seted to \code{0} in case of identification issues. 
#' @param indicator.alpha A binary, indicates whether the power's parameter \eqn{\alpha} should 
#' be estimated (1) or not (0). If \code{0}, \eqn{\alpha} will be set to \code{1} during estimation.
#' The default is 1.
#' @param frail.base Considered the heterogeneity between trial on the baseline risk (\code{1}), using 
#' the shared cluster specific frailties (\eqn{u_i}), or not (\code{0}). The default is \code{1}.
#' @param n.knots integer giving the number of knots to use. Value required in
#' the penalized likelihood estimation.  It corresponds to the (n.knots+2)
#' splines functions for the approximation of the hazard or the survival
#' functions.  We estimate I or M-splines of order 4. When the user set a
#' number of knots equals to k (n.knots=k) then the number of interior knots
#' is (k-2) and the number of splines is (k-2)+order.  Number of knots must be
#' between 4 and 20. (See \code{\link{frailtyPenal}} for more details).
#' @param nb.dataset Number of dataset to analyze. The default is \code{1}.
#' @param nbSubSimul Number of subjects.
#' @param ntrialSimul Number of trials.
#' @param LIMparam Convergence threshold of the Marquardt algorithm for the
#' parameters, \eqn{10^{-3}} by default (See \code{\link{frailtyPenal}} for more details).
#' @param LIMlogl Convergence threshold of the Marquardt algorithm for the
#' log-likelihood, \eqn{10^{-3}} by default (See \code{\link{frailtyPenal}} for more details).
#' @param LIMderiv Convergence threshold of the Marquardt algorithm for the gradient, \eqn{10^{-3}} by default 
#' (See \code{\link{frailtyPenal}} for more details).
#' @param nb.mc Number of samples considered in the Monte-Carlo integration. Required in case 
#' \code{int.method} is equals to \code{0}, \code{2} or \code{4}. A value between 100 and 300 most often gives 
#' good results. However, beyond 300, the program takes a lot of time to estimate the parameters.
#' The default is \code{300}.
#' @param nb.gh Number of nodes for the Gaussian-Hermite quadrature. It can
#' be chosen among 5, 7, 9, 12, 15, 20 and 32. The default is 32.
#' @param nb.gh2 Number of nodes for the Gauss-Hermite quadrature used to re-estimate the model, 
#' in case of non-convergence, defined as previously. The default is \code{20}.
#' @param adaptatif A binary, indicates whether the pseudo adaptive Gaussian-Hermite quadrature 
#' \code{(1)} or the classical Gaussian-Hermite quadrature \code{(0)} is used. The default is \code{0}.
#' @param int.method A numeric, indicates the integration method: \code{0} for Monte carlo, 
#' \code{1} for Gaussian-Hermite quadrature, \code{2} for a combination of both Gaussian-Hermite quadrature to 
#' integrate over the individual-level random effects and Monte carlo to integrate over the trial-level
#' random effects, \code{4} for a combination of both Monte carlo to integrate over 
#' the individual-level random effects and Gaussian-Hermite quadrature to integrate over the trial-level
#' random effects. The default is \code{2}.
#' @param nb.iterPGH Number of iterations before the re-estimation of the posterior random effects,
#' in case of the two-steps pseudo-adaptive Gaussian-hermite quadrature. If set to \code{0} there is no 
#' re-estimation". The default is \code{5}.
#' @param nb.MC.kendall Number of generated points used with the Monte-Carlo to estimate
#' integrals in the Kendall's \eqn{\tau} formulation. Beter to use at least 4000 points for
#' stable reseults. The default is \code{10000}.
#' @param nboot.kendall Number of samples considered in the parametric bootstrap to estimate the confidence
#' interval of the Kendall's \eqn{\tau}. The default is \code{1000}. 
#' @param true.init.val Numerical value. Indicates if the real parameter values 
#' \code{(1)}, or the given initial values to parameters \code{(0)} should be considered. 
#' If set to \code{2}, \eqn{\alpha} and \eqn{\gamma} are initialised using two separed shared frailty model 
#' (see \code{\link{frailtyPenal}} for more details); \eqn{\sigma^2_{v_S}}, \eqn{\sigma^2_{v_T}} and
#' \eqn{\sigma_{v_{ST}}} are fixed using the default initial values given by the user; \eqn{\zeta}, 
#' \eqn{\theta}, \eqn{\beta_S} and \eqn{\beta_T} are initialized using a classical joint 
#' frailty model, considering individual level random effects. If the joint frailty model is 
#' faced to convergence issues, \eqn{\beta_S} and \eqn{\beta_T} are initialized using 
#' two shared frailty models.  In all others scenarios, if the simplified model does not converge,
#' default given parameters values are used. Initial values for spline's associated parameters 
#' are fixed to \code{0.5}. The default for this argument is \code{0}.
#' @param theta.init Initial values for \eqn{\theta}, required if \code{true.init.val} 
#' is set to \code{0} or \code{2}. The default is \code{1}.
#' @param sigma.ss.init Initial values for \eqn{\sigma^2_{v_S}}, required if \code{true.init.val} 
#' is set to \code{0} or \code{2}. The default is \code{0.5}.
#' @param sigma.tt.init Initial values for \eqn{\sigma^2_{v_T}}, required if \code{true.init.val} 
#' is set to \code{0} or \code{2}. The default is \code{0.5}.
#' @param sigma.st.init Initial values for \eqn{\sigma_{v_{ST}}}, required if \code{true.init.val} 
#' is set to \code{0} or \code{2}. The default is \code{0.48}.
#' @param gamma.init Initial values for \eqn{\gamma}, required if \code{true.init.val} 
#' is set to \code{0} or \code{2}. The default is \code{0.5}.
#' @param alpha.init Initial values for \eqn{\alpha}, required if \code{true.init.val} 
#' is set to \code{0} or \code{2}. The default is \code{1}.
#' @param zeta.init Initial values for \eqn{\zeta}, required if \code{true.init.val} 
#' is set to \code{0} or \code{2}. The default is \code{1}.
#' @param betas.init Initial values for \eqn{\beta_S}, required if \code{true.init.val} 
#' is set to \code{0} or \code{2}. The default is \code{0.5}.
#' @param betat.init Initial values for \eqn{\beta_T}, required if \code{true.init.val} 
#' is set to \code{0} or \code{2}. The default is \code{0.5}.
#' @param random.generator Random number generator to use by the Fortran compiler, 
#' \code{1} for the intrinsec subroutine \code{Random_number} and \code{2} for the 
#' subroutine \code{uniran()}. The default is \code{1}. 
#' @param equi.subj.trial A binary, that indicates if the same proportion of subjects per trial
#' should be considered in the procces of data generation (1) or not (0). In case of 
#' different trial sizes, fill in \code{prop.subj.trial} the proportions
#' of subjects to be considered per trial. The default is \code{1}.
#' @param prop.subj.trial Vector of the proportions of subjects to consider per trial. 
#' Requires if the argument \code{equi.subj.trial} is different to \code{1}. The size of this vector is equal to the 
#' number of trials.
#' @param equi.subj.trt Indicates if the same proportion of treated subjects per trial should be
#' considered \code{(1)} or not \code{(0)}. If \code{0}, fill in \code{prop.subj.trt} 
#' the proportions of treated subjects to be considered per trial. The default is \code{1}.
#' @param prop.subj.trt Vector of the proportions of treated subjects to consider per trial. 
#' Requires if the argument \code{equi.subj.trt} is different to \code{0.5}. The size of this vector is equal to the 
#' number of trials.
#' @param theta2 True value for \eqn{\theta}. The default is \code{3.5}.
#' @param zeta True value for \eqn{\zeta} in case of simulation. The default is \code{1}.
#' @param gamma.ui True value for \eqn{\gamma} in case of simulation. The default is \code{2.5}.
#' @param alpha.ui True value for \eqn{\alpha} in case of simulation. The default is \code{1}.
#' @param betas True value for \eqn{\beta_S} in case of simulation. The default is \code{-1.25}.
#' @param betat True value for \eqn{\beta_T} in case of simulation. The default is \code{-1.25}.
#' @param lambdas Desired scale parameter for the \code{Weibull} distribution associated with the Surrogate
#' endpoint. The default is \code{1.8}.
#' @param nus Desired shape parameter for the \code{Weibull} distribution associated with the Surrogate
#' endpoint. The default is \code{0.0045}.
#' @param lambdat Desired scale parameter for the \code{Weibull} distribution associated with the True 
#' endpoint.The default is \code{3}.
#' @param nut Desired shape parameter for the \code{Weibull} distribution associated with the True endpoint.
#' The default is \code{0.0025}.
#' @param time.cens Censorship time. The default is \code{549}, for about \code{40\%} of censored 
#' subjects.
#' @param R2 Desired \eqn{R^2_{trial}}. The default is \code{0.81}.
#' @param sigma.s True value for \eqn{\sigma^2_S}. The default is \code{0.7}.
#' @param sigma.t True value for \eqn{\sigma^2_T}. The default is \code{0.7}.
# @param param.weibull A binary for the Weibull parametrization used. The default is \code{0}, as in 
# the frailtypack package. If \code{1} the function 
# \eqn{f(x)=\nu^\lambda . \lambda . x^{\lambda-1} . \exp(-(\nu x)^\lambda)} is used.
#' @param kappa.use A numeric, that indicates how to manage the smoothing parameters \code{k_1} 
#' and \code{k_2} in case of convergence issues. If it is set 
#' to \code{0}, the first smoothing parameters that allowed convergence on the first dataset is used 
#' for all simulations. if it is set to \code{1}, a smoothing parameter is estimated by cross-validation 
#' for each dataset generated. If it is set to \code{2}, the same process for chosing kappas as in case 
#' \code{1} is used, but in case of convergence issue, the first smoothing parameters that allowed 
#' convergence among the three previous that have worked is used. If it is set to \code{3}, the associated 
#' smoothing parameters are successively divided by 10, in case of convergence issues until 5 times. 
#' If it is set to \code{4}, the management of the smoothing 
#' parameters is as in case \code{2}, preceded by the successive division described in case \code{3} and 
#' by the changing of the number of nodes for the Gauss-Hermite quadrature. The default is \code{4}.
#' @param random A binary that says if we reset the random number generation with a different environment 
#' at each call \code{(1)} or not \code{(0)}. If it is set to \code{1}, we use the computer clock 
#' as seed. In the last case, it is not possible to reproduce the generated datasets". 
#' The default is \code{0}. Required if \code{random.generator} is set to 1.
#' @param random.nb.sim If \code{random} is set to \code{1}, a binary that indicates the number 
#' of generations that will be made, equal to \code{nb.dataset} in this case.
#' @param seed The seed to use for data generation. Required if \code{random} is set to \code{0}. 
#' The default is \code{0}.
#' @param init.kappa smoothing parameter used to penalized the log-likelihood. By default (init.kappa = NULL) the values used 
#' are obtain by cross-validation.
#' @param nb.decimal Number of decimal required for results presentation.
#' @param print.times a logical parameter to print estimation time. Default
#' is TRUE.
#' @param print.iter a logical parameter to print iteration process. Default
#' is FALSE.
#' 
#' @return
#' This function return an object of class jointSurroPenalSimul with elements :
#' 
#'    \item{theta2}{True value for \eqn{\theta};}
#'    \item{zeta}{true value for \eqn{\zeta};}
#'    \item{gamma.ui}{true value for \eqn{\gamma};}
#'    \item{alpha.ui}{true value for \eqn{\alpha};}
#'    \item{sigma.s}{true value for \eqn{\sigma_S};}
#'    \item{sigma.t}{true value for \eqn{\sigma_T};}
#'    \item{sigma.st}{true value for \eqn{\sigma_{ST}};}
#'    \item{betas}{true value for \eqn{\beta_S};}
#'    \item{betat}{true value for \eqn{\beta_T};}
#'    \item{R2}{true value for \eqn{R^2_{trial}};}
#'    \item{nb.subject}{total number of subjects used;}
#'    \item{nb.trials}{total number of trials used;}
#'    \item{nb.simul}{number of simulated datasets;}
#'    \item{nb.gh}{number of nodes for the Gaussian-Hermite quadrature;} 
#'    \item{nb.gh2}{number of nodes for the Gauss-Hermite quadrature used to re-estimate the model, in case of non-convergence;}
#'    \item{nb.mc}{number of samples considered in the Monte-Carlo integration;}
#'    \item{kappa.use}{a numeric, that indicates how to manage the smoothing parameters k_1 and k_2 in case of convergence issues;}
#'    \item{n.knots}{number of knots used for splines;}
#'    \item{int.method}{integration method used;}
#'    \item{n.iter}{mean number of iterations needed to converge;}
#'    \item{dataTkendall}{a matrix with \code{nb.dataset} line(s) and three columns, of the estimates of Kendall's \eqn{\tau} 
#'    and theirs confidence intervals using the parametric bootstrap. All non-convergence cases  are represented by a line of 0;}
#'    \item{dataR2boot}{a matrix with \code{nb.dataset} line(s) and three columns, of the estimates of \eqn{R^2_{trial}} 
#'    and theirs confidence intervals using the parametric bootstrap. All non-convergence cases are represented by a line of 0.}
#'    \item{dataParamEstim}{a dataframe including all estimates with the associated standard errors, for all simulation. 
#'    All non-convergence cases  are represented by a line of 0;}
#'    \item{dataHessian}{Dataframe of the variance-Covariance matrices  of the estimates for all simulations}
#'    \item{dataHessianIH}{Dataframe of the robust estimation of the variance matrices  of the estimates for all simulations}
#'    \item{datab}{Dataframe of the estimates for all simulations which rich convergence}
#'    
#'   
#' @seealso \code{\link{jointSurroPenal}}, \code{\link{summary.jointSurroPenalSimul}}, \code{\link{jointSurrSimul}}
#' 
#' @author Casimir Ledoux Sofeu \email{casimir.sofeu@u-bordeaux.fr}, \email{scl.ledoux@gmail.com} and 
#' Virginie Rondeau \email{virginie.rondeau@inserm.fr}
#' 
#' @references
#' Sofeu C.L., Emura T. and Rondeau V. (2018). One-step validation method for surrogate 
#' endpoints in multiple randomized cancer clinical trials with failure-time endpoints. 
#' \code{Under review}
#' 
#' @export
#' @importFrom doBy orderBy
#'
#' @examples
#' 
#' \dontrun{
#' # Surrogacy model evaluation performance study based on 10 generated data
#' # (Computation takes around 20 minutes using a processor including 40 
#' # cores and a read only memory of 378 Go)
#' # To realize a simulation study on 100 samples or more (as required), use 
#' # nb.dataset = 100
#' 
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
#' 
jointSurroPenalSimul = function(maxit=40, indicator.zeta = 1, indicator.alpha = 1, frail.base = 1, n.knots = 6,
                      nb.dataset = 1, nbSubSimul=1000, ntrialSimul=30, LIMparam = 0.001, LIMlogl = 0.001,
                      LIMderiv = 0.001, nb.mc = 300, nb.gh = 32, nb.gh2 = 20, adaptatif = 0, int.method = 2, 
                      nb.iterPGH = 5, nb.MC.kendall = 10000, nboot.kendall = 1000, true.init.val = 0, 
                      theta.init = 1, sigma.ss.init = 0.5, sigma.tt.init = 0.5, sigma.st.init = 0.48, 
                      gamma.init = 0.5, alpha.init = 1, zeta.init = 1, betas.init = 0.5, betat.init = 0.5,
                      random.generator = 1, equi.subj.trial = 1, prop.subj.trial = NULL, equi.subj.trt = 1,
                      prop.subj.trt = NULL, theta2 = 3.5, zeta = 1, 
                      gamma.ui = 2.5, alpha.ui = 1, betas = -1.25, betat = -1.25, lambdas = 1.8, nus = 0.0045, 
                      lambdat = 3, nut = 0.0025, time.cens = 549, R2 = 0.81, sigma.s = 0.7, sigma.t = 0.7, 
                      kappa.use = 4, random = 0, random.nb.sim = 0, seed = 0, init.kappa = NULL,
                      nb.decimal = 4, print.times = TRUE, print.iter = FALSE){
  
  data <- NULL
  scale <- 1
  one.dataset <- 0
  real.data <- 0
  gener.only <- 0
  param.weibull <- 0
  
  # ==============parameters checking======================
  if(!(indicator.zeta %in% c(0,1)) | !(indicator.alpha %in% c(0,1)) | !(frail.base %in% c(0,1))){
    stop("model options indicator.zeta, indicator.alpha and frail.base must be set to 0 or 1")
  }
  
  
  if(!(true.init.val %in% c(0,1,2))){
    stop("argument 'true.init.val' is a Numeric that requires as value: 0,1, 2. See the help ...")
  }
  
  if(!(equi.subj.trt %in% c(0,1))){
    stop("The argument 'equi.subj.trt' must be set to 0 or 1")
  }
  
  
  if(!(equi.subj.trial %in% c(0,1))){
    stop("The argument 'equi.subj.trial' must be set to 0 or 1")
  }
  
  if(!(random.generator %in% c(1,2))){
    stop("The argument 'random.generator' must be set to 1 or 2")
  }
  if(!(param.weibull %in% c(1,0))){
    stop("The argument 'param.weibull' must be set to 0 or 1")
  }
  
  # ============End parameters checking====================
  
  nsujet1 <- nbSubSimul
  ng <- nbSubSimul
  ntrials1 <- ntrialSimul
  
  
  maxiter <- maxit
  nst <- 2 # nb de fonction de risque (2 pour )
  
  # indique si l'on estime (1=oui, 0=non) covST
  indice_covST <- 1
  indice_gamma_st <- 0 #  indice_gamma_st: dit si l'on estime gamma_st_ut (1) ou non(0), pour les effets aleatoires correlees sur le risque de base, pas traite ici 
  
  if(frail.base==0) indicator.alpha <- 0 
  
  indice_a_estime <- c(indicator.zeta, indice_covST, indicator.alpha, indice_gamma_st,frail.base)
  
  if(indice_covST == 1){
    # we estimated at least 4 parameters correspondint to the covariance matrix \sigma and the variance of \omega_ij
    nb.frailty <- 4
    nparamfrail <- nb.frailty + indicator.zeta + indicator.alpha + frail.base
  } else{
    nb.frailty <- 3
    nparamfrail <- nb.frailty + indicator.zeta + indicator.alpha + frail.base
  }
  
  # parametre fonction de risque de base
  gamma1 <- 2 # paramertre de la loi gamma
  gamma2 <- 2 # paramertre de la loi gamma
  typeof <- 0 # type de function de risque  0:Splines,  1:Cpm  2:weib
  nbintervR <- 12 # nb interv surrogate
  nbintervDC <- 12 # nb interv true
  equidistant <- 0 # 1=equidist, 0=percentile
  nz <- n.knots
  param_risque_base <- c(typeof,nbintervR,nbintervDC,equidistant,nz)
  
  #nombre de variables explicatives
  ves <- 1 # nombre variables explicative surrogate
  ved <- 1 # nombre variables explicative deces/evenement terminal
  ver <- 1 # nombre total de variables explicative
  nbrevar <- c(ves,ved,ver) 
  
  # vecteur des noms de variables
  nomvarl<- "trt"
  # matrice d'indicatrice de prise en compte des variables explicatives pour le surrogate et le tru
  # filtre = vecteur associe au surrogate
  # filtre2 = vecteur associe au true
  filtre  <- 1
  filtre2 <- 1
  filtre0 <- as.matrix(data.frame(filtre,filtre2))
  
  # gestion de l'affichage a l'ecran
  flush.console()
  if (print.times){
    ptm<-proc.time()
    cat("\n")
    cat("Be patient. The program is computing ... \n")
  }
  
  
  # nombre de simulatin, 1 par default
  n_sim1 <- nb.dataset
  
  kapa <- "kappa_valid_crois.txt"
  kappa0 <- init.kappa
    donnees <- NULL
    death   <- NULL
    # simuler les donnees et preparer le fichier des kappa pour la validation croisee
    # nom du fichier pour les kappas obtenues par validation croisee
    kapa <- "kappa_valid_crois.txt"
    vect_kappa <- matrix(0,nrow = n_sim1,ncol = 2)
    for(j in 1:n_sim1){
      data.sim <- jointSurrSimul(n.obs=nbSubSimul, n.trial = ntrialSimul,cens.adm=time.cens, 
                      alpha = alpha.ui, theta = theta2, gamma = gamma.ui, zeta = zeta, sigma.s = sigma.s, 
                      sigma.t = sigma.t, rsqrt = R2, betas = betas, betat = betat, full.data = 0, 
                      random.generator = random.generator, seed = seed, nb.reject.data = j-1)
      data.sim$initTime <- 0
      donnees <- data.sim[,c("trialID","patienID","trt","initTime","timeS","statusS")]
      death   <- data.sim[,c("trialID","patienID","trt","initTime","timeT","statusT")]
      # conversion en double des jeux de donneees. je le fais separemment pour distinguer 
      # les cas ou j'aurai plus de variables explicatives pour un des jeux de donnees que pour l'autre
      for(i in 1:ncol(donnees)){
        donnees[,i] <- as.double(donnees[,i])
      }
      for(i in 1:ncol(death)){
        death[,i] <- as.double(death[,i])
      }
      
      if(print.iter) cat("+++++++++++estimation of Kappas by cross-validation +++++++++++")
      vect_kappa[j,] <- kappa_val_croisee(don_S = donnees, don_T = death, njeu = 1, n_obs = nsujet1,
                                     n_node = n.knots, adjust_S = 1, adjust_T = 1, kapp_0 = 0,
                                     print.times = print.iter)
      # I deleate the created text file
      file.remove(dir(pattern="kappa_valid_crois.txt"))
    }
   # utils::write.table(vect_kappa,"kappa_valid_crois.txt",sep=" ",row.names = F,col.names = F)

  # critere de convergence du modele on donne en entree les critere a respecter et en sortie on recupere ceux obtenue du programme
  EPS2 <- c(LIMparam, LIMlogl, LIMderiv)
  
  logNormal <- 1 #lognormal: indique si on a une distribution lognormale des effets aleatoires (1) ou Gamma (0)
  
  # Parametres d'integration
  nsim_node <- rep(NA,10)
  nsim_node[1] <- nb.mc # nombre de simulation pour l'integration par Monte carlo, vaut 0 si on ne veut pas faire du MC
  nsim_node[2] <- nb.gh # nombre de points de quadrature a utiliser (preference 5 points pour l'adaptatice et 32 poits pour la non adaptatice)
  nsim_node[3] <- adaptatif # doit-on faire de l'adaptative(1) ou de la non-adaptative(0)
  nsim_node[4] <- int.method# indique la methode d'integration 0=Monte carlo,1= quadrature, 2=quadrature individuel+MC essai, 3=Laplace, 4= monte carlo individuel + quadrature essai
  nsim_node[5] <- nparamfrail
  nsim_node[6] <- 1 # indique si lon fait de la vectorisation dans le calcul integral (1) ou non (0). rmq: la vectorisation permet de reduire le temps de calcul
  nsim_node[7] <- nb.frailty # indique le nombre d'effet aleatoire cas quadrature adaptative
  type.joint <- 1 # type de modele a estimer: 0=joint classique avec un effet aleatoire partage au niveau individuel,1=joint surrogate avec 1 frailty partage indiv et 2 frailties correles essai,
  # 2=joint surrogate sans effet aleatoire partage donc deux effets aleatoires a chaque fois"
  nsim_node[8] <- type.joint 
  nsim_node[9] <- nb.gh2 # nombre de point de quadrature a utiliser en cas de non convergence de prefenrence 7 ou 9 pour la pseudo adaptative et 32 pour la non adaptative
  nsim_node[10] <- nb.iterPGH # nombre d'itteration aubout desquelles reestimer les effects aleatoires a posteriori pour la pseude adaptative. si 0 pas de resestimation
  
  # Parametres associes au taux de kendall et au bootstrap
  meth.int.kendal <- 4
  method_int_kendal <- meth.int.kendal# methode d'integration pour le taux de kendall: 0= montecarle, 1= quadrature quaussienne classique, 2= approximation de Laplace
  N_MC_kendall <- nb.MC.kendall # nombre de boucle MC pour le calcul du taux de kendal en approximant l'integrale par montye carlo
  nboot_kendal <- nboot.kendall# nombre d'echantillon bootstrap pour le calcul de l'IC du taux de ke,ndall
  Param_kendall_boot <- c(method_int_kendal,N_MC_kendall,nboot_kendal)
  
  # vecteur contenant les taux de kendall et R2 issus du bootstrap
  fichier_kendall <- NULL
  fichier_R2 <- NULL
  
  # nom des fichiers de sortie
  param_estime <- "Parametre_estime.txt"
  param_empirique <- "Parametre_empirique.txt"
  param_empirique_NC <- "Parametre_empirique_NC.txt"
  tableau_rejet <- "tab_rejet.txt"
  NomFichier <- c(kapa,param_estime,param_empirique,param_empirique_NC,tableau_rejet)
  
  #true.init.val # dit si on initialise les parametres avec les vraies(1) valeurs en cas de 
  #simultation, ou des valeurs donnes par default(0). si 0, on estime 4 model de cox a fragilites partages pour les parametres par default
  
  # Parametres initiaux
  #theta.init		 # valeur initiale de theta2
  #sigma.ss.init	 # valeur initiale de sigma ss
  #sigma.tt.init	 # valeur initiale de sigma tt
  #sigma.st.init  # valeur initiale de sigma st
  #gamma.init  # valeur initiale de gamma.ui
  #alpha.init  # valeur initiale de alpha.ui
  #zeta.init  # valeur initiale de zeta_wij
  #betas.init  # valeur initiale de betas
  #betat.init  # valeur initiale de betat
  
  param_init <- c(theta.init,sigma.ss.init,sigma.tt.init,sigma.st.init,gamma.init,alpha.init,
                  zeta.init,betas.init,betat.init)
  
  revision_echelle <- scale # coefficient pour la division des temps de suivi. permet de reduire l'echelle des temps pour eviter les problemes numeriques en cas d'un nombre eleve de sujet pour certains cluster
  # random.generator <- # generateur des nombre aleatoire, (1) si Random_number() et (2) si uniran(). Random_number() me permet de gerer le seed
  sujet_equi <- equi.subj.trial # dit si on considere la meme proportion de sujet par essai au moment de la generation des donnee (1) ou non (0). dans ce dernier cas remplir les proportion des sujets par essai dansle fichier dedie
  
  # proportion des patients traites par esssai. si 0 alors on a des proportions variables suivant les essai, remplir le fichier dedie. si non on aura le meme proportion des traites par essai
  if(equi.subj.trt == 1) prop_trait <- 0.5
  else prop_trait <- 0
  
  # parametres de simultation, en cas de simultation
  # gamma1 <- # parametre de la loi gamma
  # gamma2 <- # parametre de la loi gamma
  # theta2 <- # variance des frailties  w_ij S
  eta <- zeta# Zeta associe a w_ij chez les deces
  # gamma.ui <- # variance des frailties associees au risque de base
  # alpha.ui <- # parametre de puissance associe a u_i chez les deces
  theta2_t <- 0.8# variance des frailties  w_ij chez les T, en cas de fragilite correles
  rsqrt_theta <- 0.8# niveau de correlation entre wij_s et wij_t
  gamma.uit <- 0.8# variance des frailties associees au risque de base sur le true
  rsqrt_gamma.ui <- 0.7# niveau de correlation entre us_i et ut_i
  # betas <- # effet fixe du traitement associe au surrogate
  # betat <- # effet fixe du traitement associe au true
  # parametres d'echelle et de forme pour la weibull
  # lambdas <- # lambdas
  # nus <- # nus
  # lambdat <- #   lambdat
  # nut <- # nut
  mode_cens <- 1 # 1= quantille et 2= date fixe
  temps_cens <- time.cens# censure fixe: temps a preciser
  cens0 <- 0.25 # si quentile, proportion des patients censures
  rsqrt <- sqrt(R2)# niveau de correlation souhaite pour les frailties niveau essai
  # sigma.s <- # variance des effest aleatoires au niveau essai en interaction avec le traitement, associee au surrogate
  # sigma.t <- # variance des effest aleatoires au niveau essai en interaction avec le traitement, associee au true
  paramSimul <- c(gamma1, gamma2, theta2, eta, gamma.ui, alpha.ui, theta2_t, rsqrt_theta, gamma.uit,
                  rsqrt_gamma.ui, betas, betat, lambdas, nus, lambdat,nut, mode_cens, temps_cens,
                  cens0, rsqrt, sigma.s, sigma.t)
  
  # Autres parametres de simulation
  weib <- 1# 0= on simule les temps par une loi exponentielle, 1= on simule par une weibull
  # param.weibull <- # parametrisation de la weibull utilisee: 0= parametrisation par default dans le programme de Virginie, 1= parametrisation a l'aide de la fonction de weibull donnee dans le cous de Pierre
  frailty_cor <- 1 # indique si l'on considere pour le modele de simulation deux effets aleatoire correles au niveau essai(=1) ou un effet aleatoire partage(=0) ou encore on simule sans effet aleatoire au niveau essai(=2, model conjoint classique)
  affiche_stat <- 0 # dit si l'on affiche les statistiques des donnees simulees(1) ou non (0)
  seed_ <- 0 # jeux de donnees a retenir pour la validation croisee
  une_donnee <- one.dataset # pour dire si on simule avec un seul jeu de donnees(1) ou pas (0). ceci pour tester le programme d'estimation
  donne_reel <- real.data # dit si 1 a la question precedente dit s'il sagit du jeux de donnees reel (1) ou non (0)
  #gener.only <- # dit si on voudrait seulement generer les donnees(1) ou generer et faire des simulation(0)
  #kappa.use <- # dit si on utilise un kappa a chaque generation de donnee (1) ou le premier kappa pour tous les jeux de donnees(0)
  decoup_simul <- 0# dans le cas ou l'on a decoupe les simulations en plusieurs paquets, donne le nombre de generation de donnees a ne pas considerer avant d'engager les simulations. ceci empeche de reproduire les meme jeux de donnees pour tous les paquets de simulation. vaut 0 si pas de decoupage pevu sinon pour chaque jeux de simulation mettre cette valeur a jour. Exp si 10 paquets de simul pour un total de 100, on affecte 0 pour le premier paquet, 10 pour le second, 20 pour le 3 ieme, ... 90 pour le 10ieme
  aleatoire <- random# dit si on reinitialise la generation des nombre aleatoire avec un environnement different a chaque appel (1) ou non(O).En cas de generation differente, on utilise l'horloge (heure) de l'ordinateur comme graine. Dans ce cas, il n'est pas possible de reproduire les donnees simulees
  nbre_sim <- random.nb.sim# dans le cas ou aleatoire=1, cette variable indique le nombre de generation qui vont etre faites
  graine <- seed # dans le cas ou l'on voudrait avoir la possibilite de reproduire les donnees generees alors on met la variable aleatoire=0 et on donne dans cette variable la graine a utiliser pour la generation
  autreParamSim <- c(weib,param.weibull,frailty_cor,affiche_stat,seed_,une_donnee,donne_reel,gener.only,
                     kappa.use,decoup_simul,aleatoire,nbre_sim,graine)
  
  # autres dichiers de sortie
  # vecteur des pametres
  nva=2 # deux parametres lies aux effets fixes du traitement
  effet <- 0
  if(typeof==0) np <- 2*(nz+2) + nva + nparamfrail
  if(typeof==1) np <- nbintervDC + nbintervR + nva + nparamfrail
  if(typeof==2) np <- 2*nst + nva + effet + nparamfrail
  
  # matrice hessienne
  H_hessOut <- matrix(0,np,np)
  
  # matrice hessienne corigee
  HIHOut <- matrix(0,np,np)
  resOut <- 0 # log-vraisamblance penalisee
  LCV <- c(0,0) # value de LCV(1) et AIC (2)
  
  # Parametre associe au fonctions de risques et survies de base et splines
  # extrait du Joint fortran
  mt11 <- 100
  mt12 <- 100
  
  if(typeof == 1){
    mt1 <- 3*nbintervR
    mt2 <- 3*nbintervDC
  }
  else{
    mt1 <- 100
    mt2 <- 100
  }
  
  # vecteur donnant la taille des tableau pour fortran
  sizeVect=c(np,mt1,mt2,mt11,mt12)
  
  x1Out <- rep(0, mt1) # temps pour la representation des fonction de risque et de la survies surrogate
  lamOut <- matrix(0, nrow = mt1, ncol= 3) # estimates du risque de base surrogate avec IC
  xSu1 <- matrix(0, mt11)
  suOut <- matrix(0, nrow = mt11 , ncol=3) # matrice des estimates de la survie a baseline(avec IC) surrogate et deces
  x2Out <- rep(0, mt2) # temps pour la representation des fonction de risque et de la survies true endpoint
  lam2Out <- matrix(0, nrow = mt2, ncol= 3) # estimates du risque de base deces avec IC
  xSu2 <- matrix(0, mt11)
  su2Out <- matrix(0, nrow = mt12 , ncol = 3) # matrice des estimates de la survie a baseline(avec IC) surrogate et deces
  
  ni <- 0 # nombre d'itteration pour la convergence
  ier <- 0 # informe sur le comportement du modele(-1 = erreur, k = perte de significativite le modele continu, 0 = pas d'erreur)
  istop <- 0 # critere d'arret: 1= le modele a converge, 2= on a attent le nombre max d'itteration, 3= echec inversion de la hessienne, 4= erreur dans les calculs 
  ziOut <- rep(0,nz+6)  # knots for baseline hazard estimated with splines
  Varcov <- matrix(0, nrow = 3, ncol = 3) # matrice de variance-covariance de (sigma_S,sigma_ST,sigmaT) obtenue par delta methode a partir de la hesienne, en raison du changement de variable au moment de l'estimation
  dataHessian <- matrix(0, nrow = np*n_sim1, ncol = np) # sauvegarde des matrices hessiennes des differentes simulations 
  dataHessianIH <- matrix(0, nrow = np*n_sim1, ncol = np) # pour la matrice hessienne corrigee (HIH)
  datab <- matrix(0, nrow = n_sim1, ncol = np) # sauvegarde des vecteurs de parametres des simulation 
  
  # proportion de sujet par essai
  if(sujet_equi==1){
    prop_i <- rep(1/ntrials1, ntrials1)
  }
  else{
    if(is.null(prop.subj.trial) | length(prop.subj.trial)!=ntrials1) stop("The proportion of subjects per trial are required")
    else prop_i <- prop.subj.trial
  }
  
  #proportion de sujet traites par essai
  if(equi.subj.trt==1){
    p <- rep(0.5, ntrials1)
  }
  else{
    if(is.null(prop.subj.trt) | length(prop.subj.trt)!=ntrials1) stop("The proportions of treated subjects per trial are required")
    else p <- prop.subj.trt
  }
  
  
  if(print.iter) 
    affiche.itter <- 1
  else
    affiche.itter <- 0
  
  param_estimes <- NULL
  
  ans <- .Fortran(C_jointsurrogate,
                  as.integer(nsujet1),
                  as.integer(ng),
                  as.integer(ntrials1),
                  as.integer(maxiter),
                  as.integer(nst),
                  as.integer(nparamfrail),
                  as.integer(indice_a_estime),
                  as.integer(param_risque_base),
                  as.integer(nbrevar),
                  as.integer(filtre0),
                  as.matrix(donnees),
                  as.matrix(death),
                  as.double(p),
                  as.double(prop_i),
                  as.integer(n_sim1),
                  EPS2 = as.double(c(LIMparam, LIMlogl, LIMderiv)),
                  as.double(kappa0),
                  as.double(vect_kappa),
                  as.integer(logNormal),
                  nsim_node = as.integer(nsim_node),
                  as.integer(Param_kendall_boot),
                  as.integer(true.init.val),
                  as.double(param_init),
                  as.double(revision_echelle),
                  as.integer(random.generator),
                  as.integer(sujet_equi),
                  as.double(prop_trait),
                  as.double(paramSimul),
                  as.integer(autreParamSim),
                  fichier_kendall = matrix (0,nrow = n_sim1, ncol = 3), # debut section des parametres de sortie
                  fichier_R2 = matrix (0,nrow = n_sim1, ncol = 3),
                  param_estimes = matrix (0,nrow = n_sim1, ncol = 24),
                  as.integer(sizeVect),
                  b = rep(0,np),
                  H_hessOut = matrix(0,np,np),
                  HIHOut = matrix(0,np,np),
                  resOut = 0,
                  LCV = c(0,0),
                  x1Out = rep(0, mt1),
                  lamOut = matrix(0, nrow = mt1, ncol= 3),
                  xSu1 = matrix(0, mt11),
                  suOut = matrix(0, nrow = mt11 , ncol=3),
                  x2Out = rep(0, mt2),
                  lam2Out = matrix(0, nrow = mt2, ncol= 3),
                  xSu2 = matrix(0, mt11),
                  su2Out = matrix(0, nrow = mt12 , ncol = 3),
                  ni= as.integer(0),
                  ier = 0,
                  istop = 0,
                  ziOut = rep(0,nz+6),
                  as.integer(affiche.itter),
                  Varcov = matrix(0, nrow = 3, ncol = 3),
                  dataHessian = matrix(0, nrow = np*n_sim1, ncol = np),
                  dataHessianIH = matrix(0, nrow = np*n_sim1, ncol = np),
                  datab = matrix(0, nrow = n_sim1, ncol = np),
                  PACKAGE="frailtypack"
  )
  
  # resultats a retourner:
  result <- NULL
  result$theta2  <- theta2
  result$zeta <- zeta
  result$gamma.ui  <- gamma.ui
  result$alpha.ui <- alpha.ui
  result$sigma.s  <- sigma.s
  result$sigma.t  <- sigma.t
  result$sigma.st <- sqrt(R2 * sigma.s * sigma.t)
  result$betas <- betas
  result$betat  <- betat
  result$R2  <- R2
  result$nb.subject <- nbSubSimul
  result$nb.trials <- ntrialSimul
  result$nb.simul <- nb.dataset
  result$nb.gh <- nb.gh  
  result$nb.gh2 <- nb.gh2 
  result$nb.mc  <- nb.mc
  result$kappa.use <- kappa.use
  result$n.knots <- n.knots
  result$int.method <- int.method
  result$n.iter <- ans$ni
  result$dataTkendall <- data.frame(ans$fichier_kendall)
  result$dataR2boot <- data.frame(ans$fichier_R2)
  result$dataParamEstim <- data.frame(ans$param_estimes)[,-c(21:23)] # on fait sauter les autres taux de kendall
  names(result$dataParamEstim) <- c("theta","SE.theta","zeta","SE.zeta","beta.S","SE.beta.S","beta.T","SE.beta_T","sigma.S",
                                    "SE.sigma.S","sigma.T","SE.sigma.T","sigma.ST","SE.sigma.ST","gamma","SE.gamma","alpha","SE.alpha",
                                    "R2trial","SE.R2trial","tau")
  names(result$dataTkendall) <- c("Ktau","inf.95%CI","sup.95%CI")
  names(result$dataR2boot) <- c("R2.boot","inf.95%CI","sup.95%CI")
  result$dataHessian <- data.frame(ans$dataHessian)
  result$datab <- data.frame(ans$datab)
  
  #if(is.na(result$n.iter)) result=NULL # model did not converged 
  
  # =====================++++++++++++++++++++++++++++++++++++++++++++++++++++
  # SI L'ON N'EST PAS EN SIMULATION, PENSER A SUPPRIMER LES FICHIERS CREES+++
  # =====================++++++++++++++++++++++++++++++++++++++++++++++++++++
  # if(nb.dataset == 1){
  #   # file.remove(c("outjoint", "OutJoint_Result_surrogate.txt", "OutJoint_simul.txt",
  #   #               "param_estime_scl.txt", "Parametre_estime.txt", "R2_bootstrap.txt",
  #   #               "Resultat_simulation.txt", "surrogate.txt", "tab_rejet.txt", 
  #   #               "Taux_kendall_bootst.txt", "true.txt", "kappa_valid_crois.txt"))
  # }
  
  # file.remove(c("true.txt", "surrogate.txt", "outjoint"))
  # try(file.remove("OutJoint_Result_surrogate.txt"))
  # try(file.remove("kappa_valid_crois.txt"))
  
  class(result) <- "jointSurroPenalSimul"
  
  # impression du temps de calcul
  if (print.times){
    cost<-(proc.time()-ptm)/60
    cat("The program took", round(cost[3],2), "minutes \n")
  }
  
  return(result)
}




