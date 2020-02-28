#'Generate survival times for two endpoints using the joint frailty-copula model for surrogacy
#'
#'Date are generated from the one-step joint frailty-copula model, under the Claton 
#'copula function (see \code{\link{jointSurroCopPenal}} for more details)
#'
#'We just considered in this generation, the Gaussian random effects. If the parameter \code{full.data} is set to 1,
#'this function  return a list containning severals parameters, including the generated random effects. 
#'The desired individual level correlation (Kendall's \eqn{\tau}) depend on the values of the copula parameter 
#'\eqn{\theta}, given that \eqn{\tau = \theta /(\theta + 2)} under the clayton copula model.
#'
#' @aliases jointSurrCopSimul
#' @param n.obs Number of considered  subjects. The default is \code{600}.
#' @param n.trial Number of considered  trials. The default is \code{30}.
#' @param prop.cens A value between \code{0} and \code{1}, \code{1-prop.cens} is the minimum proportion of 
#' people who are randomly censored. 
#' Represents the quantile to use for generating the random censorship time. In this case, the censorship 
#' time follows a uniform distribution in \code{1} and \code{(prop.cens)ieme} percentile of the 
#' generated death times. If this argument is set to \code{0}, the fix censorship is considered.
#' The default is \code{0}. 
#' @param cens.adm Censorship time. If argument \code{prop.cens} is set to \code{0}, it represents 
#' the administrative censorship time, else it represents the fix censoring time. The default is \code{549}, 
#' for about \code{40\%} of fix censored subjects.
#' @param alpha Fixed value for \eqn{\alpha}. The default is \code{1.5}.
#' @param gamma Fixed value for \eqn{\gamma}. The default is \code{2.5}.
#' @param sigma.s Fixed value for  \if{latex}{\eqn{\sigma^2_{v_S}}}
#' \if{html}{\eqn{\sigma}\out{<sup>2</sup><sub>v<sub>S</sub></sub>}}. The default is \code{0.7}.
#' @param sigma.t Fixed value for \if{latex}{\eqn{\sigma^2_{v_T}}}
#' \if{html}{\eqn{\sigma}\out{<sup>2</sup><sub>v<sub>T</sub></sub>}}. The default is \code{0.7}.
#' @param cor Desired level of correlation between \if{latex}{\eqn{v_{S_i}} and 
#' \eqn{v_{T_i}}}\if{html}{v\out{<sub>S<sub>i</sub></sub>} and v\out{<sub>T<sub>i</sub></sub>}}. 
#'  \if{latex}{\eqn{R^2_{trial} = cor^2}}
#'    \if{html}{\code{R}\out{<sup>2</sup><sub>trial</sub>} = cor \out{<sup>2</sup>}}. 
#' The default is \code{0.8}.
#' @param betas Vector of the fixed effects for \if{latex}{\eqn{\beta_S}} \if{html}{\eqn{\beta}\out{<sub>S</sub>}}.
#'  The size must be equal to \code{ver} 
#' The default is \code{c(-1.25,0.5)}.
#' @param betat Vector of the fixed effects for  \if{latex}{\eqn{\beta_T}} \if{html}{\eqn{\beta}\out{<sub>T</sub>}}.
#'  The size must be equal to \code{ver}
#' The default is \code{c(-1.25,0.5)}.
#' @param frailt.base Considered heterogeneity on the baseline risk \code{(1)} or not \code{(0)}. 
#' The default is \code{1}.
#' @param lambda.S Desired scale parameter for the \code{Weibull} distribution associated with the Surrogate
#' endpoint. The default is 1.8. 
#' @param nu.S Desired shape parameter for the \code{Weibull} distribution associated with the Surrogate
#' endpoint. The default is 0.0045. 
#' @param lambda.T Desired scale parameter for the \code{Weibull} distribution associated with the True endpoint.
#' The default is 3.
#' @param nu.T Desired shape parameter for the \code{Weibull} distribution associated with the True endpoint.
#' The default is 0.0025.
#' @param ver Number of covariates. The mandatory covariate is the treatment arm. The default is \code{2}.
#' @param typeOf Type of joint model used for data generation: 0 = classical joint model 
#' with a shared individual frailty effect (Rondeau, 2007), 1 = joint frailty-copula model with shared frailty 
#' effects \if{latex}{\eqn{u_i}} \if{html}{\code{u}\out{<sub>i</sub>}} and two correlated random effects treatment-by-trial interaction 
#' (\if{latex}{\eqn{v_{S_i}}, \eqn{v_{T_i}}}\if{html}{v\out{<sub>S<sub>i</sub></sub>}, v\out{<sub>T<sub>i</sub></sub>}}),
#'  see \code{\link{jointSurroCopPenal}}.
#' @param equi.subj.trial A binary variable that indicates if the same proportion of subjects should be included per trial (1) 
#' or not (0). If 0, the proportions of subject per trial are required with parameter \code{prop.subj.trial}.
#' @param equi.subj.trt A binary variable that indicates if the same proportion of subjects is randomized per trial (1) 
#' or not (0). If 0, the proportions of subject per trial are required with parameter \code{prop.subj.trt}.
#' @param prop.subj.trial The proportions of subjects per trial. Requires if \code{equi.subj.trial = 0}.
#' @param prop.subj.trt The proportions of randomized subject per trial. Requires if \code{equi.subj.trt = 0}.
#' @param full.data Specified if you want the function to return the full dataset (1), including the random effects, 
#' or the restictive dataset (0) with at least \code{7} columns as required for the function \code{\link{jointSurroCopPenal}}.
#' @param random.generator The random number generator used by the Fortran compiler, 
#' \code{1} for the intrinsec subroutine \code{Random_number} and \code{2} for the 
#' subroutine \code{uniran()}. The default is \code{1}. 
#' @param random A binary that says if we reset the random number generation with a different environment 
#' at each call \code{(1)} or not \code{(0)}. If it is set to \code{1}, we use the computer clock 
#' as seed. In the last case, it is not possible to reproduce the generated datasets. 
#' The default is \code{0}. Required if \code{random.generator} is set to 1.
#' @param random.nb.sim required if \code{random.generator} is set to 1, and if \code{random} is set to 1.
#' @param seed The seed to use for data (or samples) generation. Required if the argument \code{random.generator} is set to 1. 
#' Must be a positive value. If negative, the program do not account for seed. The default is \code{0}.
#' @param nb.reject.data Number of generation to reject before the considered dataset. This parameter is required
#' when data generation is for simulation. With a fixed parameter and \code{random.generator} set to 1,
#' all ganerated data are the same. By varying this parameter, different datasets are obtained during data generations. The default value is 0, 
#' in the event of one dataset.
#' @param filter.surr Vector of size the number of covariates, with the i-th element that indicates if the hazard for 
#' surrogate is adjusted on the i-th covariate (code 1) or not (code 0). By default, 2 covariates are considered.
#' @param thetacopule The desired value for the copula parameter. The default is \code{6}.
#' @param filter.true Vector defines as \code{filter.surr}, for the true endpoint. \code{filter.true} and \code{filter.surr}
#' should have the same size
#' @param covar.names Vector of the names of covariables. By default it contains "trt" for the 
#' tratment arm. Should contains the names of all covarites wished in the generated dataset.
#' @param pfs Is used to specify if the time to progression should be censored by the death time (0) or not (1). 
#' The default is 0. In the event with pfs set to 1, death is included in the surrogate endpoint as in the definition of PFS or DFS. 
# @param param.weibull A binary for the Weibull parametrization used. The default is \code{0}, as in 
# the frailtypack package. If \code{1} the function 
# \eqn{f(x)=\nu^\lambda . \lambda . x^{\lambda-1} . \exp(-(\nu x)^\lambda)} is used.

#' @return
#' This function returns if the parameter \code{full.data} is set to 0, a \code{\link{data.frame}} with columns :
#'    \item{patientID}{A numeric, that represents the patient's identifier, must be unique;}
#'    \item{trialID}{A numeric, that represents the trial in which each patient was randomized;}
#'    \item{trt}{The treatment indicator for each patient, with 1 = treated, 0 = untreated;}
#'    \item{timeS}{The follow up time associated with the surrogate endpoint;}
#'    \item{statusS}{The event indicator associated with the surrogate endpoint. Normally 
#'    0 = no event, 1 = event;}
#'    \item{timeT}{The follow up time associated with the true endpoint;}
#'    \item{statusT}{The event indicator associated with the true endpoint. Normally 
#'    0 = no event, 1 = event;}
#'and other covariates named \code{Var2, var3, ..., var[ver-1]} if \code{ver > 1}.
#' If the argument \code{full.data} is set to 1, additionnal colums corresponding to random effects 
#'\if{latex}{\eqn{u_i}} \if{html}{\code{u}\out{<sub>i</sub>}}, \if{latex}{\eqn{v_{S_i}} and 
#'\eqn{v_{T_i}}}\if{html}{\code{v}\out{<sub>S<sub>i</sub></sub>} and
#' \code{v}\out{<sub>T<sub>i</sub></sub>}} are returned. 
#'    
#'
#' @author Casimir Ledoux Sofeu \email{casimir.sofeu@u-bordeaux.fr}, \email{scl.ledoux@gmail.com} and 
#' Virginie Rondeau \email{virginie.rondeau@inserm.fr}
#' 
#' @references
#' 
#' Rondeau V., Mathoulin-Pelissier S., Jacqmin-Gadda H., Brouste V. and Soubeyran P. (2007).
#' Joint frailty models for recurring events and death using maximum penalized likelihood 
#' estimation: application on cancer events. Biostatistics 8(4), 708-721.
#' 
#' Sofeu, C. L., Emura, T., and Rondeau, V. (2020). A joint frailty-copula model for meta-analytic 
#' validation of failure time surrogate endpoints in clinical trials. \code{Under review}
#' 
#' @seealso \code{\link{jointSurrSimul}, \link{jointSurroCopPenal}}
#' @export
#'
#'
#' @examples
#' 
#' # dataset with 2 covariates and fixed censorship
#' data.sim <- jointSurrCopSimul(n.obs=600, n.trial = 30, prop.cens = 0, cens.adm=549, 
#'             alpha = 1.5, gamma = 2.5, sigma.s = 0.7, sigma.t = 0.7, 
#'             cor = 0.8, betas = c(-1.25, 0.5), betat = c(-1.25, 0.5), 
#'             full.data = 0, random.generator = 1,ver = 2, covar.names = "trt", 
#'             nb.reject.data = 0, thetacopule = 6, filter.surr = c(1,1), 
#'             filter.true = c(1,1), seed = 0)
#'             
#' #dataset with 2 covariates and random censorship
#' 
#' data.sim2 <- jointSurrCopSimul(n.obs=600, n.trial = 30, prop.cens = 0.75, 
#'             cens.adm = 549, alpha = 1.5, gamma = 2.5, sigma.s = 0.7, 
#'             sigma.t = 0.7, cor = 0.8, betas = c(-1.25, 0.5), 
#'             betat = c(-1.25, 0.5), full.data = 0, random.generator = 1,
#'             ver = 2, covar.names = "trt", nb.reject.data = 0, thetacopule = 6, 
#'             filter.surr = c(1,1), filter.true = c(1,1), seed = 0)
#' 
jointSurrCopSimul <- function(n.obs = 600, n.trial = 30, prop.cens = 0, cens.adm = 549, alpha = 1.5, gamma = 2.5, 
                           sigma.s = 0.7, sigma.t = 0.7,cor = 0.9, betas = c(-1.25, 0.5), betat = c(-1.25, 0.5), 
                           frailt.base = 1, lambda.S = 1.3, nu.S = 0.0025,lambda.T = 1.1, nu.T = 0.0025, ver = 2, typeOf = 1,
                           equi.subj.trial = 1 ,equi.subj.trt = 1, prop.subj.trial = NULL, prop.subj.trt = NULL,
                           full.data = 0, random.generator = 1, random = 0, random.nb.sim = 0, seed = 0, nb.reject.data = 0,
                           thetacopule = 6, filter.surr = c(1,1), filter.true = c(1,1), covar.names = "trt", pfs = 0){
  
  param.weibull <- 0
  theta <- 3.5
  zeta <- 1
  rsqrt <- cor
  # ==============parameters checking======================
  if(!(equi.subj.trt %in% c(0,1)) | !(equi.subj.trial %in% c(0,1))){
    stop("Model's parameters equi.subj.trt and equi.subj.trial must be set to 0 or 1")
  }
  
  if(((equi.subj.trial == 0) & is.null(prop.subj.trial)) | ((equi.subj.trt == 0) & is.null(prop.subj.trt))){
    stop("The proportions of randomized subjects per trial (or the proportions of subject per trial) are required in the variables 
          \bold{prop.subj.trt} (or \bold{prop.subj.trial}). If you want the same proportions, set the parameter \bold{equi.subj.trial} (or \bold{equi.subj.trt}) to 1
         model's parameters equi.subj.trt and equi.subj.trial must be set to 0 or 1")
  }
  
  if(is.null(filter.surr) | is.null(filter.true)){
    stop("The vectors filter.surr and filter.true must contain at least one element corresponding to the effect of the treatment")
  }
  
  if(!(length(betas) == ver) | !(length(betat)==ver)){
    stop("The vectors betas and betat must contain a number of elements corresponding to ver")
  }
  
  if(!(length(filter.surr) == ver) | !(length(filter.true)==ver)){
    stop("The vectors filter.surr and filter.true must contain a number of elements corresponding to ver")
  }
  
  if(!(length(filter.surr) == length(filter.true))){
    stop("The vectors filter.surr and filter.true should have the same size")
  }
  # ============end parameters checking====================
  
  if(length(filter.surr) > length(covar.names)){
    # si plus d'une variable explicatives avec des noms pas preciser, je les nome par var[numero], 
    covar.names = c(covar.names,paste("var",seq(2,length(filter.surr)), sep = ""))
  }
  n.col <- 13 + length(filter.surr) -1 #Number of columns of the simulated dataset. The required number is 13 when just the treatment effect is considered as covariate.
  data.sim <- NULL
  
  if(typeOf == 1){ 
    # joint surrogate model with shared frailty u_i and omega_ij
    lognormal <- 1
  }else{
    # joint classical model with shared individual frailty effect omega_ij, to take into account heterogeneity at 
    # the individual level
    lognormal <- 2
  }
    gamma1 <- 2 # paramertre de la loi gamma
    gamma2 <- 2 # paramertre de la loi gamma
    if(equi.subj.trt==1) p <- rep(0.5,n.trial)
    if(equi.subj.trt==0) p <- rep(0,n.trial)
    if(equi.subj.trial==1) {
      prop_i <- rep(1/n.trial,n.trial)
    }
    else{
      prop_i <- prop.subj.trial
    }
    
    don_simul <- as.double(matrix(0, nrow = n.obs , ncol = n.col))
    don_simulS1 <- as.double(matrix(0, nrow = n.obs , ncol = n.col))
    
    type.joint.simul = 2
    filtre <- filter.surr
    filtre2 <- filter.true
    # filtre <- matrix(filter.surr, nrow = 1, ncol = ver)
    # filtre2 <- matrix(filter.true, nrow = 1, ncol = ver)
    
    ans <- .Fortran(C_surrosim,
                    don_simul = as.double(matrix(0, nrow = n.obs , ncol = n.col)),
                    don_simulS1 = as.double(matrix(0, nrow = n.obs , ncol = n.col)),
                    as.integer(n.obs),
                    as.integer(n.col),
                    as.integer(lognormal),
                    as.integer(0),
                    vrai_theta=as.double(0),
                    as.integer(n.obs) ,
                    as.integer(ver) ,
                    as.double(alpha) ,
                    as.double(prop.cens),
                    as.double(cens.adm),
                    as.double(gamma1),
                    as.double(gamma2),
                    as.double(theta),
                    as.double(lambda.S),
                    as.double(nu.S),
                    as.double(lambda.T),
                    as.double(nu.T),
                    as.double(betas),
                    as.double(betat),
                    as.integer(n.trial),
                    as.double(rsqrt),
                    as.double(sigma.s),
                    as.double(sigma.t),
                    as.double(p),
                    as.double(prop_i),
                    as.double(gamma),
                    as.double(zeta) ,
                    as.integer(frailt.base),
                    as.integer(random.generator),
                    as.integer(random), 
                    as.integer(random.nb.sim) , 
                    as.integer(seed),
                    as.integer(nb.reject.data),
                    as.integer(param.weibull),
                    as.double(thetacopule),
                    as.integer(filtre), 
                    as.integer(filtre2),
                    as.integer(type.joint.simul),
                    as.integer(pfs),
                    PACKAGE="frailtypack"
                    )
    
    #ans$don_simul <- data.frame(ans$don_simul)
    #ans$don_simulS1 <- data.frame(ans$don_simulS1)
    
    ans$don_simul <- data.frame(matrix(ans$don_simul,nrow = n.obs , ncol = n.col))
    ans$don_simulS1 <- data.frame(matrix(ans$don_simulS1,nrow = n.obs , ncol = n.col))

    names(ans$don_simul) <- c("trt1","v_s1","v_t1","trialref1","timeS1","timeT1",
                              "timeC1","statusS1","statusT1","initTime1","Patienref1","u_i1", 
                              covar.names[-1])
    names(ans$don_simulS1) <- c("trt1","v_s1","v_t1","trialref1","timeS1","timeT1",
                                "timeC1","statusS1","statusT1","initTime1","Patienref1","u_i1", 
                                covar.names[-1])
    
    data.sim <- ans$don_simulS1[,c(4, 11, 1, 5, 8)] # donnees sans le true
    data.sim <- merge(data.sim,ans$don_simul[,c(11, 6, 9)], by="Patienref1") # on ajoute les donnees sur le True
    if(length(covar.names)>1)
      data.sim <- merge(data.sim,ans$don_simul[,c(11,12-1+seq(1,length(covar.names))[-1])], by="Patienref1") # on ajoute les donnees sur le True
    
    names(data.sim) <- c("patientID", "trialID", "trt", "timeS", "statusS", "timeT", "statusT", covar.names[-1])
    
    if(full.data == 1){
      data.comp <- merge(ans$don_simulS1[,c(11, 4, 1, 12, 2, 3, 5, 8)],
                         ans$don_simul[,c(11, 6, 9,12-1+seq(1,length(covar.names))[-1])],
                         by="Patienref1")
      
      names(data.comp) <- c("patientID", "trialID", "trt","u_i","v_Si","v_Ti", "timeS", "statusS", "timeT", "statusT", covar.names[-1])
      
      if(typeOf == 1) {
        if(length(covar.names) == 1) 
          data.comp <- data.comp[c(1:2,7:10,3:6)]
        else # ajout des autres covariables
          data.comp <- data.comp[c(1:2,7:10,3:6,(ncol(data.comp)-length(covar.names)+2):ncol(data.comp))]
        }
      if(typeOf == 0) {
        if(length(covar.names) == 1) 
          data.comp <- data.comp[c(1:2,7:10,3)]
        else # ajout des autres covariables
          data.comp <- data.comp[c(1:2,7:10,3, (ncol(data.comp)-length(covar.names)+2):ncol(data.comp))]
      }
      #return(ans)
      return(data.comp)
  }
  if(full.data == 0) {
    return(data.sim)
    #return(ans)
  }
}
