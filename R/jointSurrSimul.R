#'Generate survival times for two endpoints using the joint frailty surrogate model
#'
#'Date are generated from the one-step joint surrogate model (see \code{\link{jointSurroPenal}} for more details)
#'
#'We just considered in this generation, the Gaussian random effects. If the parameter \code{full.data} is set to 1,
#'this function  return a list containning severals parameters, including the generated random effects. 
#'the desired individual level correlation (Kendall's \eqn{\tau}) depend on the values of 
#'\eqn{\alpha}, \eqn{\theta}, \eqn{\gamma} and \eqn{\zeta}.
#'
#' @aliases jointSurrSimul
#' @param n.obs Number of considered  subjects. The default is \code{600}.
#' @param n.trial Number of considered  trials. The default is \code{30}.
#' @param cens.adm censorship time. The default is \code{549}, for about \code{40\%} of censored subjects.
#' @param alpha Fixed value for \eqn{\alpha}. The default is \code{1.5}.
#' @param theta Fixed value for \eqn{\theta}. The default is \code{3.5}.
#' @param gamma Fixed value for \eqn{\gamma}. The default is \code{2.5}.
#' @param zeta Fixed value for \eqn{\zeta}. The default is \code{1}.
#' @param sigma.s Fixed value for \if{latex}{\eqn{\sigma^2_{v_S}}}
#' \if{html}{\eqn{\sigma}\out{<sup>2</sup><sub>v<sub>S</sub></sub>}}. 
#' The default is \code{0.7}.
#' @param sigma.t Fixed value for \if{latex}{\eqn{\sigma^2_{v_T}}}
#' \if{html}{\eqn{\sigma}\out{<sup>2</sup><sub>v<sub>T</sub></sub>}}. 
#' The default is \code{0.7}.
#' @param cor Desired level of correlation between \if{latex}{\eqn{v_{S_i}} and 
#' \eqn{v_{T_i}}}\if{html}{v\out{<sub>S<sub>i</sub></sub>} and v\out{<sub>T<sub>i</sub></sub>}}. 
#'  \if{latex}{\eqn{R^2_{trial} = cor^2}}
#'    \if{html}{\code{R}\out{<sup>2</sup><sub>trial</sub>} = cor \out{<sup>2</sup>}}. 
#' The default is \code{0.8}.
#' @param betas Fixed value for \if{latex}{\eqn{\beta_S}} \if{html}{\eqn{\beta}\out{<sub>S</sub>}}.
#'  The default is \code{-1.25}.
#' @param betat Fixed value for \if{latex}{\eqn{\beta_T}} \if{html}{\eqn{\beta}\out{<sub>T</sub>}}. 
#'  The default is \code{-1.25}.
#' @param frailt.base considered the heterogeneity on the baseline risk \code{(1)} or not \code{(0)}. 
#' The default is \code{1}.
#' @param lambda.S Desired scale parameter for the \code{Weibull} distribution associated with the Surrogate
#' endpoint. The default is 1.8. 
#' @param nu.S Desired shape parameter for the \code{Weibull} distribution associated with the Surrogate
#' endpoint. The default is 0.0045. 
#' @param lambda.T Desired scale parameter for the \code{Weibull} distribution associated with the True endpoint.
#' The default is 3.
#' @param nu.T Desired shape parameter for the \code{Weibull} distribution associated with the True endpoint.
#' The default is 0.0025.
#' @param ver Number of covariates. For surrogte evaluation, we just considered one covatiate, the treatment arm
#' @param typeOf Type of joint model used for data generation: 0 = classical joint model 
#' with a shared individual frailty effect (Rondeau, 2007), 1 = joint surrogate model with shared frailty 
#' effects \if{latex}{\eqn{u_i}} \if{html}{\code{u}\out{<sub>i</sub>}} and \if{latex}{\eqn{\omega_{ij}}} 
#' \if{html}{\eqn{\omega}\out{<sub>ij</sub>}}, and two correlated random effects treatment-by-trial interaction 
#' (\if{latex}{\eqn{v_{S_i}}, \eqn{v_{T_i}}}\if{html}{v\out{<sub>S<sub>i</sub></sub>}, v\out{<sub>T<sub>i</sub></sub>}}) 
#' as described in Sofeu et al. (2018).
#' @param equi.subj.trial A binary variable that indicates if the same proportion of subjects should be included per trial (1) 
#' or not (0). If 0, the proportions of subject per trial are required in parameter \code{prop.subj.trial}.
#' @param equi.subj.trt A binary variable that indicates if the same proportion of subjects is randomized per trial (1) 
#' or not (0). If 0, the proportions of subject per trial are required in parameter \code{prop.subj.trt}.
#' @param prop.subj.trial The proportions of subjects per trial. Requires if \code{equi.subj.trial=0}.
#' @param prop.subj.trt The proportions of randomized subject per trial. Requires if \code{equi.subj.trt=0}.
#' @param full.data Specified if you want the function to return the full dataset (1), including the random effects, 
#' or the restictive dataset (0) with \code{7} columns required for the function \code{\link{jointSurroPenal}}.
#' @param random.generator Random number generator used by the Fortran compiler, 
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
#' all ganerated data are the same. By varying this parameter, different datasets are obtained during data genarations. The default value is 0, 
#' in the event of one dataset.
#' @param pfs Is used to specify if the time to progression should be censored by the death time (0) or not (1). 
#' The default is 0. In the event with pfs set to 1, death is included in the surrogate endpoint as in the definition of PFS or DFS. 
# @param param.weibull A binary for the Weibull parametrization used. The default is \code{0}, as in 
# the frailtypack package. If \code{1} the function 
# \eqn{f(x)=\nu^\lambda . \lambda . x^{\lambda-1} . \exp(-(\nu x)^\lambda)} is used.

#' @return
#' This function return if the parameter \code{full.data} is set to 0, a \code{\link{data.frame}} with columns :
#'    \item{patientID}{A numeric, that represents the patient's identifier, must be unique;}
#'    \item{trialID}{A numeric, that represents the trial in which each patient was randomized;}
#'    \item{trt}{The treatment indicator for each patient, with 1 = treated, 0 = untreated;}
#'    \item{timeS}{The follow up time associated with the surrogate endpoint;}
#'    \item{statusS}{The event indicator associated with the surrogate endpoint. Normally 
#'    0 = no event, 1 = event;}
#'    \item{timeT}{The follow up time associated with the true endpoint;}
#'    \item{statusT}{The event indicator associated with the true endpoint. Normally 
#'    0 = no event, 1 = event;}
#'If the argument \code{full.data} is set to 1, additionnal colums corresponding to random effects 
#'\if{latex}{\eqn{\omega_{ij}}} \if{html}{\eqn{\omega}\out{<sub>ij</sub>}}, 
#'\if{latex}{\eqn{u_i}} \if{html}{\code{u}\out{<sub>i</sub>}}, \if{latex}{\eqn{v_{S_i}} and \eqn{v_{T_i}}}\if{html}{v\out{<sub>S<sub>i</sub></sub>} and
#' v\out{<sub>T<sub>i</sub></sub>}} are returned. Note that
#'\if{latex}{\eqn{u_i}} \if{html}{\code{u}\out{<sub>i</sub>}}, \if{latex}{\eqn{v_{S_i}} and \eqn{v_{T_i}}}\if{html}{v\out{<sub>S<sub>i</sub></sub>} and
#' v\out{<sub>T<sub>i</sub></sub>}} are returned if \code{typeOf} is set to \code{1} 
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
#' Sofeu, C. L., Emura, T., and Rondeau, V. (2019). One-step validation method for surrogate 
#' endpoints using data from multiple randomized cancer clinical trials with failure-time endpoints. 
#' Statistics in Medicine 38, 2928-2942.
#' 
#' @seealso \code{\link{jointSurrSimul}}
#' @export
#'
#'
#' @examples
#' 
#' data.sim <- jointSurrSimul(n.obs=600, n.trial = 30,cens.adm=549.24, 
#'             alpha = 1.5, theta = 3.5, gamma = 2.5, sigma.s = 0.7, 
#'             zeta = 1, sigma.t = 0.7, cor = 0.8, betas = -1.25, 
#'             betat = -1.25, full.data = 0, random.generator = 1, 
#'             seed = 0, nb.reject.data = 0, pfs = 0)
#' 
jointSurrSimul <- function(n.obs = 600, n.trial = 30, cens.adm = 549.24, alpha = 1.5, theta = 3.5, gamma = 2.5, zeta = 1, 
                           sigma.s = 0.7, sigma.t = 0.7,cor = 0.8, betas = -1.25, betat = -1.25, frailt.base = 1,
                           lambda.S = 1.8, nu.S = 0.0045,lambda.T = 3, nu.T = 0.0025, ver = 1, typeOf = 1,
                           equi.subj.trial = 1 ,equi.subj.trt = 1, prop.subj.trial = NULL, prop.subj.trt = NULL,
                           full.data = 0, random.generator = 1, random = 0, random.nb.sim = 0, seed = 0, 
                           nb.reject.data = 0, pfs = 0){
  
  param.weibull <- 0
  n.col <- 13 #Number of columns of the simulated dataset. The required number is 13.
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
  
  # ============end parameters checking====================
  
  data.sim=NULL
  
  if(typeOf==1){ 
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
    
    # == initialisation sans utilisation ==
    thetacopule <- 0
    filtre <- 1
    filtre2 <- 1
    # == Fin initialisation sans utilisation ==
    type.joint.simul <- 1

      
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
                    as.double(0),
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
                    as.double(filtre), 
                    as.double(filtre2),
                    as.integer(type.joint.simul),
                    as.integer(pfs),
                    PACKAGE="frailtypack"
                    )
    
    #ans$don_simul <- data.frame(ans$don_simul)
    #ans$don_simulS1 <- data.frame(ans$don_simulS1)
    
    ans$don_simul <- data.frame(matrix(ans$don_simul,nrow = n.obs , ncol = n.col))
    ans$don_simulS1 <- data.frame(matrix(ans$don_simulS1,nrow = n.obs , ncol = n.col))

    names(ans$don_simul) <- c("trt1","v_s1","v_t1","trialref1","w_ij1","timeS1","timeT1",
                              "timeC1","statusS1","statusT1","initTime1","Patienref1","u_i1")
    names(ans$don_simulS1) <- c("trt1","v_s1","v_t1","trialref1","w_ij1","timeS1","timeT1",
                                "timeC1","statusS1","statusT1","initTime1","Patienref1","u_i1")
    
    data.sim <- ans$don_simulS1[,c(4, 12, 1, 6, 9)] # donnees sans le true
    data.sim <- merge(data.sim,ans$don_simul[,c(12, 7, 10)], by="Patienref1") # on ajoute les donnees sur le True
    
    names(data.sim) <- c("patientID", "trialID", "trt", "timeS", "statusS", "timeT", "statusT")
  
  if(full.data == 1){
    data.comp <- merge(ans$don_simulS1[,c(12, 4, 1, 5, 13, 2, 3, 6, 9)],
                       ans$don_simul[,c(12, 7, 10)],
                       by="Patienref1")
    names(data.comp) <- c("patientID", "trialID", "trt","w_ij","u_i","v_Si","v_Ti", "timeS", "statusS", "timeT", "statusT")
    if(typeOf == 1) data.comp <- data.comp[c(1:3,8:11,4:7)]
    if(typeOf == 0) data.comp <- data.comp[c(1:3,8:11,4)]
    return(data.comp)
  }
  if(full.data == 0) return(data.sim)
}
