
# This function aims to perform simulation studies using small packages and save interim results in the .RData files
# All the ".RData" are save in ./packageSimul"num.paquet", where "num.paquet" is a number give to the simulation design 

simulationPackages <- function(nsim = 500, nsubjet = 600, ntrial = 30, int.method = 0,
                               nb.mc = 1000, nb.gh = 20, nb.gh2 = 32, adaptatif = 0,
                               nspline = 6, kappa.use = 4, type.joint.estim = 3, 
                               typecopula = 1, type.joint.simul = 3, theta.copula = 3, 
                               time.cens = 349, prop.cens = 0, true.init.val = 1, R2 = 0.81, maxit = 40, 
                               nb.paquet = 10, num.paquet = 1, program.to.estimate = 1, 
                               list.package = NULL
){
  
  # =====description of the parameters=====
   # nsim : number of dataset to generate
   # nsubjet : number of subjects
   # ntrial : number of trials
   # nspline : number of knots for spline
   # nb.paquet : number of packages to considered
   # num.paquet : a number give to the simulation design
   # program.to.estimate : 1 by default, 2 = case with two covariates, 4 = case with init values for theta set to 1 and alpha set to 0.5
     # 5 = case with Variation of the number of subject per trial
   # list.package: list of the packages numbers to execute. If NULL, all nb.package are executed
   # More details for all the arguments can be found by displaying the help on the "jointSurroPenalSimul" function
      # i.e. help("jointSurroPenalSimul")
  # =======end description============
  
  # some initializations
  ckappa1 <- 0 
  ckappa2 <- 0
  pfs     <- 0
  
  #Management of the list of packages on which simulation studies should be performed
  
  if(is.null(list.package)){
    list.packages <- seq(1,nb.paquet)
  }else{
    list.packages <- sort(list.package)
  }
  
  # creation of the work dicrectory for the results
  current <- getwd()
  wd <- paste("./packageSimul",num.paquet,sep="")
  dir.create(wd)
  
  i <- 1
  nb.reject.data <- 0
  N.sim.data     <- 0
  
  while(i <= length(list.packages)){
    if(i %in% list.packages){
      # initialization of require arguments: number of dataset to consider per package and number of dataset to reject before simulation
      q <- floor(nsim / nb.paquet)
      modulo <- nsim %% nb.paquet
      if(i == 0){
        nb.reject.data <- 0
        if(modulo == 0){
          N.sim.data <- q
        }else{
          N.sim.data <- q+1
        }
      }else{
        k <- i - 1
        # number of  dataset to reject before simulation
        if(modulo >= k){
          nb.reject.data <- k * q + k
        }else{
          nb.reject.data <- k * q + modulo
        }
        
        # number of dataset to consider per package
        if(modulo >= i){
          N.sim.data <- q + 1
        }else{
          N.sim.data <- q
        }
      }
      
      # Simulation for the ongoing package
      if(program.to.estimate == 1){ # The default
        joint.simul2 <- jointSurroPenalSimul(nb.dataset = N.sim.data, nbSubSimul = nsubjet, ntrialSimul = ntrial,
                                             int.method = int.method, nb.mc = nb.mc, nb.gh = nb.gh, nb.gh2 = nb.gh2, 
                                             adaptatif = adaptatif, n.knots = nspline, kappa.use = kappa.use, 
                                             type.joint.estim = type.joint.estim, type.joint.simul = type.joint.simul, 
                                             time.cens = time.cens, lambdas = 1.3, nus = 0.0025, lambdat = 1.1, nut = 0.0025, 
                                             seed = 0, betas = c(-1.25), betat = c(-1.25), filter.surr = c(1),
                                             filter.true = c(1), betas.init = c(-0.25), betat.init = c(-0.25), init.kappa = NULL, 
                                             maxit = maxit, true.init.val = true.init.val, theta.copula = theta.copula, 
                                             thetacopula.init = 3, R2 = R2, typecopula = typecopula, ckappa = c(ckappa1, ckappa2), 
                                             gamma.ui = 0.8, nb.reject.data = nb.reject.data, pfs = pfs, print.iter = T)
                                                          
      }
      
      if(program.to.estimate==2){ # case with two covariates
        joint.simul2 <- jointSurroPenalSimul(nb.dataset = N.sim.data, nbSubSimul = nsubjet, ntrialSimul = ntrial,
                                             int.method = int.method, nb.mc = nb.mc, nb.gh = nb.gh, nb.gh2 = nb.gh2, 
                                             adaptatif = adaptatif, n.knots = nspline, kappa.use = kappa.use, 
                                             type.joint.estim = type.joint.estim, type.joint.simul = type.joint.simul, 
                                             time.cens = time.cens, lambdas = 1.3, nus = 0.0025, lambdat = 1.1, nut = 0.0025, 
                                             seed = 0, betas = c(-1.25, -0.5), betat = c(-1.25,-0.5), filter.surr = c(1,1),
                                             filter.true = c(1,1), betas.init = c(-0.25, -0.25), betat.init = c(-0.25, -0.25), 
                                             init.kappa = c(10000,10000), maxit = maxit, true.init.val = true.init.val, 
                                             theta.copula = theta.copula, thetacopula.init = 3, R2 = R2, typecopula = typecopula, 
                                             ckappa = c(ckappa1, ckappa2), gamma.ui = 0.8, nb.reject.data = nb.reject.data, 
                                             pfs = pfs, print.iter = F )
                                                          
      }                                                  
      
      if(program.to.estimate==4){ # case with init values for theta set to 1 and alpha set to 0.5
        joint.simul2 <- jointSurroPenalSimul(nb.dataset = N.sim.data, nbSubSimul = nsubjet, ntrialSimul = ntrial,
                                             int.method = int.method, nb.mc = nb.mc, nb.gh = nb.gh, nb.gh2 = nb.gh2, 
                                             adaptatif = adaptatif, n.knots = nspline, kappa.use = kappa.use, 
                                             type.joint.estim = type.joint.estim, type.joint.simul = type.joint.simul, 
                                             time.cens = time.cens, lambdas = 1.3, nus = 0.0025, lambdat = 1.1, 
                                             nut = 0.0025, seed = 0, betas = c(-1.25), betat = c(-1.25), filter.surr = c(1),
                                             filter.true = c(1), betas.init = c(-0.25), betat.init = c(-0.25), 
                                             init.kappa = NULL, maxit = maxit, true.init.val = true.init.val, 
                                             theta.copula = theta.copula, thetacopula.init = 1, alpha.init = 0.5, R2 = R2,
                                             typecopula = typecopula, ckappa = c(ckappa1, ckappa2), gamma.ui = 0.8, 
                                             nb.reject.data = nb.reject.data, pfs = pfs, print.iter = F)
      }
      
      if(program.to.estimate == 5){ # Variation of the number of subject per trial
        # utils::data("dataOvarian", envir = environment(), package = "frailtypack")
        # consider the proportion from the advanced avarian cancer dataset including 50 trials
        # prop = table(dataOvarian$trialID)/nrow(dataOvarian)
        prop <- c(0.229865772, 0.104865772, 0.057885906, 0.098154362, 0.017617450, 0.006711409, 
                  0.041107383, 0.033557047, 0.005872483, 0.009228188, 0.020973154, 0.008389262,
                  0.014261745, 0.020973154, 0.009228188, 0.014261745, 0.006711409, 0.005033557, 
                  0.010906040, 0.015100671, 0.005033557, 0.001677852, 0.015939597, 0.012583893,
                  0.017617450, 0.012583893, 0.004194631, 0.010067114, 0.007550336, 0.002516779,
                  0.006711409, 0.010067114, 0.001677852, 0.014261745, 0.007550336, 0.033557047,
                  0.003355705, 0.002516779, 0.026006711, 0.002516779, 0.005033557, 0.003355705, 
                  0.004194631, 0.014261745, 0.010067114, 0.005033557, 0.004194631, 0.014261745,
                  0.002516779, 0.008389262)
        joint.simul2 <- jointSurroPenalSimul(nb.dataset = N.sim.data, nbSubSimul = nsubjet, ntrialSimul = ntrial, 
                                             int.method = int.method, nb.mc = nb.mc, nb.gh = nb.gh, nb.gh2 = nb.gh2, 
                                             adaptatif = adaptatif, n.knots = nspline, kappa.use = kappa.use, 
                                             type.joint.estim = type.joint.estim, type.joint.simul = type.joint.simul, 
                                             time.cens = time.cens, lambdas = 1.3, nus = 0.0025, lambdat = 1.1, nut = 0.0025, 
                                             seed = 0, betas = c(-1.25), betat = c(-1.25), filter.surr = 1, filter.true = 1, 
                                             betas.init = c(-0.25), betat.init = c(-0.25), init.kappa = NULL, maxit = maxit,
                                             true.init.val = true.init.val, theta.copula = theta.copula, thetacopula.init = 3, 
                                             R2 = R2, typecopula = typecopula, ckappa = c(ckappa1, ckappa2), gamma.ui = 0.8, 
                                             nb.reject.data = nb.reject.data, equi.subj.trial = 0, prop.subj.trial = prop, 
                                             pfs = 0, print.iter = F)
      }

      if(program.to.estimate==6){ # Change of the weibull parameters in order to have times in year, with early censorship
        joint.simul2 <- jointSurroPenalSimul(nb.dataset = nsim, nbSubSimul=nsubjet, ntrialSimul = ntrial,
                                             int.method = int.method, nb.mc = nb.mc, nb.gh = nb.gh, nb.gh2 = nb.gh2, adaptatif = adaptatif,
                                             n.knots = nspline, kappa.use = kappa.use, type.joint.estim = type.joint.estim, print.iter = T,
                                             type.joint.simul = type.joint.simul, time.cens = 25, lambdas = 1.0, nus = 0.1185,
                                             lambdat = 0.8, nut = 0.1185, seed = 0, betas = c(-1.25), betat = c(-1.25), filter.surr = c(1),
                                             filter.true = c(1), betas.init = c(-0.25), betat.init = c(-0.25), init.kappa = NULL, maxit = maxit,
                                             true.init.val = true.init.val, theta.copula = theta.copula, thetacopula.init = 3, R2 = R2,
                                             typecopula = typecopula, ckappa = c(ckappa1, ckappa2), gamma.ui = 0.8, nb.reject.data = nb.reject.data, 
                                             pfs = pfs)
      }
      
      if(program.to.estimate==7){ # case with true alpha set to 0.1
        joint.simul2 <- jointSurroPenalSimul(nb.dataset = nsim, nbSubSimul=nsubjet, ntrialSimul = ntrial,
                                             int.method = int.method, nb.mc = nb.mc, nb.gh = nb.gh, nb.gh2 = nb.gh2, adaptatif = adaptatif,
                                             n.knots = nspline, kappa.use = kappa.use, type.joint.estim = type.joint.estim, print.iter = T,
                                             type.joint.simul = type.joint.simul, time.cens = time.cens, lambdas = 1.3, nus = 0.0025,
                                             lambdat = 1.1, nut = 0.0025, seed = 0, betas = c(-1.25), betat = c(-1.25), filter.surr = c(1),
                                             filter.true = c(1), betas.init = c(-0.25), betat.init = c(-0.25), init.kappa = NULL, maxit = maxit,
                                             true.init.val = true.init.val, theta.copula = theta.copula, thetacopula.init = 3, R2 = R2,
                                             typecopula = typecopula, ckappa = c(ckappa1, ckappa2), gamma.ui = 0.8, nb.reject.data = nb.reject.data, 
                                             pfs = pfs, alpha.ui = 0.1
        )
      }
      
      if(program.to.estimate==8){ # case with random censorship
        joint.simul2 <- jointSurroPenalSimul(nb.dataset = nsim, nbSubSimul=nsubjet, ntrialSimul = ntrial,
                                             int.method = int.method, nb.mc = nb.mc, nb.gh = nb.gh, nb.gh2 = nb.gh2, adaptatif = adaptatif,
                                             n.knots = nspline, kappa.use = kappa.use, type.joint.estim = type.joint.estim, print.iter = T,
                                             type.joint.simul = type.joint.simul, prop.cens = 1, time.cens = time.cens, lambdas = 1.3, nus = 0.0025,
                                             lambdat = 1.1, nut = 0.0025, seed = 0, betas = c(-1.25), betat = c(-1.25), filter.surr = c(1),
                                             filter.true = c(1), betas.init = c(-0.25), betat.init = c(-0.25), init.kappa = NULL, maxit = maxit,
                                             true.init.val = true.init.val, theta.copula = theta.copula, thetacopula.init = 3, R2 = R2,
                                             typecopula = typecopula, ckappa = c(ckappa1, ckappa2), gamma.ui = 0.8, nb.reject.data = nb.reject.data, 
                                             pfs = pfs, alpha.ui = 1
        )
      }
      
      if(program.to.estimate == 9){ # Variation of the number of subject per trial and random censorship
        # utils::data("dataOvarian", envir = environment(), package = "frailtypack")
        # consider the proportion from the advanced avarian cancer dataset including 50 trials
        # prop = table(dataOvarian$trialID)/nrow(dataOvarian)
        prop <- c(0.229865772, 0.104865772, 0.057885906, 0.098154362, 0.017617450, 0.006711409, 
                  0.041107383, 0.033557047, 0.005872483, 0.009228188, 0.020973154, 0.008389262,
                  0.014261745, 0.020973154, 0.009228188, 0.014261745, 0.006711409, 0.005033557, 
                  0.010906040, 0.015100671, 0.005033557, 0.001677852, 0.015939597, 0.012583893,
                  0.017617450, 0.012583893, 0.004194631, 0.010067114, 0.007550336, 0.002516779,
                  0.006711409, 0.010067114, 0.001677852, 0.014261745, 0.007550336, 0.033557047,
                  0.003355705, 0.002516779, 0.026006711, 0.002516779, 0.005033557, 0.003355705, 
                  0.004194631, 0.014261745, 0.010067114, 0.005033557, 0.004194631, 0.014261745,
                  0.002516779, 0.008389262)
        joint.simul2 <- jointSurroPenalSimul(nb.dataset = nsim, nbSubSimul = nsubjet, ntrialSimul = ntrial, 
                                             int.method = int.method, nb.mc = nb.mc, nb.gh = nb.gh, nb.gh2 = nb.gh2, 
                                             adaptatif = adaptatif, n.knots = nspline, kappa.use = kappa.use, 
                                             type.joint.estim = type.joint.estim, type.joint.simul = type.joint.simul, 
                                             time.cens = time.cens, lambdas = 1.3, nus = 0.0025, lambdat = 1.1, nut = 0.0025, 
                                             seed = 0, betas = c(-1.25), betat = c(-1.25), filter.surr = 1, filter.true = 1, 
                                             betas.init = c(-0.25), betat.init = c(-0.25), init.kappa = NULL, maxit = maxit,
                                             true.init.val = true.init.val, theta.copula = theta.copula, thetacopula.init = 3, 
                                             R2 = R2, typecopula = typecopula, ckappa = c(ckappa1, ckappa2), gamma.ui = 0.8, 
                                             nb.reject.data = nb.reject.data, equi.subj.trial = 0, prop.subj.trial = prop, 
                                             pfs = 0, print.iter = F, prop.cens = 1)
        
      }      
      
      # save of the RData
      save(joint.simul2, file = paste(wd,"/joint.simul2_", num.paquet, i, ".RData", sep = "")) 
    }
    
    i <- i + 1
  }
  
  # Merge from all packages
  joint.simul <- mergeJointSurroSimul(nb.packet = nb.paquet, envir.name = "joint.simul2_", envir.num.base = num.paquet,
                                                    wd = wd)
  setwd(current)
  
  return(joint.simul)
}