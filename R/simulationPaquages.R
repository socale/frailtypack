
# This function aims to perform simulation studies using small paquages and save interim results in the .RData files
# All the ".RData" are save in ./paquageSimul"num.paquet", where "num.paquet" is a number give to the simulation design 

simulationPaquages <- function(nsim = 500, nsubjet = 600, ntrial = 30, int.method = 0,
                               nb.mc = 1000, nb.gh = 20, nb.gh2 = 32, adaptatif = 0,
                               nspline = 6, kappa.use = 4, type.joint.estim = 3, 
                               typecopula = 1, type.joint.simul = 3, theta.copula = 3, 
                               time.cens = 349, true.init.val = 1, R2 = 0.81, maxit = 40, 
                               nb.paquet = 10, num.paquet = 1, program.to.estimate = 1, 
                               list.paquage = NULL
){
  
  # =====description of the parameters=====
   # nsim : number of dataset to generate
   # nsubjet : number of subjects
   # ntrial : number of trials
   # nspline : number of knots for spline
   # nb.paquet : number of paquages to considered
   # num.paquet : a number give to the simulation design
   # program.to.estimate : 1 by default, 2 = case with two covariates, 4 = case with init values for theta set to 1 and alpha set to 0.5
     # 5 = case with Variation of the number of subject per trial
   # list.paquage: list of the paquages numbers to execute. If NULL, all nb.paquage are executed
   # More details for all the arguments can be found by displaying the help on the "jointSurroPenalSimul" function
      # i.e. help("jointSurroPenalSimul")
  # =======end description============
  
  # some initializations
  ckappa1 <- 0 
  ckappa2 <- 0
  pfs     <- 0
  
  #Management of the list of paquages on which simulation studies should be performed
  
  if(is.null(list.paquage)){
    list.paquages <- seq(1,nb.paquet)
  }else{
    list.paquages <- sort(list.paquage)
  }
  
  # creation of the work dicrectory for the results
  current <- getwd()
  wd <- paste("./paquageSimul",num.paquet,sep="")
  dir.create(wd)
  
  i <- 1
  nb.reject.data <- 0
  N.sim.data     <- 0
  
  while(i <= length(list.paquages)){
    if(i %in% list.paquages){
      # initialization of require arguments: number of dataset to consider per paquage and number of dataset to reject before simulation
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
        
        # number of dataset to consider per paquage
        if(modulo >= i){
          N.sim.data <- q + 1
        }else{
          N.sim.data <- q
        }
      }
      
      # Simulation for the ongoing paquage
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
                                             nb.reject.data = nb.reject.data, pfs = pfs, print.iter = T)
      }
      
      if(program.to.estimate == 5){ # Variation of the number of subject per trial
        data(dataOvarian, package = "frailtypack")
        # consider the proportion from the advanced avarian cancer dataset including 50 trials
        prop = (tab1(dataOvarian$trialID,graph = F)$output.table[-51,2])/100
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

    }
    
    i <- i + 1
    
    # save of the RData
    save.image(paste("wd/joint.simul2_",numsimul, ".RData", sep = ""))
  }
  
  # Merge from all paquages
  joint.simul <- frailtypack:::mergeJointSurroSimul(nb.packet = nb.paquet, envir.name = "joint.simul2_", envir.num.base = num.paquet,
                                                    wd = wd)
  setwd(current)
  
  return(joint.simul)
}