# autres fonctions
param.empirique = function(nsim = 100, ver = 2, dec = 2, variatio.seed = 0,
                           n.obs = 600, n.trial = 30, prop.cens = 1, cens.adm=549, lambda.S = 1.3,
                           nu.S = 0.0025,lambda.T = 1.1, nu.T = 0.0025, 
                           seed = 0,alpha = 1.5, gamma = 2.5, sigma.s = 0.7, sigma.t = 0.7,
                           cor = 0.8, betas = c(-1.25, 0.5), betat = c(-1.25, 0.5), 
                           filter.surr = c(1,1), filter.true = c(1,1), frailt.base = 1, typecopula = 1,
                           thetacopule = 6, random.generator = 1, prop.i = rep(1/n.trial, n.trial),
                           nb.reject.data = 0){
  # variatio.seed : si = 1, je fais varier le seed, seulement pour des fins de verification des statistiques empirique.
  # afin de reproduire les jeux de donnees utilisees dans les simulations, il faut plutot faire varier nb.reject.data et fixer le seed, 
  # et par consequent cet argument doit avoir la valeur 0. Toutefois on retient que les stat empirique sont meilleures
  # lorque cet argumet est fixe a 0
  nbre.covar = ver
  full.data = 1
  np = 12 + nbre.covar
  d = data.frame(matrix(0, nrow = nsim, ncol = np))
  if(nbre.covar > 1) 
    names(d) = c("MuvS", "sigmaS", "MuvT", "sigmaT","SigmaST","Muui", "gamma",
                 "median.S", "median.T", "prop.S", "propT", "prop.trt",
               paste("propVar",seq(2,nbre.covar),sep=""), "Ktau")
  else
    names(d) = c("MuvS", "sigmaS", "MuvT", "sigmaT","SigmaST","Muui", "gamma", 
                 "median.S", "median.T", "prop.S", "propT", "prop.trt", "Ktau")
  
  # nombre de sujets par essai:
  n_i <- as.integer(n.obs*prop.i) # nombre de sujet par essai
  
  min_n <- min(n_i)
  max_n <- max(n_i)
  if(sum(n_i) < n.obs){
    n_i[n_i == min_n][1] <- n_i[n_i == min_n][1] + (n.obs-sum(n_i)) # on ajoute au premier essai de plus petite taille le nombre d'individu non encore affecte (1 generalement) a cause des problemes d'arrondi
  }
  if(sum(n_i) > n.obs) { 
    n_i[n_i == max_n][1] <- n_i[n_i == max_n][1]-(sum(n_i)-n.obs) # on soustrait au premier essai de plus grande taille le nombre d'individu affecte en trop (1 generalement) a cause des problemes d'arrondi
  }
  
  sujet.essai <- n_i
  pour.essai <- rep(0, n.trial)# vecteur des positions des premiers sujets par essai 
  for (j in 1: n.trial) pour.essai[j] <- (j-1)*sujet.essai[j] +1
  
  for (i in 1:nsim){
    if(variatio.seed == 1){
      seed1 <- seed + i 
      nb.reject.data1 <- nb.reject.data
    }
    else{
      seed1 <- seed 
      nb.reject.data1 <-  nb.reject.data + i -1
    }
    data.sim <- jointSurrCopSimul(n.obs = n.obs, n.trial = n.trial, prop.cens = prop.cens, cens.adm = cens.adm, 
                   lambda.S = lambda.S, nu.S = nu.S, lambda.T = lambda.T, nu.T = nu.T, full.data = full.data,
                   seed = seed1, nb.reject.data = nb.reject.data1, alpha = alpha, gamma = gamma, 
                   sigma.s = sigma.s, sigma.t = sigma.t, filter.surr = filter.surr,
                   cor = cor, betas = betas, betat = betat, filter.true= filter.true,
                   frailt.base = frailt.base, thetacopule = thetacopule, ver = ver,
                   random.generator = random.generator
                   )
                   
                                  
    d[i,1] <- mean(data.sim$v_Si[pour.essai])
    d[i,2] <- stats::var(data.sim$v_Si[pour.essai])
    d[i,3] <- mean(data.sim$v_Ti[pour.essai])
    d[i,4] <- stats::var(data.sim$v_Ti[pour.essai])
    d[i,5] <- stats::cov(data.sim$v_Si[pour.essai],data.sim$v_Ti[pour.essai])
    d[i,6] <- mean(data.sim$u_i[pour.essai])
    d[i,7] <- stats::var(data.sim$u_i[pour.essai])
    d[i,8] <- stats::median(data.sim$timeS)
    d[i,9] <- stats::median(data.sim$timeT)
    d[i,10] <- prop.table(table(data.sim$statusS))["1"]
    d[i,11] <- prop.table(table(data.sim$statusT))["1"]
    d[i,12] <- prop.table(table(data.sim$trt))["1"]
    if(nbre.covar >1)
      for(j in 1:(nbre.covar-1)){
        
        d[i,12+j] <- prop.table(table(data.sim[,10+j]))["1"]
      }
    
    ## Compute Kendall's tau at each study ##
    Tau <- numeric(n.trial)
    for(k in 1:n.trial){
      Tau[k] <- stats::cor(data.sim$timeS[(data.sim$trialID == k) & (data.sim$trt == 1)],
                 data.sim$timeT[(data.sim$trialID == k) & (data.sim$trt == 1)],method="kendall")
    }
    
    ## difference at most 0.01##
    d[i,np] <- mean(Tau)
    #print(d[i,])
    
  }
  result = data.frame(names(d))
  names(result) = "Parameters"
  if(typecopula == 1)
    result$True <- c(0, sigma.s, 0, sigma.t, round(cor * sqrt(sigma.s * sigma.t),dec), 
                     0, gamma, "-", "-", "-", "-", rep(0.5, nbre.covar), thetacopule/(thetacopule+2))
  
  if(typecopula == 2)
    result$True <- c(0, sigma.s, 0, sigma.t, round(cor * sqrt(sigma.s * sigma.t),dec), 
                     0, gamma, "-", "-", "-", "-", rep(0.5, nbre.covar), thetacopule/(thetacopule+1))
  
  result$Mean <- NA
  result$Median <- NA
  result$SD <- NA
  
  for(k in 1 : np){
    result$Mean[k] <- round(mean(stats::na.omit(d[,k])),dec)
    result$Median[k] <- round(stats::median(stats::na.omit(d[,k])),dec)
    result$SD[k] <- round(stats::sd(stats::na.omit(d[,k])),dec)
  }
  
  #prettyR::describe(d)
  return(result)
}


# calcul de l'image de l'effet du traitement sur T sachant l'effet du traitement sur S (x)
f <- function(x, object, var.used, alpha., pred.int.use){
  beta  <- object$beta.t
  alpha <- object$beta.s
  
  #===proble dans l'indexation de la matrice des coefficient. scl: 11/04/2020
  # dab   <- object$Coefficients$Estimate[nrow(object$Coefficients)-4]
  # daa   <- object$Coefficients$Estimate[nrow(object$Coefficients)-6]
  # dbb   <- object$Coefficients$Estimate[nrow(object$Coefficients)-5]
  dab   <- object$Coefficients[rownames(object$Coefficients)=="sigma_sT",1]
  daa   <- object$Coefficients[rownames(object$Coefficients)=="sigma_s",1]
  dbb   <- object$Coefficients[rownames(object$Coefficients)=="sigma_t",1]
  #====
  
  
  #x <- matrixPred$beta.S[i]
  x.     <- t(matrix(c(1, -dab/daa),1,2))
  
  # Vmu (sigma_ST, sigma_SS). on utilise la matrice obtenu par delta methode a partir de la hessienne
  Vmu   <- matrix(c(object$varcov.Sigma[3,3], object$varcov.Sigma[3,1],
                    object$varcov.Sigma[3,1], object$varcov.Sigma[1,1]),2,2)
  nparam <- nrow(object$varH)
  
  # VD (bete_T, beta_S). on utilise la hesienne directement car pas de changement de variable
  VD     <- matrix(c(object$varH[nparam,nparam], object$varH[nparam,nparam-1],
                     object$varH[nparam-1,nparam], object$varH[nparam -1,nparam - 1]),2,2)
  R2trial <- object$Coefficients[rownames(object$Coefficients)=="R2trial",1] # scl: 11/04/2020
  
  if(var.used == "error.estim") {
    if(pred.int.use == "lw"){
      return(
        sapply(x, function(x, beta, dab, daa, dbb, alpha, alpha., x., Vmu, VD, R2trial) 
          beta + (dab/daa) * (x - alpha) - qnorm(1-alpha./2) * 
            sqrt(t(x.) %*% (Vmu + (((x - alpha)/daa)**2) * VD) %*% x.
                 + dbb * (1 - R2trial)), 
          beta = beta, dab = dab, daa = daa, dbb = dbb, alpha = alpha, alpha. = alpha., 
          x. = x., Vmu = Vmu, VD = VD, R2trial = R2trial
        )
      )
    }
    else{
      return(
        sapply(x, function(x, beta, dab, daa, dbb, alpha, alpha., x., Vmu, VD, R2trial) 
          beta + (dab/daa) * (x - alpha) + qnorm(1-alpha./2) * sqrt(t(x.) %*% 
                                                                      (Vmu + (((x - alpha)/daa)**2) * VD) %*% x. + dbb * (1 - R2trial)), 
          beta = beta, dab = dab, daa = daa, dbb = dbb, alpha = alpha, alpha. = alpha., 
          x. = x., Vmu = Vmu, VD = VD, R2trial = R2trial
        )
      )
    }
  }
  else{
    if(pred.int.use == "lw"){
      return(
        sapply(x, function(x, beta, dab, daa, dbb, alpha, alpha., R2trial) 
          beta + (dab/daa) * (x - alpha) - qnorm(1-alpha./2) * 
            sqrt(dbb * (1 - R2trial)),
          beta = beta, dab = dab, daa = daa, dbb = dbb, alpha = alpha, alpha. = alpha., 
          R2trial = R2trial
        )
      )
    }
    else{
      return(
        
        sapply(x, function(x, beta, dab, daa, dbb, alpha, alpha., R2trial) 
          beta + (dab/daa) * (x - alpha) + qnorm(1-alpha./2) *
            sqrt(dbb * (1 - R2trial)),
          beta = beta, dab = dab, daa = daa, dbb = dbb, alpha = alpha, alpha. = alpha., 
          R2trial = R2trial
        )
      )
    }
  }
  
}


