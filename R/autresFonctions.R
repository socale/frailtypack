# autres fonctions
param.empirique = function(nsim = 100, nbre.covar = 2, dec = 2, variatio.seed = 1,
                           n.obs = 600, n.trial = 30, cens.adm=549, lambda.S = 1.3,
                           nu.S = 0.0025,lambda.T = 1.1, nu.T = 0.0025,full.data = 1, 
                           seed = 0,alpha = 1.5, gamma = 2.5, sigma.s = 0.7, sigma.t = 0.7,
                           rsqrt = 0.8, betas = c(-1.25, 0.5), betat = c(-1.25, 0.5), 
                           filter.surr = c(1,1), filter.true = c(1,1), frailt.base = 1,
                           thetacopule = 6, ver =2){
  # variatio.seed : si = 1, je fais varier le seed, seulement pour des fins de verification des statistiques empirique.
  # afin de reproduire les jeux de donnees utilisees dans le simulations, il faut plutot faire varier nb.reject.data et fixer le seed, 
  # et par consedent cet argument doit prendre pour valeur 0. Toutefois on retient que les stat empirique sont meilleures
  # lorque cet argumet est fixe a 0
  
  np = 11 + nbre.covar
  d = data.frame(matrix(0, nrow = nsim, ncol = np))
  if(nbre.covar > 1) 
    names(d) = c("MuvS", "sigmaS", "MuvT", "sigmaT","SigmaST","Muui", "gamma",
                 "median.S", "median.T", "prop.S", "propT", "prop.trt",
               paste("propVar",seq(2,nbre.covar),sep=""))
  else
    names(d) = c("MuvS", "sigmaS", "MuvT", "sigmaT","SigmaST","Muui", "gamma", 
                 "median.S", "median.T", "prop.S", "propT", "prop.trt")
  
  for (i in 1:nsim){
    if(variatio.seed == 1){
      seed1 <- seed + i 
      nb.reject.data1 <- 0
    }
    else{
      seed1 <- seed 
      nb.reject.data1 <-  i -1
    }
    data.sim <- jointSurrCopSimul(n.obs = n.obs, n.trial = n.trial, cens.adm = cens.adm, lambda.S = lambda.S, nu.S = nu.S, 
                                  lambda.T = lambda.T, nu.T = nu.T, full.data = full.data,
                                  seed = seed1, nb.reject.data = nb.reject.data1, alpha = alpha, gamma = gamma, 
                                  sigma.s = sigma.s, sigma.t = sigma.t, filter.surr = filter.surr,
                                  rsqrt = rsqrt, betas = betas, betat = betat, filter.true= filter.true,
                                  frailt.base = frailt.base, thetacopule = thetacopule, ver = ver)
    d[i,1] <- mean(data.sim$v_Si)
    d[i,2] <- (sd(data.sim$v_Si))**2
    d[i,3] <- mean(data.sim$v_Ti)
    d[i,4] <- (sd(data.sim$v_Ti))**2
    d[i,5] <- cov(data.sim$v_Si,data.sim$v_Ti)
    d[i,6] <- mean(data.sim$u_i)
    d[i,7] <- mean(sd(data.sim$u_i))**2
    d[i,8] <- median(data.sim$timeS)
    d[i,9] <- median(data.sim$timeT)
    d[i,10] <- prop.table(table(data.sim$statusS))["1"]
    d[i,11] <- prop.table(table(data.sim$statusT))["1"]
    d[i,12] <- prop.table(table(data.sim$trt))["1"]
    if(nbre.covar >1)
      for(j in 1:(nbre.covar-1)){
        
        d[i,12+j] <- prop.table(table(data.sim[,10+j]))["1"]
      }
  }
  result = data.frame(names(d))
  names(result) = "Parameters"
  result$True <- c(0, sigma.s, 0, sigma.t, round(sqrt(rsqrt * sigma.s * sigma.t),dec), 
                   0, gamma, "-", "-", "-", "-", rep(0.5, nbre.covar))
  result$Mean <- NA
  result$Median <- NA
  result$SD <- NA
  
  for(k in 1 : np){
    result$Mean[k] <- round(mean(d[,k]),dec)
    result$Median[k] <- round(median(d[,k]),dec)
    result$SD[k] <- round(sd(d[,k]),dec)
  }
  
  #prettyR::describe(d)
  return(result)
}


