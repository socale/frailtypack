#==================================================================================
# fonction pour l'estimation des kappa pour le risque de base par validation croisee
#==================================================================================

#Used to estimate smoothing parameters using cross-validation
#
#Function used to estimated smoothing parameters associated with the surrogate endpoint and the 
#true endpoint by cross-validation. Are used two separate cox proportional hazard models.
#This function created a file named \bold{kappa_val_croisee.txt} which contains the obtained
#smoothing parameters. This is particularly useful when doing simulations study
# 
# @aliases kappa_val_croisee
# @param don_S h
# @param don_T h
# @param njeu h
# @param n_obs h
# @param n_node h
# @param adjust_S h
# @param adjust_T h
# @param kapp_0 h
# @param kappa_1 h
# @param kappa_2 h
# @param print.times h
#
# @return Create a text file named \bold{kappa_val_croisee.txt} in the work directory containing 
# all estimated smoothing parameters. in addition, this function return an equivalent dataframe
# of two coluumns for kappa associated with the surrogate endpoint and the true endpoint.
#
# 
kappa_val_croisee <- function(don_S, don_T, njeu, n_obs, n_node = 6, adjust_S = 1, adjust_T = 1,
                           kapp_0 = 100, kappa_1 = 955000, kappa_2 = 975000, print.times = T,
                           scale = 1){
  # don_S et don_T: les "njeu" jeux de donnees pour lesquelles il faut estimer les kappa, toutes dans un seul jeu de donnees
  # njeu: represente le nombre de data a considerer dans don
  # n_obs: nombre d'observation par fichier de donnees
  # n_node: nombre de noeuds spline
  # adjust: le terme d'ajustement, on pultipli kappa par adjust
  # kapp_0: valeur attribuee a kappa s'il est inferieur a 0
  # kappa_1= valeur initiale de kappa 1
  # kappa_2= valeur initiale de kappa 2
  
  # gestion de lechelle
  done_S <- don_S
  done_T <- don_T
  done_S$timeS <- done_S$timeS/scale
  done_T$timeT <- done_T$timeT/scale
  
  kapa=matrix(0.0, nrow = njeu, ncol = 2)
  j=1
  k=1
  for(i in 1:njeu){
    j <- i*n_obs
    cox_surr <- try(frailtyPenal(Surv(timeS,statusS)~1, data = done_S[k:j,],cross.validation = T,
                              n.knots = n_node, kappa = kappa_1, print.times = print.times), silent = TRUE)
    cox_true <- try(frailtyPenal(Surv(timeT,statusT)~1, data = done_T[k:j,], cross.validation = T,
                              n.knots = n_node, kappa=kappa_2, print.times = print.times), silent = TRUE)
    
    if((class(cox_surr)=="try-error") | (class(cox_true)=="try-error")){
      if(i==1){
        #cat("probleme d'etimiation avec ce jeu de donnee, affectation du kappa par defaut")
        kapa[i,1] <- kappa_1
        kapa[i,2] <- kappa_2
      }else{
        #cat("probleme d'etimiation avec ce jeu de donnee, affectation du kappa precedent")
        kapa[i,1] <- kapa[i-1,1]
        kapa[i,2] <- kapa[i-1,2]
      }
    }else{
      kapa[i,1] <- cox_surr$kappa
      kapa[i,2] <- cox_true$kappa
    }
    
    if(!(kapp_0==0)){
      if (kapa[i,1]<1) kapa[i,1] <- kapp_0
    
      if (kapa[i,2]<1) kapa[i,2] <- kapp_0
    }
    
    k=j+1
    cat(i,"k=", kapa[i,1]*adjust_S, kapa[i,2]*adjust_T, fill=T)
  }
  kapa[,1] <- kapa[,1]*adjust_S
  kapa[,2] <- kapa[,2]*adjust_T
  utils::write.table(kapa, "kappa_valid_crois.txt", sep=" ", row.names = F, col.names = F)
  return(kapa)
}