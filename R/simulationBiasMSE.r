# Function used to summarize the simulation results. results in this case just include Biase and MSE
# param.estim is the dataset of the estimates from the simulation study
simulationBiasMSE <- function(param.estim, R2 = 0.81, ktau = 0.378, nb.simul, dec = 3){
  
  # recherche des lignes correspondantes aux simulations qui n'ont pas convergees (concerne le cas des programmes esimees par MPI-OpenMP)
  somme_row <- NULL
  for(i in 1:nrow(param.estim))
    somme_row[i] <- sum(param.estim[i,])
  
  #on les exclu du dataframe
  donnee <- param.estim[somme_row!=0,]
  
  
  d <- data.frame(parametre = c("Ktau","R2_adj","R:n(%)"))
  d$True <- c(ktau, R2, "")
  d$Mean <- NA
  d$Biais <- NA
  d$MSE <- NA
  #Ktau
  m_ktau <- mean(donnee$tau[!is.na(donnee$tau)])
  b_tau <- ktau-m_ktau # biais Ktau
  mse_tau <- mean((ktau-donnee$tau[!is.na(donnee$tau)])^2) # MSE
                
  #R2_adj
  m_r2_ad <- mean(donnee$R2trial[!is.na(donnee$R2trial)])
  b_r2_ad <- R2-m_r2_ad # biais Ktau
  mse_r2_ad <- mean((R2-donnee$R2trial[!is.na(donnee$R2trial)])^2) # MSE
                
  #NA
  n_NA <- nb.simul - nrow(donnee)
  cat(c("nb.simul", nb.simul) , "\n")
  cat(c("nrow(donnee)", nrow(donnee)), "\n")
  # affichage
  #cat(c("Ktau",round(m_ktau,3),round(b_tau,3),round(mse_tau,3)),sep='&',fill=T)
  #cat(c("R2_adj",round(m_r2_ad,3),round(b_r2_ad,3),round(mse_r2_ad,3)),sep='&',fill=T)
  #cat(c("R:n(%)","",paste(n_NA,"(",round(100*n_NA/nrow(donnee)),")",sep=""),""),sep='&',fill=T)
                
  #sauvagarde dans d
  d[1,3:5] <- c(round(m_ktau, dec),round(b_tau, dec),round(mse_tau, dec))
  d[2,3:5] <- c(round(m_r2_ad, dec),round(b_r2_ad, dec),round(mse_r2_ad, dec))
  #d[3,3:5] <- c(NA,paste(n_NA,"(",round(100*n_NA/nrow(donnee)),")",sep = ""),NA)
  d[3,3:5] <- c("",paste(n_NA,"(",round(100*n_NA/nrow(donnee)),")",sep = ""),"")

  return(d)
}