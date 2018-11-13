# calcul des taux de kendall a partir des echantillons bootstrat
kendall_bootstrap=function(fichierTau, nsim=500,nboot=1000){
  #fichierTau: fichier concerne contenant les echantillons bootstrap
  #nsim= nombre de simulation 
  #nboot= nombre d'echantillon bootstrap par simulation
  
  a=utils::read.table(fichierTau,sep="")
  if(ncol(a)==1){# on fait les calculs uniquement si on a 3 colonne dans le fichier, sinon les calculs on deja ete faits
    result=data.frame(matrix(c(NA,NA,NA),c(1,3)))
    
    for(i in 1:nsim){
      j=(i-1)*nboot+1
      d=a[c(j:(j+nboot-1)),1]
      cat(c(j,(j+nboot-1)),sep=" ",fill=T)
      result[i,1]=mean(d)
      result[i,2:3]=quantile(d,probs=c(0.025,0.975)) 
    }
    utils::write.table(result,fichierTau,sep=" ",row.names = F,col.names = F)
    return(result)
  }
  else return(a)
}