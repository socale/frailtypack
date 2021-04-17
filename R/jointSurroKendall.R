
# test de gitHub et Git IB


jointSurroKendall <- function(theta, gamma, alpha = 1, eta = 1, adaptative = 0, npg = 20, ui = 1,
                              ui.chap.Essai=rep(0,4), invBi.chol.Essai.k = matrix(rep(0,16), nrow = 4, ncol = 4)){
  
  # recherche des points et poids de quadrature
  
  # nbre_coeur <- parallel::detectCores()
  # ncores <- min(nbre_coeur-1,npg)
  # # initialisation de l'envireonnement de parallelisation
  # cl <- parallel::makeCluster(ncores)
  # # declaration des clusters precedemment initialises
  # doParallel::registerDoParallel(cl)
  # # separation des donnees en autant de groupe que de coeurs
  # vect_iter <- itertools::isplitIndices(n = npg, chunks = ncores)
  
  
  xx1 <- statmod::gauss.quad(npg,kind="hermite")$nodes
  ww1 <- statmod::gauss.quad(npg,kind="hermite")$weights * exp(xx1**2) # w = w*exp(x**2)
  
  if(ui == 1){ # model complet avec prise en compte de l'heterogeneite sur les risque de base
    xxl <- rep(0,4)
    ss <- 0
    for(ll in 1:npg){
      ss3 <- 0
      for(kk in 1:npg){
        ss2 <- 0
        for(jj in 1:npg){
          ss1 <- 0
          for(ii in 1:npg){
            xxl[1] <- xx1[ll]
            xxl[2] <- xx1[kk]        
            xxl[3] <- xx1[jj]
            xxl[4] <- xx1[ii]
        
            # changement de variable en cas de quadrature adaptative
            if(adaptative==1){ # on effectue le changement de variable
              m <- invBi.chol.Essai.k %*% xxl
              xxl <- ui.chap.Essai+sqrt(2)*m
            }
        
            auxfunca <- funcJointSurroKendall(xxl[1], xxl[2], xxl[3], xxl[4], theta, gamma, alpha, eta, 1)
            ss1 <- ss1 + ww1[ii] * auxfunca
          }
          ss2 <- ss2 + ww1[jj] * ss1
        }
        ss3 <- ss3 + ww1[kk] * ss2
      }
      ss <- ss + ww1[ll] * ss3
    }
  }
  else{ # dans ce cas, ui n'est pas pris en compte dans le modele
      ui.chap.Essai <- ui.chap.Essai[1:2]
      invBi.chol.Essai.k <- invBi.chol.Essai.k[1:2,1:2]
      xxl <- rep(0,2)
      ss <- 0
        for(jj in 1:npg){
          ss1 <- 0
          for(ii in 1:npg){
            xxl[1] <- xx1[jj]
            xxl[2] <- xx1[ii]        
            
            #changement de variable en cas de quadrature adaptative
            if(adaptative==1){ # on effectue le changement de variable
              m <- invBi.chol.Essai.k %*% xxl
              xxl <- ui.chap.Essai+sqrt(2)*m        
            }
            
            auxfunca <- funcJointSurroKendall(0,xxl[1],0, xxl[2], theta, gamma, alpha, eta, 0)
            ss1 <- ss1 + ww1[ii] * auxfunca
          }
        ss <- ss + ww1[jj] * ss1
      }
  }
  
  return(2*ss -1)
}