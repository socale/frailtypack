# synthese des resultats
synthese_result_modele_reduit=function(param_esti,ktauboot,R2boot,nb_paquet=1,ndec=3,nsim=100,
                                        theta_S=1,zeta=1,gamma_S=0.8,alpha=1,sigma_S=0.7,sigma_T=0.7,
                                        sigma_ST=0.6,beta_S=-1.25,beta_T=-1.25,R2trial=0.36,tau=0.378,
                                        n_bootstrap=1000,ick=0, R2parboot = 0, type.joint = 1){
  # wd= repertoire dans lequel se trouve les scripts R
  # rep_courant=  chemin d'access au repertoire contenant les sous dossiers des paquets de simultation
  # nb_paquet= nombre de pacquets de donnees consideres
  # ndec= nombre de parties decimale pour les parametres
  # nsim=nombre de simultations
  # param_init_R= vecteur des parametres initiaux medele reduit
  # n_bootstrap= nombre d'echantillons pour le bootstrap
  # ick= dit si on tient compte du calcul de l'IC du tau de kendall(1) ou non (0)
  # R2parboot =  dit si on calcul le taux de couverture du R2 par bootstrap (1) ou par delta-method (0)
  # type.joint = modele joint surrogate(1) ou modele joint frailty-copula(3)
  # setwd(wd)
  # source("fusion_resultats_simul.r")
  # source("kendall_bootstrap.r")
  # source("fusion_resultats_Tkandal.r")

  # synthese resultats simulations avec 500 boucles MC
  
  # estimates du model complet
    param_esti = param_esti
    
    # recherche des lignes correspondantes aux simulations qui n'ont pas convergees (concerne le cas des programmes esimees par MPI-OpenMP)
    somme_row=NULL
    for(i in 1:nrow(param_esti))
      somme_row[i]=sum(param_esti[i,])
    
    #on les exclu du dataframe
    estimates_complet2=param_esti[somme_row!=0,]
    
    ves <- length(beta_S)
    vet <- length(beta_T)
    #str(estimates_complet2)
    if(type.joint == 1){
      names(estimates_complet2)=c("Theta_S","se_theta_S","zeta","se_zeta","beta_S","se_beta_S","beta_T","se_beta_T","sigma_S",
                                  "se_sigma_S","sigma_T","se_sigma_T","sigma_ST","se_sigma_ST","gamma_S","se_gamma_S","alpha","se_alpha",
                                  "R2trial","se_R2trial","tau_00")
    }else{
      
      entete <- c("Theta_S","se_theta_S","zeta","se_zeta","beta_S","se_beta_S","beta_T","se_beta_T","sigma_S",
                  "se_sigma_S","sigma_T","se_sigma_T","sigma_ST","se_sigma_ST","gamma_S","se_gamma_S","alpha","se_alpha",
                  "R2trial","se_R2trial","tau_00","SE.KendTau")
      if(ves>1){
        for(h in 2:ves){
          entete <- c(entete, paste("beta_S_", h, sep = ""))
          entete <- c(entete, paste("se_beta_S_", h, sep = ""))
        }
      }
      if(vet>1){
        for(h in 2:vet){
          entete <- c(entete, paste("beta_T_", h, sep = ""))
          entete <- c(entete, paste("se_beta_T_", h, sep = ""))
        }
      }
      names(estimates_complet2) <- entete
    }
    
    summary(estimates_complet2)
    
    # ajout des coefficient de correlation
    #       estimates_complet2$rho_ui=estimates_complet2$gamma_ST/sqrt(estimates_complet2$gamma_S*estimates_complet2$gamma_T)
    #       estimates_complet2$rho_ui2=(estimates_complet2$rho_ui)^2
    #       estimates_complet2$rho_vi=estimates_complet2$sigma_ST/sqrt(estimates_complet2$sigma_S*estimates_complet2$sigma_T)
    #       estimates_complet2$rho_vi2=(estimates_complet2$rho_vi)^2
    #       estimates_complet2$rho_wij=estimates_complet2$Theta_ST/sqrt(estimates_complet2$Theta_S*estimates_complet2$Theta_T)
    #       estimates_complet2$rho_wij2=(estimates_complet2$rho_wij)^2
    
    #ajout des taux de couverture
    param_esti=estimates_complet2
    
    # parametres de simulation
    #param_init=c("theta_S"=1,"theta_T"=0.8,"theta_ST"=0,"gamma_S"=0.8,"gamma_T"=0.7,"gamma_ST"=0,"sigma_S"=0.7,"sigma_T"=0.7,
    #             "sigma_ST"=0,"beta_S"=-1.25,"beta_T"=-1.25,"rho_wijST"=0.8,"rho2_wijST"=0.64,"rho_vi_st"=0.6,"rho2_vi_st"=0.36,
    #             "rho_ui_st"=0.6,"rho2_ui_st"=0.36)
    #       param_init=c("theta_S"=1,"theta_T"=0.8,"theta_ST"=0,"gamma_S"=0.8,"gamma_T"=0.7,"gamma_ST"=0,"sigma_S"=0.7,"sigma_T"=0.7,
    #                    "sigma_ST"=0,"beta_S"=-1.25,"beta_T"=-1.25,"R2trial"=0.36,"rho_wijST"=0.8,"rho2_wijST"=0.64,"rho_vi_st"=0.6,
    #                    "rho2_vi_st"=0.36,"rho_ui_st"=0.6,"rho2_ui_st"=0.36,"tau_11"=0.24877,"tau_10"=0.22286,"tau_01"=0.21387,"tau_00"=0.243328)
    
    param_init=c("theta_S"=theta_S,"zeta"=zeta,"gamma_S"=gamma_S,"alpha"=alpha,"sigma_S"=sigma_S,"sigma_T"=sigma_T,
                 "sigma_ST"=sigma_ST,"beta_S"=beta_S,"beta_T"=beta_T,"R2trial"=R2trial,"tau"=tau)
    
    #       param_init["sigma_ST"]=param_init["rho_vi_st"]*sqrt(param_init["sigma_S"])*sqrt(param_init["sigma_T"])
    #       param_init["gamma_ST"]=param_init["rho_ui_st"]*sqrt(param_init["gamma_S"])*sqrt(param_init["gamma_T"])
    #       param_init["theta_ST"]=param_init["rho_wijST"]*sqrt(param_init["theta_S"])*sqrt(param_init["theta_T"])
    #       
    # taux de couverture thetaS
    param_esti$bi_se_theta=param_esti$Theta_S-1.96*param_esti$se_theta_S
    param_esti$bs_se_theta=param_esti$Theta_S+1.96*param_esti$se_theta_S
    param_esti$couverture_thetaS=ifelse((param_init["theta_S"]>=param_esti$bi_se_theta) & (param_init["theta_S"]<=param_esti$bs_se_theta),1,0)

    # taux de couverture zeta
    param_esti$bi_se_theta=param_esti$zeta-1.96*param_esti$se_zeta
    param_esti$bs_se_theta=param_esti$zeta+1.96*param_esti$se_zeta
    param_esti$couverture_zeta=ifelse((param_init["zeta"]>=param_esti$bi_se_theta) & (param_init["zeta"]<=param_esti$bs_se_theta),1,0)
    
    
    # taux de couverture sigma_S
    param_esti$bi_se_theta=param_esti$sigma_S-1.96*param_esti$se_sigma_S
    param_esti$bs_se_theta=param_esti$sigma_S+1.96*param_esti$se_sigma_S
    param_esti$couverture_sigma_S=ifelse((param_init["sigma_S"]>=param_esti$bi_se_theta) & (param_init["sigma_S"]<=param_esti$bs_se_theta),1,0)
    
    # taux de couverture sigma_T
    param_esti$bi_se_theta=param_esti$sigma_T-1.96*param_esti$se_sigma_T
    param_esti$bs_se_theta=param_esti$sigma_T+1.96*param_esti$se_sigma_T
    param_esti$couverture_sigma_T=ifelse((param_init["sigma_T"]>=param_esti$bi_se_theta) & (param_init["sigma_T"]<=param_esti$bs_se_theta),1,0)
    
    # taux de couverture sigma_ST
    param_esti$bi_se_theta=param_esti$sigma_ST-1.96*param_esti$se_sigma_ST
    param_esti$bs_se_theta=param_esti$sigma_ST+1.96*param_esti$se_sigma_ST
    param_esti$couverture_sigma_ST=ifelse((param_init["sigma_ST"]>=param_esti$bi_se_theta) & (param_init["sigma_ST"]<=param_esti$bs_se_theta),1,0)
    
    # taux de couverture gamma_S
    param_esti$bi_se_theta=param_esti$gamma_S-1.96*param_esti$se_gamma_S
    param_esti$bs_se_theta=param_esti$gamma_S+1.96*param_esti$se_gamma_S
    param_esti$couverture_gamma_S=ifelse((param_init["gamma_S"]>=param_esti$bi_se_theta) & (param_init["gamma_S"]<=param_esti$bs_se_theta),1,0)
    
    # taux de couverture theta
    param_esti$bi_se_theta=param_esti$alpha-1.96*param_esti$se_alpha
    param_esti$bs_se_theta=param_esti$alpha+1.96*param_esti$se_alpha
    param_esti$couverture_alpha=ifelse((param_init["alpha"]>=param_esti$bi_se_theta) & (param_init["alpha"]<=param_esti$bs_se_theta),1,0)
    
    # taux de couverture beta_S
    if(ves>1){
      param_esti$bi_se_theta=param_esti$beta_S-1.96*param_esti$se_beta_S
      param_esti$bs_se_theta=param_esti$beta_S+1.96*param_esti$se_beta_S
      param_esti$couverture_beta_S = ifelse((param_init["beta_S1"] >=param_esti$bi_se_theta) & (param_init["beta_S1"]<=param_esti$bs_se_theta),1,0)

      for(h in 2:ves){
        param_esti$bi_se_theta= as.numeric(as.matrix(param_esti[paste("beta_S_", h, sep = "")]-1.96*param_esti[paste("se_beta_S_", h, sep = "")]))
        param_esti$bs_se_theta= as.numeric(as.matrix(param_esti[paste("beta_S_", h, sep = "")]+1.96*param_esti[paste("se_beta_S_", h, sep = "")]))
        param_esti$couverture_beta_S_e = as.numeric(ifelse((param_init[paste("beta_S", h, sep = "")]>=param_esti$bi_se_theta) & (param_init[paste("beta_S", h, sep = "")]<=param_esti$bs_se_theta),1,0))
        names(param_esti)[ncol(param_esti)] <- paste("couverture_beta_S_", h, sep = "")
      }
    }else{
      param_esti$bi_se_theta=param_esti$beta_S-1.96*param_esti$se_beta_S
      param_esti$bs_se_theta=param_esti$beta_S+1.96*param_esti$se_beta_S
      param_esti$couverture_beta_S=ifelse((param_init["beta_S"]>=param_esti$bi_se_theta) & (param_init["beta_S"]<=param_esti$bs_se_theta),1,0)
    }
    
    # taux de couverture beta_T
    if(vet>1){
      param_esti$bi_se_theta=param_esti$beta_T-1.96*param_esti$se_beta_T
      param_esti$bs_se_theta=param_esti$beta_T+1.96*param_esti$se_beta_T
      param_esti$couverture_beta_T = ifelse((param_init["beta_T1"] >=param_esti$bi_se_theta) & (param_init["beta_T1"]<=param_esti$bs_se_theta),1,0)
      for(h in 2:vet){
        param_esti$bi_se_theta=as.numeric(as.matrix(param_esti[paste("beta_T_", h, sep = "")]-1.96*param_esti[paste("se_beta_T_", h, sep = "")]))
        param_esti$bs_se_theta=as.numeric(as.matrix(param_esti[paste("beta_T_", h, sep = "")]+1.96*param_esti[paste("se_beta_T_", h, sep = "")]))
        param_esti$temp = as.numeric(ifelse((param_init[paste("beta_T", h, sep = "")]>=param_esti$bi_se_theta) & (param_init[paste("beta_T", h, sep = "")]<=param_esti$bs_se_theta),1,0))
        names(param_esti)[ncol(param_esti)] <- paste("couverture_beta_T_", h, sep = "")
      }
    }else{
      param_esti$bi_se_theta=param_esti$beta_T-1.96*param_esti$se_beta_T
      param_esti$bs_se_theta=param_esti$beta_T+1.96*param_esti$se_beta_T
      param_esti$couverture_beta_T=ifelse((param_init["beta_T"]>=param_esti$bi_se_theta) & (param_init["beta_T"]<=param_esti$bs_se_theta),1,0)
     }


    # taux de couverture R2trial
    #param_esti$se_R2trial= param_esti$se_R2trial*sqrt(2)
    param_esti$bi_se_theta=param_esti$R2trial-1.96*param_esti$se_R2trial
    param_esti$bs_se_theta=param_esti$R2trial+1.96*param_esti$se_R2trial
    param_esti$couverture_R2trial=ifelse((param_init["R2trial"]>=param_esti$bi_se_theta) & (param_init["R2trial"]<=param_esti$bs_se_theta),1,0)
    
    if(type.joint == 3){
      # taux de couverture KTau: "tau_00","SE.KendTau"
      #param_esti$se_R2trial= param_esti$se_R2trial*sqrt(2)
      param_esti$bi_se_theta=param_esti$tau_00-1.96*param_esti$SE.KendTau
      param_esti$bs_se_theta=param_esti$tau_00+1.96*param_esti$SE.KendTau
      param_esti$couverture_KTau=ifelse((param_init["tau"]>=param_esti$bi_se_theta) & (param_init["tau"]<=param_esti$bs_se_theta),1,0)
      
      # reorganisation des donnees
      entete <- c("Theta_S","se_theta_S","couverture_thetaS","gamma_S","se_gamma_S","couverture_gamma_S",
                  "alpha","se_alpha","couverture_alpha","sigma_S","se_sigma_S","couverture_sigma_S","sigma_T","se_sigma_T","couverture_sigma_T",
                  "sigma_ST","se_sigma_ST","couverture_sigma_ST","beta_S","se_beta_S","couverture_beta_S")
      if(ves>1){
        for(h in 2:ves){
          entete <- c(entete, paste("beta_S_", h, sep = ""))
          entete <- c(entete, paste("se_beta_S_", h, sep = ""))
          entete <- c(entete, paste("couverture_beta_S_", h, sep = ""))
        }
        
      }
      entete = c(entete,"beta_T", "se_beta_T","couverture_beta_T")
      if(vet>1){
        for(h in 2:vet){
          entete <- c(entete, paste("beta_T_", h, sep = ""))
          entete <- c(entete, paste("se_beta_T_", h, sep = ""))
          entete <- c(entete, paste("couverture_beta_T_", h, sep = ""))
        }
      }
      
      entete = c(entete,"R2trial","se_R2trial","couverture_R2trial","tau_00","SE.KendTau","couverture_KTau")
      param_esti2=param_esti[,entete]
      
      param_init=c("theta_S"=theta_S,"gamma_S"=gamma_S,"alpha"=alpha,"sigma_S"=sigma_S,"sigma_T"=sigma_T,
                   "sigma_ST"=sigma_ST,"beta_S" = beta_S[1])
      if(ves>1){
        for(h in 2:ves){
          param_init <- c(param_init, "aaa" = beta_S[h])
          names(param_init)[length(param_init)] = paste("beta_S_", h, sep = "")
        }
      }
      param_init=c(param_init,"beta_T"=beta_T[1])
      
      if(vet>1){
        for(h in 2:vet){
          param_init <- c(param_init, "aaa" = beta_T[h])
          names(param_init)[length(param_init)] <- paste("beta_T_", h, sep = "")
        }
      }
      
      param_init <- c(param_init,"R2trial"=R2trial,"tau"=tau)
      
      entete <- c("theta.latex","se.theta.S","couverture.thetaS","gamma.latex","se.gamma.S","couverture.gamma.S",
              "alpha.latex","se.alpha","couverture.alpha","sigma.S.latex","se.sigma.S","couverture.sigma.S","sigma.T.latex","se.sigma.T","couverture.sigma.T",
              "sigma.ST.latex","se.sigma.ST","couverture.sigma.ST","beta.S.latex","se.beta.S.latex","couverture.beta.S")
      
      if(ves>1){
        for(h in 2:ves){
          entete <- c(entete, paste("beta.S.latex.", h, "=", beta_S[h], sep = "")
                          , paste("se.beta.S.", h, "=", beta_S[h], sep = "")
                          , paste("couverture.beta.S.", h, "=", beta_S[h], sep = ""))
        }
      }
      
      entete <- c(entete, "beta.T.latex", "se.beta.T","couverture.beta.T")
      
      if(vet>1){
        for(h in 2:vet){
          entete <- c(entete, paste("beta.T.latex.", h, "=", beta_S[h], sep = "")
                          , paste("se.beta.T.", h, "=", beta_S[h], sep = "")
                          , paste("couverture.beta.T.", h, "=", beta_S[h], sep = ""))
        }
      }
      
      entete <- c(entete, "R2trial.latex","se.R2trial","couverture.R2trial","K.tau.latex","se.KTau","couverture.KTau")
      
      names(param_esti2)= entete
      
    }else{
      # reorganisation des donnees
      param_esti2=param_esti[,c("Theta_S","se_theta_S","couverture_thetaS","zeta","se_zeta","couverture_zeta","gamma_S","se_gamma_S","couverture_gamma_S",
                                "alpha","se_alpha","couverture_alpha","sigma_S","se_sigma_S","couverture_sigma_S","sigma_T","se_sigma_T","couverture_sigma_T",
                                "sigma_ST","se_sigma_ST","couverture_sigma_ST","beta_S","se_beta_S","couverture_beta_S","beta_T",
                                "se_beta_T","couverture_beta_T","R2trial","se_R2trial","couverture_R2trial","tau_00")]
      
      names(param_esti2)=c("theta.latex","se.theta.S","couverture.thetaS","zeta.latex","se.zeta","couverture.zeta","gamma.latex","se.gamma.S","couverture.gamma.S",
                         "alpha.latex","se.alpha","couverture.alpha","sigma.S.latex","se.sigma.S","couverture.sigma.S","sigma.T.latex","se.sigma.T","couverture.sigma.T",
                         "sigma.ST.latex","se.sigma.ST","couverture.sigma.ST","beta.S.latex","se.beta.S.latex","couverture.beta.S","beta.T.latex",
                         "se.beta.T","couverture.beta.T","R2trial.latex","se.R2trial","couverture.R2trial","K.tau.latex")
    }
    
    
    # impression des resultats
    d=data.frame(matrix(NA,ncol=6,nrow=1))
    names(d)=c("Parameters","True value","Mean","Mean SE","Empirical SE","CP(%)")
    j=1
    i=1
    while(i<ncol(param_esti2)){
      #if(j==1) cat("Parameters","True value","Mean","Mean SE","Empirical Se","CP(%)",fill=T,append=T,sep=";")
      #cat(names(param_esti2)[i],round(param_init[j],2),round(mean(param_esti2[,i]),ndec),round(mean(param_esti2[,i+1]),ndec),round(sd(param_esti2[,i]),ndec),
      #    round(prop.table(table(param_esti2[,i+2])),2)["1"],fill=T,append=T,sep=";")
      d=rbind(d,(c(names(param_esti2)[i],round(as.numeric(param_init[j]),ndec),round(mean(param_esti2[,i]),ndec),round(mean(param_esti2[,i+1]),ndec),round(sd(param_esti2[,i]),ndec),
                   100*round(prop.table(table(param_esti2[,i+2])),2)["1"])))
      i=i+3
      j=j+1
    }
    
    # impression du rho
    if(type.joint==1){
      i=ncol(param_esti2)
      d=rbind(d,c(names(param_esti2)[i],round(param_init[j],ndec),round(mean(param_esti2[,i]),ndec),"-",round(sd(param_esti2[,i]),ndec),"-"))
    }
    # ajout du poucentage de rejet
    d=rbind(d,c("R : n(%)","-",paste(nsim-nrow(param_esti2),"(",round(100*(nsim-nrow(param_esti2))/nsim),")",sep=""),"-","-","-"))
   
    # Permutation de la position SE et SD
    d=d[,c("Parameters","True value","Mean","Empirical SE","Mean SE","CP(%)")]
    # ajout de l'IC du tau de kendall
    if(ick==1){
      if(nb_paquet==1){
        #b=kendall_bootstrap(fichierTau="Taux_kendall_bootst.txt",nsim=nrow(estimates_complet2),nboot=n_bootstrap)
        kendall=ktauboot
      }
      if(nb_paquet>1){ # on entrera jamais ici
        # kendall=fusion_resultats_Tkandal(rep_courant,nb_paquet)
      }
      couv_kendall = 0
      couv_R2boot = 0
      for(i in 1:nrow(kendall)){
        #j=n_bootstrap*(i-1)+1
        #dkendall_i=kendall[c(j:(j+n_bootstrap-1)),1]
        #IC=quantile(dkendall_i, probs = c(0.025,.975))
        #if(tau>=IC[1] & tau <= IC[2]) couv_kendall=couv_kendall+1
        if(tau>=kendall[i,2] & tau <= kendall[i,3]) couv_kendall=couv_kendall+1
        if(R2trial>=R2boot[i,2] & R2trial <= R2boot[i,3]) couv_R2boot=couv_R2boot+1
        
        #cat(c("suis la",i),fill=T)
      }
        
      if(type.joint == 3){ # je contrpole le calcul par bootstrat du tau de couverture du KTau
        if(R2parboot == 1) {
          d[nrow(d)-1,ncol(d)]=100*round(couv_kendall/nrow(kendall),2)
          d[nrow(d)-1,ncol(d)-1]=NA # Evidemment 
        }
      }else{
        d[nrow(d)-1,ncol(d)]=100*round(couv_kendall/nrow(kendall),2)
      }
      #cat(100*round(couv_kendall/nrow(kendall),2))
      if(R2parboot == 1) { # on presente les taux de couverture de R2 obtenus par bootstrap
        d[nrow(d)-2,ncol(d)]=100*round(couv_R2boot/nrow(R2boot),2)
        d[nrow(d)-2,ncol(d)-1]=NA # Evidemment 
      }
    }
    return(d[-1,])
}
