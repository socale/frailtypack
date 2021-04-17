# fonction permettant de fusionner les RData issus des simulations par paquet a 
# partir d'un modele joint surrogate ou joint frailty-copula et de produire le resume
# des resultats de simulation dans un object de la classe jointSurroPenalSimul. On pourra
# donc appliquer la fonction summary() sur cet objet pour avoir le resume de la fonction
# L'utilisation de cette fonction demande de nommer les fichiers RData issus des paquets de 
# simulation de la meme maniere et de les distinquer par un numero d'ordre Exemp : joint.simul2_500,
# joint.simul2_501, joint.simul2_502, .... Par la suite, il convient de ranger ces fichiers dans
# un repertioire de travail dont le nom sera par la suite fourni dans l'argument wd de la fonction

mergeJointSurroSimul = function(nb.packet = 2, envir.name = "joint.simul2_", envir.num.base = 500,
                           wd = "G:/socale/PHD-Thesis/programmes/Creation_Package/package_CRAN/Version_github/frailtypack/EspacePaquetsSimul"
                           ){
  # nb.packet = number of packets of simulation, correspond to the number of .RData files
  # envir.name = the main part of the .RData names
  # envir.num.base = the reference number of the packets
  # wd = the work directory that contains the .RData files
  
  setwd(wd)
  filename <- paste(envir.name, envir.num.base,1, ".RData", sep = "")
  load(filename)
  joint.simul <- joint.simul2
  joint.simul2 <- NULL
  if(nb.packet>1){
    for(i in 2:(nb.packet)){
      filename <- paste(envir.name, envir.num.base, i , ".RData", sep = "")
      loadwd = try(load(filename),silent=TRUE)
      if(class(loadwd)=="try-error"){
        cat(paste("packet",i,"is missing", sep = " "), fill = T)
      }
      else{
        joint.simul$dataParamEstim <- rbind(joint.simul$dataParamEstim, joint.simul2$dataParamEstim)
        joint.simul$dataTkendall <- rbind(joint.simul$dataTkendall, joint.simul2$dataTkendall)
        joint.simul$dataR2boot <- rbind(joint.simul$dataR2boot, joint.simul2$dataR2boot)
        joint.simul$nb.simul <- joint.simul$nb.simul + joint.simul2$nb.simul
        joint.simul2 <- NULL
      }
    }
  }
  return(joint.simul)
}