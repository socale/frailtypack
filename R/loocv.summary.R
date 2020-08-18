# function permettant de presenter les sorties des modeles issus du loocv. pour chaque essai,
# le modele est estime sur l'ensemble des essais exclut ce dernier

loocv.summary <- function(loocv.object, trialused, nb.parameters, names.parameters){
  n <- length(trialused)
  #n <- length(loocv.object$different.models)
  result <- data.frame(matrix(0, n, nb.parameters))
  names(result) <- names.parameters
  coef <- sapply(1:n, function(i){
    loocv.object$different.models[[i]]$Coefficient[,1]
  }
  )

  for(i in 1:n) 
    if (!is.null(coef[[i]])) 
      result[i,] <- as.numeric(coef[,i])
  result$trial <- trialused
  return(result[sapply(1:n, function(i) sum(result[i,-ncol(result)])) != 0, ])
}
