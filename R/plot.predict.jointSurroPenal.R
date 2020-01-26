plot.predict.jointSurroPenal= function(object, from = -2, to = 2, type = "Coef", var.used = "error.estim"){
  # type  = "coef" or "HR"
  beta  <- object$beta.t
  dab   <- object$Coefficients$Estimate[nrow(object$Coefficients)-4]
  daa   <- object$Coefficients$Estimate[nrow(object$Coefficients)-6]
  dbb   <- object$Coefficients$Estimate[nrow(object$Coefficients)-5]
  alpha <- object$beta.s
  # alpha0 <- matrixPred$beta.S[i]
  x     <- t(matrix(c(1, -dab/daa),1,2))
  
  Vmu   <- matrix(c(object$varcov.Sigma[3,3], object$varcov.Sigma[3,1], 
                    object$varcov.Sigma[3,1], object$varcov.Sigma[1,1]),2,2)
  nparam <- nrow(object$varH)
  
  VD     <- matrix(c(object$varH[nparam,nparam], object$varH[nparam,nparam-1],
                     object$varH[nparam-1,nparam], object$varH[nparam -1,nparam - 1]),2,2)
  R2trial <- object$Coefficients$Estimate[nrow(object$Coefficients)-1] 
  # matrixPred$beta.T.i[i] <- beta + (dab/daa) * (alpha0 - alpha)
  if(type == "Coef"){ # log HR
    curve (expr = beta + (dab/daa) * (x - alpha), from = from, to = to)
  }
  else{ # HR
    curve (expr = exp(beta + (dab/daa) * (x - alpha)), from = from, to = to)
  }
  variance.inf <- dbb * (1 - R2trial) 
  variance.N <- t(x) %*% (Vmu + (((alpha0 - alpha)/daa)**2) * VD) %*% x
  + variance.inf
  
  if(var.used == "error.estim") 
    variance <- variance.N
  else 
    variance <- variance.inf
  matrixPred$Inf.95.CI[i] <- matrixPred$beta.T.i[i] - qnorm(1-alpha./2) * sqrt(variance)
  matrixPred$Sup.95.CI[i] <- matrixPred$beta.T.i[i] + qnorm(1-alpha./2) * sqrt(variance)
}