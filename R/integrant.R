# contient les integrans

# Integrant pour le calcul du taux de kendall pour le modele conjoint de validation d'un 
# critere de substitution, en utilisant R
#
# @param w f
# @param wp f
# @param u f
# @param up f
# @param theta f
# @param gamma f
# @param alpha f
# @param eta f
# @param ui f
# 
funcJointSurroKendall <- function(u, w, up, wp, theta, gamma, alpha, eta, ui = 1){
  if(ui==1){
    num <- (exp(w + u + eta * w + alpha * u) + exp(wp + up + eta * wp + alpha * up)) *
           exp(-(wp**2)/(2 * theta)) * exp(-(up**2)/(2 * gamma)) * 
           exp(-(w**2)/(2 * theta)) * exp((-u**2)/(2 * gamma))
    
    den <- (exp(wp + up) + exp(w + u)) * (exp(eta * wp + alpha * up) + exp(eta * w + alpha * u)) *
           4 * pi**2 * theta * gamma
  }else{
    num <- (exp(w + eta * w ) + exp(wp + eta * wp )) *
      exp(-(wp**2)/(2 * theta)) * exp(-(w**2)/(2 * theta)) 
    
    den <- (exp(wp) + exp(w)) * (exp(eta * wp) + exp(eta * w)) *
      2 * pi * theta
  }
  
  return (num/den)
}