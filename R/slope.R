#' Identify variable associated with the random slope
#' 
#' This is a special function used in the context of survival additive models.
#' It identifies the variable which is in interaction with the random slope
#' (\bold{\eqn{v_i}}). Generally, this variable is the treatment variable.
#' Using \code{interaction()} in a formula implies that an additive frailty
#' model is fitted.
#' 
#' 
#' @usage slope(x)
#' @param x A factor, a character or a numerical variable
#' @return \item{x}{The variable in interaction with the random slope}
#' @note It is necessary to specify which variable is in interaction with the
#' random slope, even if only one explanatory variable is included in the
#' model.
#' @seealso \code{\link{additivePenal}}
#' @keywords misc
#' @export
#' @examples
#' 
#' 
#' \dontrun{
#' 
#' data(dataAdditive)
#' 
#' ##-- Additive with one covariate --##
#' 
#' modAdd1cov <- additivePenal(Surv(t1,t2,event)~cluster(group)+var1+
#' slope(var1),data=dataAdditive,n.knots=8,kappa=10000,hazard="Splines")
#' 
#' ##-- Additive with two covariates --##
#' 
#' set.seed(1234)
#' dataAdditive$var2 <- rbinom(nrow(dataAdditive),1,0.5)
#' 
#' modAdd2cov <- additivePenal(Surv(t1,t2,event)~cluster(group)+var1+
#' var2+slope(var1),data=dataAdditive,n.knots=8,kappa=10000,
#' hazard="Splines")
#' 
#' ##-- Additive with 2 covariates and stratification --##
#' 
#' dataAdditive$var2 <- rbinom(nrow(dataAdditive),1,0.5)
#' 
#' modAddstrat <- additivePenal(Surv(t1,t2,event)~cluster(group)+
#' strata(var2)+var1+slope(var1),data=dataAdditive,n.knots=8,
#' kappa=c(10000,10000),hazard="Splines")
#' 
#' }
#' 
#' 
"slope"<-function(x)
 {
  x
 }
