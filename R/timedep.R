#' Identify time-varying effects
#' 
#' @description{
#' This is a special function used in the context of Cox models and shared and
#' joint frailty models. It identifies time-varying effects of covariates in
#' the model. It is used in 'frailtyPenal' on the right hand side of formula or
#' of formula.terminalEvent.
#' 
#' When considering time-varying effects in a survival model, regression
#' coefficients can be modeled with a linear combination of B-splines
#' \eqn{B(t)} with coefficients \eqn{\zeta} of order \eqn{q} with \eqn{m}
#' interior knots :
#' 
#' \deqn{\beta(t)=\sum_{j=-q+1}^m\zeta_jB_{j,q}(t)}
#' 
#' You can notice that a linear combination of B-splines of order 1 without any
#' interior knots (0 interior knot) is the same as a model without time-varying
#' effect (or with constant effect over time).
#' 
#' Statistical tests (likelihood ratio tests) can be done in order to know
#' whether the time-dependent coefficients are significantly different from
#' zero or to test whether a covariate has a time-dependent effect
#' significantly different from zero or not. These tests are correct only with
#' a parametric approach yet.
#' 
#' - Proportional Hazard assumption ?
#' 
#' Time-dependency of a covariate effect can be tested. We need to estimate
#' \eqn{m+q} parameters \eqn{\zeta_j} for \eqn{j=-q+1,...,m} for a time-varying
#' coefficient. Only one (\eqn{q=1},\eqn{m=0}) parameter is estimated for a
#' constant effect. A global test is done.
#' 
#' \deqn{H_0:\beta (t)=\beta}
#' 
#' The corresponding LR statistic has a \eqn{\chi^2} distribution of degree
#' \eqn{m+q-1}.
#' 
#' - Significant association ?
#' 
#' We can also use a LR test to test whether a covariate has a significant
#' effect on the hazard function. The null hypothesis is :
#' 
#' \deqn{H_0:\beta (t)=0}
#' 
#' For that we fit a model considering the covariate with a regression
#' coefficent modeled using B-splines and a model without the covariate. Hence,
#' the LR statistic has a \eqn{\chi^2} distribution of degree \eqn{m+q}.
#' }
#' 
#' 
#' @usage timedep(x)
#' @param x A numerical or a factor variable that would have a time-varying
#' effect on the event
#' @return \item{x}{A variable identified with a time-varying effect}
#' @references Y. Mazroui, A. Mauguen, S. Mathoulin-Pelissier, G. MacGrogan, V.
#' Brouste, V. Rondeau (2013). Time-varying coefficients in a multivariate
#' frailty model: Application to breast cancer recurrences of several types and
#' death. To appear.
#' @keywords misc
#' @export
#' @examples
#' 
#' 
#' \dontrun{
#' 
#' data(readmission)
#' 
#' ###--- Shared Frailty model with time-varying effect ---###
#' 
#' sha.time <- frailtyPenal(Surv(time,event)~cluster(id)+dukes+charlson+
#' timedep(sex)+chemo,data=readmission,n.knots=8,kappa=1,
#' betaknots=3,betaorder=3)
#' 
#' #-- print results of the fit and the associated curves for the
#' #-- time-dependent effects
#' print(sha.time)
#' 
#' ###--- Joint Frailty model with time-varying effect ---###
#' 
#' joi.time <- frailtyPenal(Surv(time,event)~cluster(id)+timedep(sex)+
#' chemo+terminal(death),formula.terminalEvent=~timedep(sex)+chemo,
#' data=readmission,n.knots=8,kappa=c(1,1),betaknots=3,betaorder=3)
#' 
#' print(joi.time)
#' 
#' }
#' 
#' 
"timedep" <- function(x)
 {
	x
 }
