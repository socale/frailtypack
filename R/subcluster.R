#' Identify subclusters
#' 
#' This is a special function used in the context of survival nested or joint
#' nested models.  It identifies correlated groups of observations within other
#' groups defined by using 'cluster' function from 'survival' package, and is
#' used on the right hand side of 'frailtyPenal' formula for fitting a nested
#' or joint nested model.  Using \code{subcluster()} in a formula implies that
#' a nested or a joint nested frailty model is estimated.
#' 
#' 
#' @usage subcluster(x)
#' @param x A character, factor, or numeric variable which is supposed to
#' indicate the variable subgroup
#' @return \item{x}{A variable identified as a subcluster }
#' @seealso \code{\link{frailtyPenal}}
#' @keywords misc
#' @export
#' @examples
#' 
#' 
#' \dontrun{
#' 
#' data(dataNested)
#' modClu <- frailtyPenal(Surv(t1,t2,event)~cluster(group)+
#' subcluster(subgroup)+cov1+cov2,data=dataNested,
#' n.knots=8,kappa=c(50000,50000),hazard="Splines")
#' 
#' print(modClu)
#' 
#' #-- here is generated cluster (30 clusters)
#' readmissionNested <- transform(readmission,group=id%%30+1)
#' 
#' modJointNested_Splines <- frailtyPenal(formula = Surv(t.start, t.stop, event)
#' 	~ subcluster(id) + cluster(group) + dukes + 
#' 	terminal(death), formula.terminalEvent = ~dukes, 
#' 	data = readmissionNested, recurrentAG = TRUE, n.knots = 8, 
#' 	kappa = c(9.55e+9, 1.41e+12), initialize = TRUE)
#' 
#' }
#' 
#' 
"subcluster"<-function(x)
 {
  x
 }
