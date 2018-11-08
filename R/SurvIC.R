#' Create a survival object for interval censoring and possibly left truncated
#' data
#' 
#' This is a function used in case of interval-censoring as a response variable
#' in a model formula only for Cox proportional hazard or shared frailty model.
#' Sometimes, an unobserved event might occur in a time interval [L,U].
#' RecurrentAG argument gets invalid with the use of SurvIC. Note that this
#' function used a Kronecker product which can suffer from computation issue
#' when the number of subjects in each cluster is high. Time dependent
#' variables are not allowed.
#' 
#' Typical usages are \code{SurvIC(lower,upper,event)} or
#' \code{SurvIC(t0,lower,upper,event)}
#' 
#' @usage SurvIC(t0, lower, upper, event)
#' @param t0 Truncation time for left truncated data only. To be ignored
#' otherwise.
#' @param lower Starting time of the interval for interval-censored data. Time
#' of right-censoring instead.
#' @param upper Ending time of the interval for interval-censored data. For
#' right-censored data, lower and upper time must be equal (for numerical
#' reason).
#' @param event Status indicator 0=right-censored, 1=interval-censored
#' @export
#' @examples
#' 
#' 
#' \dontrun{
#' 
#' data(bcos)
#' bcos$event <- ifelse(bcos$left!=bcos$right,1,0)
#' 
#' ###---  Cox proportional hazard model with interval censoring ---###
#' 
#' cox.ic <- frailtyPenal(SurvIC(left,right,event)~treatment,
#' data=bcos,n.knots=8,kappa=10000)
#' 
#' ###---  Shared model with interval censoring ---###
#' 
#' bcos$group <- c(rep(1:20,4),1:14)
#' 
#' sha.ic <- frailtyPenal(SurvIC(left,right,event)~cluster(group)+
#' treatment,data=bcos,n.knots=8,kappa=10000)
#' 
#' }
#' 
#' 
"SurvIC" <- function(t0,lower,upper,event) {

# nombre d'arguments manquants
ng <- missing(t0) + missing(lower) + missing(upper) + missing(event)

# si un seul argument est fourni, c'est que c'est un lower
if (ng==3) lower <- t0
# si deux arguments fournis, c'est lower et upper
if (ng==2){
	upper <- lower
	lower <- t0
}
# si trois arguments fournis, lower upper event
if (ng==1){
	event <- upper
	upper <- lower
	lower <- t0
}
if (missing(lower)) stop("Must have a lower time argument")
if (!is.numeric(lower)) stop("Lower time variable is not numeric")
if (any(is.na(lower))) stop("There is some NA in your lower time")

if (missing(upper)) stop("Must have an upper time argument")
if (!is.numeric(upper)) stop("Upper time variable is not numeric")
if (any(is.na(upper))) stop("There is some NA in your upper time")

if (missing(event)) stop("Must have an event argument")
if (!is.numeric(event)) stop("Event variable is not numeric")
if (any(is.na(event))) stop("There is some NA in your event argument")

if (length(lower)!=length(upper)) stop("Lower and upper time are different lengths")
if (any(lower>upper)) stop("Lower time has to be less than upper time")

if (length(lower==upper) != length(event==0)) warning("There may be an error in the right censored times")

if ((ng==0) & (any(t0>lower))) warning("Be careful, some of your truncation times are higher than the lower")

# s'il manque le t0
if (ng==1){
	ss <- cbind(lower=lower,upper=upper,status=event)
	attr(ss,"type") <- "interval"
}else{
	ss <- cbind(t0=t0,lower=lower,upper=upper,status=event)
	attr(ss,"type") <- "intervaltronc"
}

if (is.R()) { class(ss) <- "SurvIC" }
else { oldClass(ss) <- "SurvIC" }
return(ss)
}
