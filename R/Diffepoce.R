#' Difference of Expected Prognostic Observed Cross-Entropy (EPOCE) estimators
#' and its 95\% tracking interval between two joint models.
#' 
#' 
#' This function computes the difference of two EPOCE estimates (CVPOL and
#' MPOL) and its 95\% tracking interval between two joint models estimated
#' using \code{frailtyPenal}, \code{longiPenal} or \code{trivPenal}. Difference
#' in CVPOL is computed when the EPOCE was previously estimated on the same
#' dataset as used for estimation (using an approximated cross-validation), and
#' difference in MPOL is computed when the EPOCE was previously estimated on an
#' external dataset.
#' 
#' 
#' From the EPOCE estimates and the individual contributions to the prognostic
#' observed log-likelihood obtained with \code{epoce} function on the same
#' dataset from two different estimated joint models, the difference of CVPOL
#' (or MPOL) and its 95\% tracking interval is computed. The 95\% tracking
#' interval is : Delta(MPOL) +/- qnorm(0.975)*sqrt(VARIANCE) for an external
#' dataset Delta(CVPOL) +/- qnorm(0.975)*sqrt(VARIANCE) for the dataset used in
#' \code{frailtyPenal}, \code{longiPenal} or \code{trivPenal} where
#' Delta(CVPOL) (or Delta(MPOL)) is the difference of CVPOL (or MPOL) of the
#' two joint models, and VARIANCE is the empirical variance of the difference
#' of individuals contributions to the prognostic observed log-likelihoods of
#' the two joint models.
#' 
#' The estimators of EPOCE from arguments epoce1 and epoce2 must have been
#' computed on the same dataset and with the pred.times.
#' 
#' @usage Diffepoce(epoce1, epoce2)
#' @param epoce1 a first object inheriting from class epoce.
#' @param epoce2 a second object inheriting from class epoce.
#' @return \item{new.data}{a boolean which is FALSE if computation is done on
#' the same data as for estimation, and TRUE otherwise} \item{pred.times}{time
#' or vector of times used in the function} \item{DEPOCE}{the difference
#' between the two MPOL or CVPOL for each time} \item{TIinf}{lower confidence
#' band for the difference} \item{TIsup}{upper confidence band for the
#' difference}
#' @references D. Commenges, B. Liquet, C. Proust-Lima (2012). Choice of
#' prognostic estimators in joint models by estimating differences of expected
#' conditional Kullback-Leibler risks. \emph{Biometrics}, \bold{68(2)},
#' 380-387.
#' @keywords misc
#' @export
#' @examples
#' 
#' 
#' \dontrun{
#' 
#' #Example for joint frailty models
#' data(readmission)
#' 
#' # first joint frailty model
#' joint1 <- frailtyPenal(Surv(t.start,t.stop,event)~ cluster(id) +
#'   dukes + charlson + sex + chemo + terminal(death),
#'   formula.terminalEvent = ~ dukes + charlson + sex + chemo ,
#'   data = readmission, n.knots = 8, kappa = c(2.11e+08,9.53e+11),
#'   recurrentAG=TRUE)
#' 
#' # second joint frailty model without dukes nor charlson as covariates
#' joint2 <- frailtyPenal(Surv(t.start,t.stop,event)~ cluster(id) +
#'   sex + chemo + terminal(death),
#'   formula.terminalEvent = ~ sex + chemo ,
#'   data = readmission, n.knots = 8, kappa = c(2.11e+08,9.53e+11),
#'   recurrentAG=TRUE)
#' 
#' temps <- c(200,500,800,1100)
#' 
#' # computation of estimators of EPOCE for the two models
#' epoce1 <- epoce(joint1,temps)
#' epoce2 <- epoce(joint2,temps)
#' 
#' # computation of the difference
#' diff <- Diffepoce(epoce1,epoce2)
#' 
#' print(diff)
#' plot(diff)
#' 
#' 
#' #Example for joint models with a biomarker
#' data(colorectal)
#' data(colorectalLongi)
#' 
#' # Survival data preparation - only terminal events 
#' colorectalSurv <- subset(colorectal, new.lesions == 0)
#' 
#' # first joint model for a biomarker and a terminal event
#' modLongi <- longiPenal(Surv(time0, time1, state) ~ age +
#' treatment + who.PS, tumor.size ~  year*treatment + age +
#' who.PS, colorectalSurv, data.Longi =colorectalLongi,
#' random = c("1", "year"),  id = "id", link = "Random-effects", 
#' left.censoring = -3.33, hazard = "Weibull", 
#' method.GH = "Pseudo-adaptive")
#' 
#' # second joint model for a biomarker, recurrent events and a terminal event
#' # (computation takes around 30 minutes)
#' modTriv <- model.weib.RE.gap <-trivPenal(Surv(gap.time, new.lesions) 
#' ~ cluster(id) + age + treatment + who.PS + prev.resection + terminal(state),
#' formula.terminalEvent =~ age + treatment + who.PS + prev.resection, 
#' tumor.size ~ year * treatment + age + who.PS, data = colorectal,
#' data.Longi = colorectalLongi, random = c("1", "year"), id = "id", 
#' link = "Random-effects", left.censoring = -3.33, recurrentAG = FALSE,
#' hazard = "Weibull", method.GH="Pseudo-adaptive", n.nodes=7)
#' 
#' time <- c(1, 1.5, 2, 2.5)
#' 
#' # computation of estimators of EPOCE for the two models
#' epoce1 <- epoce(modLongi, time)
#' # (computation takes around 10 minutes)
#' epoce2 <- epoce(modTriv, time)
#' 
#' 
#' # computation of the difference
#' diff <- Diffepoce(epoce1, epoce2)
#' 
#' print(diff)
#' plot(diff)
#' }
#' 
#' 
Diffepoce <- function(epoce1, epoce2){

	if (class(epoce1)!="epoce") stop("Diffepoce allows only arguments of classe 'epoce'")
	if (class(epoce2)!="epoce") stop("Diffepoce allows only arguments of classe 'epoce'")
	if (!(all.equal(epoce1$pred.times,epoce2$pred.times))) stop("The two epoce objects should have the same prediction times")
	if (!(all.equal(epoce1$data,epoce2$data))) stop("The two epoce objects should be derived from the same dataset")

	DEPOCE <- as.double(epoce1$cvpol)-as.double(epoce2$cvpol)
	if (epoce1$new.data==TRUE){
		DEPOCE <- as.double(epoce1$mpol)-as.double(epoce2$mpol)
	}
	
	diff <- epoce1$IndivContrib-epoce2$IndivContrib
	diff_sq <- (epoce1$IndivContrib-epoce2$IndivContrib)^2

	DMPOL <- apply(diff,FUN=sum,MARGIN=2)/apply(epoce1$IndivContrib!=0,FUN=sum,MARGIN=2)

	omega <- (apply(diff_sq,FUN=sum,MARGIN=2)/apply(epoce1$IndivContrib!=0,FUN=sum,MARGIN=2)) - (DMPOL)^2
	omega <- omega/apply(epoce1$IndivContrib!=0,FUN=sum,MARGIN=2)

	TIinf <- as.vector(DEPOCE) - qnorm(0.975) * sqrt(omega)
	TIsup <- as.vector(DEPOCE) + qnorm(0.975) * sqrt(omega)

	out <- NULL
	out$new.data <- epoce1$new.data
	out$pred.times <- epoce1$pred.times
	out$DEPOCE <- DEPOCE
	out$TIinf <- TIinf
	out$TIsup <- TIsup

	class(out) <- c("Diffepoce")
	out
}
