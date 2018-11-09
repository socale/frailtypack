#' General Frailty models: shared, joint and nested frailty models with
#' prediction
#' 
#' Frailtypack fits several classes of frailty models using a penalized
#' likelihood estimation on the hazard function but also a parametric
#' estimation.  1) A shared frailty model and Cox proportional hazard model.
#' Clustered and recurrent survival times can be studied.  2) Additive frailty
#' models for proportional hazard models with two correlated random effects
#' (intercept random effect with random slope).  3) Nested frailty models for
#' hierarchically clustered data (with 2 levels of clustering) by including two
#' iid gamma random effects.  4) Joint frailty models in the context of joint
#' modelling for recurrent events with terminal event for clustered data or
#' not. A joint frailty model for two semi-competing risks for clustered data
#' is also proposed.  5) Joint General frailty models in the context of a joint
#' modelling for recurrent events with terminal event data with two independent
#' frailty terms.  6) Joint Nested frailty models in the context of joint
#' modelling for recurrent events with terminal event, for hierarchically
#' clustered data (with two levels of clustering) by including two iid gamma
#' random effects.  7) Multivariate joint frailty models for two types of
#' recurrent events and a terminal event.  8) Joint models for longitudinal
#' data and a terminal event.  9) Trivariate joint models for longitudinal
#' data, recurrent events and a terminal event.  Prediction values are
#' available. Left truncated (not for the joint models), right-censored data,
#' interval-censored data (only for Cox proportional hazard and shared frailty
#' model) and strata are allowed. In each model, the random effects have the
#' gamma or normal distribution. Now, you can also consider time-varying effect
#' covariates in Cox, shared and joint frailty models. The package includes
#' concordance measures for Cox proportional hazards models and for shared
#' frailty models.
#' 
#' \tabular{ll}{ Package: \tab frailtypack\cr Type: \tab Package\cr Version:
#' \tab 2.8.3\cr Date: \tab 2016-01-12\cr License: \tab GPL (>= 2.0)\cr
#' LazyLoad: \tab no\cr }
#' 
#' @name frailtypack-package
#' @aliases frailtypack-package frailtypack
#' @docType package
#' @author Virginie Rondeau, Juan R. Gonzalez, Yassin Mazroui, Audrey Mauguen,
#' Agnieszka Krol, Amadou Diakite and Alexandre Laurent
#' @references V. Rondeau, Y. Mazroui and J. R. Gonzalez (2012). Frailtypack:
#' An R package for the analysis of correlated survival data with frailty
#' models using penalized likelihood estimation or parametric estimation.
#' \emph{Journal of Statistical Software} \bold{47}, 1-28.
#' 
#' Y. Mazroui, S. Mathoulin-Pelissier,P. Soubeyranb and Virginie Rondeau (2012)
#' General joint frailty model for recurrent event data with a dependent
#' terminalevent: Application to follicular lymphoma data. \emph{Statistics in
#' Medecine}, \bold{31}, 11-12, 1162-1176.
#' 
#' V. Rondeau and J. R. Gonzalez (2005). Frailtypack: A computer program for
#' the analysis of correlated failure time data using penalized likelihood
#' estimation. \emph{Computer Methods and Programs in Biomedicine} \bold{80},
#' 2, 154-164.
#' 
#' V. Rondeau, S. Michiels, B. Liquet, and J. P. Pignon (2008). Investigating
#' trial and treatment heterogeneity in an individual patient data
#' meta-analysis of survival data by mean of the maximum penalized likelihood
#' approach. \emph{Statistics in Medecine}, \bold{27}, 1894-1910.
#' 
#' V. Rondeau, S. Mathoulin-Pellissier, H. Jacqmin-Gadda, V. Brouste, P.
#' Soubeyran (2007). Joint frailty models for recurring events and death using
#' maximum penalized likelihood estimation:application on cancer events.
#' \emph{Biostatistics}, \bold{8}, 4, 708-721.
#' 
#' V. Rondeau, D. Commenges, and P. Joly (2003). Maximum penalized likelihood
#' estimation in a gamma-frailty model. \emph{Lifetime Data Analysis} \bold{9},
#' 139-153.
#' 
#' D. Marquardt (1963). An algorithm for least-squares estimation of nonlinear
#' parameters. \emph{SIAM Journal of Applied Mathematics}, 431-441.
#' 
#' V. Rondeau, L. Filleul, P. Joly (2006). Nested frailty models using maximum
#' penalized likelihood estimation. \emph{Statistics in Medecine}, \bold{25},
#' 4036-4052.
#' @useDynLib "frailtypack", .registration = TRUE, .fixes = "C_"
##' @import survival boot MASS survC1 nlme
##' @importFrom graphics abline legend lines matlines matplot par plot
##' @importFrom stats .getXlevels aggregate as.formula complete.cases
##' contrasts get_all_vars is.empty.model model.extract model.matrix 
##' pchisq pnorm qnorm quantile rgamma terms update var model.frame na.pass
##' @importFrom utils flush.console
##' @importFrom statmod gauss.quad
##' @importFrom nlme lme
#' @keywords package
#' @examples
#' 
#' 
#' \dontrun{
#' 
#' ###--- Additive model with 1 covariate ---###
#' 
#' data(dataAdditive)
#' modAdd <- additivePenal(Surv(t1,t2,event)~
#' cluster(group)+var1+slope(var1),
#' correlation=TRUE,data=dataAdditive,
#' n.knots=8,kappa=10000,hazard="Splines")
#' 
#' ###--- Joint model (recurrent and terminal events) with 2 covariates ---###
#' 
#' data(readmission)
#' modJoint.gap <- frailtyPenal(Surv(time,event)~
#' cluster(id)+sex+dukes+charlson+terminal(death),
#' formula.terminalEvent=~sex+dukes+charlson,
#' data=readmission,n.knots=10,kappa=c(100,100),
#' recurrentAG=FALSE,hazard="Splines")
#' 
#' ###--- General Joint model (recurrent and terminal events) with 2 covariates ---###
#' data(readmission)
#' modJoint.general <- frailtyPenal(Surv(time,event) ~ cluster(id) + dukes +
#' charlson + sex + chemo + terminal(death),
#' formula.terminalEvent = ~ dukes + charlson + sex + chemo,
#' data = readmission, jointGeneral = TRUE, n.knots = 8,
#' kappa = c(2.11e+08, 9.53e+11))
#' 
#' ###--- Nested model (or hierarchical model) with 2 covariates ---###
#' 
#' data(dataNested)
#' modClu <- frailtyPenal(Surv(t1,t2,event)~
#' cluster(group)+subcluster(subgroup)+cov1+cov2,
#' data=dataNested,n.knots=8,kappa=50000,hazard="Splines")
#' 
#' ###--- Joint Nested Frailty model ---###
#' 
#' #-- here is generated cluster (30 clusters)
#' readmissionNested <- transform(readmission,group=id%%30+1)
#' 
#' modJointNested_Splines <- frailtyPenal(formula = Surv(t.start, t.stop, event) 
#' ~ subcluster(id) + cluster(group) + dukes + terminal(death), 
#' formula.terminalEvent = ~dukes, data = readmissionNested, recurrentAG = TRUE, 
#' n.knots = 8, kappa = c(9.55e+9, 1.41e+12), initialize = TRUE)
#' 
#' modJointNested_Weib <- frailtyPenal(Surv(t.start,t.stop,event)~subcluster(id)
#' +cluster(group)+dukes+ terminal(death),formula.terminalEvent=~dukes, 
#' hazard = ('Weibull'), data=readmissionNested,recurrentAG=TRUE, initialize = FALSE)
#' 
#' JoiNes-GapSpline <- frailtyPenal(formula = Surv(time, event) 
#' ~ subcluster(id) + cluster(group) + dukes + terminal(death), 
#' formula.terminalEvent = ~dukes, data = readmissionNested, recurrentAG = FALSE, 
#' n.knots = 8, kappa = c(9.55e+9, 1.41e+12), initialize = TRUE,
#' init.Alpha = 1.091, Ksi = "None")
#' 
#' ###--- Semiparametric Shared model ---###
#' 
#' data(readmission)
#' sha.sp <- frailtyPenal(Surv(t.start,t.stop,event)~
#' sex+dukes+charlson+cluster(id),data=readmission,
#' n.knots=6,kappa=5000,recurrentAG=TRUE,
#' cross.validation=TRUE,hazard="Splines")
#' 
#' ###--- Parametric Shared model ---###
#' 
#' data(readmission)
#' sha.p <- frailtyPenal(Surv(t.start,t.stop,event)~
#' cluster(id)+sex+dukes+charlson,
#' data=readmission,recurrentAG=TRUE,
#' hazard="Piecewise-per",nb.int=6)
#' 
#' ###--- Joint model for longitudinal ---###
#' ###--- data and a terminal event ---###
#' 
#' data(colorectal)
#' data(colorectalLongi)
#' 
#' # Survival data preparation - only terminal events 
#' colorectalSurv <- subset(colorectal, new.lesions == 0)
#' 
#' model.weib.RE <- longiPenal(Surv(time1, state) ~ age + treatment + who.PS 
#' + prev.resection, tumor.size ~  year * treatment + age + who.PS ,
#' colorectalSurv,	data.Longi = colorectalLongi, 
#' random = c("1", "year"), id = "id", link = "Random-effects", 
#' left.censoring = -3.33, hazard = "Weibull")
#' 
#' ###--- Trivariate joint model for longitudinal ---###
#' ###--- data, recurrent and terminal events ---###
#' 
#' data(colorectal)
#' data(colorectalLongi)
#' 
#' # (computation takes around 40 minutes)
#' 
#' model.spli.RE.cal <-trivPenal(Surv(time0, time1, new.lesions) ~ cluster(id)
#' + age + treatment + who.PS +  terminal(state),
#' formula.terminalEvent =~ age + treatment + who.PS + prev.resection, 
#' tumor.size ~ year * treatment + age + who.PS, data = colorectal, 
#' data.Longi = colorectalLongi, random = c("1", "year"), id = "id", 
#' link = "Random-effects", left.censoring = -3.33, recurrentAG = TRUE,
#' n.knots = 6, kappa=c(0.01, 2), method.GH="Pseudo-adaptive",
#' n.nodes=7, init.B = c(-0.07, -0.13, -0.16, -0.17, 0.42, #recurrent events covariates
#' -0.23, -0.1, -0.09, -0.12, 0.8, -0.23, #terminal event covariates
#' 3.02, -0.30, 0.05, -0.63, -0.02, -0.29, 0.11, 0.74)) #biomarker covariates
#' 
#' }
#' 
#' 
NULL





#' Plot Method for an Additive frailty model.
#' 
#' Plots estimated baseline survival and hazard functions of an additive
#' frailty model, more generally of a class `additivePenal' object. Confidence
#' bands are allowed.
#' 
#' @name plot.additivePenal
#' @aliases plot.additivePenal lines.additivePenal
#' @usage
#' 
#' \method{plot}{additivePenal}(x, type.plot="Hazard", conf.bands=TRUE,
#' pos.legend="topright", cex.legend=0.7, main, color=2, Xlab = "Time", Ylab =
#' "Hazard function", ...)
#' @param x A fitted additive frailty model (output from calling
#' \code{additivePenal})
#' @param type.plot a character string specifying the type of curve. Possible
#' value are "Hazard", or "Survival". The default is "Hazard". Only the first
#' words are required, e.g "Haz", "Su"
#' @param conf.bands logical value. Determines whether confidence bands will be
#' plotted. The default is to do so.
#' @param pos.legend The location of the legend can be specified by setting
#' this argument to a single keyword from the list '"bottomright"', '"bottom"',
#' '"bottomleft"', '"left"', '"topleft"', '"top"', '"topright"', '"right"' and
#' '"center"'. The default is '"topright"'
#' @param cex.legend character expansion factor *relative* to current
#' 'par("cex")'. Default is 0.7
#' @param main plot title
#' @param color curve color (integer)
#' @param Xlab Label of x-axis. Default is '"Time"'
#' @param Ylab Label of y-axis. Default is '"Hazard function"'
#' @param \dots Other graphical parameters like those in
#' \code{\link{plot.frailtyPenal}}
#' @return Print a plot of HR and survival function of a class
#' \code{additivePenal} object
#' @seealso \code{\link{additivePenal}}
#' @keywords methods
#' @export
#' @examples
#' 
#' 
#' \dontrun{
#' 
#' data(dataAdditive)
#' 
#' modAdd <- additivePenal(Surv(t1,t2,event)~cluster(group)+var1+slope(var1),
#' correlation=TRUE,data=dataAdditive,n.knots=8,kappa=862,hazard="Splines")
#' 
#' #-- 'var1' is boolean as a treatment variable
#' 
#' plot(modAdd)
#' 
#' }
#' 
#' 
NULL





#' Print a short summary of results of prediction function.
#' 
#' Print a short summary of results of prediction function.
#' 
#' @name print.prediction
#' @aliases print.predFrailty print.predJoint print.predLongi
#' @usage \method{print}{predFrailty}(x, digits = 3, ...)
#' \method{print}{predJoint}(x, digits = 3, ...) \method{print}{predLongi}(x,
#' digits = 3, ...)
#' @param x An object from the 'prediction' function, objects inheriting from
#' \code{predFrailty}, \code{predJoint} and \code{predLongi} classes.
#' @param digits Number of digits to print
#' @param \dots Other unused arguments
#' @return
#' 
#' Print the probabilities estimated.
#' @seealso \code{\link{prediction}}
#' @keywords methods
NULL
