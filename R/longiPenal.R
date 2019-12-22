#' Fit a Joint Model for Longitudinal Data and a Terminal Event
#'
#' @description{
#' \if{html}{Fit a joint model for longitudinal data and a terminal event
#' using a semiparametric penalized likelihood estimation or a parametric
#' estimation on the hazard function.
#'
#' The longitudinal outcomes y\out{<sub>i</sub>}(t\out{<sub>ik</sub>})
#' (k=1,\ldots,n\out{<sub>i</sub>}, i=1,\ldots,N) for N subjects are described
#' by a linear mixed model and the risk of the terminal event is represented by
#' a proportional hazard risk model. The joint model is constructed assuming
#' that the processes are linked via a latent structure (Wulfsohn and Tsiatsis
#' 1997):
#'
#' {\figure{longimodel.png}{options: width="100\%"}}
#'
#' where \bold{\eqn{X}}\out{<sub>Li</sub>}(t) and
#' \bold{\eqn{X}}\out{<sub>Ti</sub>} are vectors of fixed effects covariates and
#' \bold{\eqn{\beta}}\out{<sub>L</sub>} and \bold{\eqn{\beta}}\out{<sub>T</sub>}
#' are the associated coefficients. Measurements errors
#' \eqn{\epsilon}\out{<sub>i</sub>}(t\out{<sub>ik</sub>}) are iid normally
#' distributed with mean 0 and variance
#' \eqn{\sigma}\out{<sub>\epsilon</sub>}\out{<sup>2</sup>}. The random
#' effects \bold{b}\out{<sub>i</sub>} =
#' (b\out{<sub>0i</sub>},\ldots,b\out{<sub>qi</sub>})\out{<sup>T</sup>}
#' \out{<span>&#126;</span>} \bold{\eqn{N}}(0,\bold{B}\out{<sub>1</sub>}) are
#' associated to covariates \bold{\eqn{Z}}\out{<sub>Li</sub>}(t) and independent
#' from the measurement error. The relationship between the two processes is
#' explained via
#' h(\bold{b}\out{<sub>i</sub>},\bold{\eqn{\beta}}\out{<sub>L</sub>},\bold{\eqn{Z}}\out{<sub>Li</sub>}(t),\bold{\eqn{X}}\out{<sub>Li</sub>}(t))
#' with coefficients \eqn{\eta}\out{<sub>T</sub>}. Two forms of the function
#' h(.) are available: the random effects \bold{b}\out{<sub>i</sub>}
#' and the current biomarker level
#' b\out{<sub>i</sub>}(t)=\bold{\eqn{X}\out{<sub>Li</sub>}}(t\out{<sub>ik</sub>})\out{<sup>T</sup>}
#' \bold{\eqn{\beta}}\out{<sub>L</sub>} +\bold{ \eqn{Z}}\out{<sub>i</sub>}(t\out{<sub>ik</sub>})\out{<sup>T</sup>} \bold{b}\out{<sub>i</sub>}.
#'
#' We consider that the longitudinal outcome can be a subject to a
#' quantification limit, i.e. some observations, below a level of detection
#' s cannot be quantified (left-censoring).
#'
#' Alternatively, a two-part model is proposed to fit a semicontinuous biomarker.
#' The two-part model decomposes the biomarker's distribution into a binary
#' outcome (zero vs. positive values) and a continuous outcome (positive values).
#' A logistic mixed effects model fits the binary outcome and a linear mixed effects
#' model fits the continuous outcome.
#' }
#' \if{latex}{Fit a joint model for longitudinal data and a terminal event using a
#' semiparametric penalized likelihood estimation or a parametric estimation on
#' the hazard function.
#'
#' The longitudinal outcomes \eqn{y_i(t_{ik})} (\eqn{k=1,\ldots,n_i},
#' \eqn{i=1,\ldots,N}) for \eqn{N} subjects are described by a linear mixed
#' model and the risk of the terminal event is represented by a proportional
#' hazard risk model. The joint model is constructed assuming that the
#' processes are linked via a latent structure (Wulfsohn and Tsiatsis 1997):
#'
#' \deqn{\left\{ \begin{array}{ll} y_{i}(t_{ik})=\bold{X}_{Li}(t_{ik})^\top
#' \bold{\beta}_L +\bold{ Z}_i(t_{ik})^\top \bold{b}_i + \epsilon_i(t_{ik}) &
#' \mbox{(Longitudinal)} \\
#' \lambda_i(t|\bold{b}_i)=\lambda_0(t)\exp(\bold{X}_{Ti}(t)\bold{\beta}_T+h(\bold{b}_i,\bold{\beta}_L,\bold{Z}_i(t),\bold{X}_{Li}(t))^\top
#' \bold{\eta}_T ) & \mbox{(Terminal)} \\ \end{array} \right. }
#'
#' where \eqn{\bold{X}_{Li}(t)} and \eqn{\bold{X}_{Ti}} are vectors of fixed
#' effects covariates and \eqn{\bold{\beta}_L} and \eqn{\bold{\beta}_T} are the
#' associated coefficients. Measurements errors \eqn{\epsilon_i(t_{ik})} are
#' iid normally distributed with mean 0 and variance \eqn{\sigma_{\epsilon}^2}.
#' The random effects \eqn{\bold{b}_i = (b_{0i},\ldots, b_{qi})^\top\sim
#' \mathcal{N}(0,\bold{B}_1)} are associated to covariates \eqn{\bold{Z}_i(t)}
#' and independent from the measurement error. The relationship between the two
#' processes is explained via
#' \eqn{h(\bold{b}_i,\bold{\beta}_L,\bold{Z}_i(t),\bold{X}_{Li}(t))} with
#' coefficients \eqn{\bold{\eta}_T}. Two forms of the function \eqn{h(\cdot)}
#' are available: the random effects \eqn{\bold{b}_i} and the current biomarker
#' level \eqn{m_i(t)=\bold{X}_{Li}(t_{ik})^\top \bold{\beta}_L +\bold{
#' Z}_i(t_{ik})^\top \bold{b}_i}.
#'
#' We consider that the longitudinal outcome can be a subject to a
#' quantification limit, i.e. some observations, below a level of detection
#' \eqn{s} cannot be quantified (left-censoring).}
#' }
#'
#' @details{ Typical usage for the joint model
#' \preformatted{longiPenal(Surv(time,event)~var1+var2, biomarker ~ var1+var2,
#' data, data.Longi, ...)}
#'
#' The method of the Gauss-Hermite quadrature for approximations of the
#' multidimensional integrals, i.e. length of \code{random} is 2, can be chosen
#' among the standard, non-adaptive, pseudo-adaptive in which the quadrature
#' points are transformed using the information from the fitted mixed-effects
#' model for the biomarker (Rizopoulos 2012) or multivariate non-adaptive
#' procedure proposed by Genz et al. 1996 and implemented in FORTRAN subroutine
#' HRMSYM.  The choice of the method is important for estimations. The standard
#' non-adaptive Gauss-Hermite quadrature (\code{"Standard"}) with a specific
#' number of points gives accurate results but can be time consuming. The
#' non-adaptive procedure (\code{"HRMSYM"}) offers advantageous computational
#' time but in case of datasets in which some individuals have few repeated
#' observations (biomarker measures or recurrent events), this method may be
#' moderately unstable.  The pseudo-adaptive quadrature uses transformed
#' quadrature points to center and scale the integrand by utilizing estimates of
#' the random effects from an appropriate linear mixed-effects model. This
#' method enables using less quadrature points while preserving the estimation
#' accuracy and thus lead to a better computational time.The Monte-Carlo method
#' is also proposed for approximations of the multidimensional integrals.
#'
#'
#' NOTE. Data frames \code{data} and \code{data.Longi} must be consistent. Names
#' and types of corresponding covariates must be the same, as well as the number
#' and identification of individuals. }
#'
#' @usage longiPenal(formula, formula.LongitudinalData, data, data.Longi,
#'   formula.Binary=FALSE, random, random.Binary=FALSE, id, intercept = TRUE,
#'   link = "Random-effects", timevar=FALSE, left.censoring =
#'   FALSE, n.knots, kappa, maxit = 350, hazard = "Splines", init.B,
#'   init.Random, init.Eta, method.GH = "Standard", seed.MC=FALSE, n.nodes, LIMparam = 1e-3,
#'   LIMlogl = 1e-3, LIMderiv = 1e-3, print.times = TRUE)
#' @param formula a formula object, with the response on the left of a
#'   \eqn{\sim} operator, and the terms on the right. The response must be a
#'   survival object as returned by the 'Surv' function like in survival
#'   package. Interactions are possible using * or :.
#' @param formula.LongitudinalData a formula object, only requires terms on the
#'   right to indicate which variables are modelling the longitudinal outcome.
#'   It must follow the standard form used for linear mixed-effects models.
#'   Interactions are possible using * or :.
#' @param formula.Binary a formula object, only requires terms on the
#'   right to indicate which variables are modelling the binary part of the
#'   two-part model fitting the longitudinal semicontinuous outcome.
#'   It must follow the standard form used for linear mixed-effects models.
#'   Interactions are possible using * or :.
#' @param data a 'data.frame' with the variables used in \code{formula}.
#' @param data.Longi a 'data.frame' with the variables used in
#'   \code{formula.LongitudinalData}.
#' @param random Names of variables for the random effects of the longitudinal
#'   outcome. Maximum 3 random effects are possible at the moment. The random
#'   intercept is chosen using \code{"1"}.
#' @param random.Binary Names of variables for the random effects of the binary
#' part of the two-part model fitting the longitudinal semicontinuous outcome.
#' The random intercept is chosen using \code{"1"}.
#' @param id Name of the variable representing the individuals.
#' @param intercept Logical value. Is the fixed intercept of the biomarker
#'   included in the mixed-effects model? The default is \code{TRUE}.
#' @param link Type of link function for the dependence between the biomarker
#'   and death: \code{"Random-effects"} for the association directly via the
#'   random effects of the biomarker, \code{"Current-level"} for the association
#'   via the true current level of the biomarker.  The option
#'   \code{"Current-level"} can be chosen only if the biomarker random effects
#'   are associated with the intercept and time (following this order).
#'   \code{"Two-part"}, this structure is only applicable with two-part models,
#'   the effect of the current probability of positive value and the effect of
#'   the expected value among positive values on the risk of event is evaluated
#'   separately. The default is \code{"Random-effects"}.
#' @param timevar Indicates the time varying variables to take into account this
#' evolution over time in the link with the survival model (useful with
#' 'Current-level' and 'Two-part' links)
#' @param left.censoring Is the biomarker left-censored below a threshold
#'   \eqn{s}? The default is \code{FALSE}, ie. no left-censoring. In case of a
#'   left-censored biomarker, this argument must be equal to the threshold
#'   \eqn{s}.
#' @param n.knots Integer giving the number of knots to use. Value required in
#'   the penalized likelihood estimation.  It corresponds to the (n.knots+2)
#'   splines functions for the approximation of the hazard or the survival
#'   functions.  We estimate I or M-splines of order 4. When the user set a
#'   number of knots equals to k (n.knots=k) then the number of interior knots
#'   is (k-2) and the number of splines is (k-2)+order.  Number of knots must be
#'   between 4 and 20. (See Note in \code{frailtyPenal} function)
#' @param kappa Positive smoothing parameter in the penalized likelihood
#'   estimation.  The coefficient kappa of the integral of the squared second
#'   derivative of hazard function in the fit (penalized log likelihood). To
#'   obtain an initial value for \code{kappa}, a solution is to fit the
#'   corresponding Cox model using cross validation (See \code{cross.validation}
#'   in function \code{frailtyPenal}).  We advise the user to identify several
#'   possible tuning parameters, note their defaults and look at the sensitivity
#'   of the results to varying them.
#' @param maxit Maximum number of iterations for the Marquardt algorithm. The
#'   default is 350.
#' @param hazard Type of hazard functions: \code{"Splines"} for semiparametric
#'   hazard functions using equidistant intervals or \code{"Splines-per"} using
#'   percentile with the penalized likelihood estimation, \code{"Weibull"} for
#'   the parametric Weibull functions. The default is \code{"Splines"}.
#' @param init.B Vector of initial values for regression coefficients. This
#'   vector should be of the same size as the whole vector of covariates with
#'   the first elements for the covariates related to the terminal event and
#'   then for the covariates related to the biomarker (interactions in the end
#'   of each component). Default is 0.5 for each.
#' @param init.Random Initial value for variance of the elements of the matrix
#'   of the distribution of the random effects. Default is 0.5 for each element.
#' @param init.Eta Initial values for regression coefficients for the link
#'   function. Default is 0.5 for each.
#' @param method.GH Method for the Gauss-Hermite quadrature: \code{"Standard"}
#'   for the standard non-adaptive Gaussian quadrature, \code{"Pseudo-adaptive"}
#'   for the pseudo-adaptive Gaussian quadrature, \code{"Monte-carlo"} for the
#'   Monte-carlo method and \code{"HRMSYM"} for the
#'   algorithm for the multivariate non-adaptive Gaussian quadrature (see
#'   Details). The default is \code{"Standard"}.
#' @param n.nodes Number of nodes for the Gauss-Hermite quadrature or the
#' Monte-carlo method. They can be chosen among 5, 7, 9, 12, 15, 20 and 32
#' for the GH quadrature and any number for the Monte-carlo method. The default is 9.
#' @param seed.MC Monte-carlo integration points selection (1=fixed, 0=random)
#' @param LIMparam Convergence threshold of the Marquardt algorithm for the
#'   parameters (see Details of \code{frailtyPenal} function), \eqn{10^{-3}} by
#'   default.
#' @param LIMlogl Convergence threshold of the Marquardt algorithm for the
#'   log-likelihood (see Details of \code{frailtyPenal} function), \eqn{10^{-3}}
#'   by default.
#' @param LIMderiv Convergence threshold of the Marquardt algorithm for the
#'   gradient (see Details of \code{frailtyPenal} function), \eqn{10^{-3}} by
#'   default.
#' @param print.times a logical parameter to print iteration process. The
#'   default is TRUE.
#' @return
#'
#' The following components are included in a 'longiPenal' object for each
#' model:
#'
#' \item{b}{The sequence of the corresponding estimation of the coefficients for
#' the hazard functions (parametric or semiparametric), the random effects
#' variances and the regression coefficients.} \item{call}{The code used for the
#' model.} \item{formula}{The formula part of the code used for the terminal
#' event part of the model.} \item{formula.LongitudinalData}{The formula part of
#' the code used for the longitudinal part of the model.} \item{formula.Binary}{The
#' formula part of the code used for the binary part of the two-part model.}
#' \item{coef}{The
#' regression coefficients (first for the terminal event and then for the
#' biomarker.} \item{groups}{The number of groups used in the fit.}
#' \item{kappa}{The value of the smoothing parameter in the penalized likelihood
#' estimation corresponding to the baseline hazard function for the terminal
#' event.} \item{logLikPenal}{The complete marginal penalized log-likelihood in
#' the semiparametric case.} \item{logLik}{The marginal log-likelihood in the
#' parametric case.} \item{n.measurements}{The number of biomarker observations
#' used in the fit.} \item{max_rep}{The maximal number of repeated measurements
#' per individual.} \item{n.deaths}{The number of events observed in the fit.}
#' \item{n.iter}{The number of iterations needed to converge.}
#' \item{n.knots}{The number of knots for estimating the baseline hazard
#' function in the penalized likelihood estimation.} \item{n.strat}{The number
#' of stratum.}
#'
#' \item{varH}{The variance matrix of all parameters (before positivity
#' constraint transformation for the variance of the measurement error, for
#' which the delta method is used).} \item{varHIH}{The robust estimation of the
#' variance matrix of all parameters.}
#'
#' \item{xD}{The vector of times where both survival and hazard function of the
#' terminal event are estimated. By default seq(0,max(time),length=99), where
#' time is the vector of survival times.} \item{lamD}{The array (dim=3) of
#' baseline hazard estimates and confidence bands (terminal event).}
#' \item{survD}{The array (dim=3) of baseline survival estimates and confidence
#' bands (terminal event).}
#' \item{median}{The value of the median survival and its confidence bands.}
#' \item{typeof}{The type of the baseline hazard
#' functions (0:"Splines", "2:Weibull").} \item{npar}{The number of parameters.}
#' \item{nvar}{The vector of number of explanatory variables for the terminal
#' event and biomarker.} \item{nvarEnd}{The number of explanatory variables for
#' the terminal event.} \item{nvarY}{The number of explanatory variables for the
#' biomarker.} \item{noVarEnd}{The indicator of absence of the explanatory
#' variables for the terminal event.} \item{noVarY}{The indicator of absence of
#' the explanatory variables for the biomarker.} \item{LCV}{The approximated
#' likelihood cross-validation criterion in the semiparametric case (with H
#' minus the converged Hessian matrix, and l(.) the full
#' log-likelihood).\deqn{LCV=\frac{1}{n}(trace(H^{-1}_{pl}H) - l(.))}}
#' \item{AIC}{The Akaike information Criterion for the parametric
#' case.\deqn{AIC=\frac{1}{n}(np - l(.))}} \item{n.knots.temp}{The initial value
#' for the number of knots.} \item{shape.weib}{The shape parameter for the
#' Weibull hazard function.} \item{scale.weib}{The scale parameter for the
#' Weibull hazard function.} \item{martingaledeath.res}{The martingale residuals
#' for each individual.} \item{conditional.res}{The conditional residuals for
#' the biomarker (subject-specific):
#' \eqn{\bold{R}_i^{(m)}=\bold{y}_i-\bold{X}_{Li}^\top\widehat{\bold{\beta}}_L-\bold{Z}_i^\top\widehat{\bold{b}}_i}.}
#' \item{marginal.res}{The marginal residuals for the biomarker (population
#' averaged):
#' \eqn{\bold{R}_i^{(c)}=\bold{y}_i-\bold{X}_{Li}^\top\widehat{\bold{\beta}}_L}.}
#' \item{marginal_chol.res}{The Cholesky marginal residuals for the biomarker:
#' \eqn{\bold{R}_i^{(m)}=\widehat{\bold{U}_i^{(m)}}\bold{R}_i^{(m)}}, where
#' \eqn{\widehat{\bold{U}_i^{(m)}}} is an upper-triangular matrix obtained by
#' the Cholesky decomposition of the variance matrix
#' \eqn{\bold{V}_{\bold{R}_i^{(m)}}=\widehat{\bold{V}_i}-\bold{X}_{Li}(\sum_{i=1}^N\bold{X}_{Li}\widehat{\bold{V}_i}^{-1}\bold{X}_{Li})^{-1}\bold{X}_{Li}^\top}.}
#' \item{conditional_st.res}{The standardized conditional residuals for the
#' biomarker.} \item{marginal_st.res}{The standardized marginal residuals for
#' the biomarker.} \item{random.effects.pred}{ The empirical Bayes predictions
#' of the random effects (ie. using conditional posterior distributions).}
#' \item{pred.y.marg}{The marginal predictions of the longitudinal outcome.}
#' \item{pred.y.cond}{The conditional (given the random effects) predictions of
#' the longitudinal outcome.} \item{lineardeath.pred}{The linear predictor for
#' the terminal part.} \item{global_chisq_d}{The vector with values of each
#' multivariate Wald test for the terminal part.} \item{dof_chisq_d}{The vector
#' with degrees of freedom for each multivariate Wald test for the terminal
#' part.} \item{global_chisq.test_d}{The binary variable equals to 0 when no
#' multivariate Wald is given, 1 otherwise (for the terminal part).}
#' \item{p.global_chisq_d}{The vector with the p_values for each global
#' multivariate Wald test for the terminal part.} \item{global_chisq}{The vector
#' with values of each multivariate Wald test for the longitudinal part.}
#' \item{dof_chisq}{The vector with degrees of freedom for each multivariate
#' Wald test for the longitudinal part.} \item{global_chisq.test}{The binary
#' variable equals to 0 when no multivariate Wald is given, 1 otherwise (for the
#' longitudinal part).} \item{p.global_chisq}{The vector with the p_values for
#' each global multivariate Wald test for the longitudinal part.}
#' \item{names.factordc}{The names of the "as.factor" variables for the terminal
#' part.} \item{names.factor}{The names of the "as.factor" variables for the
#' longitudinal part.}
#'
#' \item{intercept}{The logical value. Is the fixed intercept included in the
#' linear mixed-effects model?} \item{B1}{The variance matrix of the random
#' effects for the longitudinal outcome.} \item{ResidualSE}{The standard
#' deviation of the measurement error.} \item{eta}{The regression coefficients
#' for the link function.} \item{ne_re}{The number of random effects used in the
#' fit.} \item{names.re}{The names of variables for the random effects.}
#' \item{link}{The name of the type of the link function.}
#'
#' \item{eta_p.value}{p-values of the Wald test for the estimated regression
#' coefficients for the link function.} \item{beta_p.value}{p-values of the Wald
#' test for the estimated regression coefficients.}
#'
#' \item{leftCensoring}{The logical value. Is the longitudinal outcome
#' left-censored?} \item{leftCensoring.threshold}{For the left-censored
#' biomarker, the value of the left-censoring threshold used for the fit.}
#' \item{prop.censored}{The fraction of observations subjected to the
#' left-censoring.}
#'
#' \item{methodGH}{The method used for approximations of the
#' multidimensional integrals.}
#' \item{n.nodes}{The number of integration points.}
#' @seealso
#' \code{\link{plot.longiPenal}},\code{\link{print.longiPenal}},\code{\link{summary.longiPenal}}
#' @references
#'
#' A. Krol, A. Mauguen, Y. Mazroui, A. Laurent, S. Michiels and V. Rondeau
#' (2017). Tutorial in Joint Modeling and Prediction: A Statistical Software for
#' Correlated Longitudinal Outcomes, Recurrent Events and a Terminal Event.
#' \emph{Journal of Statistical Software} \bold{81}(3), 1-52.
#'
#' A. Krol, L. Ferrer, JP. Pignon, C. Proust-Lima, M. Ducreux, O. Bouche, S.
#' Michiels, V. Rondeau (2016). Joint Model for Left-Censored Longitudinal Data,
#' Recurrent Events and Terminal Event: Predictive Abilities of Tumor Burden for
#' Cancer Evolution with Application to the FFCD 2000-05 Trial.
#' \emph{Biometrics} \bold{72}(3) 907-16.
#'
#' D. Rizopoulos (2012). Fast fitting of joint models for longitudinal and event
#' time data using a pseudo-adaptive Gaussian quadrature rule.
#' \emph{Computational Statistics and Data Analysis} \bold{56}, 491-501.
#'
#' M.S. Wulfsohn, A.A. and Tsiatis, A. A. (1997). A joint model for survival and
#' longitudinal data measured with error. \emph{Biometrics} \bold{53}, 330-9.
#'
#' A. Genz and B. Keister (1996). Fully symmetric interpolatory rules for
#' multiple integrals over infinite regions with Gaussian weight. \emph{Journal
#' of Computational and Applied Mathematics} \bold{71}, 299-309.
#'
#' D. Rustand, L. Briollais, C. Tournigand and V. Rondeau. Two-part joint model for a
#' longitudinal semicontinuous marker and a terminal event with application
#' to metastatic colorectal cancer data. \emph{Under
#' revision}.
#' @keywords models
#' @export
#' @examples
#'
#'
#' \dontrun{
#'
#' ###--- Joint model for longitudinal data and a terminal event ---###
#'
#' data(colorectal)
#' data(colorectalLongi)
#'
#' # Survival data preparation - only terminal events
#' colorectalSurv <- subset(colorectal, new.lesions == 0)
#'
#' # Baseline hazard function approximated with splines
#' # Random effects as the link function
#'
#' model.spli.RE <- longiPenal(Surv(time1, state) ~ age + treatment + who.PS
#' + prev.resection, tumor.size ~  year * treatment + age + who.PS ,
#' data=colorectalSurv,	data.Longi = colorectalLongi, random = c("1", "year"),
#' id = "id", link = "Random-effects", left.censoring = -3.33,
#' n.knots = 7, kappa = 2)
#'
#' # Weibull baseline hazard function
#' # Current level of the biomarker as the link function
#'
#' model.weib.CL <- longiPenal(Surv(time1, state) ~ age + treatment + who.PS
#' + prev.resection, tumor.size ~  year * treatment + age + who.PS , timevar="year",
#' data=colorectalSurv, data.Longi = colorectalLongi, random = c("1", "year"),
#' id = "id", link = "Current-level", left.censoring = -3.33, hazard = "Weibull")
#'
#'
#' ###--- Two-part Joint model for semicontinuous
#' #      longitudinal data and a terminal event ---###
#'
#' data(colorectal)
#' data(colorectalLongi)
#' colorectalSurv <- subset(colorectal, new.lesions == 0)
#'
#' # Box-cox back transformation (lambda=0.3) and apply logarithm (with a 1 unit shift)
#' colorectalLongi$Yo <- (colorectalLongi$tumor.size*0.3+1)^(1/0.3)
#' colorectalLongi$Y <- log(colorectalLongi$Y+1) # log transformation with shift=1
#'
#' # Two-part joint model - random-effects association structure (~15min)
#'
#' TwoPartJoint_re <-longiPenal(Surv(time1, state)~age + treatment +
#' who.PS+ prev.resection, Y~year*treatment, formula.Binary=Y~year*treatment,
#' data = colorectalSurv, data.Longi = colorectalLongi, random = c("1"),
#' random.Binary=c("1"), id = "id", link ="Random-effects", left.censoring = F,
#' n.knots = 7, kappa = 2, hazard="Splines-per")
#'
#' print(TwoPartJoint_re)
#'
#' # Two-part joint model - current-level association structure (~15min)
#' # Simulated dataset (github.com/DenisRustand/TPJM_sim)
#' data(longDat)
#' data(survDat)
#' tte <- frailtyPenal(Surv(deathTimes, d)~trt,n.knots=5,kappa=0, data=survDat,cross.validation = T)
#' kap <- round(tte$kappa,2);kap # smoothing parameter
#'   TPJM <- longiPenal(Surv(deathTimes, d)~trt, Y~timej*trtY,
#'   data=survDat, data.Longi = longDat,
#'   random = c("1","timej"), formula.Binary=Y~timej*trtY,
#'   random.Binary=c("1"), timevar="timej", id = "id",
#'   link = "Current-level", n.knots = 5, kappa = kap,
#'   hazard="Splines-per", method.GH="Monte-carlo",
#'   n.nodes=500, seed.MC=1)
#'
#'   print(TPJM)
#' }
"longiPenal" <-
  function (formula, formula.LongitudinalData, data,  data.Longi, formula.Binary=FALSE, random, random.Binary=FALSE, id, intercept = TRUE, link="Random-effects",timevar=FALSE,left.censoring=FALSE,n.knots, kappa,
            maxit=350, hazard="Splines",   init.B,
            init.Random, init.Eta, method.GH = "Standard",seed.MC=FALSE, n.nodes, LIMparam=1e-3, LIMlogl=1e-3, LIMderiv=1e-3, print.times=TRUE)
  {

    # Ajout de la fonction minmin issue de print.survfit, permettant de calculer la mediane
    minmin <- function(y, x) {
      tolerance <- .Machine$double.eps^.5   #same as used in all.equal()
      keep <- (!is.na(y) & y <(.5 + tolerance))
      if (!any(keep)) NA
      else {
        x <- x[keep]
        y <- y[keep]
        if (abs(y[1]-.5) <tolerance  && any(y< y[1]))
          (x[1] + x[min(which(y<y[1]))])/2
        else x[1]
      }
    }
GLMlog=FALSE
MTP=FALSE
    OrderLong <- data.Longi[,id]
    OrderDat <- data[,id]

    m2 <- match.call()
    m2$formula <-  m2$data <- m2$random <- m2$random.Binary <- m2$id <- m2$link <-m2$timevar <- m2$n.knots <- m2$kappa <- m2$maxit <- m2$hazard  <- m2$init.B <- m2$LIMparam <- m2$LIMlogl <- m2$LIMderiv <- m2$print.times <- m2$left.censoring <- m2$GLMlog <- m2$MTP <- m2$init.Random <- m2$init.Eta <- m2$method.GH <- m2$seed.MC<- m2$intercept <- m2$n.nodes <- m2$... <- NULL
    Names.data.Longi <- m2$data.Longi

    # TWO-PART indicator
    TwoPart <- ifelse(formula.Binary==FALSE, FALSE, TRUE) # true if two-part activated
    # more generally, variables related to the two-part model systematically contains the
    # letter B as for 'binary part', which is the main addition compared to a standard joint model.

    # association proba + current level only with two-part
        if(link=="Two-part"){
    if (TwoPart==F) {
    stop("The link 'Two-part' can only be used with a Two-part joint model.") }
}
    #### Frailty distribution specification ####
    if (!(all(random %in% c(1,names(data.Longi))))) { stop("Random effects can be only related to variables from the longitudinal data or the intercept (1)") }
    if (!(id %in% c(names(data.Longi))) || !(id %in% c(1,names(data)))) { stop("Identification for individuals can be only related to variables from both data set") }
    if(TwoPart){
    #### Frailty distribution specification (binary part) ####
    if (!(all(random.Binary %in% c(1,names(data.Longi))))) {
    stop("Random effects (binary part) can be only related to variables from the longitudinal data or the intercept (1)") }
}
    if(MTP & (!TwoPart | !GLMlog)){
    stop("Marginal two-part model requires activation of two-part and GLMlog")}

    #### Link function specification ####
    if(!(link %in% c("Random-effects","Current-level","Two-part"))){
      stop("Only 'Random-effects','Current-level' or 'Two-part' link function can be specified in link argument.")}
        ### Time variable
    if(link %in% c("Current-level","Two-part") & FALSE %in% timevar){
      stop("You must indicate the time variable in 'timevar' argument in order to use the 'Current-level' or 'Two-part' association")
    }

    ### Left-censoring
    if(!is.null(left.censoring) && left.censoring!=FALSE){
      if(!is.numeric(left.censoring))stop("If you want to include left-censored longitudinal outcome you must give the threshold value as the argument of 'left.censoring'")
    }

    ### Intercept
    if(!is.logical(intercept))stop("The argument 'intercept' must be logical")
    ### Gauss-Hermite method

    if(all(!(c("Standard","Pseudo-adaptive","HRMSYM", "Monte-carlo") %in% method.GH))){
      stop("Only 'Standard', 'Pseudo-adaptive', 'Monte-carlo' and 'HRMSYM' method can be specified as a method for the Gaussian quadrature")
    }
    GH <- switch(method.GH,"Standard"=0,"Pseudo-adaptive"=1,"HRMSYM"=2,"Monte-carlo"=3)

  if(GH<=2){
    if(!missing(n.nodes) ){
      if(!n.nodes%in%c(5,7,9,12,15,20,32)) stop("Number of points used in the numerical integration must be chosen from following: 5, 7, 9, 12, 15, 20, 32")
      if(n.nodes%in%c(5,7,9,12,15,20,32) && GH==2) warning("Using HRMSYM algorithm the number of points cannot be chosen")
    }else{
      n.nodes <- 9
    }}else if(GH==3){
    if(!missing(n.nodes) ){
      if(n.nodes<1) stop("Number of simulations for Monte-carlo must be positive.")
    }else{
      n.nodes <- 200 #default number of simulations for monte-carlo
    }
  }

    ##### hazard specification ######
    haztemp <- hazard
    hazard <- strsplit(hazard,split="-")
    hazard <- unlist(hazard)
    if(!(length(hazard) %in% c(1,2))){stop("Please check and revise the hazard argument according to the format specified in the help.")}

    ### longueur hazard = 1
    if((all.equal(length(hazard),1)==T)==T){
      if(!(hazard %in% c("Weibull","Piecewise","Splines"))){
        stop("Only 'Weibull', 'Splines' or 'Piecewise' hazard can be specified in hazard argument.")
      }else{
        typeof <- switch(hazard,"Splines"=0,"Piecewise"=1,"Weibull"=2)
        ### Splines (equidistant par defaut)
        if (typeof == 0){

          size1 <- 100
          size2 <- 100
          equidistant <- 1
          nbintervR <- 0
          nbintervDC <- 0
        }
        ### Weibull
        if (typeof == 2){

          size1 <- 100
          size2 <- 100
          equidistant <- 2
          nbintervR <- 0
          nbintervDC <- 0
        }
        if (typeof == 1){
          stop ("The hazard argument is incorrectly specified. Type of hazard are required ('per' or 'equi'). Please refer to the help file of frailtypack.")
        }
      }
    }else{
      #### longueur hazard > 1
      if(all(!(c("Splines","Piecewise") %in% hazard))){
        stop("Only 'Splines' and 'Piecewise' hazard can be specified in hazard argument in this case")
      }
      ### Splines percentile
      if ("Splines" %in% hazard){
        typeof <- 0
        if(!(all(hazard %in% c("Splines","per")))){
          stop ("The hazard argument is incorrectly specified. Only 'per' is allowed with 'Splines'. Please refer to the help file of frailtypack.")
        }else{
          size1 <- 100
          size2 <- 100
          equidistant <- 0
          nbintervR <- 0
          nbintervDC <- 0
        }
      }
      ### Piecewise (per or equi)
      if ("Piecewise" %in% hazard){
        typeof <- 1
        if(!(all(hazard %in% c("Piecewise","per","equi")))){
          stop ("The hazard argument is incorrectly specified. Type of hazard are required ('per' or 'equi'). Please refer to the help file of frailtypack.")
        }else{
          if (!(haztemp %in% c("Piecewise-per","Piecewise-equi"))){
            stop ("The hazard argument is incorrectly specified. Please refer to the help file of frailtypack.")
          }
          equidistant <- switch(haztemp,"Piecewise-per"=0,"Piecewise-equi"=1)
        }
      }
    }


    if (missing(formula))stop("The argument formula must be specified in every model")
    if (missing(formula.LongitudinalData))stop("The argument formula.LongitudinalData must be specified in every model") #AK

    if(class(formula)!="formula")stop("The argument formula must be a formula")
    if(TwoPart) if(formula.Binary!=FALSE) if(class(formula.Binary)!="formula")stop("The argument formula.Binary must be a formula")
    if(typeof == 0){
      if (missing(n.knots))stop("number of knots are required")
      n.knots.temp <- n.knots
      if (n.knots<4) n.knots<-4
      if (n.knots>20) n.knots<-20
      if (missing(kappa))stop("smoothing parameter (kappa) is required")

    }else{
      if (!(missing(n.knots)) || !(missing(kappa)) ){
        stop("When parametric hazard function is specified, 'kappa', 'n.knots' arguments must be deleted.")
      }
      n.knots <- 0
      kappa <- 0

    }
    call <- match.call()

    m <- match.call(expand.dots = FALSE) # recupere l'instruction de l'utilisateur

m$formula.LongitudinalData <- m$formula.Binary <- m$data.Longi <- m$n.knots <- m$random <- m$random.Binary <- m$link <- m$timevar <- m$id <- m$kappa <- m$maxit <- m$hazard  <- m$init.B <- m$LIMparam <- m$LIMlogl <- m$LIMderiv <- m$left.censoring <- m$GLMlog <- m$MTP <- m$print.times <- m$init.Random <- m$init.Eta <- m$method.GH <- m$seed.MC <- m$intercept <- m$n.nodes <- NULL


    special <- c("strata", "cluster", "subcluster", "terminal","num.id","timedep")

    #========= Longitudinal Data preparation =========================
    if(TwoPart){
        data.Binary=data.Longi # all data for binary part / add TwoPart
    }

    TermsY <- if (missing(data.Longi)){
      terms(formula.LongitudinalData, special)
    }else{
      terms(formula.LongitudinalData, special, data = data.Longi)
    }

    ord <- attr(TermsY, "order") # longueur de ord=nbre de var.expli

    #       if (length(ord) & any(ord != 1))stop("Interaction terms are not valid for this function")
    #si pas vide tous si il ya au moins un qui vaut 1 on arrete


    # create a table with all measuremnts for binary outcome, and a zero-free table for positive continuous outcome
    if(TwoPart){
    biom <- which(names(data.Longi)==as.character(attr(TermsY, "variables")[[2]])) # identifying biomarker
    data.Longi=data.Longi[data.Longi[,biom]>min(data.Longi[,biom]),]  # only positive biomarker values for semi-continuous part


    OrderLong <- data.Longi[, id]
    OrderBinary <- data.Binary[, id]
    }


    m2$formula <- TermsY


    m2[[1]] <- as.name("model.frame") # m[[1]]=frailtypenal, il le remplace par model.frame en fait


    #       m <- eval(m, sys.parent()) #ici la classe de m est un data.frame donc il recupere ce qu'on lui donne en argument

    clusterY <- attr(TermsY, "specials")$cluster # (indice) nbre de var qui sont en fonction de cluster()

    #Al : tri du jeu de donnees par cluster croissant
    if (length(clusterY)){
      tempc <- untangle.specials(TermsY, "cluster", 1:10)
      ord <- attr(TermsY, "order")[tempc$terms]
      if (any(ord > 1))stop("Cluster can not be used in an interaction")
      m2 <- m2[order(m2[,tempc$vars]),] # soit que des nombres, soit des caracteres
      ordre <- as.integer(row.names(m2)) # recupere l'ordre du data set
      clusterY <- strata(m2[, tempc$vars], shortlabel = TRUE)
      uni.clusterY <- unique(clusterY)
    }



    if (NROW(m2) == 0)stop("No (non-missing) observations") #nombre ligne different de 0

    llY <- attr(TermsY, "term.labels")#liste des variables explicatives


    #=========================================================>

    name.Y <- as.character(attr(TermsY, "variables")[[2]])
    Y <- data.Longi[,which(names(data.Longi)==name.Y)]

    # on identifie les variables explicatives facteurs avec nombre de niveau plus que 2

    ind.placeY <- which(llY%in%names(which(lapply(data.Longi[,which(names(data.Longi)%in%llY)],function(x) length(levels(x)))>2)))

    defined.factor <- llY[grep("factor",llY)]

    vec.factorY.tmp <- NULL
    if(length(defined.factor)>0){
      mat.factorY.tmp <- matrix(defined.factor,ncol=1,nrow=length(defined.factor))

      # Fonction servant a prendre les termes entre "as.factor"
      vec.factorY.tmp <-apply(mat.factorY.tmp,MARGIN=1,FUN=function(x){
        if (length(grep("factor",x))>0){
          if(length(grep(":",x))>0){
            if(grep('\\(',unlist(strsplit(x,split="")))[1]<grep(":",unlist(strsplit(x,split="")))[1] && length(grep('\\(',unlist(strsplit(x,split=""))))==1){

              pos1 <- grep("r",unlist(strsplit(x,split="")))[1]+2
              pos2 <- grep(":",unlist(strsplit(x,split="")))[1]-2
              pos3 <- grep(":",unlist(strsplit(x,split="")))[1]
              pos4 <- length(unlist(strsplit(x,split="")))
              if(length(levels(as.factor(data.Longi[,which(names(data.Longi)==substr(x,start=pos1,stop=pos2))])))>2)return(paste(substr(x,start=pos1,stop=pos2),substr(x,start=pos3,stop=pos4),sep=""))
              else return(NaN)
            }else if(grep("\\(",unlist(strsplit(x,split="")))[1]>grep(":",unlist(strsplit(x,split="")))[1] && length(grep('\\(',unlist(strsplit(x,split=""))))==1){
              pos2 <- grep(":",unlist(strsplit(x,split="")))[1]
              pos3 <- grep("\\(",unlist(strsplit(x,split="")))[1]+1
              pos4 <- length(unlist(strsplit(x,split="")))-1
              if(length(levels(as.factor(data.Longi[,which(names(data.Longi)==substr(x,start=pos3,stop=pos4))])))>2)return(paste(substr(x,start=1,stop=pos2),substr(x,start=pos3,stop=pos4),sep=""))
              else return(NaN)
            }else{#both factors
              pos1 <- grep("r",unlist(strsplit(x,split="")))[1]+2
              pos2 <- grep(":",unlist(strsplit(x,split="")))[1]-2
              pos3 <- grep("\\(",unlist(strsplit(x,split="")))[2]+1
              pos4 <- length(unlist(strsplit(x,split="")))-1
              if(length(levels(as.factor(data.Longi[,which(names(data.Longi)==substr(x,start=pos1,stop=pos2))])))>2 || length(levels(as.factor(data.Longi[,which(names(data.Longi)==substr(x,start=pos3,stop=pos4))])))>2)return(paste(substr(x,start=pos1,stop=pos2),":",substr(x,start=pos3,stop=pos4),sep=""))
              else return(NaN)
            }
          }else{
            pos1 <- grep("r",unlist(strsplit(x,split="")))[1]+2
            pos2 <- length(unlist(strsplit(x,split="")))-1
            if(length(levels(as.factor(data.Longi[,which(names(data.Longi)==substr(x,start=pos1,stop=pos2))])))>2)return(substr(x,start=pos1,stop=pos2))
            else return(NaN)
          }
        }else{
          return(x)
        }})

      vec.factorY.tmp <- vec.factorY.tmp[which(vec.factorY.tmp!="NaN")]

      if(length(vec.factorY.tmp)>0){
        for(i in 1:length(vec.factorY.tmp)){
          if(length(grep(":",vec.factorY.tmp[i]))==0){
            if(length(levels(as.factor(data.Longi[,which(names(data.Longi)==vec.factorY.tmp[i])])))>2)ind.placeY <- c(ind.placeY,which(llY%in%paste("as.factor(",vec.factorY.tmp[i],")",sep="")))
          }

        }}

    }

    ind.placeY <- sort(ind.placeY)
    #=========================================================>
    # On determine le nombre de categorie pour chaque var categorielle
    strats <- attr(TermsY, "specials")$strata #nbre de var qui sont en fonction de strata()
    cluster <- attr(TermsY, "specials")$cluster #nbre de var qui sont en fonction de cluster()
    num.id <- attr(TermsY, "specials")$num.id #nbre de var qui sont en fonction de patkey()
    vartimedep <- attr(TermsY, "specials")$timedep #nbre de var en fonction de timedep()

    #booleen pour savoir si au moins une var depend du tps
    if (is.null(vartimedep)) timedepY <- 0
    else timedepY <- 1


    if (timedepY==1) stop("The option 'timedep' is not allowed in this model.")

    if(is.null(num.id)){
      joint.clust <- 1
    }else{
      joint.clust <- 0
    }

    subcluster <- attr(TermsY, "specials")$subcluster #nbre de var qui sont en fonction de subcluster()

    if (length(subcluster))stop("'subcluster' is not an allowed option")
    if (length(cluster))stop("Only the argument 'id' can represent the clusters")

    Names.cluster <- id # nom du cluster

    if (length(num.id))stop("'num.id' is not an allowed option")

    if (length(strats))stop("Stratified analysis is not an allowed option yet")




    # which_n<-which(names(data.Longi)%in%random)
    #  data.Longi[which(data.Longi[,which_n]==0),which_n]<- 0.01

    mat.factorY2 <- matrix(llY,ncol=1,nrow=length(llY))

    # Fonction servant a prendre les termes entre "as.factor"
    llY2 <-apply(mat.factorY2,MARGIN=1,FUN=function(x){
      if (length(grep("factor",x))>0  && length(grep(":",x))==0 && unlist(strsplit(x,split=""))[length(unlist(strsplit(x,split="")))]==")"){
        pos1 <- grep("r",unlist(strsplit(x,split="")))[1]+2
        pos2 <- length(unlist(strsplit(x,split="")))-1
        x<-substr(x,start=pos1,stop=pos2)
        return(paste(x,levels(as.factor(data.Longi[,which(names(data.Longi)==x)]))[2],sep=""))
      }else{
        return(x)
      }})

    # Fonction servant a prendre les termes entre "as.factor" - without the name of the level
    llY3 <-apply(mat.factorY2,MARGIN=1,FUN=function(x){
      if (length(grep("factor",x))>0  && length(grep(":",x))==0 && unlist(strsplit(x,split=""))[length(unlist(strsplit(x,split="")))]==")"){
        pos1 <- grep("r",unlist(strsplit(x,split="")))[1]+2
        pos2 <- length(unlist(strsplit(x,split="")))-1
        return(substr(x,start=pos1,stop=pos2))
      }else{
        return(x)
      }})

    llY.real.names <- llY3
    llY3 <- llY3[!llY2%in%llY]

    data.Longi <- data.Longi[order(OrderLong),]
    if(is.factor(data.Longi[,names(data.Longi)==llY.real.names[1]])){
      X_L<- as.numeric(data.Longi[,names(data.Longi)==llY.real.names[1]])-1
    }
    else X_L <- data.Longi[,names(data.Longi)==llY.real.names[1]]



    if(length(llY)>1){
      for(i in 2:length(llY.real.names)){

        if(is.factor(data.Longi[,names(data.Longi)==llY.real.names[i]])){
          X_L <- cbind(X_L,as.numeric(data.Longi[,names(data.Longi)==llY.real.names[i]])-1)
        }
        else X_L <- cbind(X_L,data.Longi[,names(data.Longi)==llY.real.names[i]])
      }}


    #X_L<- data.Longi[,names(data.Longi)%in%(llY)]

    llY.fin <- llY.real.names
    llY <- llY.real.names
    j2<-0
    j3<-0
    if(sum(ord)>length(ord)){

      for(i in 1:length(ord)){
        if(ord[i]>1){

          name_v1 <- strsplit(as.character(llY[i]),":")[[1]][1]
          name_v2 <- strsplit(as.character(llY[i]),":")[[1]][2]

          if(length(grep("factor",name_v1))>0){name_v1<-substring(name_v1,11,nchar(name_v1)-1)
          v1 <- as.factor(data.Longi[,names(data.Longi)==name_v1])}
          else{v1 <- data.Longi[,names(data.Longi)==name_v1]}
          if(length(grep("factor",name_v2))>0){name_v2<-substring(name_v2,11,nchar(name_v2)-1)
          v2 <- as.factor(data.Longi[,names(data.Longi)==name_v2])}
          else{v2 <- data.Longi[,names(data.Longi)==name_v2]}

          llY[i] <- paste(name_v1,":",name_v2,sep="")
          #   if(is.factor(v1) && length(levels(v1))>2)stop("Interactions not allowed for factors with 3 or more levels (yet)")
          #   if(is.factor(v2) && length(levels(v2))>2)stop("Interactions not allowed for factors with 3 or more levels (yet)")
          if(is.factor(v1) && !is.factor(v2)){

            dummy <- model.matrix( ~ v1 - 1)
            # if(length(levels(v1)>2))vec.factorY <- c(vec.factorY,paste(name_v1,":",name_v2,sep=""))
            for(j in 2:length(levels(v1))){
              X_L <- cbind(X_L,dummy[,j]*v2)
              if(i>1 && i<length(llY.fin))llY.fin <- c(llY.fin[1:(i-1+j-2)],paste(name_v1,".",levels(v1)[j],":",name_v2,sep=""),llY.fin[(i+1+j-2):length(llY.fin)])
              else if(i==length(llY.fin))llY.fin <- c(llY.fin[1:(i-1+j-2)],paste(name_v1,".",levels(v1)[j],":",name_v2,sep=""))
              else llY.fin <- c(paste(name_v1,".",levels(v1)[j],":",name_v2,sep=""),llY.fin[(2+j-2):length(llY.fin)])
            }

          }else if(!is.factor(v1) && is.factor(v2)){

            dummy <- model.matrix( ~ v2 - 1)
            #  if(length(levels(v2)>2))vec.factorY <- c(vec.factorY,paste(name_v1,":",name_v2,sep=""))
            for(j in 2:length(levels(v2))){

              X_L <- cbind(X_L,dummy[,j]*v1)
              ## add because bug when interaction 2 modalities X 2 modalities
              if(j3==1) j4<-1 else j4<-0
              if(j2==1) j3 <- 1 else j3 <- 0
              if(j==3) j2<-1 else j2<-0

              if(i>1 && i<length(llY.fin) && j4==0)llY.fin <- c(llY.fin[1:(i-1+j+j3-2)],paste(name_v1,":",name_v2,levels(v2)[j],sep=""),llY.fin[(i+1+j-j2-2):length(llY.fin)])
              else if(j4==1)llY.fin <- c(llY.fin[1:(i-1+j-1)],paste(name_v1,":",name_v2,levels(v2)[j],sep=""))
              else if(i==length(llY.fin))llY.fin <- c(llY.fin[1:(i-1+j-2)],paste(name_v1,":",name_v2,levels(v2)[j],sep=""))
              else llY.fin <- c(paste(name_v1,":",name_v2,levels(v2)[j],sep=""),llY.fin[(2+j-2):length(llY.fin)])
            }
          }else if(is.factor(v1) && is.factor(v2)){


            dummy1 <- model.matrix( ~ v1 - 1)
            dummy2 <- model.matrix( ~ v2 - 1)
            #   if(length(levels(v1)>2) || length(levels(v2)>2))vec.factorY <- c(vec.factorY,paste(name_v1,":",name_v2,sep=""))
            for(j in 2:length(levels(v1))){
              for(k in 2:length(levels(v2))){

                X_L <- cbind(X_L,dummy1[,j]*dummy2[,k])
                if(i>1 && i<length(llY.fin))llY.fin <- c(llY.fin[1:(i-1+j-2+k-2)],paste(name_v1,levels(v1)[j],":",name_v2,levels(v2)[k],sep=""),llY.fin[(i+1+j-2+k-2):length(llY.fin)])
                else if(i==length(llY.fin))llY.fin <- c(llY.fin[1:(i-1+j-2+k-2)],paste(name_v1,levels(v1)[j],":",name_v2,levels(v2)[k],sep=""))
                else llY.fin <- c(paste(name_v1,levels(v1)[j],":",name_v2,levels(v2)[k],sep=""),llY.fin[(2+j-2+k-2):length(llY.fin)])
              }
            }
          }else{

            X_L <- cbind(X_L,v1*v2)
          }

        }
      }
    }

    if(length(grep(":",llY))>0){
      for(i in 1:length(grep(":",llY))){
        if(length(levels(data.Longi[,which(names(data.Longi)%in%strsplit(llY[grep(":",llY)[i]],":")[[1]])[1]]))>2 || length(levels(data.Longi[,which(names(data.Longi)%in%strsplit(llY[grep(":",llY)[i]],":")[[1]])[2]]))>2){
          ind.placeY <- c(ind.placeY,grep(":",llY)[i])
          #     vec.factorY <- c(vec.factorY,llY[grep(":",llY)[i]])
        }
      }
    }

    vec.factorY <- NULL

    if(length(vec.factorY.tmp)>0)vec.factorY <- c(llY[ind.placeY],vec.factorY.tmp)
    else vec.factorY <- c(vec.factorY,llY[ind.placeY])

    vec.factorY <- unique(vec.factorY)

    mat.factorY <- matrix(vec.factorY,ncol=1,nrow=length(vec.factorY))
    # Fonction servant a prendre les termes entre "as.factor" et (AK 04/11/2015) interactions
    vec.factorY <-apply(mat.factorY,MARGIN=1,FUN=function(x){
      if (length(grep("factor",x))>0){
        if(length(grep(":",x))>0){
          if(grep('\\(',unlist(strsplit(x,split="")))[1]<grep(":",unlist(strsplit(x,split="")))[1] && length(grep('\\(',unlist(strsplit(x,split=""))))==1){

            pos1 <- grep("r",unlist(strsplit(x,split="")))[1]+2
            pos2 <- grep(":",unlist(strsplit(x,split="")))[1]-2
            pos3 <- grep(":",unlist(strsplit(x,split="")))[1]
            pos4 <- length(unlist(strsplit(x,split="")))
            return(paste(substr(x,start=pos1,stop=pos2),substr(x,start=pos3,stop=pos4),sep=""))
          }else if(grep("\\(",unlist(strsplit(x,split="")))[1]>grep(":",unlist(strsplit(x,split="")))[1] && length(grep('\\(',unlist(strsplit(x,split=""))))==1){
            pos2 <- grep(":",unlist(strsplit(x,split="")))[1]
            pos3 <- grep("\\(",unlist(strsplit(x,split="")))[1]+1
            pos4 <- length(unlist(strsplit(x,split="")))-1
            return(paste(substr(x,start=1,stop=pos2),substr(x,start=pos3,stop=pos4),sep=""))
          }else{#both factors
            pos1 <- grep("r",unlist(strsplit(x,split="")))[1]+2
            pos2 <- grep(":",unlist(strsplit(x,split="")))[1]-2
            pos3 <- grep("\\(",unlist(strsplit(x,split="")))[2]+1
            pos4 <- length(unlist(strsplit(x,split="")))-1
            return(paste(substr(x,start=pos1,stop=pos2),":",substr(x,start=pos3,stop=pos4),sep=""))
          }
        }else{
          pos1 <- grep("r",unlist(strsplit(x,split="")))[1]+2
          pos2 <- length(unlist(strsplit(x,split="")))-1
          return(substr(x,start=pos1,stop=pos2))}
      }else{
        return(x)
      }})

    for(i in 1:length(llY.fin)){

      if(sum(names(data.Longi)==llY.fin[i])>0){
        if(is.factor(data.Longi[,names(data.Longi)==llY.fin[i]]) && length(levels(data.Longi[,names(data.Longi)==llY.fin[i]]))==2){
          llY.fin[i] <- paste(llY.fin[i],levels(data.Longi[,names(data.Longi)==llY.fin[i]])[2],sep="")}
      }
    }

    #  llY <- llY.fin
    if(dim(X_L)[2]!=length(llY.fin))stop("The variables in the longitudinal part must be in the data.Longi")
    X_L <- as.data.frame(X_L)
    names(X_L) <- llY.fin

    Intercept <- rep(1,dim(X_L)[1])

    if(intercept){
      X_L <- cbind(Intercept,X_L)
      ind.placeY <- ind.placeY+1
    }

    X_Lall<- X_L
    "%+%"<- function(x,y) paste(x,y,sep="")
    if(length(vec.factorY) > 0){
      for(i in 1:length(vec.factorY)){
        if(length(grep(":",vec.factorY[i]))==0){

          factor.spot <- which(names(X_L)==vec.factorY[i])
          if(length(factor.spot)>0){
          if(factor.spot<ncol(X_L))  X_L <- cbind(X_L[1:(factor.spot-1)],model.matrix(as.formula("~"%+%0%+%"+"%+%paste(vec.factorY[i], collapse= "+")), model.frame(~.,data.Longi,na.action=na.pass))[,-1],X_L[(factor.spot+1):ncol(X_L)])
          else X_L <- cbind(X_L[1:(factor.spot-1)],model.matrix(as.formula("~"%+%0%+%"+"%+%paste(vec.factorY[i], collapse= "+")), model.frame(~.,data.Longi,na.action=na.pass))[,-1])
          }
        } }



      vect.factY<-names(X_L)[which(!(names(X_L)%in%llY))]
      if(intercept) vect.factY <- vect.factY[-1]




      #               vect.fact <-apply(matrix(vect.fact,ncol=1,nrow=length(vect.fact)),MARGIN=1,FUN=function(x){
      #               pos1 <- grep("r",unlist(strsplit(x,split="")))[1]+2
      #               pos2 <- grep(")",unlist(strsplit(x,split="")))[1]-1
      #               return(substr(x,start=pos1,stop=pos2))})

      occurY <- rep(0,length(vec.factorY))

      #         for(i in 1:length(vec.factorY)){
      #                #occur[i] <- sum(vec.factor[i] == vect.fact)
      #               occurY[i] <- length(grep(vec.factorY[i],vect.factY))
      #      }




      interaction<-as.vector(apply(matrix(vect.factY,nrow=length(vect.factY)),MARGIN=1,FUN=function(x){length(grep(":",unlist(strsplit(x,split=""))))}))
      which.interaction <- which(interaction==1)

      for(i in 1:length(vec.factorY)){

        if(length(grep(":",unlist(strsplit(vec.factorY[i],split=""))))>0){


          pos <- grep(":",unlist(strsplit(vec.factorY[i],split="")))
          length.grep <- 0
          for(j in 1:length(vect.factY)){
            if(j%in%which.interaction){

              if(length(grep(substr(vec.factorY[i],start=1,stop=pos-1),vect.factY[j]))>0 && length(grep(substr(vec.factorY[i],start=pos+1,stop=length(unlist(strsplit(vec.factorY[i],split="")))),vect.factY[j]))>0){
                length.grep <- length.grep + 1
                which <- i}
            }}
          occurY[i] <- length.grep

        }else{


          if(length(vect.factY[-which.interaction])>0){occurY[i] <- length(grep(vec.factorY[i],vect.factY[-which.interaction]))
          }else{occurY[i] <- length(grep(vec.factorY[i],vect.factY))}
        }
      }
    }

   if (ncol(X_L) == 0){
   noVarY <- 1
  }else{
   noVarY <- 0
  }
  # X_L contains all covariates in longi part with dummys for multi-level


################
################
  ### ADD TWO-PART
#========= (longitudinal) Binary Data preparation =========================
    if(TwoPart){
        TermsB <- if (missing(data.Binary)){
                terms(formula.Binary, special)
        }else{
                terms(formula.Binary, special, data = data.Binary)
        }

        ordB <- attr(TermsB, "order") # ord identifies interaction as "2"

        clusterB <- attr(TermsB, "specials")$cluster # (indice) nbre de var qui sont en fonction de cluster()

#Al : tri du jeu de donnees par cluster croissant
        if (length(clusterB)){
                tempc <- untangle.specials(TermsB, "cluster", 1:10)
                ord <- attr(TermsB, "order")[tempc$terms]
                if (any(ord > 1))stop("Cluster can not be used in an interaction")
                m2 <- m2[order(m2[,tempc$vars]),] # soit que des nombres, soit des caracteres
                ordre <- as.integer(row.names(m2)) # recupere l'ordre du data set
                clusterB <- strata(m2[, tempc$vars], shortlabel = TRUE)
          uni.clusterB <- unique(clusterB)
        }

        llB <- attr(TermsB, "term.labels")#liste des variables explicatives


#=========================================================>
  name.B <- as.character(attr(TermsB, "variables")[[2]]) # biomarker name
  Binary <- ifelse(data.Binary[,which(names(data.Binary)==name.B)]>min(data.Binary[,which(names(data.Binary)==name.B)]), 1, 0)  # BINARY values (1=positive value)
  #left.censoring=FALSE # used for threshold, now deactivate
  # We model the probability to observe a positive value
  # We assume the binary covariates are the same as those in the semi-continuous formula
  #if(setdiff(llB, llY)!=0) stop("Covariates in the binary part must be included in the semi-continuous part.")

 # on identifie les variables explicatives facteurs avec nombre de niveau plus que 2
 # identifying 'factor' covariates with more than 2 levels
    ind.placeB <- which(llB%in%names(which(lapply(data.Binary[,which(names(data.Binary)%in%llB)],function(x) length(levels(x)))>2)))

  defined.factor <- llB[grep("factor",llB)]

  vec.factorB.tmp <- NULL
  if(length(defined.factor)>0){
    mat.factorB.tmp <- matrix(defined.factor,ncol=1,nrow=length(defined.factor))

    # Fonction servant a prendre les termes entre "as.factor"
    vec.factorB.tmp <-apply(mat.factorB.tmp,MARGIN=1,FUN=function(x){
      if (length(grep("factor",x))>0){
        if(length(grep(":",x))>0){
          if(grep('\\(',unlist(strsplit(x,split="")))[1]<grep(":",unlist(strsplit(x,split="")))[1] && length(grep('\\(',unlist(strsplit(x,split=""))))==1){

            pos1 <- grep("r",unlist(strsplit(x,split="")))[1]+2
            pos2 <- grep(":",unlist(strsplit(x,split="")))[1]-2
            pos3 <- grep(":",unlist(strsplit(x,split="")))[1]
            pos4 <- length(unlist(strsplit(x,split="")))
            if(length(levels(as.factor(data.Binary[,which(names(data.Binary)==substr(x,start=pos1,stop=pos2))])))>2)return(paste(substr(x,start=pos1,stop=pos2),substr(x,start=pos3,stop=pos4),sep=""))
            else return(NaN)
          }else if(grep("\\(",unlist(strsplit(x,split="")))[1]>grep(":",unlist(strsplit(x,split="")))[1] && length(grep('\\(',unlist(strsplit(x,split=""))))==1){
            pos2 <- grep(":",unlist(strsplit(x,split="")))[1]
            pos3 <- grep("\\(",unlist(strsplit(x,split="")))[1]+1
            pos4 <- length(unlist(strsplit(x,split="")))-1
            if(length(levels(as.factor(data.Binary[,which(names(data.Binary)==substr(x,start=pos3,stop=pos4))])))>2)return(paste(substr(x,start=1,stop=pos2),substr(x,start=pos3,stop=pos4),sep=""))
            else return(NaN)
          }else{#both factors
            pos1 <- grep("r",unlist(strsplit(x,split="")))[1]+2
            pos2 <- grep(":",unlist(strsplit(x,split="")))[1]-2
            pos3 <- grep("\\(",unlist(strsplit(x,split="")))[2]+1
            pos4 <- length(unlist(strsplit(x,split="")))-1
            if(length(levels(as.factor(data.Binary[,which(names(data.Binary)==substr(x,start=pos1,stop=pos2))])))>2 || length(levels(as.factor(data.Binary[,which(names(data.Binary)==substr(x,start=pos3,stop=pos4))])))>2)return(paste(substr(x,start=pos1,stop=pos2),":",substr(x,start=pos3,stop=pos4),sep=""))
            else return(NaN)
          }
        }else{
          pos1 <- grep("r",unlist(strsplit(x,split="")))[1]+2
          pos2 <- length(unlist(strsplit(x,split="")))-1
          if(length(levels(as.factor(data.Binary[,which(names(data.Binary)==substr(x,start=pos1,stop=pos2))])))>2)return(substr(x,start=pos1,stop=pos2))
          else return(NaN)
        }
      }else{
        return(x)
      }})

  vec.factorB.tmp <- vec.factorB.tmp[which(vec.factorB.tmp!="NaN")]

  if(length(vec.factorB.tmp)>0){
    for(i in 1:length(vec.factorB.tmp)){
        if(length(grep(":",vec.factorB.tmp[i]))==0){
          if(length(levels(as.factor(data.Binary[,which(names(data.Binary)==vec.factorB.tmp[i])])))>2)ind.placeB <- c(ind.placeB,which(llB%in%paste("as.factor(",vec.factorB.tmp[i],")",sep="")))
    }

    }}

}

ind.placeB <- sort(ind.placeB)



#=========================================================>
# On determine le nombre de categorie pour chaque var categorielle
        strats <- attr(TermsB, "specials")$strata #nbre de var qui sont en fonction de strata()
        cluster <- attr(TermsB, "specials")$cluster #nbre de var qui sont en fonction de cluster()
        num.id <- attr(TermsB, "specials")$num.id #nbre de var qui sont en fonction de patkey()
        vartimedep <- attr(TermsB, "specials")$timedep #nbre de var en fonction de timedep()

        #booleen pour savoir si au moins une var depend du tps
        if (is.null(vartimedep)) timedepB <- 0
        else timedepB <- 1


        if (timedepB==1) stop("The option 'timedep' is not allowed in this model.")
        subcluster <- attr(TermsB, "specials")$subcluster #nbre de var qui sont en fonction de subcluster()

  if (length(subcluster))stop("'subcluster' is not an allowed option")
        if (length(cluster))stop("Only the argument 'id' can represent the clusters")

                Names.cluster <- id # nom du cluster

  if (length(num.id))stop("'num.id' is not an allowed option")

        if (length(strats))stop("Stratified analysis is not an allowed option yet")





mat.factorB2 <- matrix(llB,ncol=1,nrow=length(llB))
llB2 <-apply(mat.factorB2,MARGIN=1,FUN=function(x){
    if (length(grep("factor",x))>0  && length(grep(":",x))==0 && unlist(strsplit(x,split=""))[length(unlist(strsplit(x,split="")))]==")"){
    pos1 <- grep("r",unlist(strsplit(x,split="")))[1]+2
    pos2 <- length(unlist(strsplit(x,split="")))-1
    x<-substr(x,start=pos1,stop=pos2)
    return(paste(x,levels(as.factor(data.Binary[,which(names(data.Binary)==x)]))[2],sep=""))
    }else{
    return(x)
}})

llB3 <-apply(mat.factorB2,MARGIN=1,FUN=function(x){
    if (length(grep("factor",x))>0  && length(grep(":",x))==0 && unlist(strsplit(x,split=""))[length(unlist(strsplit(x,split="")))]==")"){
    pos1 <- grep("r",unlist(strsplit(x,split="")))[1]+2
    pos2 <- length(unlist(strsplit(x,split="")))-1
    return(substr(x,start=pos1,stop=pos2))
    }else{
    return(x)
}})
llB.real.names <- llB3
llB3 <- llB3[!llB2%in%llB]

data.Binary <- data.Binary[order(OrderBinary),]

if(is.factor(data.Binary[,names(data.Binary)==llB.real.names[1]])){
    X_B<- as.numeric(data.Binary[,names(data.Binary)==llB.real.names[1]])-1
}
else X_B <- data.Binary[,names(data.Binary)==llB.real.names[1]]

if(length(llB)>1){ # number of covariates in binary formula (not including outcome)
  for(i in 2:length(llB.real.names)){
    if(is.factor(data.Binary[,names(data.Binary)==llB.real.names[i]])){
        X_B <- cbind(X_B,as.numeric(data.Binary[,names(data.Binary)==llB.real.names[i]])-1)
    }
    else X_B <- cbind(X_B,data.Binary[,names(data.Binary)==llB.real.names[i]])
}}
llB.fin <- llB.real.names
llB <- llB.real.names


if(sum(ordB)>length(ordB)){ # checking for interactions

for(i in 1:length(ordB)){
if(ordB[i]>1){

  name_v1 <- strsplit(as.character(llB[i]),":")[[1]][1] # first term of interaction
  name_v2 <- strsplit(as.character(llB[i]),":")[[1]][2] # second term

  if(length(grep("factor",name_v1))>0){name_v1<-substring(name_v1,11,nchar(name_v1)-1)
                                     v1 <- as.factor(data.Binary[,names(data.Binary)==name_v1])}
  else{v1 <- data.Binary[,names(data.Binary)==name_v1]} # vector of values for first term of interaction
  if(length(grep("factor",name_v2))>0){name_v2<-substring(name_v2,11,nchar(name_v2)-1)
                                     v2 <- as.factor(data.Binary[,names(data.Binary)==name_v2])}
  else{v2 <- data.Binary[,names(data.Binary)==name_v2]} # vector of values for second term

  llB[i] <- paste(name_v1,":",name_v2,sep="")
#   if(is.factor(v1) && length(levels(v1))>2)stop("Interactions not allowed for factors with 3 or more levels (yet)")
#   if(is.factor(v2) && length(levels(v2))>2)stop("Interactions not allowed for factors with 3 or more levels (yet)")
if(is.factor(v1) && !is.factor(v2)){

 dummy <- model.matrix( ~ v1 - 1)
# if(length(levels(v1)>2))vec.factorB <- c(vec.factorB,paste(name_v1,":",name_v2,sep=""))
 for(j in 2:length(levels(v1))){
   X_B <- cbind(X_B,dummy[,j]*v2)
   if(i>1 && i<length(llB.fin))llB.fin <- c(llB.fin[1:(i-1+j-2)],paste(name_v1,".",levels(v1)[j],":",name_v2,sep=""),llB.fin[(i+1+j-2):length(llB.fin)])
   else if(i==length(llB.fin))llB.fin <- c(llB.fin[1:(i-1+j-2)],paste(name_v1,".",levels(v1)[j],":",name_v2,sep=""))
   else llB.fin <- c(paste(name_v1,".",levels(v1)[j],":",name_v2,sep=""),llB.fin[(2+j-2):length(llB.fin)])
 }

}else if(!is.factor(v1) && is.factor(v2)){

  dummy <- model.matrix( ~ v2 - 1)
  for(j in 2:length(levels(v2))){

    X_B <- cbind(X_B,dummy[,j]*v1)

    if(i>1 && i<length(llB.fin))llB.fin <- c(llB.fin[1:(i-1+j-2)],paste(name_v1,":",name_v2,levels(v2)[j],sep=""),llB.fin[(i+1+j-2):length(llB.fin)])
    else if(i==length(llB.fin))llB.fin <- c(llB.fin[1:(i-1+j-2)],paste(name_v1,":",name_v2,levels(v2)[j],sep=""))
    else llB.fin <- c(paste(name_v1,":",name_v2,levels(v2)[j],sep=""),llB.fin[(2+j-2):length(llB.fin)])
    }
 }else if(is.factor(v1) && is.factor(v2)){

   dummy1 <- model.matrix( ~ v1 - 1)
   dummy2 <- model.matrix( ~ v2 - 1)
   for(j in 2:length(levels(v1))){
     for(k in 2:length(levels(v2))){

       X_B <- cbind(X_B,dummy1[,j]*dummy2[,k])
       if(i>1 && i<length(llB.fin))llB.fin <- c(llB.fin[1:(i-1+j-2+k-2)],paste(name_v1,levels(v1)[j],":",name_v2,levels(v2)[k],sep=""),llB.fin[(i+1+j-2+k-2):length(llB.fin)])
       else if(i==length(llB.fin))llB.fin <- c(llB.fin[1:(i-1+j-2+k-2)],paste(name_v1,levels(v1)[j],":",name_v2,levels(v2)[k],sep=""))
       else llB.fin <- c(paste(name_v1,levels(v1)[j],":",name_v2,levels(v2)[k],sep=""),llB.fin[(2+j-2+k-2):length(llB.fin)])
       }
   }
 }else{

 X_B <- cbind(X_B,v1*v2)
}

} # X_B => dataset with all binary covariates (starting with intercept) / +interaction terms (still need to create dummys for multilevels)
}
}



if(length(grep(":",llB))>0){
  for(i in 1:length(grep(":",llB))){
    if(length(levels(data.Binary[,which(names(data.Binary)%in%strsplit(llB[grep(":",llB)[i]],":")[[1]])[1]]))>2 || length(levels(data.Binary[,which(names(data.Binary)%in%strsplit(llB[grep(":",llB)[i]],":")[[1]])[2]]))>2){
      ind.placeB <- c(ind.placeB,grep(":",llB)[i])
 #     vec.factorB <- c(vec.factorB,llB[grep(":",llB)[i]])
    }
  }
}

  # Creating dummys for multi-level covariates
vec.factorB <- NULL

if(length(vec.factorB.tmp)>0)vec.factorB <- c(llB[ind.placeB],vec.factorB.tmp)
else vec.factorB <- c(vec.factorB,llB[ind.placeB])

vec.factorB <- unique(vec.factorB)

mat.factorB <- matrix(vec.factorB,ncol=1,nrow=length(vec.factorB))
# Fonction servant a prendre les termes entre "as.factor" et (AK 04/11/2015) interactions
vec.factorB <-apply(mat.factorB,MARGIN=1,FUN=function(x){
  if (length(grep("factor",x))>0){
    if(length(grep(":",x))>0){
      if(grep('\\(',unlist(strsplit(x,split="")))[1]<grep(":",unlist(strsplit(x,split="")))[1] && length(grep('\\(',unlist(strsplit(x,split=""))))==1){

        pos1 <- grep("r",unlist(strsplit(x,split="")))[1]+2
        pos2 <- grep(":",unlist(strsplit(x,split="")))[1]-2
        pos3 <- grep(":",unlist(strsplit(x,split="")))[1]
        pos4 <- length(unlist(strsplit(x,split="")))
        return(paste(substr(x,start=pos1,stop=pos2),substr(x,start=pos3,stop=pos4),sep=""))
      }else if(grep("\\(",unlist(strsplit(x,split="")))[1]>grep(":",unlist(strsplit(x,split="")))[1] && length(grep('\\(',unlist(strsplit(x,split=""))))==1){
        pos2 <- grep(":",unlist(strsplit(x,split="")))[1]
        pos3 <- grep("\\(",unlist(strsplit(x,split="")))[1]+1
        pos4 <- length(unlist(strsplit(x,split="")))-1
        return(paste(substr(x,start=1,stop=pos2),substr(x,start=pos3,stop=pos4),sep=""))
      }else{#both factors
        pos1 <- grep("r",unlist(strsplit(x,split="")))[1]+2
        pos2 <- grep(":",unlist(strsplit(x,split="")))[1]-2
        pos3 <- grep("\\(",unlist(strsplit(x,split="")))[2]+1
        pos4 <- length(unlist(strsplit(x,split="")))-1
        return(paste(substr(x,start=pos1,stop=pos2),":",substr(x,start=pos3,stop=pos4),sep=""))
      }
    }else{
      pos1 <- grep("r",unlist(strsplit(x,split="")))[1]+2
      pos2 <- length(unlist(strsplit(x,split="")))-1
      return(substr(x,start=pos1,stop=pos2))}
  }else{
    return(x)
  }})

for(i in 1:length(llB.fin)){

  if(sum(names(data.Binary)==llB.fin[i])>0){
  if(is.factor(data.Binary[,names(data.Binary)==llB.fin[i]]) && length(levels(data.Binary[,names(data.Binary)==llB.fin[i]]))==2){
    llB.fin[i] <- paste(llB.fin[i],levels(data.Binary[,names(data.Binary)==llB.fin[i]])[2],sep="")
  }
    #else if(is.factor(data.Binary[,names(data.Binary)==llB.fin[i]]) && length(levels(data.Binary[,names(data.Binary)==llB.fin[i]]))==3){
    #llB.fin[j] <- paste(llB.fin[i],levels(data.Binary[,names(data.Binary)==llB.fin[i]])[2],sep="")
    #llB.fin[j+1] <- paste(llB.fin[i],levels(data.Binary[,names(data.Binary)==llB.fin[i]])[3],sep="")
    #j=j+2
    #}
  }
}

#  llB <- llB.fin
 # if(dim(X_B)[2]!=length(llB.fin))stop("The variables in the longitudinal part must be in the data.Binary")
   X_B <- as.data.frame(X_B)
   names(X_B) <- llB.fin

Intercept.Binary <- rep(1,dim(X_B)[1])

  if(intercept){
    X_B <- cbind(Intercept.Binary,X_B)
    ind.placeB <- ind.placeB+1
  }

  X_Ball<- X_B
  "%+%"<- function(x,y) paste(x,y,sep="")
        if(length(vec.factorB) > 0){
          for(i in 1:length(vec.factorB)){
          if(length(grep(":",vec.factorB[i]))==0){

          factor.spot <- which(names(X_B)==vec.factorB[i])

              if(factor.spot<ncol(X_B))  X_B <- cbind(X_B[1:(factor.spot-1)],model.matrix(as.formula("~"%+%0%+%"+"%+%paste(vec.factorB[i], collapse= "+")), model.frame(~.,data.Binary,na.action=na.pass))[,-1],X_B[(factor.spot+1):ncol(X_B)])
     else X_B <- cbind(X_B[1:(factor.spot-1)],model.matrix(as.formula("~"%+%0%+%"+"%+%paste(vec.factorB[i], collapse= "+")), model.frame(~.,data.Binary,na.action=na.pass))[,-1])

         } }



  vect.factB<-names(X_B)[which(!(names(X_B)%in%llB))]
  if(intercept) vect.factB <- vect.factB[-1]




#               vect.fact <-apply(matrix(vect.fact,ncol=1,nrow=length(vect.fact)),MARGIN=1,FUN=function(x){
#               pos1 <- grep("r",unlist(strsplit(x,split="")))[1]+2
#               pos2 <- grep(")",unlist(strsplit(x,split="")))[1]-1
#               return(substr(x,start=pos1,stop=pos2))})

                occurB <- rep(0,length(vec.factorB))

       #         for(i in 1:length(vec.factorB)){
        #                #occur[i] <- sum(vec.factor[i] == vect.fact)
         #               occurB[i] <- length(grep(vec.factorB[i],vect.factB))
          #      }




  interaction<-as.vector(apply(matrix(vect.factB,nrow=length(vect.factB)),MARGIN=1,FUN=function(x){length(grep(":",unlist(strsplit(x,split=""))))}))
  which.interaction <- which(interaction==1)

  for(i in 1:length(vec.factorB)){

  if(length(grep(":",unlist(strsplit(vec.factorB[i],split=""))))>0){


    pos <- grep(":",unlist(strsplit(vec.factorB[i],split="")))
    length.grep <- 0
    for(j in 1:length(vect.factB)){
      if(j%in%which.interaction){

        if(length(grep(substr(vec.factorB[i],start=1,stop=pos-1),vect.factB[j]))>0 && length(grep(substr(vec.factorB[i],start=pos+1,stop=length(unlist(strsplit(vec.factorB[i],split="")))),vect.factB[j]))>0){
          length.grep <- length.grep + 1
          which <- i}
      }}
    occurB[i] <- length.grep

  }else{


    if(length(vect.factB[-which.interaction])>0){occurB[i] <- length(grep(vec.factorB[i],vect.factB[-which.interaction]))
    }else{occurB[i] <- length(grep(vec.factorB[i],vect.factB))}
  }
}
}

  if (ncol(X_B) == 0){
   noVarB <- 1
  }else{
   noVarB <- 0
  }
  }

  if(!exists("noVarB")) noVarB <- 1
  # X_B contains all covariates in longi part with dummys for multi-level
################
################


    #=========================================================>

    clusterY <- data.Longi[,which(colnames(data.Longi)==id)]
if(TwoPart) clusterB <- data.Binary[,which(colnames(data.Binary)==id)]
if(TwoPart) max_repB <- max(table(clusterB))

    max_rep <- max(table(clusterY))
    uni.cluster<-as.factor(unique(clusterY))

    if(is.null(id)) stop("grouping variable is needed")

    if(is.null(random))     stop("variable for random effects is needed")

    if(TwoPart) if(is.null(random.Binary)) stop("variable for binary part random effects is needed")
    if(length(uni.cluster)==1){
      stop("grouping variable must have more than 1 level")
    }

    if (length(subcluster))stop("'Subcluster' is not allowed")

    if (typeof==0 && missing(kappa)) stop("smoothing parameter (kappa) is required")

    if ((typeof==0) & (length(kappa)!=1)) stop("wrong length of argument 'kappa'")

    #newTerm vaut Terms - les variables dont les position sont dans drop

    #========================================>

    nvarY<-ncol(X_L) #nvar==1 correspond a 2 situations:

    # au cas ou on a aucune var explicative dans la partie rec, mais X=0
    # cas ou on a 1seul var explicative, ici X est en general different de 0

    #  varnotdepY <- colnames(X_L)[-grep("timedep",colnames(X_L))]
    #  vardepY <- colnames(X_L)[grep("timedep",colnames(X_L))]
    #  vardepY <- apply(matrix(vardepY,ncol=1,nrow=length(vardepY)),1,timedep.names)

    #    if (length(intersect(varnotdepY,vardepY)) != 0) {
    #           stop("A variable is both used as a constant and time-varying effect covariate")
    #   }

    #   nvartimedepY <- length(vardepY)

    #    filtretpsY <- rep(0,nvarY)
    #    filtretpsY[grep("timedep",colnames(X_L))] <- 1

    varY <- as.matrix(sapply(X_L, as.numeric))

    nsujety<-nrow(X_L)

    if(TwoPart) nvarB<-ncol(X_B)
    if(TwoPart) varB <- as.matrix(sapply(X_B, as.numeric))
    if(TwoPart) nsujetB <- nrow(X_B)

    #=======================================>
    #======= Construction du vecteur des indicatrice
    if(length(vec.factorY) > 0){
      #               ind.place <- ind.place -1
      k <- 0
      for(i in 1:length(vec.factorY)){
        ind.placeY[i] <- ind.placeY[i]+k
        k <- k + occurY[i]-1
      }
    }

    if(TwoPart){
        if(length(vec.factorB) > 0){
        kB <- 0
        for(i in 1:length(vec.factorB)){
            ind.placeB[i] <- ind.placeB[i]+kB
            kB <- kB + occurB[i]-1
        }
    }}
    # Random effects

    if(link=="Random-effects") link0 <- 1
    if(link=="Current-level") link0 <- 2
       if(link=="Two-part") link0 <- 3

    if(TwoPart){
      nREY <- length(random)
      nREB <- length(random.Binary)
      nRE <- nREY+nREB
    }else{
      nRE <- length(random)
    }
    ne_re <- nRE*(nRE+1)/2


    matzy <- NULL
    names.matzy <- NULL
    if(1%in%random){
      names.matzy<-c("Intercept",random[-which(random==1)])
    }else{
      names.matzy<-random
    }

    if(TwoPart){ # vector of random effects names for binary part
    matzB <- NULL
    names.matzB <- NULL


      if(1%in%random.Binary){
        names.matzB<-c("Intercept.Binary",random.Binary[-which(random.Binary==1)])
      }else{
        names.matzB<-random.Binary
      }
    }

    matzy <- data.matrix(X_L[,which(names(X_L)%in%names.matzy)])

    if(TwoPart){
        matzB <- data.matrix(X_B[,which(names(X_B)%in%names.matzB)])
    }

    if(!intercept && 1%in%random) matzy <- as.matrix(cbind(rep(1,nsujety),matzy))


    # next control is not relevant with multiple functions of time, removed for now
    #if(dim(matzy)[2]>=3 && link0 == 2)stop("The current-level link can be chosen only if the biomarker random effects are associated with the intercept and time.")

    if(link0==1)netadc <- ncol(matzy)
    if(link0==1)if(TwoPart) netadc <- netadc+ncol(matzB) #add TwoPart
    if(link0==2)netadc <- 1
    if(link0==3)netadc <- 2


    #== Left-censoring ==
    cag <- c(0,0)

    if(!is.null(left.censoring) && is.numeric(left.censoring)){
      if(TwoPart) stop("No left-censoring if Two-Part model is activated.")
      if(left.censoring<min(Y))stop("The threshold for the left censoring cannot be smaller than the minimal value of the longitudinal outcome")
      cag[1] <- 1
      cag[2] <- left.censoring
      n.censored <- length(which(Y<=left.censoring))
      prop.censored <- n.censored/nsujety
    }



    #============= pseudo-adaptive Gauss Hermite ==============
    #m <- lme(measuret ~ time+interact+treatment, data = data, random = ~ 1| idd)
    if(method.GH=="Pseudo-adaptive"){
    if(TwoPart) stop("Pseudo-adaptive not defined yet for Two-Part models.")
      if(length(random)>2){
        random_lme <- random[2]
        for(i in 3:length(random)){
          random_lme <- paste(random_lme, "+", random[i])
        }}else{random_lme <- random}
      inn <-paste("pdSymm(form=~",random_lme,")",sep="")
      rand <- list(eval(parse(text=inn)))
      names(rand) <- id

      m_lme<-lme(formula.LongitudinalData,data = data.Longi, random = rand, control=lmeControl(opt='optim'))

      b_lme <-as.matrix(ranef(m_lme))

      formula_lme <- formula(m_lme$modelStruct$reStruct[[1]])
      model_lme <- model.frame(terms(formula_lme), data = data.Longi)
      Terms_lme <- attr(model_lme, "terms")
      Z_lme <- model.matrix(formula_lme, model_lme)
      #id <- as.vector(unclass(m2$groups[[1]]))
      # Cholesky matrices for GH
      #B_lme <- lapply(pdMatrix(m_lme$modelStruct$reStruct), "*",  m_lme$sigma^2)[[1]]
      B_lme <-pdMatrix(m_lme$modelStruct$reStruct)[[1]]

      Bi_tmp  <- vector("list", length(uni.cluster) )
      invBi_chol <- matrix(rep(0,ne_re*length(uni.cluster)),nrow=length(uni.cluster) )

      invB_lme  <- solve(B_lme )
      invB_chol<- chol(B_lme)
      invB_cholDet <- det(invB_chol)
      for (i in 1:length(uni.cluster )) {
        Zi_lme <- Z_lme[clusterY  == i, , drop = FALSE]

        Bi_tmp[[i]] <- chol(solve(crossprod(Zi_lme) / m_lme$sigma^2 + invB_lme))
        for(j in 1:nRE){
          for(k in 1:nRE){
            if (k>=j)invBi_chol[i,j+k*(k-1)/2] <-  Bi_tmp[[i]][j,k]
            #  else     inv.chol.VC(j,k)=matv[k+j*(j-1)/2]
          }}
      }
      #Vs[[i]] <- solve(crossprod(Z.i) / m2$sigma^2 + inv.VC)
      #}

      #invBi_chol <- lapply(Bi_lme, function (x) solve(chol(solve(x))))
      invBi_cholDet <- sapply(Bi_tmp,  det)
    }else{

      b_lme <- matrix(rep(0,length(uni.cluster)*nRE),ncol=nRE,nrow=length(uni.cluster))
      invBi_cholDet <-  matrix(rep(0,length(uni.cluster)),ncol=1,nrow=length(uni.cluster))
      invBi_chol <- matrix(rep(0,length(uni.cluster)*ne_re),ncol=ne_re,nrow=length(uni.cluster))

    }

    #================== Survival Data ===================

    Terms <- if (missing(data)){
      terms(formula, special)
    }else{
      terms(formula, special, data = data)
    }

    ord <- attr(Terms, "order") # longueur de ord=nbre de var.expli

    # if (length(ord) & any(ord != 1))stop("Interaction terms are not valid for this function")
    #si pas vide tous si il ya au moins un qui vaut 1 on arrete

    m$formula <- Terms

    m[[1]] <- as.name("model.frame") # m[[1]]=frailtypenal, il le remplace par model.frame en fait

    #model.frame(formula = Surv(time, event) ~ cluster(id) + as.factor(dukes) +
    #as.factor(charlson) + sex + chemo + terminal(death), data = readmission)

    m <- eval(m, sys.parent()) #ici la classe de m est un data.frame donc il recupere ce qu'on lui donne en argument

    cluster <- id # (indice) nbre de var qui sont en fonction de cluster()

    subcluster <- attr(Terms, "specials")$subcluster #nbre de var qui sont en fonction de subcluster()

    # booleen pour voir si l'objet Y est reponse avant tri des donnees Surv ou SurvIC
    classofT <- attr(model.extract(m, "response"),"class")
    # attention le package pec rajoute un element dans l'attribut "class" des objets de survie
    if (length(classofT)>1) classofT <- classofT[2]

    typeofT <- attr(model.extract(m, "response"),"type") # type de reponse : interval etc..

    # verification de la structure nested si besoin
    if (length(subcluster))stop("subcluster can not be used in the model")

    # tri par ordre croissant de subcluster a l'interieur des clusters
    # ordre <- as.integer(row.names(m))

    if (NROW(m) == 0)stop("No (non-missing) observations") #nombre ligne different de 0

    T <- model.extract(m, "response") # objet de type Surv =Time

    if (classofT != "Surv") stop("Response must be a survival object")


    llT <- attr(Terms, "term.labels")#liste des variables explicatives

    #=========================================================>

    mt <- attr(m, "terms") #m devient de class "formula" et "terms"

    X_T <- if (!is.empty.model(mt))model.matrix(mt, m) #idem que mt sauf que ici les factor sont divise en plusieurs variables

    ind.placeT <- unique(attr(X_T,"assign")[duplicated(attr(X_T,"assign"))]) ### unique : changement au 25/09/2014


    vec.factorT <- NULL
    vec.factorT <- c(vec.factorT,llT[ind.placeT])

    #=========================================================>
    # On determine le nombre de categorie pour chaque var categorielle
    strats <- attr(Terms, "specials")$strata #nbre de var qui sont en fonction de strata()
    cluster <- attr(Terms, "specials")$cluster #nbre de var qui sont en fonction de cluster()
    num.id <- attr(Terms, "specials")$num.id #nbre de var qui sont en fonction de patkey()
    vartimedep <- attr(Terms, "specials")$timedep #nbre de var en fonction de timedep()

    #booleen pour savoir si au moins une var depend du tps
    if (is.null(vartimedep)) timedepT <- 0
    else timedepT <- 1

    if (timedepT==1) stop("The option 'timedep' is not allowed in this model.")

    if(length(num.id)) stop("num.id function is not allowed")

    subcluster <- attr(Terms, "specials")$subcluster #nbre de var qui sont en fonction de subcluster()

    if (length(subcluster))stop("'subcluster' is not allowed")
    if (length(cluster))stop("Clusters are defined by the argument 'id'")

    if (length(strats))stop("Stratified analysis is not allowed")


    mat.factorT <- matrix(vec.factorT,ncol=1,nrow=length(vec.factorT))

    # Fonction servant a prendre les termes entre "as.factor"
    vec.factorT <-apply(mat.factorT,MARGIN=1,FUN=function(x){
      if (length(grep("factor",x))>0){
        if(length(grep(":",x))>0){
          if(grep('\\(',unlist(strsplit(x,split="")))[1]<grep(":",unlist(strsplit(x,split="")))[1] && length(grep('\\(',unlist(strsplit(x,split=""))))==1){

            pos1 <- grep("r",unlist(strsplit(x,split="")))[1]+2
            pos2 <- grep(":",unlist(strsplit(x,split="")))[1]-2
            pos3 <- grep(":",unlist(strsplit(x,split="")))[1]
            pos4 <- length(unlist(strsplit(x,split="")))
            return(paste(substr(x,start=pos1,stop=pos2),substr(x,start=pos3,stop=pos4),sep=""))
          }else if(grep("\\(",unlist(strsplit(x,split="")))[1]>grep(":",unlist(strsplit(x,split="")))[1] && length(grep('\\(',unlist(strsplit(x,split=""))))==1){
            pos2 <- grep(":",unlist(strsplit(x,split="")))[1]
            pos3 <- grep("\\(",unlist(strsplit(x,split="")))[1]+1
            pos4 <- length(unlist(strsplit(x,split="")))-1
            return(paste(substr(x,start=1,stop=pos2),substr(x,start=pos3,stop=pos4),sep=""))
          }else{#both factors
            pos1 <- grep("r",unlist(strsplit(x,split="")))[1]+2
            pos2 <- grep(":",unlist(strsplit(x,split="")))[1]-2
            pos3 <- grep("\\(",unlist(strsplit(x,split="")))[2]+1
            pos4 <- length(unlist(strsplit(x,split="")))-1
            return(paste(substr(x,start=pos1,stop=pos2),":",substr(x,start=pos3,stop=pos4),sep=""))
          }
        }else{
          pos1 <- grep("r",unlist(strsplit(x,split="")))[1]+2
          pos2 <- length(unlist(strsplit(x,split="")))-1
          return(substr(x,start=pos1,stop=pos2))}
      }else{
        return(x)
      }})

    # On determine le nombre de categorie pour chaque var categorielle
    if(length(vec.factorT) > 0){
      vect.factT <- attr(X_T,"dimnames")[[2]]
      vect.factT <- vect.factT[grep(paste(vec.factorT,collapse="|"),vect.factT)]

      occurT <- rep(0,length(vec.factorT))
      #    for(i in 1:length(vec.factordc)){
      #      occurdc[i] <- length(grep(vec.factordc[i],vect.factdc))
      #    }
      #  }
      interactionT<-as.vector(apply(matrix(vect.factT,nrow=length(vect.factT)),MARGIN=1,FUN=function(x){length(grep(":",unlist(strsplit(x,split=""))))}))
      which.interactionT <- which(interactionT==1)
      for(i in 1:length(vec.factorT)){

        if(length(grep(":",unlist(strsplit(vec.factorT[i],split=""))))>0){
          pos <- grep(":",unlist(strsplit(vec.factorT[i],split="")))
          length.grep <- 0
          for(j in 1:length(vect.factT)){
            if(j%in%which.interactionT){
              if(length(grep(substr(vec.factorT[i],start=1,stop=pos-1),vect.factT[j]))>0 && length(grep(substr(vec.factorT[i],start=pos+1,stop=length(unlist(strsplit(vec.factorT[i],split="")))),vect.factT[j]))>0){

                length.grep <- length.grep + 1
                which <- i}
            }}
          occurT[i] <- length.grep

        }else{


          if(length(vect.factT[-which.interactionT])>0){occurT[i] <- length(grep(vec.factorT[i],vect.factT[-which.interactionT]))
          }else{occurT[i] <- length(grep(vec.factorT[i],vect.factT))}
        }
      }
    }



    #=========================================================>

    dropx <- NULL


    clusterT <- 1:nrow(m) #nrow(data) # valeurs inutiles pour un modele de Cox
    uni.clusterT <- 1:nrow(m) #nrow(data)

    if(length(uni.cluster)==1)stop("grouping variable must have more than 1 level")
    if (length(subcluster))stop("subcluster can not be used in the model")
    if (length(strats))stop("Stratified analysis not allowed")


    #type <- attr(Y, "type")
    type <- typeofT

    if (type != "right" && type != "counting" && type != "interval" && type != "intervaltronc") { # Cox supporte desormais la censure par intervalle
      stop(paste("Cox model doesn't support \"", type, "\" survival data", sep = ""))
    }

    #       if ((type == "interval" || type == "interval2" || type == "intervaltronc") && intcens == FALSE) { # rajout
    #               stop("You are trying to do interval censoring without intcens = TRUE")
    #       }



    #newTerm vaut Terms - les variables dont les position sont dans drop

    X_T <- model.matrix(Terms, m)

    assign <- lapply(attrassign(X_T, Terms)[-1], function(x) x - 1)
    Xlevels <- .getXlevels(Terms, m)
    contr.save <- attr(X_T, 'contrasts')

    # assigne donne la position pour chaque variables
    #ncol(X) : nombre de variable sans sans les fonction speciaux comme terminal()...+id
    if(length(vec.factorT) > 0){
      #========================================>
      position <- unlist(assign,use.names=F)
    }

    #========================================>

    if (ncol(X_T) == 1){
      X_T<-X_T-1
      noVarT <- 1
    }else{
      X_T <- X_T[, -1, drop = FALSE]
      noVarT <- 0
    }
    # on enleve ensuite la premiere colonne correspondant a id
    nvarT<-ncol(X_T) #nvar==1 correspond a 2 situations:


    # au cas ou on a aucune var explicative dans la partie rec, mais X=0
    # cas ou on a 1seul var explicative, ici X est en general different de 0

    #  varnotdepT <- colnames(X_T)[-grep("timedep",colnames(X_T))]
    #  vardepT <- colnames(X_T)[grep("timedep",colnames(X_T))]
    #  vardepT <- apply(matrix(vardepT,ncol=1,nrow=length(vardepT)),1,timedep.names)

    # if (length(intersect(varnotdepT,vardepT)) != 0) {
    #    stop("A variable is both used as a constant and time-varying effect covariate")
    # }

    #  nvartimedepT <- length(vardepT)

    #  filtretpsT <- rep(0,nvarT)
    #  filtretpsT[grep("timedep",colnames(X_T))] <- 1

    varT<-matrix(c(X_T),nrow=nrow(X_T),ncol=nvarT) #matrix sans id et sans partie ex terminal(death)
    varT <- varT[order(OrderDat),]
    T <- T[order(OrderDat)]  # add myriam 24/04/17
    ng<-nrow(X_T)

    #add Alexandre 04/06/2012
    #lire les donnees differemment si censure par intervalle

    if (type=="right"){
      tt0dc <- rep(0,ng)
      tt1dc <- T[,1]
      cens <- T[,2]
    } else {
      tt0dc <- T[,1]
      tt1dc <- T[,2]
      cens <- T[,3]
    }                   # attention ne pas mettre de 0 sinon en cas de left trunc probleme dans la logV


    if (min(cens)==0) cens.data<-1
    if (min(cens)==1 && max(cens)==1) cens.data<-0

    AG<-0

    #=======================================>
    #======= Construction du vecteur des indicatrice
    if(length(vec.factorT) > 0){
      #             ind.place <- ind.place -1
      k <- 0
      for(i in 1:length(vec.factorT)){
        ind.placeT[i] <- ind.placeT[i]+k
        k <- k + occurT[i]-1
      }
    }

    #
    # Begin JOINT MODEL
    #

    # Preparing data ...
    if(TwoPart){
      nvar = nvarB + nvarY + nvarT
    }else{
      nvar = nvarY + nvarT # total number of fixed parameters
    }
      Y <- Y[order(OrderLong)] # id ordering
      if(TwoPart) Binary <- Binary[order(OrderBinary)] # id ordering binary outcome

      if ((typeof == 0) | (typeof == 2)) indic.nb.int <- 0 # 0 is splines ; 2 is weibull

      if (sum(as.double(varT))==0) nvarT <- 0
      if (sum(as.double(varY))==0) nvarY <- 0
    if(TwoPart) if (sum(as.double(varB))==0) nvarB <- 0

# np is total number of parameters (including knots locations for splines baseline hazard (or weibull parameters);
# fixed parameters asociated to covariates; error term; variance and covariance of random effects and association parameter(s))
    np <- switch(as.character(typeof),
                 "0"=((as.integer(n.knots) + 2) + as.integer(nvar) + 1 + ne_re + netadc  ),
                 "2"=(2 + nvar + 1  + ne_re + netadc ))

    # traitement de l'initialisation du Beta rentre par l'utilisateur

    Beta <- rep(0.5,np)
    if (!missing(init.B)) {
      if (length(init.B) != nvar) stop("Wrong number of regression coefficients in init.B")
      #  if (timedep) stop("You can hardly know initial regression coefficient while time-varying effect")
      Beta <- c(rep(0.5,np-nvar),init.B)
    }

    if (!missing(init.Random)) {
      if (length(init.Random)!=ne_re) stop("init.Random must be of length that corresponds to the number of elements to estimate of the B1 matrix")
      Beta[(np-nvar-ne_re+1):(np-nvar)] <- init.Random
    }
    if (!missing(init.Eta)) {
      if (length(init.Eta)!=netadc) stop("init.Eta must be of length that corresponds to the dimension of the link function")
      Beta[(np -nvar-1-ne_re-netadc+1):(np-nvar-1-ne_re)] <- init.Eta
    }

    xSuT <- matrix(0,nrow=100,ncol=1)
    if (typeof==0){
      mt1 <- size1
    }else{
      mt1 <- 100
    }
    size2 <- mt1


if(link0 %in% c(2,3)){
# Current-level association structure:
# we need to evaluate the biomarker value at multiple time-points to approximate the cumulative hazard
# time-interactions are particularly tricky to handle (especially in case of non-linear time trend)
      # position of time and interactions for current-level association

  # if interaction terms are not included, get them from survival data:
  columnsT <- colnames(X_T)
for(i in 1:length(timevar)){
    interact <- 0
    columns <- names(X_L)

    if(timevar[i] %in% columns){
      if(length(grep(":", columns))>0){
        interactions <- grep(":", columns)
        timevariable <- grep(timevar[i], columns)
        interact <- intersect(interactions, timevariable)
      }
      }


    count=0
    if(i==1){
      count2=0
      positionVarTime <- NULL
      numInterac <- 0
    }
    if(interact!=0){ # continuous
    if(length(interact)==1){ # one time-interaction term in the model
      if(interact!=0){

      name_1 <- strsplit(as.character(columns[interact]),":")[[1]][1]
      name_2 <- strsplit(as.character(columns[interact]),":")[[1]][2]

      # time-interaction terms / current-level association
      # PositionVarT contains: position of variable and time-variable for interaction and position of interaction
      # (first continuous, and then binary part if two-part model)
      if(timevar[i]%in%columns & count<2){
        # save positions of interaction with time for current-level association
        positionVarTime[count2+1] <- ifelse(length(which(columns==ifelse(name_1==timevar[i], name_2, name_1)))==0,which(columnsT==ifelse(name_1==timevar[i], name_2, name_1))+100,which(columns==ifelse(name_1==timevar[i], name_2, name_1)))
        positionVarTime[count2+2] <- ifelse(length(which(columns==ifelse(name_1==timevar[i], name_1, name_2)))==0,which(columnsT==ifelse(name_1==timevar[i], name_1, name_2))+100,which(columns==ifelse(name_1==timevar[i], name_1, name_2)))
        positionVarTime[count2+3] <- interact
        positionVarTime[count2+4] <- ifelse(timevar[i]=="f1", 1, ifelse(timevar[i]=="f2", 2,0))
        count=count+1
        count2=count2+4
        numInterac=numInterac+count
      }
      }
    }else if(length(interact)>1){ #
      for(j in 1:length(interact)){

      name_1 <- strsplit(as.character(columns[interact[j]]),":")[[1]][1]
      name_2 <- strsplit(as.character(columns[interact[j]]),":")[[1]][2]


      if(timevar[i]%in%columns & count<2){
        # save positions of interaction with time for current-level association
        positionVarTime[count2+1] <- ifelse(length(which(columns==ifelse(name_1==timevar[i], name_2, name_1)))==0,which(columnsT==ifelse(name_1==timevar[i], name_2, name_1))+100,which(columns==ifelse(name_1==timevar[i], name_2, name_1)))
        positionVarTime[count2+2] <- ifelse(length(which(columns==ifelse(name_1==timevar[i], name_1, name_2)))==0,which(columnsT==ifelse(name_1==timevar[i], name_1, name_2))+100,which(columns==ifelse(name_1==timevar[i], name_1, name_2)))
        positionVarTime[count2+3] <- interact[j]
        positionVarTime[count2+4] <- ifelse(timevar[i]=="f1", 1, ifelse(timevar[i]=="f2", 2,0))
        count=count+1
        count2=count2+4
      }}
      numInterac=numInterac+count

    }}else{
      if(is.null(positionVarTime)){
    positionVarTime=0
    numInterac=0
      }

      if(length(grep(timevar[i],columns))!=0){ # no interaction, only time
        positionVarTime[count2+1] <- 0
        positionVarTime[count2+2] <- grep(timevar[i],columns)
        positionVarTime[count2+3] <- 0
        positionVarTime[count2+4] <- ifelse(timevar[i]=="f1", 1, ifelse(timevar[i]=="f2", 2,0))
        count=count+1
        count2=count2+4
        numInterac=numInterac+count

      }
    }
}
for(i in 1:length(timevar)){
  interactB<-0
  if(TwoPart) columnsB <- names(X_B)

  if(TwoPart){
    if(timevar[i] %in% columnsB){
      if(length(grep(":", columnsB))>0){
        interactionsB <- grep(":", columnsB)
        timevariableB <- grep(timevar[i], columnsB)
        interactB <- intersect(interactionsB, timevariableB)
      }
    }}

    count=0
    if(TwoPart){
        if(interactB!=0){  # binary
if(i==1){
    numInteracB <- 0
}
        if(length(interactB)==1){
          if(interactB!=0){

      name_1B <- strsplit(as.character(columnsB[interactB]),":")[[1]][1]
      name_2B <- strsplit(as.character(columnsB[interactB]),":")[[1]][2]

      # interaction terms / current-level association
      # PositionVarT contains: position of variable and time-variable for interaction and position of interaction (first continuous, and then binary)
      if(timevar[i]%in%columnsB & count<2){
        # save positions of interaction with time for current-level association
        positionVarTime[count2+1] <- ifelse(length(which(columnsB==ifelse(name_1B==timevar[i], name_2B, name_1B)))==0,which(columnsT==ifelse(name_1B==timevar[i], name_2B, name_1B))+100,which(columnsB==ifelse(name_1B==timevar[i], name_2B, name_1B)))
        positionVarTime[count2+2] <- ifelse(length(which(columnsB==ifelse(name_1B==timevar[i], name_1B, name_2B)))==0,which(columnsT==ifelse(name_1B==timevar[i], name_1B, name_2B))+100,which(columnsB==ifelse(name_1B==timevar[i], name_1B, name_2B)))
        positionVarTime[count2+3] <- interactB
        positionVarTime[count2+4] <- ifelse(timevar[i]=="f1", 1, ifelse(timevar[i]=="f2", 2,0))
        count=count+1
        count2=count2+4
        numInteracB=numInteracB+count
      }
          }
    }else if (length(interactB)>1){
      for(j in 1:length(interactB)){

      name_1B <- strsplit(as.character(columnsB[interactB[j]]),":")[[1]][1]
      name_2B <- strsplit(as.character(columnsB[interactB[j]]),":")[[1]][2]

      # interaction terms / current-level association
      if(timevar[i]%in%columnsB & count<2){
        # save positions of interaction with time for current-level association
        positionVarTime[count2+1] <- ifelse(length(which(columnsB==ifelse(name_1B==timevar[i], name_2B, name_1B)))==0,which(columnsT==ifelse(name_1B==timevar[i], name_2B, name_1B))+100,which(columnsB==ifelse(name_1B==timevar[i], name_2B, name_1B)))
        positionVarTime[count2+2] <- ifelse(length(which(columnsB==ifelse(name_1B==timevar[i], name_1B, name_2B)))==0,which(columnsT==ifelse(name_1B==timevar[i], name_1B, name_2B))+100,which(columnsB==ifelse(name_1B==timevar[i], name_1B, name_2B)))
        positionVarTime[count2+3] <- interactB[j]
        positionVarTime[count2+4] <- ifelse(timevar[i]=="f1", 1, ifelse(timevar[i]=="f2", 2,0))
        count=count+1
        count2=count2+4
      }}
      numInteracB=numInteracB+count

    }
        }else{
          if(!exists("numInteracB")){
            numInteracB=0
          }
          if(length(grep(timevar[i],columnsB))!=0){ # no interaction, only time
            positionVarTime[count2+1] <- 0
            positionVarTime[count2+2] <- grep(timevar[i],columnsB)
            positionVarTime[count2+3] <- 0
            positionVarTime[count2+4] <- ifelse(timevar[i]=="f1", 1, ifelse(timevar[i]=="f2", 2,0))
            count=count+1
            count2=count2+4
            numInteracB=numInteracB+count

          }
        }
    }else{
    numInteracB=0
    }
}
}else{
  numInteracB=0
  numInterac=1
  positionVarTime=c(404,0,0,0)
}

seed.MC=ifelse(seed.MC==F,0,seed.MC) # seed for Monte-carlo (0 if random / unspecified)

    if(GLMlog==F) GLMloglink=0 else GLMloglink=1
    if(MTP==F) Mtwopart=0 else Mtwopart=1

  if(!TwoPart){ # initialize TwoPart variables if not activated to avoid memory allocation problems
    Binary <- rep(0, length(nsujety))
    nsujetB=0
    clusterB <- 0
    matzB <- matrix(as.double(0),nrow=1,ncol=1)
    nvarB <- 0
    varB <- matrix(as.double(0),nrow=1,ncol=1)
    nREB <- 0
    noVarB <- 1
    numInteracB=0
  }
    flush.console()
    if (print.times){
      ptm<-proc.time()
      cat("\n")
      cat("Be patient. The program is computing ... \n")
    }

    # call joint_longi.f90
        ans <- .Fortran(C_joint_longi,
			VectNsujet = as.integer(c(1,nsujety, nsujetB)),

            ngnzag=as.integer(c(ng, n.knots, 1, seed.MC)),

			k0 = as.double(c(0,kappa)), # joint avec generalisation de strate
			tt00 = as.double(0),
			tt10 = as.double(0),
		    ic0 = as.integer(0),
		    groupe0 = as.integer(0),
			tt0dc0 = as.double(tt0dc),
			tt1dc0 = as.double(tt1dc),
			icdc0 = as.integer(cens),
		    link0 = as.integer(c(link0,0)),
		    yy0 = as.double(Y),
            bb0 = as.double(Binary),
		    groupey0 = as.integer(clusterY),
            groupeB0 = as.integer(clusterB),
		    Vectnb0 = as.integer(c(nRE, nREB)),
		    matzy0 =as.double(matzy),
            matzB0 =as.double(matzB),
		    cag0 = as.double(cag),
            VectNvar=as.integer(c(1, nvarT, nvarY, nvarB)),
            vax0 = matrix(as.double(0),nrow=1,ncol=1),

			vaxdc0 = as.double(varT),
			vaxy0 = as.double(varY),
            vaxB0 = as.double(varB),
			noVar = as.integer(c(0,noVarT,noVarY, noVarB)),

			maxit0 = as.integer(maxit),
			np=as.integer(np),
			neta0 = as.integer(c(netadc,0)),
			b=as.double(Beta),
			H=as.double(matrix(0,nrow=np,ncol=np)),
			HIH=as.double(matrix(0,nrow=np,ncol=np)),

			loglik=as.double(0),
			LCV=as.double(rep(0,2)),
			xR=as.double(matrix(0,nrow=1,ncol=1)),
			lamR=as.double(matrix(0,nrow=1,ncol=3)),
			xSuR=as.double(array(0,dim=100)),
			survR=as.double(array(0,dim=1)),
			xD=as.double(rep(0,100)),
			lamD=as.double(matrix(0,nrow=size1,ncol=3)),
			xSuD=as.double(xSuT),
			survD=as.double(matrix(0,nrow=size2,ncol=3)),
			as.integer(typeof),
			as.integer(equidistant),
			as.integer(c(1,size1,1,mt1)),###
			counts=as.integer(c(0,0,0)),
			ier_istop=as.integer(c(0,0)),
			paraweib=as.double(rep(0,4)),
			MartinGale=as.double(matrix(0,nrow=ng,ncol=3+nRE)),###
			ResLongi = as.double(matrix(0,nrow=nsujety,ncol=4)),
			Pred_y  = as.double(matrix(0,nrow=nsujety,ncol=2)),
            GLMlog0 = as.integer(c(GLMloglink,Mtwopart)), # glm with log link + marginal two-part

			positionVarTime = as.integer(positionVarTime),
			numInterac = as.integer(c(numInterac, numInteracB)),

			linear.pred=as.double(rep(0,ng)),
			lineardc.pred=as.double(rep(0,as.integer(ng))),
			zi=as.double(rep(0,(n.knots+6))),

			paratps=as.integer(c(0,0,0)),#for future developments
			as.integer(c(0,0,0)),#for future developments
			BetaTpsMat=as.double(matrix(0,nrow=101,ncol=1+4*0)), #for future developments
			BetaTpsMatDc=as.double(matrix(0,nrow=101,ncol=1+4*0)),#for future developments
			BetaTpsMatY = as.double(matrix(0,nrow=101,ncol=1+4*0)),#for future developments
			EPS=as.double(c(LIMparam,LIMlogl,LIMderiv)),
			GH = c(as.integer(GH),as.integer(n.nodes)),
			paGH = data.matrix(cbind(b_lme,invBi_cholDet,as.data.frame(invBi_chol)))
			)#,
    #PACKAGE = "frailtypack") #62 arguments

    MartinGale <- matrix(ans$MartinGale,nrow=ng,ncol=3+nRE)
    Residuals <- matrix(ans$ResLongi,nrow=nsujety,ncol=4)

    if (ans$ier_istop[2] == 4){
      warning("Problem in the loglikelihood computation. The program stopped abnormally. Please verify your dataset. \n")
    }

    if (ans$ier_istop[2] == 2){
      warning("Model did not converge.")
    }
    if (ans$ier_istop[2] == 3){
      warning("Matrix non-positive definite.")
    }


    #AD:
    if (noVarT==1 & noVarY==1) nvar<-0
    #AD:

    np <- ans$np
    fit <- NULL
    fit$b <- ans$b
    fit$na.action <- attr(m, "na.action")
    fit$call <- call
    fit$groups <- ng



    fit$n.deaths <- ans$counts[3]
    fit$n.measurements <- nsujety
    fit$n.measurementsB <- nsujetB # add TwoPart

    if(as.character(typeof)=="0"){
      fit$logLikPenal <- ans$loglik
    }else{
      fit$logLik <- ans$loglik
    }
    #AD:
    fit$LCV <- ans$LCV[1]
    fit$AIC <- ans$LCV[2]
    #random effects
    fit$B1 <- matrix(rep(0,nRE^2),ncol=nRE,nrow=nRE)
    Ut <- matrix(rep(0,nRE^2),ncol=nRE,nrow=nRE)
    Utt <- matrix(rep(0,nRE^2),ncol=nRE,nrow=nRE)
    for(j in 1:nRE){
      for(k in 1:j){
        Ut[j,k]=ans$b[np - nvar - ne_re + k + j*(j-1)/2]
        Utt[k,j]=ans$b[np -nvar - ne_re + k + j*(j-1)/2]
      }}
    fit$B1 <- Ut%*%Utt

    fit$ResidualSE <- sqrt(ans$b[(np  - nvar - ne_re )]^2)
    fit$eta <- ans$b[(np  - nvar - 1 - ne_re - netadc + 1):(np  - nvar - 1 -ne_re)]

    fit$npar <- np

    #AD:
    if ((noVarT==1 & noVarY==1)) {
      fit$coef <- NULL
    }
    else
    {
      fit$coef <- ans$b[(np - nvar + 1):np]
      noms <- c(factor.names(colnames(X_T)),factor.names(colnames(X_L)))
      if(TwoPart) noms <-  c(factor.names(colnames(X_T)),factor.names(colnames(X_L)),factor.names(colnames(X_B))) # add TwoPart
      #  if (timedep == 1){
      #          while (length(grep("timedep",noms))!=0){
      #                  pos <- grep("timedep",noms)[1]
      #                  noms <- noms[-pos]
      #                  fit$coef <- fit$coef[-(pos:(pos-1))]
      #          }
      # }
      names(fit$coef) <- noms


    }

  if(TwoPart==1){
  fit$names.re <- c(names.matzy, names.matzB) # add TwoPart
  }else{
  fit$names.re <- names.matzy
  }


    temp1 <- matrix(ans$H, nrow = np, ncol = np)
    temp2 <- matrix(ans$HIH, nrow = np, ncol = np)

    varH.eta <- temp1[(np  - nvar - 1 - ne_re - netadc +1):(np  - nvar - 1 - ne_re),
                      (np - nvar - 1 - ne_re - netadc +1):(np  - nvar - 1 - ne_re )]

    if(netadc>1)fit$se.eta <- sqrt(diag(varH.eta))
    if(netadc==1)fit$se.eta <- sqrt(varH.eta)

    fit$eta_p.value <- 1 - pchisq((fit$eta/fit$se.eta)^2,1)


    fit$se.ResidualSE <- sqrt(temp1[(np  - nvar - ne_re ),(np  - nvar - ne_re)])
    fit$varHtotal <- temp1
    fit$varHIHtotal <- temp2

    fit$varH <- temp1[(np  - nvar +1):np, (np  - nvar +1 ):np]
    fit$varHIH <- temp2[(np - nvar +1):np, (np - nvar +1):np]

    if(TwoPart==1){
  noms <- c("MeasurementError","B1",factor.names(colnames(X_T)),factor.names(colnames(X_L)),factor.names(colnames(X_B))) # add TwoPart
  }else{
  noms <- c("MeasurementError","B1",factor.names(colnames(X_T)),factor.names(colnames(X_L)))
  }

    #     if (timedep == 1){ # on enleve les variances des parametres des B-splines
    #              while (length(grep("timedep",noms))!=0){
    #                      pos <- grep("timedep",noms)[1]
    #                      noms <- noms[-pos]
    #                      fit$varH <- fit$varH[-(pos:(pos-1)),-(pos:(pos-1))]
    #                      fit$varHIH <- fit$varHIH[-(pos:(pos-1)),-(pos:(pos-1))]
    #              }
    #      }
    fit$nvar<-c(nvarT,nvarY, nvarB) # modif TwoPart
    fit$formula <- formula #formula(Terms)
    fit$formula.LongitudinalData <- formula.LongitudinalData #formula(TermsY)
    fit$formula.Binary <- formula.Binary # add TwoPart
    fit$xD <- matrix(ans$xD, nrow = size1, ncol = 1)

    fit$lamD <- array(ans$lamD, dim = c(size1,3,1))
    fit$xSuD <- matrix(ans$xSuD, nrow = 100, ncol = 1)
    fit$survD <- array(ans$survD, dim = c(size2,3,1))

    fit$link <- link
    fit$type <- type
    fit$n.strat <- 1
    fit$n.iter <- ans$counts[1]
    fit$typeof <- typeof
    if (typeof == 0){
      fit$n.knots<-n.knots
      fit$kappa <- ans$k0[2]
      fit$n.knots.temp <- n.knots.temp
      fit$zi <- ans$zi
    }

    #AD:

    median <- NULL
    for (i in (1:fit$n.strat)) median[i] <- ifelse(typeof==0, minmin(fit$survD[,1,i],fit$xD), minmin(fit$survD[,1,i],fit$xSuD))
    lower <- NULL
    for (i in (1:fit$n.strat)) lower[i] <- ifelse(typeof==0, minmin(fit$survD[,2,i],fit$xD), minmin(fit$survD[,2,i],fit$xSuD))
    upper <- NULL
    for (i in (1:fit$n.strat)) upper[i] <- ifelse(typeof==0, minmin(fit$survD[,3,i],fit$xD), minmin(fit$survD[,3,i],fit$xSuD))
    fit$median <- cbind(lower,median,upper)

    fit$noVarEnd <- noVarT
    fit$noVarY <- noVarY
    fit$noVarB <- noVarB # add TwoPart


    fit$nvarEnd <- nvarT
    fit$nvarY <- nvarY
    fit$nvarB <- nvarB # add TwoPart
    fit$istop <- ans$ier_istop[2]

    fit$shape.weib <- ans$paraweib[2]#ans$shape.weib
    fit$scale.weib <- ans$paraweib[4]#ans$scale.weib

    if(ans$cag0[1]==1)fit$leftCensoring <- TRUE
    if(ans$cag0[1]==0)fit$leftCensoring <- FALSE

    if(fit$leftCensoring){fit$leftCensoring.threshold <-ans$cag0[2]
    fit$prop.censored <- prop.censored}
    #AD:

    # verif que les martingales ont ete bien calculees
    msg <- "Problem in the estimation of the random effects (perhaps high number of events in some clusters)"

    #       if (any(MartinGale[,2]==0)){

    #               fit$martingaledeath.res <- msg


    #               fit$linear.pred <- msg

    #       }else{

    #fit$martingale.res <- MartinGale[,1]#ans$martingale.res
    fit$martingaledeath.res <- MartinGale[,2]#ans$martingaledc.res
    fit$conditional.res <- Residuals[,1]
    fit$marginal.res <- Residuals[,3]
    fit$marginal_chol.res <- Residuals[,4]

    fit$conditional_st.res <- Residuals[,2]
    fit$marginal_st.res <- Residuals[,3]/fit$ResidualSE


    fit$random.effects.pred <- MartinGale[,3:(3+nRE-1)]#ans$frailty.pred
    fit$pred.y.marg <- matrix(ans$Pred_y,ncol=2)[,2]
    fit$pred.y.cond <- matrix(ans$Pred_y,ncol=2)[,1]
    fit$lineardeath.pred <- ans$lineardc.pred
    #       }



    #    if (joint.clust==0){
    #        fit$kendall <- matrix(ans$kendall,nrow=4,ncol=2)
    #    }

    #  fit$BetaTpsMat <- matrix(ans$BetaTpsMat,nrow=101,ncol=1+4*nvartimedepT)
    #  fit$BetaTpsMatDc <- matrix(ans$BetaTpsMatDc,nrow=101,ncol=1+4*nvartimedepY)
    #  fit$nvartimedep <- c(nvartimedepT,nvartimedepY)

    #  fit$Names.vardepT <- vardepT
    #  fit$Names.vardepY <- vardepY

    fit$EPS <- ans$EPS

    fit$ne_re <- nRE
    fit$netadc<-netadc



#================================> For the Binary
#========================= Test de Wald
if(TwoPart){
        if ((length(vec.factorB) > 0)){
                Beta <- ans$b[(np-nvar + 1):np]
                VarBeta <- fit$varH
                nfactorB <- length(vec.factorB)
                p.waldB <- rep(0,nfactorB)

                if(fit$istop == 1) fit$global_chisq_B <- waldtest(N=nvarB,nfact=nfactorB,place=ind.placeB,
                modality=occurB,b=Beta,Varb=VarBeta,Lfirts=(nvarT+nvarY),Ntot=nvar)# modif TwoPart
                else fit$global_chisq_B <- 0

                fit$dof_chisq_B <- occurB
                fit$global_chisq.test_B <- 1
# Calcul de pvalue globale
                for(i in 1:length(vec.factorB)){
                        p.waldB[i] <- signif(1 - pchisq(fit$global_chisq_B[i], occurB[i]), 3)
                }
                fit$p.global_chisq_B <- p.waldB
                fit$names.factor_B <- vec.factorB
        }else{
                fit$global_chisq.test_B <- 0

        }
}
    #================================> For the longitudinal
    #========================= Test de Wald

    if ((length(vec.factorY) > 0)){
      Beta <- ans$b[(np-nvar + 1):np]
      VarBeta <- fit$varH
      nfactor <- length(vec.factorY)
      p.wald <- rep(0,nfactor)

    if(fit$istop == 1) fit$global_chisq <- waldtest(N=nvarY,nfact=nfactor,place=ind.placeY,
    modality=occurY,b=Beta,Varb=VarBeta,Lfirts=nvarT, Llast=nvarB,Ntot=nvar)# modif TwoPart
    else fit$global_chisq <- 0

      fit$dof_chisq <- occurY
      fit$global_chisq.test <- 1
      # Calcul de pvalue globale
      for(i in 1:length(vec.factorY)){
        p.wald[i] <- signif(1 - pchisq(fit$global_chisq[i], occurY[i]), 3)
      }
      fit$p.global_chisq <- p.wald
      fit$names.factor <- vec.factorY
    }else{
      fit$global_chisq.test <- 0

    }

    #================================> For the death
    #========================= Test de Wald


    if ((length(vec.factorT) > 0) ){
      Beta <- ans$b[(np - nvar + 1):(np)]

      # if(npbetatps>0){VarBeta <- diag(diag(fit$varH)[-c((np-npbetatps):np)])
      # }else{
      VarBeta <- fit$varH
      #}
      nfactor <- length(vec.factorT)
      p.waldT <- rep(0,nfactor)
      fit$global_chisq_d <- waldtest(N=nvarT,nfact=nfactor,place=ind.placeT,
      modality=occurT,b=Beta,Varb=VarBeta,Llast=(nvarY+nvarB),Ntot=nvar) # modif TwoPart
      fit$dof_chisq_d <- occurT
      fit$global_chisq.test_d <- 1
      # Calcul de pvalue globale
      for(i in 1:length(vec.factorT)){
        p.waldT[i] <- signif(1 - pchisq(fit$global_chisq_d[i], occurT[i]), 3)
      }
      fit$p.global_chisq_d <- p.waldT
      fit$names.factordc <- vec.factorT
    }else{
      fit$global_chisq.test_d <- 0
    }

    if (!is.null(fit$coef)){
      if(nvar != 1){
        seH <- sqrt(diag(fit$varH))
      }else{
        seH <- sqrt(fit$varH)
      }
      fit$beta_p.value <- 1 - pchisq((fit$coef/seH)^2, 1)
    }

    if(intercept)fit$intercept <- TRUE
    else fit$intercept <- FALSE
    fit$Frailty <- FALSE
    fit$max_rep <- max_rep
    fit$joint.clust <- 1
    fit$methodGH <- method.GH
    fit$n.nodes <- n.nodes
    fit$TwoPart <- TwoPart # add TwoPart

    if(GLMloglink==1) fit$GLMlog <- TRUE else fit$GLMlog <- FALSE
     if(MTP==1) fit$MTP <- TRUE else fit$MTP <- FALSE
   class(fit) <- "longiPenal"
    if (print.times){
      cost<-proc.time()-ptm
      cat("The program took", round(cost[3],2), "seconds \n")
    }
    fit

  }

