#' Fit a Non-Linear Trivariate Joint Model for Recurrent Events and a Terminal
#' Event with a Biomarker Described with an ODE Population Model
#' 
#' @description{
#' \if{html}{Fit a non-linear trivariate joint model for a longitudinal biomarker,
#' recurrent events and a terminal event using a semiparametric penalized
#' likelihood estimation or a parametric estimation on the hazard functions.
#' 
#' The values y\out{<sub>i</sub>}(t) (i=1,...,N) for N subjects represent
#' the individual evolution of the biomarker e.g. tumor size expressed as the
#' sum of the longest diameters (SLD) of target lesions. The dynamics of the
#' biomarker are described by an ordinary differential equation (ODE) that
#' includes the effect of the natural net growth and the treatment effect:
#' 
#' {\figure{trivNLmodel1.png}{options: width="100\%"}}
#' 
#' The model includes the following parameters (using the interpretation of
#' tumor dynamics): exp(\eqn{K}\out{<sub>G,0</sub>}) the constant tumor growth rate,
#' exp(\eqn{K}\out{<sub>D,0</sub>}) the drug-induced tumor decline rate, \eqn{\lambda}
#' resistance effect to drug (exponential tumor decay change with time),
#' exp(y\out{<sub>0</sub>}) the initial level of the biomarker and d\out{<sub>i</sub>} is the
#' treatment concentration (e.g dose). The random effects \bold{b}\out{<sub>i</sub>}\out{<sup>T</sup>}
#' = (b\out{<sub>y0,i</sub>},b\out{<sub>G,i</sub>},b\out{<sub>D,i</sub>},b\out{<sub>\lambda,i</sub>})\out{<sup>T</sup>} are gaussian variables
#' with a diagonal covariance matrix \bold{B}\out{<sub>i</sub>}. In the trivariate model
#' we use the analytical solution of the equation with the population-based
#' approach of the non-linear mixed effects model. We can also assume a
#' transformation for the observations of the biomarker (one parameter Box-Cox
#' transformation) and we include a gaussian measurement error, for individual
#' i and observation k (k=1,...,n\out{<sub>i</sub>}),
#' \eqn{\epsilon}\out{<sub>ik</sub>} \out{&#126;} \bold{\eqn{N}}(0,\eqn{\sigma}\out{<sub>\epsilon</sub>}\out{<sup>2</sup>}).
#' 
#' The risks of the recurrent (r\out{<sub>ij</sub>}(.) the risk of the j\out{<sup>th</sup>}
#' event of the individual i) and terminal events (\eqn{\lambda}\out{<sub>i</sub>} the
#' risk of the event of the individual i) are represented by proportional
#' hazard risk models. The joint model is constructed assuming that the
#' processes are linked via a latent structure and includes the non-linear
#' mixed effects model for the longitudinal data:
#' 
#' {\figure{trivNLmodel2.png}{options: width="100\%"}}
#' 
#' where \bold{\eqn{X}}\out{<sub>G,i</sub>}(t), \bold{\eqn{X}}\out{<sub>D,i</sub>}(t),
#' \bold{\eqn{X}}\out{<sub>R,ij</sub>}(t) and \bold{\eqn{X}}\out{<sub>T,i</sub>}(t) are vectors of possible
#' time-varying fixed effects covariates and \bold{\eqn{\beta}}\out{<sub>G</sub>},
#' \bold{\eqn{\beta}}\out{<sub>D</sub>}, \bold{\eqn{\beta}}\out{<sub>R</sub>} and \bold{\eqn{\beta}}\out{<sub>T</sub>} are the
#' associated coefficients. The random effects \bold{b}\out{<sub>i</sub>} are independent
#' from the measurement error. The relationship between the biomarker and
#' recurrent events is explained via g(y\out{<sub>i</sub>}(t)) with coefficients
#' \bold{\eqn{\eta}}\out{<sub>R</sub>} and between the biomarker and terminal event is
#' explained via h(y\out{<sub>i</sub>}(t)) with coefficients \bold{\eqn{\eta}}\out{<sub>T</sub>}.
#' Currently, only one form of the functions g(.) and h(.)
#' is available: the random effects \bold{b}\out{<sub>i</sub>}. The frailty term
#' v\out{<sub>i</sub>} is gaussian with mean 0 and variance \eqn{\sigma}\out{<sub>v</sub>}. Together with
#' \bold{b}\out{<sub>i</sub>} constitutes the random effects of the model: 
#' 
#' {\figure{trivNLmodel3.png}{options: width="100\%"}}
#' 
#' Any combination of the random effects \bold{b}\out{<sub>i</sub>}, e.g.
#' \bold{b}\out{<sub>i</sub>}=b\out{<sub>y0,i</sub>} or \bold{b}\out{<sub>i</sub>} =
#' \{b\out{<sub>G,i</sub>},b\out{<sub>D,i</sub>},b\out{<sub>\lambda,i</sub>}\} can be chosen for the model.
#' 
#' We consider that the longitudinal outcome can be a subject to a
#' quantification limit, i.e. some observations, below a level of detection
#' s cannot be quantified (left-censoring).
#' }
#' \if{latex}{Fit a non-linear trivariate joint model for a longitudinal biomarker,
#' recurrent events and a terminal event using a semiparametric penalized
#' likelihood estimation or a parametric estimation on the hazard functions.
#' 
#' The values \eqn{y_i(t)} (\eqn{i=1,\ldots,N}) for \eqn{N} subjects represent
#' the individual evolution of the biomarker e.g. tumor size expressed as the
#' sum of the longest diameters (SLD) of target lesions. The dynamics of the
#' biomarker are described by an ordinary differential equation (ODE) that
#' includes the effect of the natural net growth and the treatment effect:
#' 
#' \deqn{\left\{ \begin{array}{ll}
#' \frac{dy_{i}(t)}{dt}&=\exp(K_{G,0}+b_{G,i}+\bold{X}_{G,i}(t)^\top\bold{\beta}_G)y_{i}(t)\\
#' &-d_i\exp(K_{D,0}+b_{D,i}-t\times\exp(\lambda+b_{\lambda,i})+\bold{X}_{D,i}(t)^\top\bold{\beta}_D)y_{i}(t)\\
#' y_{i}(0)&=\exp(y_0+b_{y_0,i}) \\ \end{array}, \right. }
#' 
#' The model includes the following parameters (using the interpretation of
#' tumor dynamics): \eqn{\exp(K_{G,0})} the constant tumor growth rate,
#' \eqn{\exp(K_{D,0})} the drug-induced tumor decline rate, \eqn{\lambda}
#' resistance effect to drug (exponential tumor decay change with time),
#' \eqn{\exp(y_0)} the initial level of the biomarker and \eqn{d_i} is the
#' treatment concentration (e.g dose). The random effects \eqn{\bold{b}_i^\top
#' = (b_{y_0,i},b_{G,i},b_{D,i},b_{\lambda,i})^\top} are gaussian variables
#' with a diagonal covariance matrix \eqn{\bold{B}_1}. In the trivariate model
#' we use the analytical solution of the equation with the population-based
#' approach of the non-linear mixed effects model. We can also assume a
#' transformation for the observations of the biomarker (one parameter Box-Cox
#' transformation) and we include a gaussian measurement error, for individual
#' \eqn{i} and observation \eqn{k} (\eqn{k=1,\ldots,n_i}),
#' \eqn{\epsilon_{ik}\sim N(0,\sigma_{\epsilon}^2)}.
#' 
#' The risks of the recurrent (\eqn{r_{ij}(\cdot)} the risk of the \eqn{j}-th
#' event of the individual \eqn{i}) and terminal events (\eqn{\lambda_{i}} the
#' risk of the event of the individual \eqn{i}) are represented by proportional
#' hazard risk models. The joint model is constructed assuming that the
#' processes are linked via a latent structure and includes the non-linear
#' mixed effects model for the longitudinal data:
#' 
#' \deqn{\left\{ \begin{array}{ll}
#' y(t_{ik})&=\exp[y_0+b_{y_0,i}+t_{ik}\times\exp(K_{G,0}+b_{G,i}+\bold{X}_{G,i}(t)^\top\bold\beta_{G})\\
#' &+d_i\times
#' \exp(K_{D,0}+b_{D,i}-\lambda-b_{\lambda,i}+\bold{X}_{D,i}(t)^\top\bold{\beta}_{D})\\
#' &\times(\exp(-\exp(\lambda+b_{\lambda,i})t_{ik})-1)]+\epsilon_{ik}\\
#' r_{ij}(t|\bold{b}_i)&=r_0(t)\exp(v_i+\bold{X}_{Rij}(t)^\top\bold{\beta}_R+g(y_i(t))^\top
#' \bold{\eta}_R ) \\ \lambda_i(t|\bold{b}_i)&=\lambda_0(t)\exp(\alpha
#' v_i+\bold{X}_{Ti}(t)^\top\bold{\beta}_T+h(y_i(t))^\top \bold{\eta}_T ) \\
#' \end{array} \right. }
#' 
#' where \eqn{\bold{X}_{G,i}(t)}, \eqn{\bold{X}_{D,i}(t)},
#' \eqn{\bold{X}_{R,ij}(t)} and \eqn{\bold{X}_{T,i}(t)} are vectors of possible
#' time-varying fixed effects covariates and \eqn{\bold{\beta}_G},
#' \eqn{\bold{\beta}_D}, \eqn{\bold{\beta}_R} and \eqn{\bold{\beta}_T} are the
#' associated coefficients. The random effects \eqn{\bold{b}_i} are independent
#' from the measurement error. The relationship between the biomarker and
#' recurrent events is explained via \eqn{g(y_i(t))} with coefficients
#' \eqn{\bold{\eta}_R} and between the biomarker and terminal event is
#' explained via \eqn{h(y_i(t))} with coefficients \eqn{\bold{\eta}_T}.
#' Currently, only one form of the functions \eqn{g(\cdot)} and \eqn{h(\cdot)}
#' is available: the random effects \eqn{\bold{b}_i}. The frailty term
#' \eqn{v_i} is gaussian with mean 0 and variance \eqn{\sigma_v}. Together with
#' \eqn{\bold{b}_i} constitutes the random effects of the model: \deqn{
#' \bold{u}_i=\left(\begin{array}{c} \bold{b}_{i}\\v_i \\ \end{array}\right)
#' \sim \mathcal{N}\left(\bold{0}, \left(\begin{array} {cc} \bold{B}_1&\bold{0}
#' \\ \bold{0} & \sigma_v^{2}\\\end{array}\right)\right), }
#' 
#' Any combination of the random effects \eqn{\bold{b}_i}, e.g.
#' \eqn{\bold{b}_i=b_{y_0,i}} or \eqn{\bold{b}_i =
#' \{b_{G,i},b_{D,i},b_{\lambda,i}\}} can be chosen for the model.
#' 
#' We consider that the longitudinal outcome can be a subject to a
#' quantification limit, i.e. some observations, below a level of detection
#' \eqn{s} cannot be quantified (left-censoring).}
#' }
#' @details{
#' 
#' Typical usage for the joint model
#' \preformatted{trivPenalNL(Surv(time,event)~cluster(id) + var1 + var2 +
#' terminal(death), formula.terminalEvent =~ var1 + var3, biomarker =
#' "biomarker.name", dose = "dose.name", time.biomarker = "time", formula.KG ~
#' var1, formula.KD ~ var2, data, data.Longi, ...)}
#' 
#' The method of the Gauss-Hermite quadrature for approximations of the
#' multidimensional integrals, i.e. length of \code{random} more than 2, can be
#' chosen among the standard (non-adaptive) and pseudo-adaptive in which the
#' quadrature points are transformed using the information from the fitted
#' mixed-effects model for the biomarker (Rizopoulos 2012) or multivariate
#' non-adaptive procedure proposed by Genz et al. 1996 and implemented in
#' FORTRAN subroutine HRMSYM.  The choice of the method is important for
#' estimations. The standard non-adaptive Gauss-Hermite quadrature
#' (\code{"Standard"}) with a specific number of points gives accurate results
#' but can be time consuming.  The pseudo-adaptive quadrature uses transformed
#' quadrature points to center and scale the integrand by utilizing estimates
#' of the random effects from an appropriate non-linear mixed-effects model
#' (this transformation does not include the frailty in the trivariate model,
#' for which the standard method, with 20 quadrature points, is used). This
#' method enables using less quadrature points while preserving the estimation
#' accuracy and thus lead to a better computational time.
#' 
#' NOTE. Data frames \code{data} and \code{data.Longi} must be consistent.
#' Names and types of corresponding covariates must be the same, as well as the
#' number and identification of individuals.
#' }
#' 
#' @usage
#' 
#' trivPenalNL(formula, formula.terminalEvent, biomarker, formula.KG,
#' formula.KD, dose, time.biomarker, data, data.Longi, random, id, link =
#' "Random-effects", BoxCox = FALSE, left.censoring = FALSE, recurrentAG =
#' FALSE, n.knots, kappa, maxit = 300, hazard = "Splines", init.B, init.Random,
#' init.Eta, init.Alpha, init.Biomarker, method.GH = "Standard", init.GH =
#' FALSE, n.nodes, LIMparam = 1e-3, LIMlogl = 1e-3, LIMderiv = 1e-3,
#' print.times = TRUE)
#' @param formula a formula object, with the response on the left of a
#' \eqn{\sim} operator, and the terms on the right. The response must be a
#' survival object as returned by the 'Surv' function like in survival package.
#' Interactions are possible using * or :.
#' @param formula.terminalEvent A formula object, only requires terms on the
#' right to indicate which variables are modelling the terminal event.
#' Interactions are possible using * or :.
#' @param biomarker Name of the variable representing the longitudinal
#' biomarker.
#' @param formula.KG A formula object, only requires terms on the right to
#' indicate which covariates related to the biomarker growth are included in
#' the longitudinal sub-model.  It must follow the standard form used for
#' linear mixed-effects models. Interactions are possible using * or :.
#' @param formula.KD A formula object, only requires terms on the right to
#' indicate which covariates related to the biomarker drug-induced decline are
#' included in the longitudinal sub-model.  It must follow the standard form
#' used for linear mixed-effects models. Interactions are possible using * or
#' :.
#' @param dose Name of the variable representing the drug concentration
#' indicator.
#' @param time.biomarker Name of the variable of times of biomarker
#' measurements.
#' @param data A 'data.frame' with the variables used in \code{formula}.
#' @param data.Longi A 'data.frame' with the variables used in
#' \code{formula.KG}, \code{formula.KD}, \code{biomarker}, \code{dose},
#' \code{time.biomarker} and \code{id}.
#' @param random Names of parameters for which the random effects are included
#' in the mixed model. The names must be chosen among \code{"y0"}, \code{"KG"},
#' \code{"KD"} and \code{"lambda"}. Any combination of the mentioned names is
#' allowed.
#' @param id Name of the variable representing the individuals.
#' @param link Type of link functions for the dependence between the biomarker
#' and death and between the biomarker and the recurrent events: only
#' \code{"Random-effects"} for the association directly via the random effects
#' of the biomarker is allowed for the moment (option for a future extension).
#' @param BoxCox Should the Box-Cox transformation be used for the longitudinal
#' biomarker? If there is no transformation, the argument must be equal to
#' \code{FALSE}, otherwise the of the transformation parameter must be given,
#' then the transformed values are \eqn{y^*=(y^{\xi}-1)/\xi}, where \eqn{\xi}
#' is the Box-Cox parameter.
#' @param left.censoring Is the biomarker left-censored below a threshold
#' \eqn{s}? If there is no left-censoring, the argument must be equal to
#' \code{FALSE}, otherwise the value of the threshold must be given.
#' @param recurrentAG Logical value.  Is Andersen-Gill model fitted?  If so
#' indicates that recurrent event times with the counting process approach of
#' Andersen and Gill is used. This formulation can be used for dealing with
#' time-dependent covariates.  The default is FALSE.
#' @param n.knots Integer giving the number of knots to use. Value required in
#' the penalized likelihood estimation.  It corresponds to the (n.knots+2)
#' splines functions for the approximation of the hazard or the survival
#' functions.  We estimate I or M-splines of order 4. When the user set a
#' number of knots equals to k (n.knots=k) then the number of interior knots is
#' (k-2) and the number of splines is (k-2)+order.  Number of knots must be
#' between 4 and 20. (See Note in \code{frailtyPenal} function)
#' @param kappa Positive smoothing parameters in the penalized likelihood
#' estimation.  The coefficient kappa of the integral of the squared second
#' derivative of hazard function in the fit (penalized log likelihood). To
#' obtain an initial value for \code{kappa}, a solution is to fit the
#' corresponding Cox model using cross validation (See \code{cross.validation}
#' in function \code{frailtyPenal}).  We advise the user to identify several
#' possible tuning parameters, note their defaults and look at the sensitivity
#' of the results to varying them.
#' @param maxit Maximum number of iterations for the Marquardt algorithm.
#' Default is 300
#' @param hazard Type of hazard functions: \code{"Splines"} for semiparametric
#' hazard functions using equidistant intervals or \code{"Splines-per"} using
#' percentile with the penalized likelihood estimation, \code{"Weibull"} for
#' the parametric Weibull functions. The default is \code{"Splines"}.
#' @param init.B Vector of initial values for regression coefficients. This
#' vector should be of the same size as the whole vector of covariates with the
#' first elements for the covariates related to the recurrent events, then to
#' the terminal event and then to the biomarker (interactions in the end of
#' each component). Default is 0.5 for each.
#' @param init.Random Initial value for variance of the elements of the matrix
#' of the distribution of the random effects.
#' @param init.Eta Initial values for regression coefficients for the link
#' functions, first for the recurrent events (\eqn{\bold{\eta}_R}) and for the
#' terminal event (\eqn{\bold{\eta}_T}).
#' @param init.Alpha Initial value for parameter alpha
#' @param init.Biomarker Initial values for biomarker parameters: \eqn{y_0},
#' \eqn{K_{G,0}}, \eqn{K_{D,0}} and \eqn{\lambda} (using this order).
#' @param method.GH Method for the Gauss-Hermite quadrature: \code{"Standard"}
#' for the standard non-adaptive Gaussian quadrature and
#' \code{"Pseudo-adaptive"} for the pseudo-adaptive Gaussian quadrature (see
#' Details). The default is \code{"Standard"}.  When the option
#' \code{"Pseudo-adaptive"} is chosen, then a univariate model (non-linear
#' mixed model for the biomarker) is fitted in order to obtain the estimations
#' of the random effects \eqn{\bold{b}_i}.
#' @param init.GH Only when the opiton \code{"Pseudo-adaptive"} of the argument
#' \code{method.GH} is chosen. If \code{TRUE}, the estimations of the biomarker
#' parameters (\eqn{y_0}, \eqn{K_{G,0}}, \eqn{K_{D,0}} and \eqn{\lambda}),
#' \eqn{\sigma_{\epsilon}}, \eqn{\bold{\beta}_G} and \eqn{\bold{\beta}_D} from
#' the univariate mixed model are used as the initial values of the parameters
#' related to the biomarker.
#' @param n.nodes Number of nodes for the Gauss-Hermite quadrature (from 5 to
#' 32). The default is 9.
#' @param LIMparam Convergence threshold of the Marquardt algorithm for the
#' parameters (see Details), \eqn{10^{-3}} by default.
#' @param LIMlogl Convergence threshold of the Marquardt algorithm for the
#' log-likelihood (see Details), \eqn{10^{-3}} by default.
#' @param LIMderiv Convergence threshold of the Marquardt algorithm for the
#' gradient (see Details), \eqn{10^{-3}} by default.
#' @param print.times a logical parameter to print iteration process. Default
#' is TRUE.
#' @return
#' 
#' The following components are included in a 'trivPenalNL' object for each
#' model:
#' 
#' \item{b}{The sequence of the corresponding estimation of the coefficients
#' for the hazard functions (parametric or semiparametric), the random effects
#' variances and the regression coefficients.} \item{call}{The code used for
#' the model.} \item{formula}{The formula part of the code used for the
#' recurrent event part of the model.} \item{formula.terminalEvent}{The formula
#' part of the code used for the terminal event part of the model.}
#' \item{formula.KG}{The formula part of the code used for the longitudinal
#' part of the model, for the biomarker growth dynamics.} \item{formula.KD}{The
#' formula part of the code used for the longitudinal part of the model, for
#' the biomarker decline dynamics.} \item{coef}{The regression coefficients
#' (first for the recurrent events, then for the terminal event, then for the
#' biomarker growth and then for the biomarker decline.} \item{groups}{The
#' number of groups used in the fit.} \item{kappa}{The values of the smoothing
#' parameters in the penalized likelihood estimation corresponding to the
#' baseline hazard functions for the recurrent and terminal events.}
#' \item{logLikPenal}{The complete marginal penalized log-likelihood in the
#' semiparametric case.} \item{logLik}{The marginal log-likelihood in the
#' parametric case.} \item{n.measurements}{The number of biomarker observations
#' used in the fit.} \item{max_rep}{The maximal number of repeated measurements
#' per individual.} \item{n}{The number of observations in 'data' (recurrent
#' and terminal events) used in the fit.} \item{n.events}{The number of
#' recurrent events observed in the fit.} \item{n.deaths}{The number of
#' terminal events observed in the fit.} \item{n.iter}{The number of iterations
#' needed to converge.} \item{n.knots}{The number of knots for estimating the
#' baseline hazard function in the penalized likelihood estimation.}
#' \item{n.strat}{The number of stratum.}
#' 
#' \item{varH}{The variance matrix of all parameters (before positivity
#' constraint transformation for the variance of the measurement error, for
#' which the delta method is used).} \item{varHIH}{The robust estimation of the
#' variance matrix of all parameters.}
#' 
#' \item{xR}{The vector of times where both survival and hazard function of the
#' recurrent events are estimated. By default seq(0,max(time),length=99), where
#' time is the vector of survival times.} \item{lamR}{The array (dim=3) of
#' baseline hazard estimates and confidence bands (recurrent events).}
#' \item{survR}{The array (dim=3) of baseline survival estimates and confidence
#' bands (recurrent events).}
#' 
#' \item{xD}{The vector of times where both survival and hazard function of the
#' terminal event are estimated. By default seq(0,max(time),length=99), where
#' time is the vector of survival times.} \item{lamD}{The array (dim=3) of
#' baseline hazard estimates and confidence bands.} \item{survD}{The array
#' (dim=3) of baseline survival estimates and confidence bands.}
#' \item{medianR}{The value of the median survival and its confidence bands for the recurrent event.}
#' \item{medianD}{The value of the median survival and its confidence bands for the terminal event.}
#' \item{typeof}{The type of the baseline hazard function (0:"Splines",
#' "2:Weibull").} \item{npar}{The number of parameters.} \item{nvar}{The vector
#' of number of explanatory variables for the recurrent events, terminal event,
#' biomarker growth and biomarker decline.} \item{nvarRec}{The number of
#' explanatory variables for the recurrent events.} \item{nvarEnd}{The number
#' of explanatory variables for the terminal event.} \item{nvarKG}{The number
#' of explanatory variables for the biomarker growth.} \item{nvarKD}{The number
#' of explanatory variables for the biomarker decline.} \item{noVarRec}{The
#' indicator of absence of the explanatory variables for the recurrent events.}
#' \item{noVarEnd}{The indicator of absence of the explanatory variables for
#' the terminal event.} \item{noVarKG}{The indicator of absence of the
#' explanatory variables for the biomarker growth.} \item{noVarKD}{The
#' indicator of absence of the explanatory variables for the biomarker
#' decline.} \item{LCV}{The approximated likelihood cross-validation criterion
#' in the semiparametric case (with H minus the converged Hessian matrix, and
#' l(.) the full log-likelihood).\deqn{LCV=\frac{1}{n}(trace(H^{-1}_{pl}H) -
#' l(.))}} \item{AIC}{The Akaike information Criterion for the parametric
#' case.\deqn{AIC=\frac{1}{n}(np - l(.))}} \item{n.knots.temp}{The initial
#' value for the number of knots.} \item{shape.weib}{The shape parameter for
#' the Weibull hazard functions (the first element for the recurrences and the
#' second one for the terminal event).} \item{scale.weib}{The scale parameter
#' for the Weibull hazard functions (the first element for the recurrences and
#' the second one for the terminal event).}
#' 
#' \item{random.effects.pred}{ The empirical Bayes predictions of the random
#' effects (ie. using conditional posterior distributions).}
#' \item{global_chisq.testR}{The binary variable equals to 0 when no
#' multivariate Wald is given, 1 otherwise (for the recurrent part).}
#' \item{global_chisq.testT}{The binary variable equals to 0 when no
#' multivariate Wald is given, 1 otherwise (for the terminal part).}
#' \item{global_chisq.testKG}{The binary variable equals to 0 when no
#' multivariate Wald is given, 1 otherwise (for the biomarker growth).}
#' \item{global_chisq.testKD}{The binary variable equals to 0 when no
#' multivariate Wald is given, 1 otherwise (for the biomarker decline).}
#' 
#' \item{AG}{The logical value. Is Andersen-Gill model fitted? }
#' 
#' \item{B1}{The variance matrix of the random effects for the longitudinal
#' outcome.} \item{sigma2}{The variance of the frailty term (\eqn{\sigma_v}).}
#' \item{alpha}{The coefficient \eqn{\alpha} associated with the frailty
#' parameter in the terminal hazard function.} \item{ResidualSE}{The variance
#' of the measurement error.} \item{etaR}{The regression coefficients for the
#' link function \eqn{g(\cdot)}.} \item{etaT}{The regression coefficients for
#' the link function \eqn{h(\cdot)}.} \item{ne_re}{The number of random effects
#' b used in the fit.} \item{names.re}{The names of variables for the random
#' effects \eqn{\bold{b}_i}.} \item{link}{The name of the type of the link
#' functions.}
#' 
#' \item{leftCensoring}{The logical value. Is the longitudinal outcome
#' left-censored?} \item{leftCensoring.threshold}{For the left-censored
#' biomarker, the value of the left-censoring threshold used for the fit.}
#' \item{prop.censored}{The fraction of observations subjected to the
#' left-censoring.}
#' 
#' \item{methodGH}{The Gaussian quadrature method used in the fit.}
#' \item{n.nodes}{The number of nodes used for the Gaussian quadrature in the
#' fit.}
#' 
#' \item{K_G0}{Value of the estimate of the biomarker growth parameter.}
#' \item{K_D0}{Value of the estimate of the biomarker decay parameter.}
#' \item{lambda}{Value of the estimate of the biomarker resistance to drug.}
#' \item{y_0}{Value of the estimate of the biomarker intial level.}
#' 
#' \item{biomarker}{Name of the variable associated with the biomarker in the
#' data.} \item{time.biomarker}{Name of the variable associated with the time
#' of measurements of the biomarker in the data.} \item{dose}{Name of the
#' variable associated with the drug concentration in the data.}
#' 
#' \item{BoxCox}{The logical value. Is the BoxCox transformation applied for
#' the biomarker?} \item{BoxCox_parameter}{The value of the BoxCox
#' transformation parameter.}
#' 
#' \item{alpha_p.value}{p-value of the Wald test for the estimated coefficient
#' \eqn{\alpha}.} \item{sigma2_p.value}{p-value of the Wald test for the
#' estimated variance of the frailty term (\eqn{\sigma_v}).}
#' \item{etaR_p.value}{p-values of the Wald test for the estimated regression
#' coefficients for the link function \eqn{g(\cdot)}.}
#' \item{etaT_p.value}{p-values of the Wald test for the estimated regression
#' coefficients for the link function \eqn{h(\cdot)}.}
#' \item{y_0_p.value}{p-value of the Wald test for the estimated biomarker
#' intial level.} \item{K_G0_p.value}{p-value of the Wald test for the
#' estimated biomarker growth parameter.} \item{K_D0_p.value}{p-value of the
#' Wald test for the estimated biomarker decay parameter.}
#' \item{lambda_p.value}{p-value of the Wald test for the estimated biomarker
#' resistance to drug.} \item{beta_p.value}{p-values of the Wald test for the
#' estimated regression coefficients.}
#' @note It is recommended to initialize the parameter values using the results
#' from a corresponding reduced model (\code{frailtyPenal} for the recurrent
#' and terminal part). See example.
#' 
#' Estimations of models with more than three random effects can be very long.
#' @seealso
#' \code{\link{plot.trivPenalNL}},\code{\link{print.trivPenalNL}},\code{\link{summary.trivPenalNL}}
#' @references A. Krol, C. Tournigand, S. Michiels and V. RondeauS (2018).
#' Multivariate joint frailty model for the analysis of nonlinear tumor
#' kinetics and dynamic predictions of death. \emph{Statistics in Medicine}.
#' 
#' A. Krol, L. Ferrer, JP. Pignon, C. Proust-Lima, M. Ducreux, O. Bouche, S.
#' Michiels, V. Rondeau (2016). Joint Model for Left-Censored Longitudinal
#' Data, Recurrent Events and Terminal Event: Predictive Abilities of Tumor
#' Burden for Cancer Evolution with Application to the FFCD 2000-05 Trial.
#' \emph{Biometrics} \bold{72}(3) 907-16.
#' 
#' D. Rizopoulos (2012). Fast fitting of joint models for longitudinal and
#' event time data using a pseudo-adaptive Gaussian quadrature rule.
#' \emph{Computational Statistics and Data Analysis} \bold{56}, 491-501.
#' 
#' L. Claret, P. Girard, P.M. Hoff, E. Van Cutsem, K.P. Zuideveld, K. Jorga, J.
#' Fagerberg, R Bruno (2009). Model-based prediction of phase III overall
#' survival in colorectal cancer on the basis of phase II tumor dynamics.
#' \emph{Journal of Clinical Oncology} \bold{27}(25), 4103-8.
#' @keywords models
#' @export
#' @examples
#' 
#' 
#' \dontrun{
#' 
#' ###--- Non-linear trivariate joint model for longitudinal data, ---###
#' ###--- recurrent events and a terminal event ---###
#' 
#' data(colorectal)
#' data(colorectalLongi)
#' 
#' # No information on dose - creation of a dummy variable 
#' colorectalLongi$dose <- 1
#' 
#' 
#' # Parameters initialisation - estimation of a simplified model
#' # with two random effects (a frailty term and a random effect 
#' # related to biomarker growth (KG))
#' initial.model <- trivPenalNL(Surv(time0, time1, new.lesions) ~ cluster(id)
#'  + age + treatment + terminal(state), formula.terminalEvent =~ age + treatment, 
#'  biomarker = "tumor.size", formula.KG ~ 1, formula.KD ~ treatment, dose = "dose",
#'  time.biomarker = "year", data = colorectal, data.Longi =colorectalLongi, 
#'  random = "KG", id = "id", recurrentAG = TRUE, n.knots = 5, kappa = c(0.01, 2),
#'  method.GH = "Pseudo-adaptive")
#' 
#' 
#' # Trivariate joint model with initial values for parameters
#' # (computation takes around 40 minutes)
#' 
#' model <- trivPenalNL(Surv(time0, time1, new.lesions) ~ cluster(id) + age + treatment
#'  + terminal(state), formula.terminalEvent =~ age + treatment, biomarker = "tumor.size",
#'  formula.KG ~ 1, formula.KD ~ treatment, dose = "dose", time.biomarker = "year", 
#'  data = colorectal, data.Longi =colorectalLongi, random = c("y0", "KG"), id = "id", 
#'  init.B = c(-0.22, -0.16, -0.35, -0.19, 0.04, -0.41, 0.23), init.Alpha = 1.86,
#'  init.Eta = c(0.5, 0.57, 0.5, 2.34), init.Biomarker = c(1.24, 0.81, 1.07, -1.53),
#'  recurrentAG = TRUE, n.knots = 5, kappa = c(0.01, 2), method.GH = "Pseudo-adaptive")
#' 
#' }
#' 
#' 
"trivPenalNL" <-
  function (formula, formula.terminalEvent, biomarker, formula.KG, formula.KD, dose, time.biomarker, data,  data.Longi, random, id, 
            link="Random-effects", BoxCox = FALSE,
            left.censoring=FALSE, recurrentAG=FALSE, n.knots, kappa,
            maxit=300, hazard="Splines", init.B,
            init.Random, init.Eta, init.Alpha, init.Biomarker,
            method.GH = "Standard", init.GH = FALSE, n.nodes, LIMparam=1e-3, LIMlogl=1e-3, LIMderiv=1e-3, print.times=TRUE)
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
    
    m3 <- match.call() # longitudinal (KG)
    m3$formula <- m3$formula.terminalEvent <- m3$biomarker <- m3$formula.KD <- m3$dose <- m3$data <- m3$recurrentAG <- m3$random <- m3$id <- m3$link <- m3$n.knots <- m3$kappa <- m3$maxit <- m3$hazard <- m3$init.B <- m3$LIMparam <- m3$LIMlogl <- m3$LIMderiv <- m3$print.times <- m3$left.censoring <- m3$init.Random <- m3$init.Eta <- m3$init.Alpha <- m3$method.GH <- m3$n.nodes <- m3$time.biomarker <- m3$init.GH <- m3$BoxCox <- m3$init.Biomarker <- m3$... <- NULL
    Names.data.Longi <- m3$data.Longi
    
    m4 <- match.call() # longitudinal (KD)
    m4$formula <- m4$formula.terminalEvent <- m4$biomarker <- m4$formula.KG <- m4$dose <- m4$data <- m4$recurrentAG <- m4$random <- m4$id <- m4$link <- m4$n.knots <- m4$kappa <- m4$maxit <- m4$hazard <- m4$init.B <- m4$LIMparam <- m4$LIMlogl <- m4$LIMderiv <- m4$print.times <- m4$left.censoring <- m4$init.Random <- m4$init.Eta <- m4$init.Alpha <- m4$method.GH <- m4$n.nodes <- m4$time.biomarker <- m4$init.GH <- m4$BoxCox <-  m4$init.Biomarker <- m4$... <- NULL
    
    
    m2 <- match.call() #terminal
    m2$formula <- m2$formula.terminalEvent <- m2$formula.LongitudinalData <- m2$formula.KG <- m2$formula.KD <- m2$biomarker <- m2$dose <- m2$data.Longi <- m2$recurrentAG <- m2$random <- m2$id <- m2$link <- m2$n.knots <- m2$kappa <- m2$maxit <- m2$hazard  <-  m2$init.B <- m2$LIMparam <- m2$LIMlogl <- m2$LIMderiv <- m2$print.times <- m2$left.censoring <- m2$init.Random <- m2$init.Eta <- m2$init.Alpha <- m2$method.GH <- m2$n.nodes <- m2$time.biomarker <- m2$init.GH <- m2$BoxCox <-  m2$init.Biomarker <- m2$... <- NULL
    Names.data.Terminal <- m2$data
    
    
    #### Frailty distribution specification ####
    if (!(all(random %in% c("y0","KG","KD","lambda")))) { stop("Random effects can be only related to variables from the longitudinal data or the intercept (1)") }
    if (!(id %in% c(names(data.Longi))) || !(id %in% c(1,names(data)))) { stop("Identification for individuals can be only related to variables from both data set") }
    
    if(length(random) == 1){
      random.which <- switch(random, "y0" = 1, "KG" = 2, "KD" = 3, "lambda" = 4)
    }else if(length(random) == 2){
      if(all(random%in%c("y0", "KG"))&all(c("y0", "KG")%in%random)){
        random.which <- 5;random <- c("y0", "KG")} #to maintain good order of random effects names
      else if(all(random%in%c("y0", "KD"))&all(c("y0", "KD")%in%random)){
        random.which <- 6; random <- c("y0", "KD")}
      else if(all(random%in%c("y0", "lambda"))&all(c("y0", "lambda")%in%random)){
        random.which <- 7; random <- c("y0", "lambda")}
      else if(all(random%in%c("KG", "KD"))&all(c("KG", "KD")%in%random)){
        random.which <- 8; random <- c("KG", "KD")}
      else if(all(random%in%c("KG", "lambda"))&all(c("KG", "lambda")%in%random)){
        random.which <- 9; random <- c("KG", "lambda")}
      else if(all(random%in%c("KD", "lambda"))&all(c("KD", "lambda")%in%random)){
        random.which <- 10; random <- c("KD", "lambda")}
      
      
    }else if(length(random) == 3){
      if(all(random%in%c("y0", "KG", "KD"))&all(c("y0", "KG", "KD")%in%random)){
        random.which <- 11; random <- c("y0", "KG", "KD")}
      else if(all(random%in%c("y0", "KG", "lambda"))&all(c("y0", "KG", "lambda")%in%random)){
        random.which <- 12; random <- c("y0", "KG", "lambda")}
      else if(all(random%in%c("y0", "KD", "lambda"))&all(c("y0", "KD", "lambda")%in%random)){
        random.which <- 13; random <- c("y0", "KD", "lambda")}
      else if(all(random%in%c("KG", "KD", "lambda"))&all(c("KG", "KD", "lambda")%in%random)){
        random.which <- 14; random <- c("KG", "KD", "lambda")}
      
      
    }else{random.which <- 15; random <- c("y0", "KG", "KD", "lambda")}
    
    
    #### Link function specification ####
    # if(!(link %in% c("Random-effects","Current-level"))){
    #    stop("Only 'Random-effects' or 'Current-level' link function can be specified in link argument.")}
    
    if(link=="Current-level"){
      stop("The link function 'Current-level' is not available yet.")}
    
    if(!(link=="Random-effects")){
      stop("Only 'Random-effects' link function can be specified in link argument.")}
    
    
    ### Left-censoring
    if(!is.null(left.censoring) && left.censoring!=FALSE){
      if(!is.numeric(left.censoring))stop("If you want to include left-censored longitudinal outcome you must give the threshold value as the argument of 'left.censoring'")
    }
    ### Biomarker, time and dose
    if(!(biomarker %in% names(data.Longi)))stop("The biomarker must be included in the data.")
    if(!(time.biomarker %in% names(data.Longi)))stop("The time of biomarker measurements must be included in the data.")
    if(!(dose %in% names(data.Longi)))stop("The dose must be included in the data.")
    ### Gauss-Hermite method
    if(all(!(c("Standard","Pseudo-adaptive") %in% method.GH))){
      stop("Only 'Standard' and 'Pseudo-adaptive' can be specified as a method for the Gaussian quadrature")
    }
    if(!init.GH%in%c(TRUE,FALSE)){
      stop("Only logical value for the option init.GH")
    }
    if(init.GH == TRUE && method.GH == "Standard"){
      stop("Initialization is only for the pseudo-adaptive Gaussian quadrature")
    }
    GH <- switch(method.GH,"Standard"=0,"Pseudo-adaptive"=1)
    
    if(!missing(n.nodes) ){
      if(n.nodes>32 || n.nodes<5 ) stop("Number of points used in the numerical integration must be chosen from the interval [5, 32]")
    }else{
      n.nodes <- 9
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
          
        }
        ### Weibull
        if (typeof == 2){
          
          size1 <- 100
          size2 <- 100
          equidistant <- 2
          
        }
        
      }
    }else{
      #### longueur hazard > 1
      if(all(!(c("Splines") %in% hazard))){
        stop("Only 'Splines'  hazard can be specified in hazard argument in this case")
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
          
        }
      }
    }
    
    
    if (missing(formula))stop("The argument formula must be specified in every model")
    if (missing(formula.KG))stop("The argument formula.KG must be specified in every model. If no covariates should be included, make it equal to 1") #AK
    if (missing(formula.KD))stop("The argument formula.KD must be specified in every model. If no covariates should be included, make it equal to 1") #AK
    if (missing(formula.terminalEvent))stop("The argument formula.terminalEvent must be specified in every model") #AK
    
    if(class(formula)!="formula")stop("The argument formula must be a formula")
    
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
    
    m <- match.call(expand.dots = FALSE) # recurrent events
    
    m$formula.KG <- m$formula.KD <- m$formula.terminalEvent <- m$biomarker <- m$dose <- m$recurrentAG <- m$data.Longi <- m$n.knots <- m$random <- m$link  <- m$id <- m$kappa <- m$maxit <- m$hazard   <-  m$init.B <- m$LIMparam <- m$LIMlogl <- m$LIMderiv <- m$left.censoring <- m$print.times <- m$init.Random <- m$init.Eta <- m$init.Alpha <- m$method.GH <- m$n.nodes <- m$time.biomarker <- m$init.GH <- m$BoxCox <-  m$init.Biomarker <- m$... <- NULL
    
    
    special <- c("strata", "cluster", "subcluster", "terminal","num.id","timedep")
    
    Terms <- if (missing(data)){
      terms(formula, special)
    }else{
      terms(formula, special, data = data)
    }
    
    ord <- attr(Terms, "order") # longueur de ord=nbre de var.expli
    
    #  if (length(ord) & any(ord != 1))stop("Interaction terms are not valid for this function")
    #si pas vide tous si il ya au moins un qui vaut 1 on arrete
    
    m$formula <- Terms
    
    
    
    m[[1]] <- as.name("model.frame") # m[[1]]=frailtypenal, il le remplace par model.frame en fait
    
    #model.frame(formula = Surv(time, event) ~ cluster(id) + as.factor(dukes) +
    #as.factor(charlson) + sex + chemo + terminal(death), data = readmission)
    
    m <- eval(m, sys.parent()) #ici la classe de m est un data.frame donc il recupere ce qu'on lui donne en argument
    
    cluster <- attr(Terms, "specials")$cluster # (indice) nbre de var qui sont en fonction de cluster()
    
    
    subcluster <- attr(Terms, "specials")$subcluster #nbre de var qui sont en fonction de subcluster()
    
    # booleen pour voir si l'objet Y est reponse avant tri des donnees Surv ou SurvIC
    classofR <- attr(model.extract(m, "response"),"class")
    # attention le package pec rajoute un element dans l'attribut "class" des objets de survie
    if (length(classofR)>1) classofR <- classofR[2]
    
    typeofR <- attr(model.extract(m, "response"),"type") # type de reponse : interval etc..
    
    #Al : tri du jeu de donnees par cluster croissant
    if (length(cluster)){
      tempc <- untangle.specials(Terms, "cluster", 1:10)
      ord <- attr(Terms, "order")[tempc$terms]
      #  if (any(ord > 1))stop("Cluster can not be used in an interaction")
      m <- m[order(m[,tempc$vars]),] # soit que des nombres, soit des caracteres
      ordre <- as.integer(row.names(m)) # recupere l'ordre du data set
      cluster <- strata(m[, tempc$vars], shortlabel = TRUE)
      uni.cluster <- unique(cluster)
    }
    
    # verification de la sutructure nested si besoin
    if (length(subcluster)) stop("'subcluster' is not an allowed option")
    
    
    if (NROW(m) == 0)stop("No (non-missing) observations") #nombre ligne different de 0
    
    R <- model.extract(m, "response") # objet de type Surv =Time
    
    
    if (classofR != "Surv") stop("Response in the 'formula' must be a survival object")
    
    
    ll <- attr(Terms, "term.labels")#liste des variables explicatives
    
    #cluster(id) as.factor(dukes) as.factor(charlson) sex chemo terminal(death)
    
    #=========================================================>
    
    mt <- attr(m, "terms") #m devient de class "formula" et "terms"
    
    X <- if (!is.empty.model(mt))model.matrix(mt, m) #idem que mt sauf que ici les factor sont divise en plusieurs variables
    
    ind.place <- unique(attr(X,"assign")[duplicated(attr(X,"assign"))]) ### unique : changement au 25/09/2014
    
    
    vec.factor <- NULL
    vec.factor <- c(vec.factor,ll[ind.place])
    
    #=========================================================>
    # On determine le nombre de categorie pour chaque var categorielle
    strats <- attr(Terms, "specials")$strata #nbre de var qui sont en fonction de strata()
    cluster <- attr(Terms, "specials")$cluster #nbre de var qui sont en fonction de cluster()
    num.id <- attr(Terms, "specials")$num.id #nbre de var qui sont en fonction de patkey()
    vartimedep <- attr(Terms, "specials")$timedep #nbre de var en fonction de timedep()
    
    #booleen pour savoir si au moins une var depend du tps
    if (is.null(vartimedep)) timedep <- 0
    else timedep <- 1
    
    if (timedep==1) stop("The option 'timedep' is not allowed in this model.")
    
    if(!is.null(num.id)) stop("num.id is not an allowed option")
    
    
    subcluster <- attr(Terms, "specials")$subcluster #nbre de var qui sont en fonction de subcluster()
    
    if (length(subcluster))stop("'subcluster' is not an allowed option")
    if (length(strats))stop("Stratified analysis is not allowed")
    
    if (length(subcluster)){
      ll <- ll[-grep("subcluster",ll)]
    }
    if (length(cluster)){
      ll_tmp <- ll[grep("cluster",ll)]
      ll <- ll[-grep("cluster",ll)]
      
      pos1 <- grep("r",unlist(strsplit(ll_tmp,split="")))[1]+2
      pos2 <- length(unlist(strsplit(ll_tmp,split="")))-1
      Names.cluster <- substr(ll_tmp,start=pos1,stop=pos2) # nom du cluster
    }
    if (length(strats)){
      ll <- ll[-grep("strata",ll)]
    }
    
    #   plus besoin de as.factor() pour afficher le test de Wald global
    if (length(grep("strata",vec.factor))) vec.factor <- vec.factor[-grep("strata",vec.factor)]
    if (length(grep("cluster",vec.factor))) vec.factor <- vec.factor[-grep("cluster",vec.factor)]
    if (length(grep("subcluster",vec.factor))) vec.factor <- vec.factor[-grep("subcluster",vec.factor)]
    if (length(grep("num.id",vec.factor))) vec.factor <- vec.factor[-grep("num.id",vec.factor)]
    
    
    
    mat.factor <- matrix(vec.factor,ncol=1,nrow=length(vec.factor))
    # Fonction servant a prendre les termes entre "as.factor" et (AK 04/11/2015) interactions
    vec.factor <-apply(mat.factor,MARGIN=1,FUN=function(x){
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
    
    
    if(length(grep("terminal",ll))>0){ind.place <- grep(paste(vec.factor,collapse="|"),ll[-grep("terminal",ll)])
    }else{ind.place <- grep(paste(vec.factor,collapse="|"),ll)}
    
    if(length(vec.factor) > 0){
      vect.fact <- attr(X,"dimnames")[[2]]
      
      #vect.fact <- vect.fact[grep("factor",vect.fact)]
      vect.fact <- vect.fact[grep(paste(vec.factor,collapse="|"),vect.fact)]
      
      #   	vect.fact <-apply(matrix(vect.fact,ncol=1,nrow=length(vect.fact)),MARGIN=1,FUN=function(x){
      # 		pos1 <- grep("r",unlist(strsplit(x,split="")))[1]+2
      # 		pos2 <- grep(")",unlist(strsplit(x,split="")))[1]-1
      # 		return(substr(x,start=pos1,stop=pos2))})
      occur <- rep(0,length(vec.factor))
      
      interaction<-as.vector(apply(matrix(vect.fact,nrow=length(vect.fact)),MARGIN=1,FUN=function(x){length(grep(":",unlist(strsplit(x,split=""))))}))
      which.interaction <- which(interaction==1)
      
      for(i in 1:length(vec.factor)){
        
        if(length(grep(":",unlist(strsplit(vec.factor[i],split=""))))>0){
          
          
          pos <- grep(":",unlist(strsplit(vec.factor[i],split="")))
          length.grep <- 0
          for(j in 1:length(vect.fact)){
            if(j%in%which.interaction){
              
              if(length(grep(substr(vec.factor[i],start=1,stop=pos-1),vect.fact[j]))>0 && length(grep(substr(vec.factor[i],start=pos+1,stop=length(unlist(strsplit(vec.factor[i],split="")))),vect.fact[j]))>0){
                length.grep <- length.grep + 1
                which <- i}
            }}
          occur[i] <- length.grep
          
        }else{
          
          
          if(length(vect.fact[-which.interaction])>0){occur[i] <- length(grep(vec.factor[i],vect.fact[-which.interaction]))
          }else{occur[i] <- length(grep(vec.factor[i],vect.fact))}
        }
      }
    }
    
    
    
    #=========================================================>
    
    terminalEvent <- attr(Terms, "specials")$terminal #nbre de var qui sont en fonction de terminal()
    
    dropx <- NULL
    
    if (length(cluster) ){
      tempc <- untangle.specials(Terms, "cluster", 1:10)
      ord <- attr(Terms, "order")[tempc$terms]
      if (any(ord > 1))stop("Cluster can not be used in an interaction")
      
      cluster <- strata(m[, tempc$vars], shortlabel = TRUE)
      
      dropx <- tempc$terms
      uni.cluster<-unique(cluster)
      
    }else if (!length(cluster)){
      
      stop("grouping variable is needed in the formula")
      
    }
    
    if (length(num.id))stop("'num.id' is not an allowed option")
    
    if(length(uni.cluster)==1){
      stop("grouping variable must have more than 1 level")
    }
    
    
    if (length(subcluster))stop("subcluster is not an allowed option")
    
    
    
    if (length(strats))stop("Stratified analysis is not allowed")
    
    if (typeof==0 && (length(kappa)!=2)) stop("wrong length of argument 'kappa'. The length should be 2")
    
    #AD: indicator of terminal()
    ind.terminal <- length(terminalEvent)
    
    #AD:
    if (length(terminalEvent)){
      
      tempterm <- untangle.specials(Terms, "terminal", 1:10)
      #ici on comme terme tempterm$vars qui est le nom dans l'appel(ex;"terminal(death)"
      #et tempterm$terms qui est la position de la variable dans l'appel, ici elle vient a la position 6
      
      ord <- attr(Terms, "order")[tempterm$terms] # ord[6]=1 ici dans notre exemple
      
      if (any(ord > 1))stop("Terminal can not be used in an interaction")
      dropx <- c(dropx,tempterm$terms) # vecteur de position
      terminal <- strata(m[, tempterm$vars], shortlabel = TRUE)
      terminal <- as.numeric(as.character(terminal))
      
    }
    
    #type <- attr(Y, "type")
    type <- typeofR
    
    if (type != "right" && type != "counting") { # Cox supporte desormais la censure par intervalle
      stop(paste("Cox model doesn't support \"", type, "\" survival data", sep = ""))
    }
    
    if (type != "counting" && recurrentAG) {
      stop("recurrentAG needs counting process formulation")
    }
    
    if (length(dropx)){
      newTerms <- Terms[-dropx]
    }else{
      newTerms <- Terms
    }
    
    #newTerm vaut Terms - les variables dont les position sont dans drop
    
    X <- model.matrix(newTerms, m)
    
    assign <- lapply(attrassign(X, newTerms)[-1], function(x) x - 1)
    Xlevels <- .getXlevels(newTerms, m)
    contr.save <- attr(X, 'contrasts')
    
    
    # assigne donne la position pour chaque variables
    #ncol(X) : nombre de variable sans sans les fonction speciaux comme terminal()...+id
    if(length(vec.factor) > 0){
      #========================================>
      position <- unlist(assign,use.names=F)
    }
    
    #========================================>
    
    if (ncol(X) == 1){
      X<-X-1
      noVarR <- 1
    }else{
      X <- X[, -1, drop = FALSE]
      noVarR <- 0
    }
    # on enleve ensuite la premiere colonne correspondant a id
    
    
    nvar<-ncol(X) #nvar==1 correspond a 2 situations:
    
    # au cas ou on a aucune var explicative dans la partie rec, mais X=0
    # cas ou on a 1seul var explicative, ici X est en general different de 0
    
    #varnotdep <- colnames(X)[-grep("timedep",colnames(X))]
    #vardep <- colnames(X)[grep("timedep",colnames(X))]
    #vardep <- apply(matrix(vardep,ncol=1,nrow=length(vardep)),1,timedep.names)
    
    #if (length(intersect(varnotdep,vardep)) != 0) {
    #  stop("A variable is both used as a constant and time-varying effect covariate")
    #}
    
    #nvartimedep <- length(vardep)
    
    #filtretps <- rep(0,nvar)
    #filtretps[grep("timedep",colnames(X))] <- 1
    
    var<-matrix(c(X),nrow=nrow(X),ncol=nvar) #matrix sans id et sans partie ex terminal(death)
    
    nsujet<-nrow(X)
    
    if (type=="right"){
      tt0 <- rep(0,nsujet)
      tt1 <- R[,1]
      cens <- R[,2]
    } else {
      tt0 <- R[,1]
      tt1 <- R[,2]
      cens <- R[,3]
    }
    
    if (min(cens)==0) cens.data<-1
    if (min(cens)==1 && max(cens)==1) cens.data<-0
    
    AG<-ifelse(recurrentAG,1,0)
    
    
    #=======================================>
    #======= Construction du vecteur des indicatrice
    if(length(vec.factor) > 0){
      #         ind.place <- ind.place -1
      k <- 0
      for(i in 1:length(vec.factor)){
        ind.place[i] <- ind.place[i]+k
        k <- k + occur[i]-1
      }
    }
    
    
    
    #========= Longitudinal Data preparation =========================
    
    Y <- data.Longi[,which(names(data.Longi)==biomarker)]
    
    if(!is.null(formula.KG[3]) && formula.KG[3] != "1()"){
      
      TermsKG <- if (missing(data.Longi)){
        terms(formula.KG, special)
      }else{
        terms(formula.KG, special, data = data.Longi)
      }
      
      ord <- attr(TermsKG, "order") # longueur de ord=nbre de var.expli
      
      #si pas vide tous si il ya au moins un qui vaut 1 on arrete
      
      m3$formula.KG <- TermsKG
      
      
      m3[[1]] <- as.name("model.frame") # m[[1]]=frailtypenal, il le remplace par model.frame en fait
      
      
      if (NROW(m3) == 0)stop("No (non-missing) observations") #nombre ligne different de 0
      
      llKG <- attr(TermsKG, "term.labels")#liste des variables explicatives
      
      
      #=========================================================>
      data.Longi <- data.Longi[order(data.Longi[,id]),]
      name.KG <- as.character(attr(TermsKG, "variables")[[2]])
      KG <- data.Longi[,which(names(data.Longi)==name.KG)]
      
      
      # on identifie les variables explicatives facteurs avec nombre de niveau plus que 2
      
      # print(as.data.frame(data.Longi[,which(names(data.Longi)%in%llKG)]))
      
      
      # which(llKG%in%names(which(lapply(data.Longi[,which(names(data.Longi)%in%llKG)],function(x) length(levels(x)))>2)))
      ind.placeKG <- which(lapply(as.data.frame(data.Longi[,which(names(data.Longi)%in%llKG)]),function(x) length(levels(x)))>2)
      
      defined.factor <- llKG[grep("factor",llKG)]
      
      vec.factorKG.tmp <- NULL
      if(length(defined.factor)>0){
        mat.factorKG.tmp <- matrix(defined.factor,ncol=1,nrow=length(defined.factor))
        
        # Fonction servant a prendre les termes entre "as.factor"
        vec.factorKG.tmp <-apply(mat.factorKG.tmp,MARGIN=1,FUN=function(x){
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
        
        vec.factorKG.tmp <- vec.factorKG.tmp[which(vec.factorKG.tmp!="NaN")]
        
        if(length(vec.factorKG.tmp)>0){
          for(i in 1:length(vec.factorKG.tmp)){
            if(length(grep(":",vec.factorKG.tmp[i]))==0){
              if(length(levels(as.factor(data.Longi[,which(names(data.Longi)==vec.factorKG.tmp[i])])))>2)ind.placeKG <- c(ind.placeKG,which(llKG%in%paste("as.factor(",vec.factorKG.tmp[i],")",sep="")))
            }
            
          }}
        
      }
      
      ind.placeKG <- sort(ind.placeKG)
      
      #=========================================================>
      # On determine le nombre de categorie pour chaque var categorielle
      stratsKG <- attr(TermsKG, "specials")$strata #nbre de var qui sont en fonction de strata()
      clusterKG <- attr(TermsKG, "specials")$cluster #nbre de var qui sont en fonction de cluster()
      num.idKG <- attr(TermsKG, "specials")$num.id #nbre de var qui sont en fonction de patkey()
      vartimedepKG <- attr(TermsKG, "specials")$timedep #nbre de var en fonction de timedep()
      
      #booleen pour savoir si au moins une var depend du tps
      if (!is.null(vartimedepKG)) stop("The option 'timedep' is not allowed in this  model.")
      if (!is.null(stratsKG)) stop("Stratified analysis is not allowed in this model.")
      if (!is.null(clusterKG)) stop("Only the argument 'id' can represent the clusters.")
      if (!is.null(num.idKG))stop("'num.id' for the interval censoring is an allowed option")
      
      
      
      Names.clusterY <- id # nom du cluster
      
      
      # which_n<-which(names(data.Longi)%in%random)
      #  data.Longi[which(data.Longi[,which_n]==0),which_n]<- 0.01
      
      mat.factorKG2 <- matrix(llKG,ncol=1,nrow=length(llKG))
      
      # Fonction servant a prendre les termes entre "as.factor"
      llKG2 <-apply(mat.factorKG2,MARGIN=1,FUN=function(x){
        if (length(grep("factor",x))>0  && length(grep(":",x))==0 && unlist(strsplit(x,split=""))[length(unlist(strsplit(x,split="")))]==")"){
          pos1 <- grep("r",unlist(strsplit(x,split="")))[1]+2
          pos2 <- length(unlist(strsplit(x,split="")))-1
          x<-substr(x,start=pos1,stop=pos2)
          return(paste(x,levels(as.factor(data.Longi[,which(names(data.Longi)==x)]))[2],sep=""))
        }else{
          return(x)
        }})
      
      # Fonction servant a prendre les termes entre "as.factor" - without the name of the level
      llKG3 <-apply(mat.factorKG2,MARGIN=1,FUN=function(x){
        if (length(grep("factor",x))>0  && length(grep(":",x))==0 && unlist(strsplit(x,split=""))[length(unlist(strsplit(x,split="")))]==")"){
          pos1 <- grep("r",unlist(strsplit(x,split="")))[1]+2
          pos2 <- length(unlist(strsplit(x,split="")))-1
          return(substr(x,start=pos1,stop=pos2))
        }else{
          return(x)
        }})
      
      llKG.real.names <- llKG3  
      llKG3 <- llKG3[!llKG2%in%llKG]
      
      
      
      if(is.factor(data.Longi[,names(data.Longi)==llKG.real.names[1]]))X_KG<- as.numeric(data.Longi[,names(data.Longi)==llKG.real.names[1]])-1
      else X_KG<- data.Longi[,names(data.Longi)==llKG.real.names[1]]
      
      if(length(llKG) == 1)X_KG <- as.data.frame(X_KG)
      
      
      if(length(llKG)>1){
        
        
        for(i in 2:length(llKG.real.names)){
          
          if(is.factor(data.Longi[,names(data.Longi)==llKG.real.names[i]]))X_KG<- cbind(X_KG,as.numeric(data.Longi[,names(data.Longi)==llKG.real.names[i]])-1)
          else X_KG<- cbind(X_KG,data.Longi[,names(data.Longi)==llKG.real.names[i]])
        }
      }
      #X_KG<- data.Longi[,names(data.Longi)%in%(llY)]
      
      llKG.fin <- llKG.real.names
      llKG <- llKG.real.names
      
      if(sum(ord)>length(ord)){
        
        for(i in 1:length(ord)){
          if(ord[i]>1){
            
            name_v1 <- strsplit(as.character(llKG[i]),":")[[1]][1]
            name_v2 <- strsplit(as.character(llKG[i]),":")[[1]][2]
            
            if(length(grep("factor",name_v1))>0){name_v1<-substring(name_v1,11,nchar(name_v1)-1)
            v1 <- as.factor(data.Longi[,names(data.Longi)==name_v1])}
            else{v1 <- data.Longi[,names(data.Longi)==name_v1]}
            if(length(grep("factor",name_v2))>0){name_v2<-substring(name_v2,11,nchar(name_v2)-1)
            v2 <- as.factor(data.Longi[,names(data.Longi)==name_v2])}
            else{v2 <- data.Longi[,names(data.Longi)==name_v2]}
            
            llKG[i] <- paste(name_v1,":",name_v2,sep="")
            #   if(is.factor(v1) && length(levels(v1))>2)stop("Interactions not allowed for factors with 3 or more levels (yet)")
            #   if(is.factor(v2) && length(levels(v2))>2)stop("Interactions not allowed for factors with 3 or more levels (yet)")
            if(is.factor(v1) && !is.factor(v2)){
              
              dummy <- model.matrix( ~ v1 - 1)
              # if(length(levels(v1)>2))vec.factorY <- c(vec.factorY,paste(name_v1,":",name_v2,sep=""))
              for(j in 2:length(levels(v1))){
                X_KG <- cbind(X_KG,dummy[,j]*v2)
                if(i>1 && i<length(llKG.fin))llKG.fin <- c(llKG.fin[1:(i-1+j-2)],paste(name_v1,".",levels(v1)[j],":",name_v2,sep=""),llKG.fin[(i+1+j-2):length(llKG.fin)])
                else if(i==length(llKG.fin))llKG.fin <- c(llKG.fin[1:(i-1+j-2)],paste(name_v1,".",levels(v1)[j],":",name_v2,sep=""))
                else llKG.fin <- c(paste(name_v1,".",levels(v1)[j],":",name_v2,sep=""),llKG.fin[(2+j-2):length(llKG.fin)])
              }
              
            }else if(!is.factor(v1) && is.factor(v2)){
              
              dummy <- model.matrix( ~ v2 - 1)
              #  if(length(levels(v2)>2))vec.factorY <- c(vec.factorY,paste(name_v1,":",name_v2,sep=""))
              for(j in 2:length(levels(v2))){
                
                X_KG <- cbind(X_KG,dummy[,j]*v1)
                
                if(i>1 && i<length(llKG.fin))llKG.fin <- c(llKG.fin[1:(i-1+j-2)],paste(name_v1,":",name_v2,levels(v2)[j],sep=""),llKG.fin[(i+1+j-2):length(llKG.fin)])
                else if(i==length(llKG.fin))llKG.fin <- c(llKG.fin[1:(i-1+j-2)],paste(name_v1,":",name_v2,levels(v2)[j],sep=""))
                else llKG.fin <- c(paste(name_v1,":",name_v2,levels(v2)[j],sep=""),llKG.fin[(2+j-2):length(llKG.fin)])
              }
            }else if(is.factor(v1) && is.factor(v2)){
              
              
              dummy1 <- model.matrix( ~ v1 - 1)
              dummy2 <- model.matrix( ~ v2 - 1)
              #   if(length(levels(v1)>2) || length(levels(v2)>2))vec.factorY <- c(vec.factorY,paste(name_v1,":",name_v2,sep=""))
              for(j in 2:length(levels(v1))){
                for(k in 2:length(levels(v2))){
                  
                  X_KG <- cbind(X_KG,dummy1[,j]*dummy2[,k])
                  if(i>1 && i<length(llKG.fin))llKG.fin <- c(llKG.fin[1:(i-1+j-2+k-2)],paste(name_v1,levels(v1)[j],":",name_v2,levels(v2)[k],sep=""),llKG.fin[(i+1+j-2+k-2):length(llKG.fin)])
                  else if(i==length(llKG.fin))llKG.fin <- c(llKG.fin[1:(i-1+j-2+k-2)],paste(name_v1,levels(v1)[j],":",name_v2,levels(v2)[k],sep=""))
                  else llKG.fin <- c(paste(name_v1,levels(v1)[j],":",name_v2,levels(v2)[k],sep=""),llKG.fin[(2+j-2+k-2):length(llKG.fin)])
                }
              } 
            }else{
              
              X_KG <- cbind(X_KG,v1*v2)
            }
            
          }
        }
      }
      
      
      if(length(grep(":",llKG))>0){
        for(i in 1:length(grep(":",llKG))){
          if(length(levels(data.Longi[,which(names(data.Longi)%in%strsplit(llKG[grep(":",llKG)[i]],":")[[1]])[1]]))>2 || length(levels(data.Longi[,which(names(data.Longi)%in%strsplit(llKG[grep(":",llKG)[i]],":")[[1]])[2]]))>2){
            ind.placeKG <- c(ind.placeKG,grep(":",llKG)[i])
            #     vec.factorY <- c(vec.factorY,llY[grep(":",llY)[i]])
          }
        }
      }
      vec.factorKG <- NULL
      
      if(length(vec.factorKG.tmp)>0)vec.factorKG <- c(llKG[ind.placeKG],vec.factorKG.tmp)
      else vec.factorKG <- c(vec.factorKG,llKG[ind.placeKG])
      
      vec.factorKG <- unique(vec.factorKG)
      
      
      
      mat.factorKG <- matrix(vec.factorKG,ncol=1,nrow=length(vec.factorKG))
      # Fonction servant a prendre les termes entre "as.factor" et (AK 04/11/2015) interactions
      vec.factorKG <-apply(mat.factorKG,MARGIN=1,FUN=function(x){
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
      
      
      for(i in 1:length(llKG.fin)){
        
        if(sum(names(data.Longi)==llKG.fin[i])>0){
          if(is.factor(data.Longi[,names(data.Longi)==llKG.fin[i]]) && length(levels(data.Longi[,names(data.Longi)==llKG.fin[i]]))==2){
            llKG.fin[i] <- paste(llKG.fin[i],levels(data.Longi[,names(data.Longi)==llKG.fin[i]])[2],sep="")}
        }
      }
      
      #  llY <- llY.fin
      X_KG <- as.data.frame(X_KG)
      if(dim(X_KG)[2]!=length(llKG.fin))stop("The variables in the longitudinal part must be in the data.Longi")
      
      names(X_KG) <- llKG.fin
      
      
      X_KGall<- X_KG
      "%+%"<- function(x,y) paste(x,y,sep="")
      if(length(vec.factorKG) > 0){
        for(i in 1:length(vec.factorKG)){
          if(length(grep(":",vec.factorKG[i]))==0){
            
            factor.spot <- which(names(X_KG)==vec.factorKG[i])
            
            if(factor.spot<ncol(X_KG))  X_KG <- cbind(X_KG[1:(factor.spot-1)],model.matrix(as.formula("~"%+%0%+%"+"%+%paste(vec.factorKG[i], collapse= "+")), model.frame(~.,data.Longi,na.action=na.pass))[,-1],X_KG[(factor.spot+1):ncol(X_KG)])
            else X_KG <- cbind(X_KG[1:(factor.spot-1)],model.matrix(as.formula("~"%+%0%+%"+"%+%paste(vec.factorKG[i], collapse= "+")), model.frame(~.,data.Longi,na.action=na.pass))[,-1])
            
          } }
        
        
        
        
        
        vect.factKG<-names(X_KG)[which(!(names(X_KG)%in%llKG))]
        
        
        #               vect.fact <-apply(matrix(vect.fact,ncol=1,nrow=length(vect.fact)),MARGIN=1,FUN=function(x){
        #               pos1 <- grep("r",unlist(strsplit(x,split="")))[1]+2
        #               pos2 <- grep(")",unlist(strsplit(x,split="")))[1]-1
        #               return(substr(x,start=pos1,stop=pos2))})
        
        occurKG <- rep(0,length(vec.factorKG))
        
        #         for(i in 1:length(vec.factorY)){
        #                #occur[i] <- sum(vec.factor[i] == vect.fact)
        #               occurY[i] <- length(grep(vec.factorY[i],vect.factY))
        #      }
        
        
        
        interaction<-as.vector(apply(matrix(vect.factKG,nrow=length(vect.factKG)),MARGIN=1,FUN=function(x){length(grep(":",unlist(strsplit(x,split=""))))}))
        which.interaction <- which(interaction==1)
        
        for(i in 1:length(vec.factorKG)){
          
          if(length(grep(":",unlist(strsplit(vec.factorKG[i],split=""))))>0){
            
            
            pos <- grep(":",unlist(strsplit(vec.factorKG[i],split="")))
            length.grep <- 0
            for(j in 1:length(vect.factKG)){
              if(j%in%which.interaction){
                
                if(length(grep(substr(vec.factorKG[i],start=1,stop=pos-1),vect.factKG[j]))>0 && length(grep(substr(vec.factorKG[i],start=pos+1,stop=length(unlist(strsplit(vec.factorKG[i],split="")))),vect.factKG[j]))>0){
                  length.grep <- length.grep + 1
                  which <- i}
              }}
            occurKG[i] <- length.grep
            
          }else{
            
            
            if(length(vect.factKG[-which.interaction])>0){occurKG[i] <- length(grep(vec.factorKG[i],vect.factKG[-which.interaction]))
            }else{occurKG[i] <- length(grep(vec.factorKG[i],vect.factKG))}
          }
        }
      }
      
      if (ncol(X_KG) == 0){
        noVarKG <- 1
      }else{
        noVarKG <- 0
      }
      nvarKG<-ncol(X_KG) #nvar==1 correspond a 2 situations:
      varKG <- as.matrix(sapply(X_KG, as.numeric))
      
      nsujety<-nrow(X_KG)
      
    }
    
    
    #=========================================================>
    
    clusterY <- data.Longi[,which(colnames(data.Longi)==id)]
    
    max_rep <- max(table(clusterY))
    uni.clusterY<-as.factor(unique(clusterY))
    
    
    if(is.null(id))     stop("grouping variable is needed")
    
    if(is.null(random)) stop("variable for random effects is needed")
    
    
    if(length(uni.clusterY)==1){
      stop("grouping variable must have more than 1 level in the longitudinal part")
    }
    
    
    if (length(subcluster))stop("'Subcluster' is not allowed")
    
    
    
    if (typeof==0 && missing(kappa)) stop("smoothing parameter (kappa) is required")
    
    
    
    
    #newTerm vaut Terms - les variables dont les position sont dans drop
    
    
    
    
    
    
    #========================================>
    
    
    
    # au cas ou on a aucune var explicative dans la partie rec, mais X=0
    # cas ou on a 1seul var explicative, ici X est en general different de 0
    
    #varnotdepY <- colnames(X_KG)[-grep("timedep",colnames(X_KG))]
    #vardepY <- colnames(X_KG)[grep("timedep",colnames(X_KG))]
    #vardepY <- apply(matrix(vardepY,ncol=1,nrow=length(vardepY)),1,timedep.names)
    
    #if (length(intersect(varnotdepY,vardepY)) != 0) {
    #  stop("A variable is both used as a constant and time-varying effect covariate")
    #}
    
    #nvartimedepY <- length(vardepY)
    
    #filtretpsY <- rep(0,nvarY)
    #filtretpsY[grep("timedep",colnames(X_KG))] <- 1
    
    
    
    if(!is.null(formula.KG[3]) && formula.KG[3] != "1()"){
      
      
      
      #=======================================>
      #======= Construction du vecteur des indicatrice
      if(length(vec.factorKG) > 0){
        #         ind.place <- ind.place -1
        k <- 0
        for(i in 1:length(vec.factorKG)){
          ind.placeKG[i] <- ind.placeKG[i]+k
          k <- k + occurKG[i]-1
        }
      }
      
      
    }else{
      noVarKG <- 1 
      nvarKG <- 0 
      varKG <- c()#rep(0, dim(data.Longi)[1])
      X_KG <- c()
      nsujety <- length(Y)
      vec.factorKG <- NULL
    }
    
    #========= Longitudinal Data preparation - KD =========================
    
    if(!is.null(formula.KD[3]) && formula.KD[3] != "1()"){
      
      TermsKD <- if (missing(data.Longi)){
        terms(formula.KD, special)
      }else{
        terms(formula.KD, special, data = data.Longi)
      }
      
      ord <- attr(TermsKD, "order") # longueur de ord=nbre de var.expli
      
      #si pas vide tous si il ya au moins un qui vaut 1 on arrete
      
      m4$formula.KD <- TermsKD
      
      
      m4[[1]] <- as.name("model.frame") # m[[1]]=frailtypenal, il le remplace par model.frame en fait
      
      
      if (NROW(m4) == 0)stop("No (non-missing) observations") #nombre ligne different de 0
      
      llKD <- attr(TermsKD, "term.labels")#liste des variables explicatives
      
      
      #=========================================================>
      
      name.KD <- as.character(attr(TermsKD, "variables")[[2]])
      KD <- data.Longi[,which(names(data.Longi)==name.KD)]
      
      
      # on identifie les variables explicatives facteurs avec nombre de niveau plus que 2
      
      # print(as.data.frame(data.Longi[,which(names(data.Longi)%in%llKG)]))
      
      
      # which(llKG%in%names(which(lapply(data.Longi[,which(names(data.Longi)%in%llKG)],function(x) length(levels(x)))>2)))
      ind.placeKD <- which(lapply(as.data.frame(data.Longi[,which(names(data.Longi)%in%llKD)]),function(x) length(levels(x)))>2)
      
      defined.factor <- llKD[grep("factor",llKD)]
      
      vec.factorKD.tmp <- NULL
      if(length(defined.factor)>0){
        mat.factorKD.tmp <- matrix(defined.factor,ncol=1,nrow=length(defined.factor))
        
        # Fonction servant a prendre les termes entre "as.factor"
        vec.factorKD.tmp <-apply(mat.factorKG.tmp,MARGIN=1,FUN=function(x){
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
        
        vec.factorKD.tmp <- vec.factorKD.tmp[which(vec.factorKD.tmp!="NaN")]
        
        if(length(vec.factorKD.tmp)>0){
          for(i in 1:length(vec.factorKD.tmp)){
            if(length(grep(":",vec.factorKD.tmp[i]))==0){
              if(length(levels(as.factor(data.Longi[,which(names(data.Longi)==vec.factorKD.tmp[i])])))>2)ind.placeKD <- c(ind.placeKD,which(llKD%in%paste("as.factor(",vec.factorKD.tmp[i],")",sep="")))
            }
            
          }}
        
      }
      
      ind.placeKD <- sort(ind.placeKD)
      
      #=========================================================>
      # On determine le nombre de categorie pour chaque var categorielle
      stratsKD <- attr(TermsKD, "specials")$strata #nbre de var qui sont en fonction de strata()
      clusterKD <- attr(TermsKD, "specials")$cluster #nbre de var qui sont en fonction de cluster()
      num.idKD <- attr(TermsKD, "specials")$num.id #nbre de var qui sont en fonction de patkey()
      vartimedepKD <- attr(TermsKD, "specials")$timedep #nbre de var en fonction de timedep()
      
      #booleen pour savoir si au moins une var depend du tps
      if (!is.null(vartimedepKD)) stop("The option 'timedep' is not allowed in this  model.")
      if (!is.null(stratsKD)) stop("Stratified analysis is not allowed in this model.")
      if (!is.null(clusterKD)) stop("Only the argument 'id' can represent the clusters.")
      if (!is.null(num.idKD))stop("'num.id' for the interval censoring is an allowed option")
      
      
      
      # which_n<-which(names(data.Longi)%in%random)
      #  data.Longi[which(data.Longi[,which_n]==0),which_n]<- 0.01
      
      mat.factorKD2 <- matrix(llKD,ncol=1,nrow=length(llKD))
      
      # Fonction servant a prendre les termes entre "as.factor"
      llKD2 <-apply(mat.factorKD2,MARGIN=1,FUN=function(x){
        if (length(grep("factor",x))>0  && length(grep(":",x))==0 && unlist(strsplit(x,split=""))[length(unlist(strsplit(x,split="")))]==")"){
          pos1 <- grep("r",unlist(strsplit(x,split="")))[1]+2
          pos2 <- length(unlist(strsplit(x,split="")))-1
          x<-substr(x,start=pos1,stop=pos2)
          return(paste(x,levels(as.factor(data.Longi[,which(names(data.Longi)==x)]))[2],sep=""))
        }else{
          return(x)
        }})
      
      # Fonction servant a prendre les termes entre "as.factor" - without the name of the level
      llKD3 <-apply(mat.factorKD2,MARGIN=1,FUN=function(x){
        if (length(grep("factor",x))>0  && length(grep(":",x))==0 && unlist(strsplit(x,split=""))[length(unlist(strsplit(x,split="")))]==")"){
          pos1 <- grep("r",unlist(strsplit(x,split="")))[1]+2
          pos2 <- length(unlist(strsplit(x,split="")))-1
          return(substr(x,start=pos1,stop=pos2))
        }else{
          return(x)
        }})
      
      llKD.real.names <- llKD3  
      llKD3 <- llKD3[!llKD2%in%llKD]
      
      
      if(is.factor(data.Longi[,names(data.Longi)==llKD.real.names[1]]))X_KD<- as.numeric(data.Longi[,names(data.Longi)==llKD.real.names[1]])-1
      else X_KD<- data.Longi[,names(data.Longi)==llKD.real.names[1]]
      
      
      
      
      if(length(llKD)>1){
        for(i in 2:length(llKD.real.names)){
          
          if(is.factor(data.Longi[,names(data.Longi)==llKD.real.names[i]]))X_KD<- cbind(X_KD,as.numeric(data.Longi[,names(data.Longi)==llKD.real.names[i]])-1)
          else X_KD<- cbind(X_KD,data.Longi[,names(data.Longi)==llKD.real.names[i]])
        }}
      
      #X_KD<- data.Longi[,names(data.Longi)%in%(llY)]
      
      llKD.fin <- llKD.real.names
      llKD <- llKD.real.names
      
      if(sum(ord)>length(ord)){
        
        for(i in 1:length(ord)){
          if(ord[i]>1){
            
            name_v1 <- strsplit(as.character(llKD[i]),":")[[1]][1]
            name_v2 <- strsplit(as.character(llKD[i]),":")[[1]][2]
            
            if(length(grep("factor",name_v1))>0){name_v1<-substring(name_v1,11,nchar(name_v1)-1)
            v1 <- as.factor(data.Longi[,names(data.Longi)==name_v1])}
            else{v1 <- data.Longi[,names(data.Longi)==name_v1]}
            if(length(grep("factor",name_v2))>0){name_v2<-substring(name_v2,11,nchar(name_v2)-1)
            v2 <- as.factor(data.Longi[,names(data.Longi)==name_v2])}
            else{v2 <- data.Longi[,names(data.Longi)==name_v2]}
            
            llKD[i] <- paste(name_v1,":",name_v2,sep="")
            #   if(is.factor(v1) && length(levels(v1))>2)stop("Interactions not allowed for factors with 3 or more levels (yet)")
            #   if(is.factor(v2) && length(levels(v2))>2)stop("Interactions not allowed for factors with 3 or more levels (yet)")
            if(is.factor(v1) && !is.factor(v2)){
              
              dummy <- model.matrix( ~ v1 - 1)
              # if(length(levels(v1)>2))vec.factorY <- c(vec.factorY,paste(name_v1,":",name_v2,sep=""))
              for(j in 2:length(levels(v1))){
                X_KD <- cbind(X_KD,dummy[,j]*v2)
                if(i>1 && i<length(llKD.fin))llKD.fin <- c(llKD.fin[1:(i-1+j-2)],paste(name_v1,".",levels(v1)[j],":",name_v2,sep=""),llKD.fin[(i+1+j-2):length(llKD.fin)])
                else if(i==length(llKD.fin))llKD.fin <- c(llKD.fin[1:(i-1+j-2)],paste(name_v1,".",levels(v1)[j],":",name_v2,sep=""))
                else llKD.fin <- c(paste(name_v1,".",levels(v1)[j],":",name_v2,sep=""),llKD.fin[(2+j-2):length(llKD.fin)])
              }
              
            }else if(!is.factor(v1) && is.factor(v2)){
              
              dummy <- model.matrix( ~ v2 - 1)
              #  if(length(levels(v2)>2))vec.factorY <- c(vec.factorY,paste(name_v1,":",name_v2,sep=""))
              for(j in 2:length(levels(v2))){
                
                X_KD <- cbind(X_KD,dummy[,j]*v1)
                
                if(i>1 && i<length(llKD.fin))llKD.fin <- c(llKD.fin[1:(i-1+j-2)],paste(name_v1,":",name_v2,levels(v2)[j],sep=""),llKD.fin[(i+1+j-2):length(llKD.fin)])
                else if(i==length(llKD.fin))llKD.fin <- c(llKD.fin[1:(i-1+j-2)],paste(name_v1,":",name_v2,levels(v2)[j],sep=""))
                else llKD.fin <- c(paste(name_v1,":",name_v2,levels(v2)[j],sep=""),llKD.fin[(2+j-2):length(llKD.fin)])
              }
            }else if(is.factor(v1) && is.factor(v2)){
              
              
              dummy1 <- model.matrix( ~ v1 - 1)
              dummy2 <- model.matrix( ~ v2 - 1)
              #   if(length(levels(v1)>2) || length(levels(v2)>2))vec.factorY <- c(vec.factorY,paste(name_v1,":",name_v2,sep=""))
              for(j in 2:length(levels(v1))){
                for(k in 2:length(levels(v2))){
                  
                  X_KD <- cbind(X_KD,dummy1[,j]*dummy2[,k])
                  if(i>1 && i<length(llKD.fin))llKD.fin <- c(llKD.fin[1:(i-1+j-2+k-2)],paste(name_v1,levels(v1)[j],":",name_v2,levels(v2)[k],sep=""),llKG.fin[(i+1+j-2+k-2):length(llKG.fin)])
                  else if(i==length(llKD.fin))llKD.fin <- c(llKD.fin[1:(i-1+j-2+k-2)],paste(name_v1,levels(v1)[j],":",name_v2,levels(v2)[k],sep=""))
                  else llKD.fin <- c(paste(name_v1,levels(v1)[j],":",name_v2,levels(v2)[k],sep=""),llKD.fin[(2+j-2+k-2):length(llKD.fin)])
                }
              } 
            }else{
              
              X_KD <- cbind(X_KD,v1*v2)
            }
            
          }
        }
      }
      
      
      
      if(length(grep(":",llKD))>0){
        for(i in 1:length(grep(":",llKD))){
          if(length(levels(data.Longi[,which(names(data.Longi)%in%strsplit(llKD[grep(":",llKD)[i]],":")[[1]])[1]]))>2 || length(levels(data.Longi[,which(names(data.Longi)%in%strsplit(llKD[grep(":",llKD)[i]],":")[[1]])[2]]))>2){
            ind.placeKD <- c(ind.placeKD,grep(":",llKD)[i])
            #     vec.factorY <- c(vec.factorY,llY[grep(":",llY)[i]])
          }
        }
      }
      
      vec.factorKD <- NULL
      
      if(length(vec.factorKD.tmp)>0)vec.factorKD <- c(llKD[ind.placeKD],vec.factorKD.tmp)
      else vec.factorKD <- c(vec.factorKD,llKD[ind.placeKD])
      
      vec.factorKD <- unique(vec.factorKD)
      
      
      
      mat.factorKD <- matrix(vec.factorKD,ncol=1,nrow=length(vec.factorKD))
      # Fonction servant a prendre les termes entre "as.factor" et (AK 04/11/2015) interactions
      vec.factorKD <-apply(mat.factorKD,MARGIN=1,FUN=function(x){
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
      
      
      for(i in 1:length(llKD.fin)){
        
        if(sum(names(data.Longi)==llKD.fin[i])>0){
          if(is.factor(data.Longi[,names(data.Longi)==llKD.fin[i]]) && length(levels(data.Longi[,names(data.Longi)==llKD.fin[i]]))==2){
            llKD.fin[i] <- paste(llKD.fin[i],levels(data.Longi[,names(data.Longi)==llKD.fin[i]])[2],sep="")}
        }
      }
      
      #  llY <- llY.fin
      X_KD <- as.data.frame(X_KD)
      if(dim(X_KD)[2]!=length(llKD.fin))stop("The variables in the longitudinal part must be in the data.Longi")
      
      names(X_KD) <- llKD.fin
      
      
      X_KDall<- X_KD
      "%+%"<- function(x,y) paste(x,y,sep="")
      if(length(vec.factorKD) > 0){
        for(i in 1:length(vec.factorKD)){
          if(length(grep(":",vec.factorKD[i]))==0){
            
            factor.spot <- which(names(X_KD)==vec.factorKD[i])
            
            if(factor.spot<ncol(X_KD))  X_KD <- cbind(X_KD[1:(factor.spot-1)],model.matrix(as.formula("~"%+%0%+%"+"%+%paste(vec.factorKD[i], collapse= "+")), model.frame(~.,data.Longi,na.action=na.pass))[,-1],X_KD[(factor.spot+1):ncol(X_KD)])
            else X_KD <- cbind(X_KD[1:(factor.spot-1)],model.matrix(as.formula("~"%+%0%+%"+"%+%paste(vec.factorKD[i], collapse= "+")), model.frame(~.,data.Longi,na.action=na.pass))[,-1])
            
          } }
        
        
        
        
        
        vect.factKD<-names(X_KD)[which(!(names(X_KD)%in%llKD))]
        
        
        #               vect.fact <-apply(matrix(vect.fact,ncol=1,nrow=length(vect.fact)),MARGIN=1,FUN=function(x){
        #               pos1 <- grep("r",unlist(strsplit(x,split="")))[1]+2
        #               pos2 <- grep(")",unlist(strsplit(x,split="")))[1]-1
        #               return(substr(x,start=pos1,stop=pos2))})
        
        occurKD <- rep(0,length(vec.factorKD))
        
        #         for(i in 1:length(vec.factorY)){
        #                #occur[i] <- sum(vec.factor[i] == vect.fact)
        #               occurY[i] <- length(grep(vec.factorY[i],vect.factY))
        #      }
        
        
        
        interaction<-as.vector(apply(matrix(vect.factKD,nrow=length(vect.factKD)),MARGIN=1,FUN=function(x){length(grep(":",unlist(strsplit(x,split=""))))}))
        which.interaction <- which(interaction==1)
        
        for(i in 1:length(vec.factorKD)){
          
          if(length(grep(":",unlist(strsplit(vec.factorKD[i],split=""))))>0){
            
            
            pos <- grep(":",unlist(strsplit(vec.factorKD[i],split="")))
            length.grep <- 0
            for(j in 1:length(vect.factKD)){
              if(j%in%which.interaction){
                
                if(length(grep(substr(vec.factorKD[i],start=1,stop=pos-1),vect.factKD[j]))>0 && length(grep(substr(vec.factorKD[i],start=pos+1,stop=length(unlist(strsplit(vec.factorKD[i],split="")))),vect.factKD[j]))>0){
                  length.grep <- length.grep + 1
                  which <- i}
              }}
            occurKD[i] <- length.grep
            
          }else{
            
            
            if(length(vect.factKD[-which.interaction])>0){occurKD[i] <- length(grep(vec.factorKD[i],vect.factKD[-which.interaction]))
            }else{occurKD[i] <- length(grep(vec.factorKD[i],vect.factKD))}
          }
        }
      }
      
      
      if (ncol(X_KD) == 0){
        noVarKD <- 1
      }else{
        noVarKD <- 0
      }
      
      #========================================>
      
      
      nvarKD<-ncol(X_KD) #nvar==1 correspond a 2 situations:
      
      # au cas ou on a aucune var explicative dans la partie rec, mais X=0
      # cas ou on a 1seul var explicative, ici X est en general different de 0
      
      #varnotdepY <- colnames(X_KD)[-grep("timedep",colnames(X_KD))]
      #vardepY <- colnames(X_KD)[grep("timedep",colnames(X_KD))]
      #vardepY <- apply(matrix(vardepY,ncol=1,nrow=length(vardepY)),1,timedep.names)
      
      #if (length(intersect(varnotdepY,vardepY)) != 0) {
      #  stop("A variable is both used as a constant and time-varying effect covariate")
      #}
      
      #nvartimedepY <- length(vardepY)
      
      #filtretpsY <- rep(0,nvarY)
      #filtretpsY[grep("timedep",colnames(X_KD))] <- 1
      
      varKD <- as.matrix(sapply(X_KD, as.numeric))
      
      nsujety<-nrow(X_KD)
      
      
      
      
      #=======================================>
      #======= Construction du vecteur des indicatrice
      if(length(vec.factorKD) > 0){
        #         ind.place <- ind.place -1
        k <- 0
        for(i in 1:length(vec.factorKD)){
          ind.placeKD[i] <- ind.placeKD[i]+k
          k <- k + occurKD[i]-1
        }
      }
      
    }else{
      noVarKD <- 1 
      nvarKD <- 0 
      varKD <- c()#rep(0, dim(data.Longi)[1])
      vec.factorKD <- NULL
    }
    
    
    
    
    
    # Random effects
    
    if(link=="Random-effects") link0 <- 1
    if(link=="Current-level") link0 <- 2
    
    nRE <- length(random)
    
    ne_re <- nRE#*(nRE+1)/2
    
    
    matzy <- NULL       # here matzy is for time and dose
    matzy <- cbind(data.Longi[,which(colnames(data.Longi)==time.biomarker)],data.Longi[,which(colnames(data.Longi)==dose)])
    
    
    matzy <- as.matrix(matzy)
    
    if(link0==1){netadc <- nRE
    netar <- nRE}
    if(link0==2){netadc <- 1
    netar <- 1}
    
    
    #== Left-censoring ==
    cag <- c(0,0)
    if(!is.null(left.censoring) && left.censoring!=FALSE){
      if(left.censoring<min(Y))stop("The threshold for the left censoring cannot be smaller than the minimal value of the longitudinal outcome")
      cag[1] <- 1
      cag[2] <- left.censoring
      n.censored <- length(which(Y<=left.censoring))
      prop.censored <- n.censored/nsujety
    }
    
    #== Box Cox transformation of the biomarker ==
    boxcox <- c(0,1)
    if(!is.null(BoxCox) && BoxCox!=FALSE){
      boxcox[1] <- 1
      boxcox[2] <- BoxCox
    }
    #============= pseudo-adaptive Gauss Hermite ==============
    #m <- lme(measuret ~ time+interact+treatment, data = data, random = ~ 1| idd)
    random2 <- c()
    for(i in 1:length(random)){
      random2 <- paste(random2,random[i],sep="+")
    }
    random2 <- substring(random2, 2)
    inn <-paste("pdSymm(form=~",random2,")",sep="")
    
    rand <- list(eval(parse(text=inn)))
    names(rand) <- id
    
    
    
    
    
    #==============================================
    #=== preparation survival data ===============
    
    if (!recurrentAG)
    {
      
      tt1T<-aggregate(tt1,by=list(cluster),FUN=sum)[,2]
      tt0T<-rep(0,length(tt1T))
      clusterT <- 0
      ligneT <- 0
      tempT <- 0
      
    }else{
      tt1T<-aggregate(tt1,by=list(cluster),FUN=function(x) x[length(x)])[,2]
      tt0T<-rep(0,length(tt1T))
      clusterT <- 0
      ligneT <- 0
      tempT <- 0
      
    }
    
    TermsT <- terms(formula.terminalEvent, special, data = data)
    
    
    #AD:
    if (!missing(formula.terminalEvent)){
      ord2 <- attr(TermsT, "order")
      
      #   if (length(ord2) & any(ord2 != 1)){
      #     stop("Interaction terms are not valid for terminal event formula")
      #   }
    }
    
    #AD:Joint model needs "terminal()"
    if (!ind.terminal)stop(" Joint frailty model miss specified ")
    
    terminalEvent<-aggregate(terminal,by=list(cluster),FUN=function(x) x[length(x)])[,2]
    # terminalEvent might be 0-1
    if (!all(terminalEvent%in%c(1,0)))  stop("terminal must contain a variable coded 0-1 and a non-factor variable")
    
    m2$formula <- TermsT
    
    m2[[1]] <- as.name("model.frame")
    m2 <- eval(m2, sys.parent()) #ici il prend les colonne associe au var.exp, si rien il prend tout
    
    
    match.noNA<-dimnames(m2)[[1]]%in%dimnames(m)[[1]]#masque logique pour ne pas avoir de NA
    
    m2<-m2[match.noNA, ,drop=FALSE]#m2 inchanger si pas de NA
    
    if (!missing(formula.terminalEvent))newTermsT<-TermsT
    
    #=========================================================>
    if (!missing(formula.terminalEvent)){
      X_T <- model.matrix(newTermsT, m2)
      
      llT <- attr(newTermsT,"term.labels")
      #ind.placedc <- grep("factor",lldc)
      ind.placeT <- unique(attr(X_T,"assign")[duplicated(attr(X_T,"assign"))])#changement unique le 26/09/2014
      
      vec.factorT <- NULL
      vec.factorT <- c(vec.factorT,llT[ind.placeT])
      
      mat.factorT <- matrix(vec.factorT,ncol=1,nrow=length(vec.factorT))
      # Fonction servant a prendre les termes entre "as.factor"
      vec.factorT <-apply(mat.factorT,MARGIN=1,FUN=function(x){
        if (length(grep("factor",x))>0){
          if(length(grep(":",x))>0){
            if(grep('\\(',unlist(strsplit(x,split="")))[1]<grep(":",unlist(strsplit(x,split="")))[1]  && length(grep('\\(',unlist(strsplit(x,split=""))))==1){
              
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
      assign <- lapply(attrassign(X_T, newTermsT)[-1], function(x) x - 1)
      XlevelsT <- .getXlevels(newTermsT, m2)
      contr.saveT <- attr(X_T, 'contrasts')
      #========================================>
      if(length(vec.factorT) > 0){
        positionT <- unlist(assign,use.names=F)
      }
      #========================================>
      
      if (ncol(X_T) == 1)
      {
        X_T<-X_T-1
        noVarT <- 1
      }else{
        X_T <- X_T[, -1, drop = FALSE]
        noVarT <- 0
      }
      
      nvarT <- ncol(X_T)
      
      vartimedepT <- attr(TermsT, "specials")$timedep #nbre de var en fonction de timedep()
      
      # verifier qu'il y ait du timedep dans la deuxieme formule
      if (!is.null(vartimedepT)) timedepT <- 1
      else timedepT <- 0
      if (timedepT==1) stop("The option 'timedep' is not allowed in this model.")
      
      #varnotdepT <- colnames(X_T)[-grep("timedep",colnames(X_T))]
      #vardepT <- colnames(X_T)[grep("timedep",colnames(X_T))]
      #vardepT <- apply(matrix(vardepT,ncol=1,nrow=length(vardepT)),1,timedep.names)
      
      #if (length(intersect(varnotdepT,vardepT)) != 0) {
      #  stop("A variable is both used as a constant and time-varying effect covariate in the formula of terminal event")
      #}
      
      #nvartimedepT <- length(vardepT)
      
      #filtretpsT <- rep(0,nvarT)
      #filtretpsT[grep("timedep",colnames(X_T))] <- 1
      X_T <- X_T[order(data[,colnames(data)==id]),]
      if(nvarT>1)varT.temp<-matrix(c(X_T),nrow=nrow(X_T),ncol=nvarT)
      else varT.temp <- X_T
      
      
      if(is.null(nrow(m2)))
      {
        if (length(m2) != nrow(m)){
          stop(" There are missing values in the covariates modelling the terminal event. \n Prepare data only with complete cases")
        }
      }else{
        
        if (nrow(m2) != nrow(m)){
          stop(" There are missing values in the covariates modelling the terminal event. \n Prepare data only with complete cases")
        }
        
      }
      
      
      if (!is.null(ncol(varT.temp))){
        varT<-aggregate(varT.temp[,1],by=list(cluster), FUN=function(x) x[length(x)])[,2]
        
        if (ncol(varT.temp)>1){
          
          for (i in 2:ncol(varT.temp)){
            varT.i<-aggregate(varT.temp[,i],by=list(cluster), FUN=function(x) x[length(x)])[,2]
            varT<-cbind(varT,varT.i)
          }
        }
      }else{
        varT<-aggregate(varT.temp,by=list(cluster), FUN=function(x) x[length(x)])[,2]
      }
    }else{
      noVarT <- 1
      varT<-0
    }
    
    #    cluster <- aggregate(cluster,by=list(cluster),FUN=function(x) x[1])[,2]
    nvarR<-nvar
    
    if (!missing(formula.terminalEvent)){
      #=======================================>
      #======= Construction du vecteur des indicatrice
      
      if(length(vec.factorT) > 0){
        k <- 0
        for(i in 1:length(vec.factorT)){
          ind.placeT[i] <- ind.placeT[i]+k
          k <- k + occurT[i]-1
        }
      }
      
      #==================================
      
      if(is.null(nrow(varT))){
        nvarEnd<-1
      }else{
        nvarEnd<-ncol(varT)
      }
    }else{
      nvarEnd<-0
    }
    
    if (sum(as.double(var))==0) nvarR <- 0
    if ( sum(as.double(varT))==0) nvarEnd <- 0
    
    ng<-length(uni.cluster)
    
    # ... end preparing data
    
    #
    # Begin JOINT MODEL
    #
    
    
    nvar = nvarR + nvarKG + nvarKD + nvarT
    
    if ((typeof == 0) | (typeof == 2)) indic.nb.int <- 0
    
    if (sum(as.double(varT))==0) nvarT <- 0
    if (sum(as.double(varKG))==0) nvarKG <- 0
    
    if (sum(as.double(varKD))==0) nvarKD <- 0
    
    #if (timedep==0){
    #   npbetatps1 <- 0
    #  npbetatps2 <- 0
    #   npbetatps3 <- 0
    #}else{
    ##  npbetatps1 <- (betaknots+betaorder-1)*nvartimedep
    #  npbetatps2 <- (betaknots+betaorder-1)*nvartimedepT
    #  npbetatps3 <- (betaknots+betaorder-1)*nvartimedepY
    #}
    #npbetatps <- npbetatps1 + npbetatps2 + npbetatps3
    
    effet <- 1
    indic.alpha <- 1
    np <- switch(as.character(typeof),
                 "0"=(2*(as.integer(n.knots) + 2) + as.integer(nvar) + 1 + ne_re + netadc  + netar +   indic.alpha + effet + 4),
                 "2"=(2*2 + nvar + 1  + ne_re + netadc  + netar + indic.alpha + effet + 4))
    
    if (all(all.equal(as.numeric(cens),terminal)==T)){
      stop("'Recurrent event' variable and 'Terminal event' variable need to be different")
    }
    # traitement de l'initialisation du Beta rentre par l'utilisateur
    
    Beta <- rep(0.5,np)
    if (!missing(init.B)) {
      if (length(init.B) != nvar) stop("Wrong number of regression coefficients in init.B")
      # if (timedep) stop("You can hardly know initial regression coefficient while time-varying effect")
      Beta <- c(rep(0.5,np-nvar),init.B)
    }
    
    if (!missing(init.Random)) {
      if (length(init.Random)!=ne_re+1) stop("init.Random must be of length that corresponds to the number of elements to estimate of the B1 matrix")
      Beta[(np -nvar-ne_re+1):(np-nvar)] <- init.Random[-length(init.Random)]
      Beta[np -nvar-1-ne_re-netadc - netar - indic.alpha] <- init.Random[length(init.Random)]
    }
    if (!missing(init.Eta)) {
      
      if (length(init.Eta)!=netadc+netar) stop("init.Eta must be of length that corresponds to the dimension of the link function")
      Beta[(np -nvar-1-ne_re-netadc - netar +1):(np-nvar-1-ne_re)] <- init.Eta
    }
    
    if (!missing(init.Alpha)) {
      if (!is.numeric(init.Alpha)) stop("init.Alpha must be numeric")
      Beta[np -nvar-1-ne_re-netadc - netar ] <- init.Alpha
    }
    
    if (!missing(init.Biomarker)) {
      if (length(init.Biomarker)!=4) stop("init.Biomarker must be a vector of length 4 for the parameters of the biomarker: y_0, K_G0, K_D0 and lambda")
      Beta[np-nvar-ne_re-1-netadc - netar-effet-indic.alpha] <- init.Biomarker[1]   #y0
      Beta[np-nvar-ne_re-1-netadc - netar-effet-indic.alpha-3]<- init.Biomarker[2]  #K_G0
      Beta[np-nvar-ne_re-1-netadc - netar-effet-indic.alpha-2]<- init.Biomarker[3]  #K_D0
      Beta[np-nvar-ne_re-1-netadc - netar-effet-indic.alpha-1]<- init.Biomarker[4]  #lambda
    }
    
    xSuT <- matrix(0,nrow=100,ncol=1)
    xSuR <- matrix(0,nrow=100,ncol=1)
    if (typeof==0){
      mt1 <- size1
      mt2 <- size2
    }else{
      mt1 <- 100
      mt2 <- 100
    }
    
    RE_which_all <- c("y0", "KG", "KD", "lambda")
    RE_which <- which(RE_which_all%in%random)
    
    initialize <- 1
    
    npinit <- switch(as.character(typeof),
                     "0"=((n.knots + 2) + nvarR + effet),
                     "2"=(2 + nvarR + effet))
    
    lappend <- function (lst, ...){
      lst <- c(lst, list(...))
      return(lst)
    }
    nodes <- gauss.quad(n.nodes,kind="hermite")$nodes
    weights <- gauss.quad(n.nodes,kind="hermite")$weights*exp(nodes^2)
    
    nodes2 <- gauss.quad(20,kind="hermite")$nodes
    weights2 <- gauss.quad(20,kind="hermite")$weights*exp(nodes2^2)
    
    tmp <- rep(list(nodes), nRE)
    tmp <- lappend(tmp, nodes2)
    
    nodes <- as.data.frame(do.call(expand.grid,tmp))
    
    tmp <- rep(list(weights), nRE)
    tmp <- lappend(tmp, weights2)
    weights <- as.data.frame(do.call(expand.grid,tmp))
    
    nodes <- sapply(nodes, as.double)
    weights <- sapply(weights, as.double)
    
    paGH <- as.double(matrix(0,nrow = ng, ncol = 2*nRE + 3+((nRE+1)*(nRE+1)-1)/2))
    
    flush.console()
    if (print.times){
      ptm<-proc.time()
      cat("\n")
      cat("Be patient. The program is computing ... \n")
    }
    
    if(GH == 1){
      tmp2 <- rep(list(nodes2), nRE)
      
      nodes2 <- as.data.frame(do.call(expand.grid,tmp2))
      
      tmp2 <- rep(list(weights2), nRE)
      weights2 <- as.data.frame(do.call(expand.grid,tmp2))
      
      nodes2 <- sapply(nodes2, as.double)
      weights2 <- sapply(weights2, as.double)
      
      nvar_uni <- nvarKG+nvarKD
      np_uni <- nvarKG+nvarKD+  nRE +1+4
      b_uni <- rep(0.5, np_uni)
      
      if (!missing(init.B)) {
        b_uni <- c(rep(0.5,np_uni-nvar_uni),init.B)
      }
      
      if (!missing(init.Random)) {
        b_uni[(np_uni -nvar_uni-nRE+1):(np_uni-nvar_uni)] <- init.Random[-length(init.Random)]
      }
      
      
      if (!missing(init.Biomarker)) {
        b_uni[np_uni-nvar_uni-nRE-1] <- init.Biomarker[1]   #y0
        b_uni[np_uni-nvar_uni-nRE-1-3]<- init.Biomarker[2]  #K_G0
        b_uni[np_uni-nvar_uni-nRE-1-2]<- init.Biomarker[3]  #K_D0
        b_uni[np_uni-nvar_uni-nRE-1-1]<- init.Biomarker[4]  #lambda
      }
      
      #  b_uni[	1	]=	1.6778988440443428 
      #  b_uni[	2	]=	  0.77513711109230388   
      #  b_uni[	3	]= -3.4656893965811508 
      #  b_uni[	4	]=	4.2968974314271859 
      #  b_uni[	5	]=	3.1040023897610336 
      #  b_uni[	6	]=	0.54246286557317647  
      #  b_uni[	7	]=	 5.2583198545407918E-002
      
      
      ans <- .Fortran(C_longiuninl,
                      as.integer(nsujety),
                      as.integer(ng),
                      yy0 = as.double(Y),
                      groupey0 = as.integer(clusterY),
                      nb0 = as.integer(nRE),
                      which_random0 = as.integer(random.which),
                      box_cox = as.double(boxcox),
                      matzy0 = as.double(matzy),
                      cag0 = as.double(cag),
                      nva30 = as.integer(nvarKG),
                      nva40 = as.integer(nvarKD),
                      vaxy0 = as.double(cbind(varKG,varKD)),
                      noVar = as.integer(c(noVarKG,noVarKD)),
                      as.integer(maxit),
                      as.integer(np_uni),
                      b_init = as.double(b_uni),
                      H_uni = as.double(matrix(0,nrow=np_uni,ncol=np_uni)),
                      HIH_uni = as.double(matrix(0,nrow=np_uni,ncol=np_uni)),
                      loglik=as.double(0),
                      LCV=as.double(rep(0,2)),
                      counts=as.integer(c(0,0,0)),
                      ier_istop=as.integer(c(0,0)),
                      EPS=as.double(c(LIMparam,LIMlogl,LIMderiv)),
                      GH = as.double(c(0,as.integer(20))),
                      paGH = as.double(matrix(0, nrow = ng,ncol = 2*nRE+1+ (nRE*(nRE-1))/2)),
                      b_pred = as.double(matrix(0,nrow = ng, ncol = 2*nRE + 1+(nRE*(nRE-1))/2)),
                      as.double(weights2),
                      as.double(nodes2),
                      as.integer(20**nRE),
                      as.integer(1),
                      k0_dummy = as.double(c(1,1)))
      
      b_init <- ans$b_init
      paGH <- ans$b_pred
      
      
      if(nvarKG >= 1)Beta[(np-nvar+nvarR+nvarT+1):(np-nvar+nvarR+nvarT+nvarKG)] <- b_init[(np_uni-nvar_uni+1):(np_uni-nvar_uni+nvarKG)]
      if(nvarKD >= 1)Beta[(np-nvar+nvarR+nvarT+nvarKG+1):(np-nvar+nvarR+nvarT+nvarKG+nvarKD)] <- b_init[(np_uni-nvar_uni+nvarKG+1):(np_uni-nvar_uni+nvarKG+nvarKD)]
      
      Beta[np-nvar-nRE+1]=	b_init[np_uni-nvar_uni-nRE+1]
      if(nRE >= 2)Beta[np-nvar-nRE+2]= b_init[np_uni-nvar_uni-nRE+2]
      if(nRE >= 3)Beta[np-nvar-nRE+3] = b_init[np_uni-nvar_uni-nRE+3] 
      if(nRE == 4)Beta[np-nvar-nRE+4] = b_init[np_uni-nvar_uni-nRE+4] 
      
      
      Beta[np-nvar-nRE] =  b_init[np_uni-nvar_uni-nRE] 
      
      
      Beta[np-nvar-nRE-1-netadc - netar-indic.alpha-effet-3] <- b_init[np_uni-nvar_uni-nRE-1-3] 
      Beta[np-nvar-nRE-1-netadc - netar-indic.alpha-effet-2] <- b_init[np_uni-nvar_uni-nRE-1-2]
      Beta[np-nvar-nRE-1-netadc - netar-indic.alpha-effet-1] <- b_init[np_uni-nvar_uni-nRE-1-1] 
      
      Beta[np-nvar-nRE-1-netadc - netar-indic.alpha-effet] <- b_init[np_uni-nvar_uni-nRE-1] 
      
    }
    
    
    
    ans <- .Fortran(C_jointlonginl,
                    as.integer(nsujet),
                    as.integer(nsujety),
                    as.integer(ng),
                    as.integer(n.knots),
                    k0=as.double(kappa), # joint avec generalisation de strate
                    as.double(tt0),
                    as.double(tt1),
                    as.integer(cens),
                    as.integer(cluster),
                    as.double(tt0T),
                    as.double(tt1T),
                    as.integer(terminalEvent),
                    link0 = as.integer(c(link0,link0)),
                    yy0 = as.double(Y),
                    groupey0 = as.integer(clusterY),
                    nb0 = as.integer(nRE),
                    which_random0 = as.integer(random.which),
                    box_cox = as.double(boxcox ),
                    matzy0 =as.double(matzy),
                    cag0 = as.double(cag),
                    as.integer(nvarR),
                    as.double(var),
                    as.integer(nvarT),
                    as.double(varT),
                    nva30 = as.integer(nvarKG),
                    nva40 = as.integer(nvarKD),
                    vaxy0 = as.double(cbind(varKG,varKD)),
                    noVar = as.integer(c(noVarR,noVarT,noVarKG,noVarKD)),
                    ag0 = as.integer(AG),
                    as.integer(maxit),
                    np=as.integer(np),
                    neta0 = as.integer(c(netadc,netar)),
                    b=as.double(Beta),
                    H=as.double(matrix(0,nrow=np,ncol=np)),
                    HIH=as.double(matrix(0,nrow=np,ncol=np)),
                    
                    loglik=as.double(0),
                    LCV=as.double(rep(0,2)),
                    xR=as.double(matrix(0,nrow=size1,ncol=1)),
                    lamR=as.double(matrix(0,nrow=size1,ncol=3)),
                    xSuR=as.double(xSuR),
                    survR=as.double(matrix(0,nrow=mt1,ncol=3)),
                    xD=as.double(rep(0,size2)),
                    lamD=as.double(matrix(0,nrow=size2,ncol=3)),
                    xSuD=as.double(xSuT),
                    survD=as.double(matrix(0,nrow=mt2,ncol=3)),
                    as.integer(typeof),
                    as.integer(equidistant),
                    as.integer(c(size1,size2,mt1,mt2)),###
                    counts=as.integer(c(0,0,0)),
                    ier_istop=as.integer(c(0,0)),
                    paraweib=as.double(rep(0,4)),
                    #  MartinGale=as.double(matrix(0,nrow=ng,ncol=5)),###
                    #   ResLongi = as.double(matrix(0,nrow=nsujety,ncol=4)),
                    #    Pred_y  = as.double(matrix(0,nrow=nsujety,ncol=2)),
                    
                    zi=as.double(rep(0,(n.knots+6))),
                    
                    
                    EPS=as.double(c(LIMparam,LIMlogl,LIMderiv)),
                    GH = c(as.integer(GH),as.integer(n.nodes), as.integer(init.GH)),
                    paGH = as.double(matrix(paGH,nrow = ng, ncol = 2*nRE + 1+nRE*(nRE-1)/2)),
                    b_pred = as.double(matrix(0,nrow = ng, ncol = 2*(nRE+1) + 1+((nRE+1)*nRE)/2)),
                    effet = as.integer(1),
                    indicalpha = as.integer(1),
                    weights =weights,
                    nodes = nodes,
                    as.integer(n.nodes^nRE*20),
                    as.integer(RE_which))
    
    
    
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
    if (noVarT==1 & noVarKG==1 & noVarKD == 1 & noVarR==1) nvar<-0
    #AD:
    
    np <- ans$np
    fit <- NULL
    fit$b <- ans$b
    
    fit$na.action <- attr(m, "na.action")
    fit$call <- call
    fit$n <- nsujet
    fit$groups <- ng
    
    fit$n.events <- ans$counts[2]
    fit$n.deaths <- ans$counts[3]
    fit$n.measurements <- nsujety
    
    if(as.character(typeof)=="0"){
      fit$logLikPenal <- ans$loglik
    }else{
      fit$logLik <- ans$loglik
    }
    #AD:
    fit$LCV <- ans$LCV[1]
    fit$AIC <- ans$LCV[2]
    
    fit$b_pred <- ans$b_pred
    random.effects.pred<- matrix(ans$b_pred, nrow = ng)[,1:(nRE+1)]
    colnames(random.effects.pred) <- c(random,"frailty")
    fit$random.effects.pred <- matrix(ans$b_pred, nrow = ng)[,1:(nRE+1)]
    fit$weights <- ans$weights
    fit$nodes <- ans$nodes
    
    fit$sigma2 <- ans$b[np-nvar-1-ne_re-netadc - netar-indic.alpha]^2
    fit$alpha <- ans$b[np -nvar-1-ne_re-netadc - netar ]
    #random effects
    fit$B1 <- matrix(rep(0,nRE^2),ncol=nRE,nrow=nRE)
    for(j in 1:nRE){
      fit$B1[j,j] <- ans$b[np - nvar - nRE + j]^2
    }
    
    
    fit$ResidualSE <- ans$b[(np  - nvar - ne_re)]^2
    fit$etaR <- ans$b[(np  - nvar - 1 - ne_re - netadc - netar + 1):(np  - nvar - 1 -ne_re - netadc)]
    fit$etaT <- ans$b[(np - nvar - 1 - ne_re - netadc + 1):(np - nvar - 1 -ne_re)]
    
    
    
    
    fit$npar <- np
    
    #AD:
    if ((noVarT==1 & noVarKG==1 & noVarKD == 1)) {
      fit$coef <- NULL
    }
    else
    {
      fit$coef <- ans$b[(np - nvar + 1):np]
      if(is.null(X_KG))names_KG <- c()
      else names_KG <- factor.names(colnames(X_KG))
      if(is.null(X_KD))names_KD <- c()
      else names_KD <- factor.names(colnames(X_KD))
      if(is.null(colnames(X_T)))names_T <- llT
      else names_T<-factor.names(colnames(X_T))
      noms <- c(factor.names(colnames(X)),names_T,names_KG,names_KD)
      #  if (timedep == 1){
      #    while (length(grep("timedep",noms))!=0){
      #      pos <- grep("timedep",noms)[1]
      #     noms <- noms[-pos]
      #      fit$coef <- fit$coef[-(pos:(pos+betaorder+betaknots-1))]
      #    }
      #  }
      names(fit$coef) <- noms
      
      
    }
    
    fit$names.re <- random
    
    temp1 <- matrix(ans$H, nrow = np, ncol = np)
    temp2 <- matrix(ans$HIH, nrow = np, ncol = np)
    
    varH.etaR <- temp1[(np  - nvar - 1 - ne_re - netadc - netar +1):(np  - nvar - 1 - ne_re - netar),
                       (np  - nvar - 1 - ne_re - netadc - netar +1):(np  - nvar - 1 - ne_re - netar )]
    if(netar>1)fit$se.etaR <- sqrt(diag(varH.etaR))
    if(netar==1)fit$se.etaR <- sqrt(varH.etaR)
    varH.etaT<- temp1[(np  - nvar - 1 - ne_re - netadc +1):(np  - nvar - 1 - ne_re),
                      (np  - nvar - 1 - ne_re - netadc +1):(np  - nvar - 1 - ne_re )]
    if(netadc>1)fit$se.etaT <- sqrt(diag(varH.etaT))
    if(netadc==1)fit$se.etaT <- sqrt(varH.etaT)
    fit$se.ResidualSE <- sqrt(temp1[(np  - nvar - ne_re),(np  - nvar - ne_re)])
    fit$varHtotal <- temp1
    fit$varHIHtotal <- temp2
    
    
    fit$varH <- temp1[(np  - nvar +1):np, (np  - nvar +1 ):np]
    fit$varHIH <- temp2[(np  - nvar +1):np, (np  - nvar +1):np]
    noms <- c("alpha","Eta","MeasurementError","B1",factor.names(colnames(X)),names_T,names_KG, names_KD)
    
    
    fit$K_G0 <- ans$b[(np - nvar - ne_re - 1 - netadc - netar-2-3)]
    fit$K_D0 <- ans$b[(np - nvar - ne_re - 1 - netadc - netar-2-2)]
    fit$lambda <- ans$b[(np - nvar - ne_re - 1 - netadc - netar-2-1)]
    fit$y_0 <- ans$b[(np - nvar - ne_re - 1 - netadc - netar-2)]
    
    fit$se.K_G0 <- sqrt(temp1[(np - nvar - ne_re - 1 - netadc - netar-2-3),(np - nvar - ne_re - 1 - netadc - netar-2-3)])
    fit$se.K_D0 <- sqrt(temp1[(np - nvar - ne_re - 1 - netadc - netar-2-2),(np - nvar - ne_re - 1 - netadc - netar-2-2)])
    fit$se.lambda <- sqrt(temp1[(np - nvar - ne_re - 1 - netadc - netar-2-1),(np - nvar - ne_re - 1 - netadc - netar-2-1)])
    fit$se.y_0 <- sqrt(temp1[(np - nvar - ne_re - 1 - netadc - netar-2),(np - nvar - ne_re - 1 - netadc - netar-2)])
    
    fit$alpha_p.value <- 1 - pchisq((fit$alpha/sqrt(diag(fit$varH))[2])^2,1)
    seH.frail <- sqrt(((2 * (fit$sigma2^0.5))^2) *diag(fit$varH)[1]) # delta methode
    fit$sigma2_p.value <- 1 - pnorm(fit$sigma2/seH.frail)
    
    
    if(netadc>1)fit$se.etaT <- sqrt(diag(varH.etaT))
    if(netadc==1)fit$se.etaT <- sqrt(varH.etaT)
    
    if(netar>1)fit$se.etaR <- sqrt(diag(varH.etaR))
    if(netar==1)fit$se.etaR <- sqrt(varH.etaR)
    
    fit$etaR_p.value <- 1 - pchisq((fit$etaR/fit$se.etaR)^2, 1)
    fit$etaT_p.value <- 1 - pchisq((fit$etaT/fit$se.etaT)^2, 1)
    
    fit$y_0_p.value <- 1 - pchisq((fit$y_0/fit$se.y_0)^2, 1) 
    fit$K_G0_p.value <- 1 - pchisq((fit$K_G0/fit$se.K_G0)^2, 1)
    fit$K_D0_p.value <- 1 - pchisq((fit$K_D0/fit$se.K_D0)^2, 1)
    fit$lambda_p.value <- 1 - pchisq((fit$lambda/fit$se.lambda)^2, 1)
    #if (timedep == 1){ # on enleve les variances des parametres des B-splines
    #  while (length(grep("timedep",noms))!=0){
    #    pos <- grep("timedep",noms)[1]
    #    noms <- noms[-pos]
    #    fit$varH <- fit$varH[-(pos:(pos+betaorder+betaknots-1)),-(pos:(pos+betaorder+betaknots-1))]
    #    fit$varHIH <- fit$varHIH[-(pos:(pos+betaorder+betaknots-1)),-(pos:(pos+betaorder+betaknots-1))]
    #  }
    #}
    fit$nvar<-c(nvarR,nvarT,nvarKG,nvarKD)
    #fit$nvarnotdep<-c(nvarR-nvartimedep,nvarT-nvartimedepT,nvarY-nvartimedepY)
    fit$formula <- formula #formula(Terms)
    fit$formula.KG <- formula.KG 
    fit$formula.KD <- formula.KD 
    fit$biomarker <- biomarker
    fit$formula.terminalEvent <- formula.terminalEvent
    
    fit$xR <- matrix(ans$xR, nrow = size1, ncol = 1)
    fit$lamR <- array(ans$lamR, dim = c(size1,3,1))
    fit$xSuR <- matrix(ans$xSuR, nrow = 100, ncol =1)
    fit$survR <- array(ans$survR, dim = c(mt1,3,1))
    
    
    fit$xD <-matrix(ans$xD, nrow = size2, ncol = 1)
    fit$lamD <- matrix(ans$lamD, nrow = size2, ncol = 3)
    fit$xSuD <- matrix(ans$xSuD, nrow = 100, ncol =1)
    fit$survD <- matrix(ans$survD, nrow = mt2, ncol = 3)
    
    fit$link <- link
    fit$type <- type
    fit$n.strat <- 1
    fit$n.iter <- ans$counts[1]
    fit$typeof <- typeof
    if (typeof == 0){
      fit$n.knots<-n.knots
      fit$kappa <- ans$k0#[2]
      fit$n.knots.temp <- n.knots.temp
      fit$zi <- ans$zi
    }
    
    medianR <- NULL
    for (i in (1:fit$n.strat)) medianR[i] <- ifelse(typeof==0, minmin(fit$survR[,1,i],fit$xR), minmin(fit$survR[,1,i],fit$xSuR))
    lowerR <- NULL
    for (i in (1:fit$n.strat)) lowerR[i] <- ifelse(typeof==0, minmin(fit$survR[,2,i],fit$xR), minmin(fit$survR[,2,i],fit$xSuR))
    upperR <- NULL
    for (i in (1:fit$n.strat)) upperR[i] <- ifelse(typeof==0, minmin(fit$survR[,3,i],fit$xR), minmin(fit$survR[,3,i],fit$xSuR))
    fit$medianR <- cbind(lowerR,medianR,upperR)
    
    medianD <- ifelse(typeof==0, minmin(fit$survD[,1],fit$xD), minmin(fit$survD[,1],fit$xSuD))
    lowerD <- ifelse(typeof==0, minmin(fit$survD[,2],fit$xD), minmin(fit$survD[,2],fit$xSuD))
    upperD <- ifelse(typeof==0, minmin(fit$survD[,3],fit$xD), minmin(fit$survD[,3],fit$xSuD))
    fit$medianD <- cbind(lowerD,medianD,upperD)
    
    #AD:
    fit$noVarRec <- noVarR
    fit$noVarEnd <- noVarT
    fit$noVarKG <- noVarKG
    fit$noVarKD <- noVarKD
    
    
    fit$nvarRec <- nvarR
    fit$nvarEnd <- nvarT
    fit$nvarKG <- nvarKG
    fit$nvarKD <- nvarKD
    fit$istop <- ans$ier_istop[2]
    
    fit$shape.weib <- ans$paraweib[1:2]#ans$shape.weib
    fit$scale.weib <- ans$paraweib[3:4]#ans$scale.weib
    
    if(ans$cag0[1]==1)fit$leftCensoring <- TRUE
    if(ans$cag0[1]==0)fit$leftCensoring <- FALSE
    
    if(ans$box_cox[1] == 1){
      fit$BoxCox <- TRUE
      fit$BoxCox_parameter <- ans$box_cox[2]}
    if(ans$box_cox[1] == 0)fit$BoxCox <- FALSE
    
    fit$random.which <- random.which
    
    if(fit$leftCensoring){fit$leftCensoring.threshold <-ans$cag0[2]
    fit$prop.censored <- prop.censored}
    #AD:
    
    # verif que les martingales ont ete bien calculees
    msg <- "Problem in the estimation of the random effects (perhaps high number of events in some clusters)"
    
    
    
    
    fit$AG <- recurrentAG
    
    # fit$BetaTpsMatTY <- matrix(ans$BetaTpsMatY,nrow=101,ncol=1+4*nvartimedepY)
    #  fit$BetaTpsMatR <- matrix(ans$BetaTpsMat,nrow=101,ncol=1+4*nvartimedep)
    #  fit$BetaTpsMatT <- matrix(ans$BetaTpsMatDc,nrow=101,ncol=1+4*nvartimedepT)
    #  fit$nvartimedep <- c(nvartimedep,nvartimedepT,nvartimedepY)
    
    #  fit$Names.vardepR <- vardep
    #  fit$Names.vardepT <- vardepT
    #  fit$Names.vardepY <- vardepY
    
    fit$EPS <- ans$EPS
    
    fit$netar<-netar
    fit$netadc<-netadc
    fit$ne_re <- nRE
    
    
    #================================> For the reccurrent
    #========================= Test de Wald
    
    if ((length(vec.factor) > 0) ){
      Beta <- ans$b[(np - nvar + 1):np]
      VarBeta <- fit$varH#[-c(1:(1+indic.alpha+netar+netadc+1+ne_re))])
      
      nfactor <- length(vec.factor)
      p.wald <- rep(0,nfactor)
      ntot <- nvarR + nvarT + nvarKG + nvarKD
      fit$global_chisqR <- waldtest(N=nvarR,nfact=nfactor,place=ind.place,modality=occur,b=Beta,Varb=VarBeta,Llast=nvarT+nvarKG+nvarKD,Ntot=ntot)
      
      
      fit$dof_chisqR <- occur
      fit$global_chisq.testR <- 1
      # Calcul de pvalue globale
      for(i in 1:length(vec.factor)){
        p.wald[i] <- signif(1 - pchisq(fit$global_chisqR[i], occur[i]), 3)
      }
      fit$p.global_chisqR <- p.wald
      fit$names.factorR <- vec.factor
    }else{
      fit$global_chisq.testR <- 0
    }
    
    #================================> For the longitudinal
    #========================= Test de Wald
    if ((length(vec.factorKG) > 0) ){
      
      Beta <- ans$b[(np - nvar + 1):np]
      
      VarBeta <- fit$varH#[-c(1:(1+indic.alpha+netar+netadc+1+ne_re))])
      
      
      nfactor <- length(vec.factorKG)
      p.wald <- rep(0,nfactor)
      ntot <- nvarR + nvarT + nvarKG + nvarKD
      fit$global_chisqKG <- waldtest(N=nvarKG,nfact=nfactor,place=ind.placeKG,modality=occurKG,b=Beta,Varb=VarBeta,Lfirts=nvarT+nvarR,Llast = nvarKD, Ntot=ntot)
      fit$dof_chisqKG <- occurKG
      fit$global_chisq.testKG <- 1
      # Calcul de pvalue globale
      for(i in 1:length(vec.factorKG)){
        p.wald[i] <- signif(1 - pchisq(fit$global_chisqKG[i], occurKG[i]), 3)
      }
      fit$p.global_chisqKG <- p.wald
      fit$names.factorKG <- vec.factorKG
    }else{
      fit$global_chisq.testKG <- 0
      
    }
    
    if ((length(vec.factorKD) > 0) ){
      
      Beta <- ans$b[(np - nvar + 1):np ]
      
      VarBeta <- fit$varH#[-c(1:(1+indic.alpha+netar+netadc+1+ne_re))])
      
      
      nfactor <- length(vec.factorKD)
      p.wald <- rep(0,nfactor)
      ntot <- nvarR + nvarT + nvarKG + nvarKD
      fit$global_chisqKD <- waldtest(N=nvarKD,nfact=nfactor,place=ind.placeKD,modality=occurKD,b=Beta,Varb=VarBeta,Lfirts=nvarT+nvarR + nvarKG, Ntot=ntot)
      
      fit$dof_chisqKD <- occurKD
      fit$global_chisq.testKD <- 1
      # Calcul de pvalue globale
      for(i in 1:length(vec.factorKD)){
        p.wald[i] <- signif(1 - pchisq(fit$global_chisqKD[i], occurKD[i]), 3)
      }
      fit$p.global_chisqKD <- p.wald
      fit$names.factorKD <- vec.factorKD
    }else{
      fit$global_chisq.testKD <- 0
      
    }
    
    #================================> For the death
    #========================= Test de Wald
    
    if ((length(vec.factorT) > 0) ){
      Beta <- ans$b[(np - nvar  + 1):np]
      
      VarBeta <- fit$varH#[-c(1:(1+indic.alpha+netar+netadc+1+ne_re))])
      
      nfactor <- length(vec.factorT)
      p.wald <- rep(0,nfactor)
      ntot <- nvarR + nvarT + nvarKG + nvarKD
      # print("deces")
      fit$global_chisqT <- waldtest(N=nvarT,nfact=nfactor,place=ind.placeT,modality=occurT,b=Beta,Varb=VarBeta,Llast=nvarKG+nvarKD,Lfirts=nvarR,Ntot=ntot)
      fit$dof_chisqT <- occurT
      fit$global_chisq.testT <- 1
      # Calcul de pvalue globale
      
      for(i in 1:length(vec.factorT)){
        p.wald[i] <- signif(1 - pchisq(fit$global_chisqT[i], occurT[i]), 3)
      }
      
      fit$p.global_chisqT <- p.wald
      fit$names.factorT<- vec.factorT
    }else{
      fit$global_chisq.testT <- 0
    }
    if (!is.null(fit$coef)){
      if(nvar != 1){
        seH <- sqrt(diag(fit$varH))
      }else{
        seH <- sqrt(fit$varH)
      }
      fit$beta_p.value <- 1 - pchisq((fit$coef/seH)^2, 1)
    }
    fit$max_rep <- max_rep
    fit$joint.clust <- 1
    fit$Frailty <- FALSE
    fit$methodGH <- method.GH
    fit$n.nodes <- n.nodes
    fit$time.biomarker <- time.biomarker
    fit$dose <- dose
    #fit$paGH <- data.matrix(cbind(b_lme,invBi_cholDet,as.data.frame(invBi_chol)))
    class(fit) <- "trivPenalNL"
    
    if (print.times){
      cost<-proc.time()-ptm
      cat("The program took", round(cost[3],2), "seconds \n")
    }
    fit
    
  }

