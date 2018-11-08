
\name{longiPenal}
\alias{longiPenal}


\title{ Fit a Joint Model for Longitudinal Data and a Terminal Event}

\description{

Fit a joint model for longitudinal data and a terminal event using a semiparametric penalized likelihood estimation or a parametric estimation on the hazard function.

The longitudinal outcomes \eqn{y_i(t_{ik})} (\eqn{k=1,\ldots,n_i}, \eqn{i=1,\ldots,N}) for \eqn{N} subjects are described by a linear mixed model and the risk of the terminal event is represented by a proportional hazard risk model. The joint model is constructed assuming that the processes are linked via a latent structure (Wulfsohn and Tsiatsis 1997):

\deqn{\left\{
\begin{array}{ll}
y_{i}(t_{ik})=\bold{X}_{Li}(t_{ik})^\top \bold{\beta}_L +\bold{ Z}_i(t_{ik})^\top \bold{b}_i + \epsilon_i(t_{ik}) & \mbox{(Longitudinal)} \\
\lambda_i(t|\bold{b}_i)=\lambda_0(t)\exp(\bold{X}_{Ti}(t)\bold{\beta}_T+h(\bold{b}_i,\bold{\beta}_L,\bold{Z}_i(t),\bold{X}_{Li}(t))^\top \bold{\eta}_T ) & \mbox{(Terminal)} \\
 \end{array}
 \right.
}

where \eqn{\bold{X}_{Li}(t)} and \eqn{\bold{X}_{Ti}} are vectors of fixed effects covariates and \eqn{\bold{\beta}_L} and \eqn{\bold{\beta}_T} are
the associated coefficients. Measurements errors \eqn{\epsilon_i(t_{ik})} are iid normally distributed with mean 0 and
variance \eqn{\sigma_{\epsilon}^2}. The random effects \eqn{\bold{b}_i = (b_{0i},\ldots, b_{qi})^\top\sim \mathcal{N}(0,\bold{B}_1)} are associated to covariates \eqn{\bold{Z}_i(t)} and independent from the measurement error. The
relationship between the two processes is explained via \eqn{h(\bold{b}_i,\bold{\beta}_L,\bold{Z}_i(t),\bold{X}_{Li}(t))} with coefficients \eqn{\bold{\eta}_T}. Two forms of the function \eqn{h(\cdot)} are available:
the random effects \eqn{\bold{b}_i} and the current biomarker level \eqn{m_i(t)=\bold{X}_{Li}(t_{ik})^\top \bold{\beta}_L +\bold{ Z}_i(t_{ik})^\top \bold{b}_i}. 

We consider that the longitudinal outcome can be a subject to a quantification limit, i.e.
some observations, below a level of detection \eqn{s} cannot be quantified (left-censoring). 
}

\usage{
longiPenal(formula, formula.LongitudinalData, data,  data.Longi,
      random, id, intercept = TRUE, link = "Random-effects",
      left.censoring = FALSE, n.knots, kappa, maxit = 350, 
      hazard = "Splines", init.B, init.Random, init.Eta, 
      method.GH = "Standard", n.nodes, LIMparam = 1e-3, 
      LIMlogl = 1e-3, LIMderiv = 1e-3, print.times = TRUE)
    }

\arguments{

  \item{formula}{ a formula object, with the response on the left of a \eqn{\sim}
          operator, and the terms on the right. The response must be a
          survival object as returned by the 'Surv' function like in survival package.  Interactions are possible using * or :.}

  \item{formula.LongitudinalData}{a formula object, only requires terms on the right to indicate which variables are modelling the longitudinal outcome. It must follow the standard form used for linear mixed-effects models. Interactions are possible using * or :.}

  \item{data}{a 'data.frame' with the variables used in \code{formula}.}
 \item{data.Longi}{a 'data.frame' with the variables used in \code{formula.LongitudinalData}.}

  \item{random}{ Names of variables for the random effects of the longitudinal outcome. Maximum 3 random effects are possible at the moment. The random intercept is chosen using \code{"1"}.}

  \item{id}{Name of the variable representing the individuals.}
  \item{intercept}{Logical value. Is the fixed intercept of the biomarker included in the mixed-effects model? The default is \code{TRUE}.}
        \item{link}{Type of link function for the dependence between the biomarker and death: \code{"Random-effects"} for the association directly via the random effects of the biomarker, \code{"Current-level"} for the association via the true current level of the biomarker.  The option \code{"Current-level"} can be chosen only if the biomarker random effects are associated with the intercept and time (following this order). The default is \code{"Random-effects"}. }
        \item{left.censoring}{Is the biomarker left-censored below a threshold \eqn{s}? The default is \code{FALSE}, ie. no left-censoring. In case of a left-censored biomarker, this argument must be equal to the threshold \eqn{s}.}
  \item{n.knots}{Integer giving the number of knots to use. Value required in the penalized likelihood estimation.
                It corresponds to the (n.knots+2) splines functions for the approximation of the hazard or the survival functions.
                We estimate I or M-splines of order 4. When the user set a number of knots equals to k (n.knots=k) then the number of interior knots is (k-2) and the number of splines is (k-2)+order.
                Number of knots must be between 4 and 20. (See Note in \code{frailtyPenal} function)}

  \item{kappa}{Positive smoothing parameter in the penalized likelihood estimation. 
          The coefficient kappa of the integral of the squared second derivative of hazard function in the fit (penalized log likelihood). To obtain an initial
          value for \code{kappa}, a solution is to fit the corresponding Cox model using cross validation (See \code{cross.validation} in function \code{frailtyPenal}).
          We advise the user to identify several possible tuning parameters, note their defaults and look at the sensitivity
          of the results to varying them. }

  \item{maxit}{Maximum number of iterations for the Marquardt algorithm. The default is 350.}

  \item{hazard}{Type of hazard functions: \code{"Splines"} for semiparametric hazard functions using equidistant intervals or \code{"Splines-per"} using percentile with the penalized likelihood estimation,
           \code{"Weibull"} for the parametric Weibull functions. The default is \code{"Splines"}.}

	
  \item{init.B}{Vector of initial values for regression coefficients. This vector should be of the same size as the whole vector of covariates with the first elements for the covariates related to the terminal event and then for the covariates related to the biomarker (interactions in the end of each component). Default is 0.5 for each.}

  \item{init.Random}{Initial value for variance of the elements of the matrix of the distribution of the random effects. Default is 0.5 for each element.}

  \item{init.Eta}{Initial values for regression coefficients for the link function. Default is 0.5 for each.}
\item{method.GH}{Method for the Gauss-Hermite quadrature: \code{"Standard"} for the standard non-adaptive Gaussian quadrature, \code{"Pseudo-adaptive"} for the pseudo-adaptive Gaussian quadrature and \code{"HRMSYM"} for the algorithm for the multivariate non-adaptive Gaussian quadrature (see Details). The default is \code{"Standard"}.}
\item{n.nodes}{Number of nodes for the Gauss-Hermite quadrature. They can be chosen among 5, 7, 9, 12, 15, 20 and 32. The default is 9.}
    \item{LIMparam}{Convergence threshold of the Marquardt algorithm for the parameters (see Details of \code{frailtyPenal} function), \eqn{10^{-3}} by default.}

  \item{LIMlogl}{Convergence threshold of the Marquardt algorithm for the log-likelihood (see Details of \code{frailtyPenal} function), \eqn{10^{-3}} by default.}

  \item{LIMderiv}{Convergence threshold of the Marquardt algorithm for the gradient (see Details of \code{frailtyPenal} function), \eqn{10^{-3}} by default.}

  \item{print.times}{a logical parameter to print iteration process. The default is TRUE.}
}

\details{ Typical usage for the joint model
\preformatted{longiPenal(Surv(time,event)~var1+var2, biomarker ~ var1+var2,
   data, data.Longi, ...)}

The method of the Gauss-Hermite quadrature for approximations of the multidimensional integrals, i.e. length of \code{random} is 2, can be chosen among the standard, non-adaptive, pseudo-adaptive in which the quadrature points are transformed using the information from the fitted mixed-effects model for the biomarker (Rizopoulos 2012) or multivariate non-adaptive procedure proposed by Genz et al. 1996 and implemented in FORTRAN subroutine HRMSYM. 
The choice of the method is important for estimations. The standard non-adaptive Gauss-Hermite quadrature (\code{"Standard"}) with a specific number of points gives accurate results but can be time consuming. The non-adaptive procedure (\code{"HRMSYM"}) offers advantageous computational time but in case of datasets in which some individuals have few repeated observations (biomarker measures or recurrent events), this method may be moderately unstable. 
The pseudo-adaptive quadrature uses transformed quadrature points to center and scale the integrand by utilizing estimates of the random effects from an appropriate linear mixed-effects model. This method enables using less quadrature
points while preserving the estimation accuracy and thus lead to a better computational time.

NOTE. Data frames \code{data} and \code{data.Longi} must be consistent. Names and types of corresponding covariates must be the same, as well as the number and identification of individuals. 
}

\value{

The following components are included in a 'longiPenal' object for each model:


\item{b}{The sequence of the corresponding estimation of the coefficients for the hazard functions (parametric or semiparametric), the random effects variances and the regression coefficients.}
\item{call}{The code used for the model.}
\item{formula}{The formula part of the code used for the terminal event part of the model.}
\item{formula.LongitudinalData}{The formula part of the code used for the longitudinal part of the model.}
\item{coef}{The regression coefficients (first for the terminal event and then for the biomarker.}
\item{groups}{The number of groups used in the fit.}
\item{kappa}{The value of the smoothing parameter  in the penalized likelihood estimation corresponding to the baseline hazard function for the terminal event.}
\item{logLikPenal}{The complete marginal penalized log-likelihood in the semiparametric case.}
\item{logLik}{The marginal log-likelihood in the parametric case.}
\item{n.measurements}{The number of biomarker observations used in the fit.}
\item{max_rep}{The maximal number of repeated measurements per individual.}
\item{n.deaths}{The number of events observed in the fit.}
\item{n.iter}{The number of iterations needed to converge.}
\item{n.knots}{The number of knots for estimating the baseline hazard function in the penalized likelihood estimation.}
\item{n.strat}{The number of stratum.}

\item{varH}{The variance matrix of all parameters (before positivity constraint transformation for the variance of the measurement error, for which the delta method is used).}
\item{varHIH}{The robust estimation of the variance matrix of all parameters.}


\item{xD}{The vector of times where both survival and hazard function of the terminal event are estimated. By default seq(0,max(time),length=99), where time is the vector of survival times.}
\item{lamD}{The array (dim=3) of baseline hazard estimates and confidence bands (terminal event).}
\item{survD}{The array (dim=3) of baseline survival estimates and confidence bands (terminal event).}
\item{typeof}{The type of the baseline hazard functions (0:"Splines", "2:Weibull").}
\item{npar}{The number of parameters.}
\item{nvar}{The vector of number of explanatory variables for the terminal event and biomarker.}
\item{nvarEnd}{The number of explanatory variables for the terminal event.}
\item{nvarY}{The number of explanatory variables for the biomarker.}
\item{noVarEnd}{The indicator of absence of the explanatory variables for the terminal event.}
\item{noVarY}{The indicator of absence of the explanatory variables for the biomarker.}
\item{LCV}{The approximated likelihood cross-validation criterion in the semiparametric case (with H minus the converged Hessian matrix, and l(.) the full log-likelihood).\deqn{LCV=\frac{1}{n}(trace(H^{-1}_{pl}H) - l(.))}}
\item{AIC}{The Akaike information Criterion for the parametric case.\deqn{AIC=\frac{1}{n}(np - l(.))}}
\item{n.knots.temp}{The initial value for the number of knots.}
\item{shape.weib}{The shape parameter for the Weibull hazard function.}
\item{scale.weib}{The scale parameter for the Weibull hazard function.}
\item{martingaledeath.res}{The martingale residuals for each individual.}
\item{conditional.res}{The conditional residuals for the biomarker (subject-specific): \eqn{\bold{R}_i^{(m)}=\bold{y}_i-\bold{X}_{Li}^\top\widehat{\bold{\beta}}_L-\bold{Z}_i^\top\widehat{\bold{b}}_i}.}
\item{marginal.res}{The marginal residuals for the biomarker (population averaged): \eqn{\bold{R}_i^{(c)}=\bold{y}_i-\bold{X}_{Li}^\top\widehat{\bold{\beta}}_L}.}
\item{marginal_chol.res}{The Cholesky marginal residuals for the biomarker: \eqn{\bold{R}_i^{(m)}=\widehat{\bold{U}_i^{(m)}}\bold{R}_i^{(m)}}, where \eqn{\widehat{\bold{U}_i^{(m)}}} is an upper-triangular matrix obtained by the Cholesky decomposition of the variance
matrix \eqn{\bold{V}_{\bold{R}_i^{(m)}}=\widehat{\bold{V}_i}-\bold{X}_{Li}(\sum_{i=1}^N\bold{X}_{Li}\widehat{\bold{V}_i}^{-1}\bold{X}_{Li})^{-1}\bold{X}_{Li}^\top}.}
\item{conditional_st.res}{The standardized conditional residuals for the biomarker.}
\item{marginal_st.res}{The standardized marginal residuals for the biomarker.}
\item{random.effects.pred}{ The empirical Bayes predictions of the random effects (ie. using conditional posterior distributions).}
\item{pred.y.marg}{The marginal predictions of the longitudinal outcome.}
\item{pred.y.cond}{The conditional (given the random effects) predictions of the longitudinal outcome.}
\item{lineardeath.pred}{The linear predictor for the terminal part.}
\item{global_chisq_d}{The vector with values of each multivariate Wald test for the terminal part.}
\item{dof_chisq_d}{The vector with degrees of freedom for each multivariate Wald test for the terminal part.}
\item{global_chisq.test_d}{The binary variable equals to 0 when no multivariate Wald is given, 1 otherwise (for the terminal part).}
\item{p.global_chisq_d}{The vector with the p_values for each global multivariate Wald test for the terminal part.}
\item{global_chisq}{The vector with values of each multivariate Wald test for the longitudinal part.}
\item{dof_chisq}{The vector with degrees of freedom for each multivariate Wald test for the longitudinal part.}
\item{global_chisq.test}{The binary variable equals to 0 when no multivariate Wald is given, 1 otherwise (for the longitudinal part).}
\item{p.global_chisq}{The vector with the p_values for each global multivariate Wald test for the longitudinal part.}
\item{names.factordc}{The names of the "as.factor" variables for the terminal part.}
\item{names.factor}{The names of the "as.factor" variables for the longitudinal part.}


\item{intercept}{The logical value. Is the fixed intercept included in the linear mixed-effects model?}
\item{B1}{The variance matrix of the random effects for the longitudinal outcome.}
\item{ResidualSE}{The standard deviation of the measurement error.}
\item{eta}{The regression coefficients for the link function.}
\item{ne_re}{The number of random effects used in the fit.}
\item{names.re}{The names of variables for the random effects.}
\item{link}{The name of the type of the link function.}

\item{eta_p.value}{p-values of the Wald test for the estimated regression coefficients for the link function.}
\item{beta_p.value}{p-values of the Wald test for the estimated regression coefficients.}

\item{leftCensoring}{The logical value. Is the longitudinal outcome left-censored?}
\item{leftCensoring.threshold}{For the left-censored biomarker, the value of the left-censoring threshold used for the fit.} 
\item{prop.censored}{The fraction of observations subjected to the left-censoring.}


 \item{methodGH}{The Guassian quadrature method used in the fit.}
\item{n.nodes}{The number of nodes used for the Gaussian quadrature in the fit.}
  

}



\references{

A. Krol, A. Mauguen, Y. Mazroui, A. Laurent, S. Michiels and V. Rondeau (2017). Tutorial in Joint Modeling and Prediction: A Statistical Software for Correlated Longitudinal Outcomes, 
Recurrent Events and a Terminal Event. \emph{Journal of Statistical Software} \bold{81}(3), 1-52.

A. Krol, L. Ferrer, JP. Pignon, C. Proust-Lima, M. Ducreux, O. Bouche, S. Michiels, V. Rondeau (2016). Joint Model for Left-Censored Longitudinal Data, Recurrent Events and Terminal
Event: Predictive Abilities of Tumor Burden for Cancer Evolution with Application to the
FFCD 2000-05 Trial. \emph{Biometrics} \bold{72}(3) 907-16.


D. Rizopoulos (2012). Fast fitting of joint models for longitudinal and event time data
using a pseudo-adaptive Gaussian quadrature rule. \emph{Computational Statistics and Data
Analysis} \bold{56}, 491-501.

M.S. Wulfsohn, A.A. and Tsiatis, A. A. (1997). A joint model for survival and longitudinal data
measured with error. \emph{Biometrics} \bold{53}, 330-9.

A. Genz and B. Keister (1996). Fully symmetric interpolatory rules for multiple integrals
over infinite regions with Gaussian weight. \emph{Journal of Computational and Applied
Mathematics} \bold{71}, 299-309.


}

\seealso{
\code{\link{plot.longiPenal}},\code{\link{print.longiPenal}},\code{\link{summary.longiPenal}}
}


\examples{

\dontrun{

###--- Joint model for longitudinal data and a terminal event ---###

data(colorectal)
data(colorectalLongi)

# Survival data preparation - only terminal events 
colorectalSurv <- subset(colorectal, new.lesions == 0)

# Baseline hazard function approximated with splines
# Random effects as the link function

model.spli.RE <- longiPenal(Surv(time1, state) ~ age + treatment + who.PS 
+ prev.resection, tumor.size ~  year * treatment + age + who.PS ,
colorectalSurv,	data.Longi = colorectalLongi, random = c("1", "year"),
id = "id", link = "Random-effects", left.censoring = -3.33, 
n.knots = 7, kappa = 2)

# Weibull baseline hazard function
# Current level of the biomarker as the link function

model.weib.CL <- longiPenal(Surv(time1, state) ~ age + treatment + who.PS
+ prev.resection, tumor.size ~  year * treatment + age + who.PS , 
colorectalSurv, data.Longi = colorectalLongi, random = c("1", "year"),
id = "id", link = "Current-level", left.censoring = -3.33, hazard = "Weibull")
}

}

\keyword{models}
