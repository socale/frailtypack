
#' Fit the one-step Joint surrogate model for evaluating a canditate surrogate endpoint
#'
#'@description{
#' \if{html}{\bold{Joint Frailty Surrogate model definition} 
#'
#' Fit the one-step Joint surrogate model for the evaluation of a canditate surrogate endpoint, 
#' with different integration methods on the random effects, using a semiparametric penalized 
#' likelihood estimation. This approach extends that of Burzykowski \code{et al.} (2001) by 
#' including in the same joint frailty model the individual-level and the trial-level random effects.
#'  
#' For the j\out{<sup>th</sup>} subject (j=1,...,n\out{<sub>i</sub>}) of the i\out{<sup>th</sup>} 
#' trial (i=1,...,G), the joint surrogate model is defined as follows:
#' 
#' {\figure{surromodel1.png}{options: width="100\%"}}
#' 
#' where,
#' \eqn{\omega}\out{<sub>ij</sub>} \out{&#126;} \eqn{N}(0,\eqn{\theta}), u\out{<sub>i</sub>} \out{&#126;} \eqn{N}(0,\eqn{\gamma}), \eqn{\omega}\out{<sub>i</sub>} \out{&#8869;} u\out{<sub>i</sub>},
#' u\out{<sub>i</sub>} \out{&#8869;} v\out{<sub>S<sub>i</sub></sub>}, u\out{<sub>i</sub>} \out{&#8869;} v\out{<sub>T<sub>i</sub></sub>}
#' 
#' and 
#' (v\out{<sub>S<sub>i</sub></sub>},v\out{<sub>T<sub>i</sub></sub>})\out{<sup>T</sup>} \out{&#126;} \eqn{N}(0,\eqn{\Sigma}\out{<sub>v</sub>})
#' 
#' with
#' 
#' {\figure{surromodel2.png}{options: width="100\%"}}
#' 
#' In this model, \eqn{\lambda}\out{<sub>0s</sub>}(t) is the baseline hazard function associated with the 
#' surrogate endpoint and \eqn{\beta}\out{<sub>S</sub>} the fixed treatment effect (or log-hazard ratio); 
#' \eqn{\lambda}\out{<sub>0T</sub>}(t) is the baseline hazard function associated with the true endpoint 
#' and \eqn{\beta}\out{<sub>T</sub>} the fixed treatment effect. \eqn{\omega}\out{<sub>ij</sub>} is a shared individual-level frailty that serve to take into account the 
#' heterogeneity in the data at the individual level; u\out{<sub>i</sub>} is a shared frailty effect associated 
#' with the baseline hazard function that serve to take into account the heterogeneity between trials 
#' of the baseline hazard function, associated with the fact that we have several trials in this 
#' meta-analytical design. The power parameters \eqn{\zeta} and \eqn{\alpha} distinguish 
#' both individual and trial-level heterogeneities between the surrogate and the true endpoint. 
#' v\out{<sub>S<sub>i</sub></sub>} and v\out{<sub>T<sub>i</sub></sub>} are two correlated random effects treatment-by-trial interactions. 
#' \eqn{Z}\out{<sub>ij1</sub>} represents the treatment arm to which the patient has been randomized.
#' 
#' \bold{Surrogacy evaluation}
#'
#' We proposed new definitions of Kendall's \eqn{\tau} and coefficient of determination as 
#' individual-level and trial-level association measurements, to evaluate a candidate 
#' surrogate endpoint (Sofeu \emph{et al.}, 2018). The formulations are given below.
#' 
#' \bold{Individual-level surrogacy}
#' 
#' To measure the strength of association between \eqn{S}\out{<sub>ij</sub>} and \eqn{T}\out{<sub>ij</sub>} after 
#' adjusting the marginal distributions for the trial and the treatment effects, as show in 
#' Sofeu \emph{et al.}(2018), we use the Kendall's \eqn{\tau} define by :
#' 
#' {\figure{surromodel3.png}{options: width="100\%"}}
#'        
#'  
#'  where \eqn{\theta}, \eqn{\zeta}, \eqn{\alpha} and \eqn{\gamma} are estimated using the joint surrogate model
#'  defined previously. Kendall's \eqn{\tau} is the difference between the probability of 
#'  concordance and the probability of discordance of two realizations of \eqn{S}\out{<sub>ij</sub>} and \eqn{T}\out{<sub>ij</sub>}. 
#'  It belongs to the interval [-1,1] and assumes a zero value when \eqn{S}\out{<sub>ij</sub>} and \eqn{T}\out{<sub>ij</sub>} are 
#'  independent. We estimate Kendall's \eqn{\tau} using Monte-Carlo or Gaussian Hermite
#'  quadrature integration methods. Its confidence interval is estimated using parametric 
#'  bootstrap
#'  
#'  \bold{Trial-level surrogacy}
#'  
#'  The key motivation for validating a surrogate endpoint is to be able to predict the effect 
#'  of treatment on the true endpoint, based on the observed effect of treatment on the 
#'  surrogate endpoint. As shown by Buyse \emph{et al.} (2000), the coefficenient of 
#'  determination obtains from the covariance matrix \eqn{\Sigma}\out{<sub>v</sub>} of the random effects 
#'  treatment-by-trial interaction can be used to evaluate underlined prediction, and 
#'  therefore as surrogacy evaluation measurement at trial-level. It is defined by: 
#'  
#'  {\figure{surromodel4.png}{options: width="100\%"}}
#'  
#'  The SEs of \eqn{R}\out{<sub>trial</sub>}\out{<sup>2</sup>} is calculated using the Delta-method. We also propose 
#'  \eqn{R}\out{<sub>trial</sub>}\out{<sup>2</sup>} and 95\% CI computed using the parametric bootstrap. The use of delta-method 
#'  can lead to confidence limits violating the [0,1], as noted by 
#'  (Burzykowski \emph{et al.}, 2001). However, using other methods would not significantly alter
#'  the findings of the surrogacy assessment 
#'  }
#'  \if{latex}{\bold{Joint Frailty Surrogate model definition} 
#'
#' Fit the one-step Joint surrogate model for the evaluation of a canditate surrogate endpoint, 
#' with different integration methods on the random effects, using a semiparametric penalized 
#' likelihood estimation. This approach extends that of Burzykowski \code{et al.} (2001) by 
#' including in the same joint frailty model the individual-level and the trial-level random effects.
#'  
#' For the \eqn{j^{th}} subject (\code{j=1,\ldots,}\eqn{n_i}) of the \eqn{i^{th}} 
#' trial \code{i} (\code{i=1,\ldots,}\eqn{G}), the joint surrogate model is defined as follows:
#' 
#' \deqn{\left\{ \begin{array}{ll} \lambda_{S,ij}(t|\omega_{ij},u_i,v_{S_i},{Z_{ij1}}) &= \lambda_{0S}(t) \exp(\omega_{ij}+u_i+
#' v_{S_i} Z_{ij1} +  \beta_SZ_{ij1}) \\
#' \lambda_{T,ij}(t|\omega_{ij},u_i,v_{T_i},{Z_{ij1}}) & = \lambda_{0T}(t) \exp({\zeta}\omega_{ij}+
#' \alpha u_i+v_{T_i} Z_{ij1} + \beta_TZ_{ij1}) \\ \end{array} \right. }
#' 
#' where,
#' 
#' \deqn{ \omega_{ij} \sim N (0,\theta), u_i \sim N (0,\gamma), \omega_{ij} \perp u_i, u_i \perp v_{S_i},
#'  u_i \perp v_{T_i} }
#' 
#' and 
#' \eqn{(v_{S_i},v_{T_i})^{T}\sim\mathcal{N}\left({{0}},\Sigma_{v}\right)}, with
#' \deqn{\Sigma_{v}=\left(
#'                       \begin{array}{cc} 
#'                          \sigma^2_{v_S}  &  \sigma_{v_{ST}} \\ 
#'                          \sigma_{v_{ST}} &  \sigma^2_{v_T}
#'                       \end{array}
#'                  \right)}
#' 
#' In this model, \eqn{\lambda_{0S}(t)} is the baseline hazard function associated with the 
#' surrogate endpoint and \eqn{\beta_S} the fixed treatment effect (or log-hazard ratio); 
#' \eqn{\lambda_{0T}(t)} is the baseline hazard function associated with the true endpoint 
#' and \eqn{\beta_T} the fixed treatment effect. \eqn{\omega_{ij}} is a shared individual-level frailty that serve to take into account the 
#' heterogeneity in the data at the individual level; \eqn{u_i} is a shared frailty effect associated 
#' with the baseline hazard function that serve to take into account the heterogeneity between trials 
#' of the baseline hazard function, associated with the fact that we have several trials in this 
#' meta-analytical design. The power parameters \eqn{\zeta} and \eqn{\alpha} distinguish 
#' both individual and trial-level heterogeneities between the surrogate and the true endpoint. 
#' \eqn{v_{S_i}} and \eqn{v_{T_i}} are two correlated random effects treatment-by-trial interactions. 
#' \eqn{Z_{ij1}} represents the treatment arm to which the patient has been randomized.
#' 
#' \bold{Surrogacy evaluation}
#'
#' We proposed new definitions of Kendall's \eqn{\tau} and coefficient of determination as 
#' individual-level and trial-level association measurements, to evaluate a candidate 
#' surrogate endpoint (Sofeu \emph{et al.}, 2018). The formulations are given below.
#' 
#' \bold{Individual-level surrogacy}
#' 
#' To measure the strength of association between \eqn{S_{ij}} and \eqn{T_{ij}} after 
#' adjusting the marginal distributions for the trial and the treatment effects, as show in 
#' Sofeu \emph{et al.}(2018), we use the Kendall's \eqn{\tau} define by :
#' 
#' \deqn{  \begin{array}{ll}
#'   \tau & = 2\int_{u_{i}}\int_{\omega_{ij}}
#'     \int_{u_{i'}}\int_{\omega_{i'j'}}\{ \\
#'     & \frac{\exp(\omega_{ij}+u_{i}+\zeta \omega_{ij}+\alpha u_{i})+\exp(\omega_{i'j'}+u_{i'} + 
#'     \zeta \omega_{i'j'}+\alpha u_{i'})}{( \exp(\omega_{i'j'}+u_{i'})+ \exp(\omega_{ij}+u_{i})) 
#'     (\exp(\zeta \omega_{i'j'}+\alpha u_{i'}) + \exp(\zeta \omega_{ij}+\alpha u_{i}))} \\
#'    & \frac{1}{\sqrt{2\pi \theta}}\exp\left[-\frac{1}{2}\frac{\omega^2_{i'j'}}{\theta}\right]
#'      \frac{1}{\sqrt{2\pi\gamma}}\exp\left[-\frac{1}{2}\frac{u^2_{i'}}{\gamma}\right]
#'      d\omega_{i'j'}du_{i'} \\
#'    & \frac{1}{\sqrt{2\pi\theta}}\exp\left[-\frac{1}{2}\frac{\omega^2_{ij}}{\theta}\right] \frac{1}{\sqrt{2\pi\gamma}}\exp\left[-\frac{1}{2}\frac{u^2_{i}}{\gamma}\right]
#'      d\omega_{ij}du_{i}\} 
#'      -1
#' \\ \end{array} }
#'        
#'  
#'  where \eqn{\theta, \zeta, \alpha} and \eqn{\gamma} are estimated using the joint surrogate model
#'  defined previously. Kendall's \eqn{\tau} is the difference between the probability of 
#'  concordance and the probability of discordance of two realizations of \eqn{S_{ij}} and \eqn{T_{ij}}. 
#'  It belongs to the interval [-1,1] and assumes a zero value when \eqn{S_{ij}} and \eqn{T_{ij}} are 
#'  independent. We estimate Kendall's \eqn{\tau} using Monte-Carlo or Gaussian Hermite
#'  quadrature integration methods. Its confidence interval is estimated using parametric 
#'  bootstrap
#'  
#'  \bold{Trial-level surrogacy}
#'  
#'  The key motivation for validating a surrogate endpoint is to be able to predict the effect 
#'  of treatment on the true endpoint, based on the observed effect of treatment on the 
#'  surrogate endpoint. As shown by Buyse \emph{et al.} (2000), the coefficenient of 
#'  determination obtains from the covariance matrix \eqn{\Sigma_v} of the random effects 
#'  treatment-by-trial interaction can be used to evaluate underlined prediction, and 
#'  therefore as surrogacy evaluation measurement at trial-level. It is defined by: 
#'  
#'  \deqn{ R^2_{trial}=\frac{\sigma^2_{v_{ST}}}{\sigma^2_{v_S}\sigma^2_{v_T}}
#'  }
#'  
#'  The SEs of \eqn{R^2_{trial}} is calculated using the Delta-method. We also propose 
#'  \eqn{R^2_{trial}} and 95\% CI computed using the parametric bootstrap. The use of delta-method 
#'  can lead to confidence limits violating the [0,1], as noted by 
#'  (Burzykowski \emph{et al.}, 2001). However, using other methods would not significantly alter
#'  the findings of the surrogacy assessment 
#'  }
#' 
#' }
#' 
#' @details{
#' The estimated parameter are obtained using the robust Marquardt algorithm
#' (Marquardt, 1963) which is a combination between a Newton-Raphson algorithm
#' and a steepest descent algorithm. The iterations are stopped when the
#' difference between two consecutive log-likelihoods was small
#' (<
#' \if{html}{10\out{<sup>-3</sup>}} \if{latex}{\eqn{10^{-3}}}), the estimated 
#' coefficients were stable (consecutive
#' values (< \if{html}{10\out{<sup>-3</sup>}} \if{latex}{\eqn{10^{-3}}})), and the gradient 
#' small enough (< \if{html}{10\out{<sup>-3</sup>}} \if{latex}{\eqn{10^{-3}}}), by default.
#' Cubic M-splines of order 4 are used for the hazard function, and I-splines (integrated M-splines) are
#' used for the cumulative hazard function.
#' 
#' The inverse of the Hessian matrix is the variance estimator and to deal
#' with the positivity constraint of the variance component and the spline
#' coefficients, a squared transformation is used and the standard errors are
#' computed by the \eqn{\Delta}-method (Knight & Xekalaki, 2000). The smooth
#' parameter can be chosen by maximizing a likelihood cross validation
#' criterion (Joly and other, 1998). 
#' 
#' We proposed based on the joint surrogate model a new definition 
#' of the Kendall's \eqn{\tau}. Moreover, distinct numerical integration methods are available to approximate the 
#' integrals in the marginal log-likelihood.
#' 
#' \bold{Non-convergence case management procedure}
#' 
#' Special attention must be given to initializing model parameters, the choice of the number of 
#' spline knots, the smoothing parameters and the number of quadrature points to solve convergence 
#' issues. We first initialized parameters using the user's desired strategy, as specified 
#' by the option \code{true.init.val}. When numerical or convergence problems are encountered, 
#' with \code{kappa.use} set to \code{4}, the model is fitted again using a combination of the following strategies: 
#' vary the number of quadrature point (\code{nb.gh} to \code{nb.gh2} or \code{nb.gh2} to \code{nb.gh})
#' in the event of the use of the Gaussian Hermite quadrature integration (see \code{int.method}); 
#' divided or multiplied the smoothing parameters (\if{latex}{\code{k_1}} \if{html}{k\out{<sub>1</sub>}}, 
#' \if{latex}{\code{k_2}} \if{html}{k\out{<sub>2</sub>}}) by 10 or 100 according to 
#' their preceding values, or used parameter vectors obtained during the last iteration (with a 
#' modification of the number of quadrature points and smoothing parameters). Using this strategy, 
#' we usually obtained during simulation the rejection rate less than 3\%. A sensitivity analysis 
#' was conducted without this strategy, and similar results were obtained on the converged samples, 
#' with about a 23\% rejection rate. 
#' }
#' 
#' @aliases jointSurroPenal
#' @usage 
#' jointSurroPenal(data, maxit=40, indicator.zeta = 1, 
#'    indicator.alpha = 1, frail.base = 1, n.knots = 6, 
#'    LIMparam = 0.001, LIMlogl = 0.001, LIMderiv = 0.001, 
#'    nb.mc = 300, nb.gh = 32, nb.gh2 = 20, adaptatif = 0, 
#'    int.method = 2, nb.iterPGH = 5, nb.MC.kendall = 10000, 
#'    nboot.kendall = 1000, true.init.val = 0, 
#'    theta.init = 1, sigma.ss.init = 0.5, sigma.tt.init = 0.5, 
#'    sigma.st.init = 0.48, gamma.init = 0.5, alpha.init = 1, 
#'    zeta.init = 1, betas.init = 0.5, betat.init = 0.5, scale = 1, 
#'    random.generator = 1, kappa.use = 4, random = 0, 
#'    random.nb.sim = 0, seed = 0, init.kappa = NULL, ckappa = c(0,0), 
#'    nb.decimal = 4, print.times = TRUE, print.iter=FALSE)
#'
#' @param data A \code{\link{data.frame}} containing at least seven variables entitled: 
#'    \itemize{
#'    \item{\code{patientID:} A numeric, that represents the patient's identifier and must be unique;}
#'    \item{\code{trialID:} A numeric, that represents the trial in which each patient was randomized;}
#'    \item{\code{timeS:} The follow-up time associated with the surrogate endpoint;}
#'    \item{\code{statusS:} The event indicator associated with the surrogate endpoint. Normally 
#'    0 = no event, 1 = event;}
#'    \item{\code{timeT:} The follow-up time associated with the true endpoint;}
#'    \item{\code{statusT:} The event indicator associated with the true endpoint. Normally 
#'    0 = no event, 1 = event;}
#'    \item{\code{trt:} The treatment indicator for each patient, with 1 = treated, 0 = untreated.}
#'    }
#' @param maxit maximum number of iterations for the Marquardt algorithm.
#' The default being \code{40}. 
#' @param indicator.zeta A binary, indicates whether the power's parameter \eqn{\zeta} should 
#' be estimated (1) or not (0). If \code{0}, \eqn{\zeta} will be set to \code{1} during estimation. 
#' The default is \code{1}. This parameter can be seted to \code{0} in the event of convergence and 
#' identification issues. 
#' @param indicator.alpha A binary, indicating whether the power's parameter \eqn{\alpha} should 
#' be estimated (1) or not (0). If \code{0}, \eqn{\alpha} will be set to \code{1} during estimation.
#' The default is 1.
#' @param frail.base A binary, indicating whether the heterogeneity between trial on the baseline risk 
#' is  considered (\code{1}) or not (\code{0}), using 
#' the shared cluster specific frailties \if{html}{u\out{<sub>i</sub>}} \if{latex}{(\eqn{u_i})}. The default is \code{1}.
#' @param n.knots integer giving the number of knots to use. Value required in
#' the penalized likelihood estimation.  It corresponds to the (n.knots+2)
#' splines functions for the approximation of the hazard or the survival
#' functions.  We estimate I or M-splines of order 4. When the user set a
#' number of knots equals to k (n.knots=k) then the number of interior knots
#' is (k-2) and the number of splines is (k-2)+order.  Number of knots must be
#' between 4 and 20. (See \code{\link{frailtyPenal}} for more details).
#' @param LIMparam Convergence threshold of the Marquardt algorithm for the
#' parameters, \if{html}{10\out{<sup>-3</sup>}} \if{latex}{\eqn{10^{-3}}} by default (See \code{\link{frailtyPenal}} for more details).
#' @param LIMlogl Convergence threshold of the Marquardt algorithm for the
#' log-likelihood, \if{html}{10\out{<sup>-3</sup>}} \if{latex}{\eqn{10^{-3}}} by default (See \code{\link{frailtyPenal}} for more details).
#' @param LIMderiv Convergence threshold of the Marquardt algorithm for the gradient, \if{html}{10\out{<sup>-3</sup>}} \if{latex}{\eqn{10^{-3}}} by default 
#' (See \code{\link{frailtyPenal}} for more details).
#' @param nb.mc Number of samples considered in the Monte-Carlo integration. Required in the event 
#' \code{int.method} is equals to \code{0}, \code{2} or \code{4}. A value between 100 and 300 most often gives 
#' good results. However, beyond 300, the program takes a lot of time to estimate the parameters.
#' The default is \code{300}.
#' @param nb.gh Number of nodes for the Gaussian-Hermite quadrature. It can
#' be chosen among 5, 7, 9, 12, 15, 20 and 32. The default is 32.
#' @param nb.gh2 Number of nodes for the Gauss-Hermite quadrature used to re-estimate the model, 
#' in the event of non-convergence, defined as previously. The default is \code{20}.
#' @param adaptatif A binary, indicates whether the pseudo adaptive Gaussian-Hermite quadrature \code{(1)} or the classical
#' Gaussian-Hermite quadrature \code{(0)} is used. The default is \code{0}.
#' @param int.method A numeric, indicates the integration method: \code{0} for Monte carlo, 
#' \code{1} for Gaussian-Hermite quadrature, \code{2} for a combination of both Gaussian-Hermite quadrature to 
#' integrate over the individual-level random effects and Monte carlo to integrate over the trial-level
#' random effects, \code{4} for a combination of both Monte carlo to integrate over 
#' the individual-level random effects and Gaussian-Hermite quadrature to integrate over the trial-level
#' random effects. The default is \code{2}.
#' @param nb.iterPGH Number of iterations before the re-estimation of the posterior random effects,
#' in the event of the two-steps pseudo-adaptive Gaussian-hermite quadrature. If set to \code{0} there is no 
#' re-estimation". The default is \code{5}.
#' @param nb.MC.kendall Number of generated points used with the Monte-Carlo to estimate
#' integrals in the Kendall's \eqn{\tau} formulation. Beter to use at least 4000 points for
#' stable reseults. The default is \code{10000}.
#' @param nboot.kendall Number of samples considered in the parametric bootstrap to estimate the confidence
#' interval of the Kendall's \eqn{\tau}. The default is \code{1000}. 
#' @param true.init.val Numerical value. Indicates if the given initial values to parameters \code{(0)} should be considered. 
#' If set to \code{2}, \eqn{\alpha} and \eqn{\gamma} are initialised using two separed shared frailty model 
#' (see \code{\link{frailtyPenal}} for more details); \if{html}{\eqn{\sigma}\out{<sup>2</sup><sub>v<sub>S</sub></sub>}, 
#' \eqn{\sigma}\out{<sup>2</sup><sub>v<sub>T</sub></sub>} and
#' \eqn{\sigma}\out{<sub>v<sub>ST</sub></sub>}} 
#' \if{latex}{\eqn{\sigma^2_{v_S}}, \eqn{\sigma^2_{v_T}} and
#' \eqn{\sigma_{v_{ST}}}} are fixed by the user or the default values; \eqn{\zeta}, 
#' \eqn{\theta}, \if{html}{\eqn{\beta}\out{<sub>S</sub>} and \eqn{\beta}\out{<sub>T</sub>}} 
#' \if{latex}{\eqn{\beta_S} and \eqn{\beta_T}} are initialized using a classical joint 
#' frailty model, considering individual level random effects. If the joint frailty model is 
#' faced to convergence issues, \if{html}{\eqn{\beta}\out{<sub>S</sub>} and \eqn{\beta}\out{<sub>T</sub>}} 
#' \if{latex}{\eqn{\beta_S} and \eqn{\beta_T}} are initialized using 
#' two shared frailty models.  In all other scenarios, if the simplified model does not converge,
#' default given parameters values are used. Initial values for spline's associated parameters 
#' are fixed to \code{0.5}. The default for this argument is \code{0}.
#' @param theta.init Initial values for \eqn{\theta}, required if \code{true.init.val} 
#' is set to \code{0} or \code{2}. The default is \code{1}.
#' @param sigma.ss.init Initial values for \if{latex}{\eqn{\sigma^2_{v_S}}} 
#' \if{html}{\eqn{\sigma}\out{<sup>2</sup><sub>v<sub>S</sub></sub>}}, required if \code{true.init.val} 
#' is set to \code{0} or \code{2}. The default is \code{0.5}.
#' @param sigma.tt.init Initial values for \if{latex}{\eqn{\sigma^2_{v_T}}} 
#' \if{html}{\eqn{\sigma}\out{<sup>2</sup><sub>v<sub>T</sub></sub>}}, required if \code{true.init.val} 
#' is set to \code{0} or \code{2}. The default is \code{0.5}.
#' @param sigma.st.init Initial values for\if{latex}{\eqn{\sigma_{v_{ST}}}} 
#' \if{html}{\eqn{\sigma}\out{<sub>v<sub>ST</sub></sub>}}, required if \code{true.init.val} 
#' is set to \code{0} or \code{2}. The default is \code{0.48}.
#' @param gamma.init Initial values for \eqn{\gamma}, required if \code{true.init.val} 
#' is set to \code{0} or \code{2}. The default is \code{0.5}.
#' @param alpha.init Initial values for \eqn{\alpha}, required if \code{true.init.val} 
#' is set to \code{0} or \code{2}. The default is \code{1}.
#' @param zeta.init Initial values for \eqn{\zeta}, required if \code{true.init.val} 
#' is set to \code{0} or \code{2}. The default is \code{1}.
#' @param betas.init Initial values for \if{latex}{\eqn{\beta_S}} \if{html}{\eqn{\beta}\out{<sub>S</sub>}}, required if \code{true.init.val} 
#' is set to \code{0} or \code{2}. The default is \code{0.5}.
#' @param betat.init Initial values for \if{latex}{\eqn{\beta_T}} \if{html}{\eqn{\beta}\out{<sub>T</sub>}}, required if \code{true.init.val} 
#' is set to \code{0} or \code{2}. The default is \code{0.5}.
#' @param scale A numeric that allows to rescale (multiplication) the survival times, to avoid numerical 
#' problems in the event of some convergence issues. If no change is needed the argument is set to 1, the default value. 
#' eg: code{1/365} aims to convert days to years ".
#' @param random.generator Random number generator used by the Fortran compiler, 
#' \code{1} for the intrinsec subroutine \code{Random_number} and \code{2} for the 
#' subroutine \code{uniran()}. The default is \code{1}. in the event of convergence problem 
#' with \code{int.method} set to \code{0}, \code{2} or \code{4}, that requires  
#' integration by Monte-Carlo, user could change the random numbers generator.
#' @param kappa.use A numeric, that indicates how to manage the smoothing parameters \if{latex}{\code{k_1}} \if{html}{k\out{<sub>1</sub>}}  
#' and \if{latex}{\code{k_2}} \if{html}{k\out{<sub>2</sub>}}  in the event of convergence issues. If it is set to \code{1}, 
#' the given smoothing parameters or those obtained by cross-validation are used. 
#' If it is set to \code{3}, the associated smoothing parameters are successively divided by 10, 
#' in the event of convergence issues until 5 times. If it is set to \code{4}, the management of the
#' smoothing parameter is as in the event \code{1}, follows by the successive division as described 
#' in the event \code{3} and preceded by the changing of the number of nodes for the Gauss-Hermite quadrature. 
#' The default is \code{4}.
#' @param random A binary that says if we reset the random number generation with a different environment 
#' at each call \code{(1)} or not \code{(0)}. If it is set to \code{1}, we use the computer clock 
#' as seed. In the last case, it is not possible to reproduce the generated datasets. 
#' The default is \code{0}. Required if \code{random.generator} is set to 1.
#' @param random.nb.sim If \code{random} is set to \code{1}, a binary that indicates the number 
#' of generations that will be made.
#' @param seed The seed to use for data (or samples) generation. required if \code{random} is set to \code{0}. 
#' The default is \code{0}.
#' @param init.kappa smoothing parameter used to penalized the log-likelihood. By default (init.kappa = NULL) the values used 
#' are obtain by cross-validation.
#' @param ckappa Vector of two fixed values to add to the smoothing parameters. By default it is set to (0,0). this argument allows
#' to well manage the smoothing parameters in the event of convergence issues.
#' @param nb.decimal Number of decimal required for results presentation.
#' @param print.times a logical parameter to print estimation time. Default
#' is TRUE.
#' @param print.iter a logical parameter to print iteration process. Default
#' is FALSE.
#' 
#' @return
#' This function return an object of class jointSurroPenal with elements :
#' 
#'    \item{EPS}{A vector containing the obtained convergence thresholds with the Marquardt algorithm,  
#'     for the parameters, the log-likelihood and for the gradient;}
#'    \item{b}{\if{latex}{A vector containing estimates for the splines parameter's; 
#'    the power's parameter \eqn{\zeta} (if \code{indicator.zeta} is set to \code{1}),
#'     the standard error of the shared individual-level frailty \eqn{\omega_{ij}} (\eqn{\theta}),elements of the
#'     lower triangular matrix (L) from the Cholesky decomposition such that \eqn{\Sigma = LL^T}, with \eqn{\Sigma} 
#'     the covariance of the random effects \eqn{(v_{S_i},v_{T_i})}; the coefficient \eqn{\alpha} 
#'     (if \code{indicator.alpha} is set to \code{1}); the satandard error of the random effect \eqn{u_i}; and 
#'     the regression coefficients \eqn{\beta_S} and \eqn{\beta_T};}
#'     \if{html}{A vector containing estimates for the splines parameter's; 
#'     the power's parameter \eqn{\zeta} (if \code{indicator.zeta} is set to \code{1}),
#'     the standard error of the shared individual-level frailty \eqn{\omega}\out{<sub>ij</sub>} (\eqn{\theta}),elements of the
#'     lower triangular matrix (L) from the Cholesky decomposition such that \eqn{\Sigma} = LL\out{<sup>T</sup>}, with \eqn{\Sigma} 
#'     the covariance of the random effects (\out{v<sub>S<sub>i</sub></sub>},\out{v<sub>T<sub>i</sub></sub>}); 
#'     the coefficient \eqn{\alpha} (if \code{indicator.alpha} is set to \code{1}); the satandard error 
#'     of the random effect \code{u}\out{<sub>i</sub>}and the regression coefficients \eqn{\beta}\out{<sub>S</sub>} 
#'     and \eqn{\beta}\out{<sub>T</sub>};}
#'     }
#'     \item{varH}{The variance matrix of all parameters in \code{b} (before positivity constraint transformation 
#'    for the variance of the measurement error, for which the delta method is used);}
#'    \item{varHIH}{The robust estimation of the variance matrix of all parameters in \code{b};}
#'    \item{loglikPenal}{The complete marginal penalized log-likelihood;}
#'    \item{LCV}{the approximated likelihood cross-validation criterion in the semiparametric case (with \code{H}
#'     minus the converged Hessian matrix, and \code{l(.)} the full log-likelihood).
#'    \if{html}{
#'     {\figure{lcv.png}{options: width="50\%"}}}
#'     \if{latex}{\deqn{LCV = \frac{1}{n}(trace(H^{-1}_{pl}H) - l(.))};}}
#'    \item{xS}{vector of times for surrogate endpoint where both survival and hazard function are estimated. 
#'    By default seq(0,max(time),length=99), where time is the vector of survival times;}
#'    \item{lamS}{array (dim = 3) of hazard estimates and confidence bands, for surrogate endpoint;}
#'    \item{survS}{array (dim = 3) of baseline survival estimates and confidence bands, for surrogate endpoint;}
#'    \item{xT}{vector of times for true endpoint where both survival and hazard function are estimated. 
#'    By default seq(0, max(time), length = 99), where time is the vector of survival times;}
#'    \item{lamT}{array (dim = 3) of hazard estimates and confidence bands, for true endpoint;}
#'    \item{survT}{array (dim = 3) of baseline survival estimates and confidence bands, for true endpoint;}
#'    \item{n.iter}{number of iterations needed to converge;}
#'    \item{theta}{Estimate for \eqn{\theta};}
#'    \item{gamma}{Estimate for \eqn{\gamma};}
#'    \item{alpha}{Estimate for \eqn{\alpha};}
#'    \item{zeta}{Estimate for \eqn{\zeta};}
#'    \item{sigma.s}{Estimate for \if{latex}{\eqn{\sigma^2_{v_S}}}\if{html}{\eqn{\sigma}\out{<sup>2</sup><sub>v<sub>S</sub></sub>}};}
#'    \item{sigma.t}{Estimate for \if{latex}{\eqn{\sigma^2_{v_T}}}\if{html}{\eqn{\sigma}\out{<sup>2</sup><sub>v<sub>T</sub></sub>}};}
#'    \item{sigma.st}{Estimate for \if{latex}{\eqn{\sigma_{v_{ST}}}} \if{html}{\eqn{\sigma}\out{<sub>v<sub>ST</sub></sub>}};}
#'    \item{beta.s}{Estimate for \if{latex}{\eqn{\beta_S}} \if{html}{\eqn{\beta}\out{<sub>S</sub>}};}
#'    \item{beta.t}{Estimate for \if{latex}{\eqn{\beta_T}} \if{html}{\eqn{\beta}\out{<sub>T</sub>}};}
#'    \item{ui}{A binary, that indicates if the heterogeneity between trial on the baseline risk 
#'    has been Considered (\code{1}), using the shared cluster specific frailties \if{latex}{(\eqn{u_i})}
#'    \if{html}{(\code{u}\out{<sub>i</sub>})}, 
#'    or not (\code{0});}
#'    \item{ktau}{The Kendall's \eqn{\tau} with the correspondant 95  \eqn{\%} CI computed using the parametric bootstrap;}
#'    \item{R2.boot}{The \if{latex}{\eqn{R^2_{trial}}}
#'    \if{html}{\code{R}\out{<sup>2</sup><sub>trial</sub>}} with the correspondant 95 \eqn{\%} CI computed using the parametric bootstrap;}
#'    \item{Coefficients}{The estimates with the corresponding standard errors and the 95 \eqn{\%} CI}
#'    \item{kappa}{Positive smoothing parameters used for convergence. These values could be different to initial 
#'    values if \code{kappa.use} is set to \code{3} or \code{4};}
#'    \item{scale}{The value used to rescale the survival times}
#'    \item{data}{The dataset used in the model}
#'    \item{varcov.Sigma}{covariance matrix of \if{latex}{(\eqn{\sigma^2_{v_S}},\eqn{\sigma^2_{v_T}},\eqn{\sigma_{v_{ST}}})}
#'    \if{html}{the estimates of (\eqn{\sigma}\out{<sup>2</sup><sub>v<sub>S</sub></sub>},\eqn{\sigma}\out{<sup>2</sup><sub>v<sub>T</sub></sub>}, 
#'    \eqn{\sigma}\out{<sub>v<sub>ST</sub></sub>})}
#'    obtained from the delta-method}
#'    \item{parameter}{list of all arguments used in the model}
#'
#' 
#' @seealso \code{\link{jointSurrSimul}}, \code{\link{summary.jointSurroPenal}}, \code{\link{jointSurroPenalSimul}}
#' 
#' @author Casimir Ledoux Sofeu \email{casimir.sofeu@u-bordeaux.fr}, \email{scl.ledoux@gmail.com} and 
#' Virginie Rondeau \email{virginie.rondeau@inserm.fr}
#' 
#' @references
#' 
#' Burzykowski, T., Molenberghs, G., Buyse, M., Geys, H., and Renard, D. (2001). Validation
#' of surrogate end points in multiple randomized clinical trials with failure time end points. 
#' Journal of the Royal Statistical Society: Series C (Applied Statistics) 50, 405-422.
#' 
#' Buyse, M., Molenberghs, G., Burzykowski, T., Renard, D., and Geys, H. (2000). The validation
#' of surrogate endpoints in meta-analyses of randomized experiments. Biostatistics 1, 49-67
#'
#' Sofeu, C. L., Emura, T., and Rondeau, V. (2019). One-step validation method for surrogate 
#' endpoints using data from multiple randomized cancer clinical trials with failure-time endpoints. 
#' Statistics in Medicine 38, 2928-2942.
#' @export
#' @importFrom doBy orderBy
#'
#' @examples
#' 
#' # Generation of data to use 
#' data.sim <- jointSurrSimul(n.obs=600, n.trial = 30,cens.adm=549.24, 
#'          alpha = 1.5, theta = 3.5, gamma = 2.5, zeta = 1, sigma.s = 0.7, 
#'          sigma.t = 0.7, cor = 0.8, betas = -1.25, betat = -1.25, 
#'          full.data = 0, random.generator = 1, seed = 0, nb.reject.data = 0)
#' 
#' \dontrun{
#' #Surrogacy evaluation based on ganerated data with a combination of Monte Carlo 
#' #and classical Gaussian Hermite integration.*
#' # (Computation takes around 5 minutes)
#' 
#' joint.surro.sim.MCGH <- jointSurroPenal(data = data.sim, int.method = 2, 
#'                    nb.mc = 300, nb.gh = 20)
#'                    
#' #Surrogacy evaluation based on ganerated data with a combination of Monte Carlo 
#' # and Pseudo-adaptive Gaussian Hermite integration.
#' # (Computation takes around 4 minutes)
#' 
#' joint.surro.sim.MCPGH <- jointSurroPenal(data = data.sim, int.method = 2, 
#'                    nb.mc = 300, nb.gh = 20, adaptatif = 1)
#'                    
#' # Results
#' summary(joint.surro.sim.MCGH)
#' summary(joint.surro.sim.MCPGH)
#' 
#' # Data from the advanced ovarian cancer randomized clinical trials.
#' # Joint surrogate model with \eqn{\zeta} fixed to 1, 8 nodes spline 
#' # and the rescaled survival time. 
#' 
#' data(dataOvarian)
#' # (Computation takes around 20 minutes)
#'  
#' joint.surro.ovar <- jointSurroPenal(data = dataOvarian, n.knots = 8, 
#'                 init.kappa = c(2000,1000), indicator.alpha = 0, nb.mc = 200, 
#'                 scale = 1/365)
#' # results
#' summary(joint.surro.ovar)
#' 
#' # data from the adjuvant chemotherapy and resectable gastric cancer 
#' # meta-analyses :
#' # Joint surrogate model with initial values for the parameters and the 
#' # smoothing parameters, and sample for the Monte-Carlo integration
#' # generated by the subroutine \code{uniran}.
#' # (Computation takes around 14 minutes)
#' 
#' data(gastadj)
#' joint.surro.gast <- jointSurroPenal(data = gastadj, nb.mc = 100, nb.gh = 20, 
#'                 indicator.zeta = 0, indicator.alpha = 0, n.knots = 10, 
#'                 random.generator = 2, init.kappa = c(367700100,10025184521))
#'
#' # results
#' summary(joint.surro.gast)
#' 
#' }
#' 
jointSurroPenal = function(data, maxit = 40, indicator.zeta = 1, indicator.alpha = 1, frail.base = 1, 
                      n.knots = 6, LIMparam = 0.001, LIMlogl = 0.001, LIMderiv = 0.001, nb.mc = 300, 
                      nb.gh = 32, nb.gh2 = 20, adaptatif = 0, int.method = 2, nb.iterPGH = 5, 
                      nb.MC.kendall = 10000, nboot.kendall = 1000, true.init.val = 0, theta.init = 1, 
                      sigma.ss.init = 0.5, sigma.tt.init = 0.5, sigma.st.init = 0.48, gamma.init = 0.5, 
                      alpha.init = 1, zeta.init = 1, betas.init = 0.5, betat.init = 0.5, scale = 1, 
                      random.generator = 1, kappa.use = 4, random = 0, random.nb.sim = 0, seed = 0, 
                      init.kappa = NULL, ckappa = c(0,0), nb.decimal = 4, print.times = TRUE, 
                      print.iter = FALSE){
  
 # The initial followup time. The default value is 0
  data$initTime <- 0 
  pfs <- 1 # pfs : used to specified if the time to progression should be censored by the death time (0) or not (1). The default is 1. In this case, death is included in the surrogate endpoint. 
  
 # list of models parameters:
  parameter <- c(maxit = maxit,indicator.zeta = indicator.zeta, indicator.alpha = indicator.alpha,
                 frail.base = frail.base, n.knots = n.knots, LIMparam = LIMparam, LIMlogl = LIMlogl, 
                 LIMderiv = LIMderiv, nb.mc = nb.mc, nb.gh = nb.gh, nb.gh2 = nb.gh2, adaptatif = adaptatif, 
                 int.method = int.method, nb.iterPGH = nb.iterPGH, nb.MC.kendall = nb.MC.kendall, 
                 nboot.kendall = nboot.kendall, true.init.val = true.init.val, theta.init = theta.init, 
                 sigma.ss.init = sigma.ss.init, sigma.tt.init = sigma.tt.init, sigma.st.init = sigma.st.init, 
                 gamma.init = gamma.init, alpha.init = alpha.init, zeta.init = zeta.init, betas.init = betas.init, 
                 betat.init = betat.init, scale = scale, random.generator = random.generator, kappa.use = kappa.use, 
                 random = random, random.nb.sim = random.nb.sim, seed = seed, init.kappa = init.kappa, 
                 nb.decimal = nb.decimal, print.times = print.times, print.iter = print.iter)
  
 # some initializations: for all these parameters, refers to the function jointSurroPenalSimul for help (or descriptions)
  nb.dataset <- 1
  ntrialSimul <- 30
  nbSubSimul <- 1000
  equi.subject <- 1
  prop.rand <- 0.5
  theta2 <- 3.5
  zeta <- 1
  gamma.ui <- 2.5
  alpha.ui <- 1
  betas <- -1.25
  betat <- -1.25
  lambdas <- 1.8
  nus <- 0.0045
  lambdat <- 3
  nut <- 0.0025
  time.cens <- 549
  R2 <- 0.81
  sigma.s <- 0.7
  sigma.t <- 0.7
  param.weibull <- 0
  one.dataset <- 1
  real.data <- 1
  gener.only <- 0
  theta.copule <- 0.5
  
  # end initialization
  
  # ==============parameters checking======================
  if(!(indicator.zeta %in% c(0,1)) | !(indicator.alpha %in% c(0,1)) | !(frail.base %in% c(0,1))){
    stop("model options indicator.zeta, indicator.alpha and frail.base must be set to 0 or 1")
  }
  
  if(!(int.method %in% c(0, 1, 2, 4))){
    stop("The integration method should be specifized by the code: 0, 1, 2 or 4")
  }
  
  if(is.null(data) & nb.dataset == 1){
    stop("data is set to NULL. to do simulation with only one dataset, use the appropriated function for data generation.")
  }
  
  if(!is.null(data) & nb.dataset > 1){
    stop("to do simulation, you don't need to upload the datafram. just set the option data to NULL, or set nb.dataset to 1 if not")
  }
  
  if(!(true.init.val %in% c(0,2))){
    stop("argument 'true.init.val' is a Numeric that requires as value: 0, 2. See the help ...")
  }
  
  if(!(prop.rand %in% c(0,0.5))){
    stop("The argument 'prop.rand' must be set to 0 or 0.5")
  }
  
  if(!(equi.subject %in% c(0,1))){
    stop("The argument 'equi.subject' must be set to 0 or 1")
  }
  
  if(!(random.generator %in% c(1,2))){
    stop("The argument 'random.generator' must be set to 1 or 2")
  }
  if(!(param.weibull %in% c(1,0))){
    stop("The argument 'param.weibull' must be set to 0 or 1")
  }
  
  
  
  
  # ============End parameters checking====================
  
  # ================ data checking=======================
  if(!is.null(data)){
    dataUse <- data
    # dataset's names control
    varStatus=(c("initTime","timeS","statusS","timeT","statusT","trialID","patientID","trt") %in% names(data))
    if(F %in% varStatus){
      stop("Control the names of your variables. They must contain at leat 7 variables named: timeS, statusS, timeT, statusT, trialID, patientID and trt. see the help on this function")
    }
    
    # traitement des donnees
    if(max(table(data$patientID)) > 1){
      stop("Control your dataset. You probably have a duplicate on individual (patientID variable)")
    }
    
    if(!is.numeric(data$timeS)|!is.numeric(data$timeT)|!is.numeric(data$trialID)){
      stop("The variables timeS, timeT and trialID must be numeric") 
    }
    
    if(F %in% c(levels(as.factor(as.character(data$statusS))) %in% c(0,1),levels(as.factor(as.character(data$statusT))) %in% c(0,1),levels(as.factor(as.character(data$trt))) %in% c(0,1))){
      stop("The variables statusS, statusT and trt must be coded 0 or 1")
    }
    if(T %in% c((data$timeT - data$initTime) <= 0, (data$timeS - data$initTime) <= 0)){
      stop("Controll the follow up times of your sujects. the is at leat one subjects with intTime > (timeS or timeT) ")
    }
    
    # if the initTime is not equal to 0 for all individual, we reset the follow up times
    if(F %in% (dataUse$initTime == 0)){
      dataUse$timeT <- dataUse$timeT - dataUse$initTime
      dataUse$timeS <- dataUse$timeS - dataUse$initTime
      dataUse$initTime <- 0
    }
    
    if(!is.null(init.kappa)){
      if(length(init.kappa)!=2){
        stop("init.kappa must have dimension 2 else compute the kappas by cross-validation, using the default null value of 'init.kappa' ")
      }
    }
      kappa0 <- init.kappa
  }
  # ====================== end data checking ==============================
  
  if(!is.null(data)){
    # legere mise en forme des donnees ======Ne marche pas ces traitements, donc je commentaire
    
    # tri du dataframe suivant les essais
    dataUse <- orderBy(~trialID, dataUse)
    
    a <- table(dataUse$trialID)
    trial = rep(1,a[1])
    for(i in 2 : length(a)){
      trial <- c(trial, rep(i, a[i]))
    }
    dataUse$trialID <- trial
    dataUse$patientID <- 1:(nrow(dataUse))
    
    nsujet1 <- nrow(dataUse)
    ng <- nrow(dataUse)
    ntrials1 <- length(table(dataUse$trialID))
  }else{
    nsujet1 <- nbSubSimul
    ng <- nbSubSimul
    ntrials1 <- ntrialSimul
  }
  
  
  maxiter <- maxit
  nst <- 2 # nb de fonction de risque (2 pour )
  
  # indique si l'on estime (1=oui, 0=non) covST
  indice_covST <- 1
  indice_gamma_st <- 0 #  indice_gamma_st: dit si l'on estime gamma_st_ut (1) ou non(0), pour les effets aleatoires correlees sur le risque de base, pas traite ici 
  
  if(frail.base==0) indicator.alpha <- 0 
  
  indice_a_estime <- c(indicator.zeta, indice_covST, indicator.alpha, indice_gamma_st,frail.base)
  
  if(indice_covST == 1){
    # we estimated at least 4 parameters correspondint to the covariance matrix \sigma and the variance of \omega_ij
    nb.frailty <- 4
    nparamfrail <- nb.frailty + indicator.zeta + indicator.alpha + frail.base
  } else{
    nb.frailty <- 3
    nparamfrail <- nb.frailty + indicator.zeta + indicator.alpha + frail.base
  }
  
  # parametre fonction de risque de base
  gamma1 <- 2 # paramertre de la loi gamma
  gamma2 <- 2 # paramertre de la loi gamma
  typeof <- 0 # type de function de risque  0:Splines,  1:Cpm  2:weib
  nbintervR <- 12 # nb interv surrogate
  nbintervDC <- 12 # nb interv true
  equidistant <- 0 # 1=equidist, 0=percentile
  nz <- n.knots
  param_risque_base <- c(typeof,nbintervR,nbintervDC,equidistant,nz)
  
  #nombre de variables explicatives
  ves <- 1 # nombre variables explicative surrogate
  ved <- 1 # nombre variables explicative deces/evenement terminal
  ver <- 1 # nombre total de variables explicative
  nbrevar <- c(ves,ved,ver) 
  
  # vecteur des noms de variables
  nomvarl<- "trt"
  # matrice d'indicatrice de prise en compte des variables explicatives pour le surrogate et le tru
  # filtre = vecteur associe au surrogate
  # filtre2 = vecteur associe au true
  filtre  <- 1
  filtre2 <- 1
  filtre0 <- as.matrix(data.frame(filtre,filtre2))
  
  # gestion de l'affichage a l'ecran
  flush.console()
  if (print.times){
    ptm<-proc.time()
    cat("\n")
    cat("Be patient. The program is computing ... \n")
  }
  
  # nombre de simulatin, 1 par default
  n_sim1 <- nb.dataset
  
  kapa <- "kappa_valid_crois.txt"
  kappa0 <- init.kappa
  if(nb.dataset == 1){
    # jeux de donnees (6 colonnes): donnees pour surrogate et death pour true
    donnees <- dataUse[,c("trialID","patientID","trt","initTime","timeS","statusS")]
    death   <- dataUse[,c("trialID","patientID","trt","initTime","timeT","statusT")]
    # conversion en double des jeux de donneees. je le fais separemment pour distinguer 
    # les cas ou j'aurai plus de variables explicatives pour un des jeux de donnees que pour l'autre
    for(i in 1:ncol(donnees)){
      donnees[,i] <- as.double(donnees[,i])
    }
    for(i in 1:ncol(death)){
      death[,i] <- as.double(death[,i])
    }
    if(is.null(kappa0)){
      if(print.iter) cat("+++++++++++estimation of Kappas by cross-validation +++++++++++")
      # kappas obtenus par validation croisee correspondant sur le jeu de donnees reelles
      #kappa0 <- frailtypack:::kappa_val_croisee(don_S=donnees,don_T=death,njeu=1,n_obs=nsujet1,n_node=n.knots,adjust_S=1,adjust_T=1,kapp_0 = 0)
      kappa0 <- kappa_val_croisee(don_S=donnees,don_T=death,njeu=1,n_obs=nsujet1,n_node=n.knots,adjust_S=1,adjust_T=1,kapp_0 = 0, print.times = F, scale = scale)
      # I deleate the created text file
      file.remove(dir(pattern="kappa_valid_crois.txt"))
    } 
    # else, we use the kappa give by the user
  }else{
    donnees <- NULL
    death   <- NULL
    # simuler les donnees et preparer le fichier des kappa pour la validation croisee
    # nom du fichier pour les kappas obtenues par validation croisee
    kapa <- "kappa_valid_crois.txt"
    
    
  }
  
  # critere de convergence du modele on donne en entree les critere a respecter et en sortie on recupere ceux obtenue du programme
  EPS2 <- c(LIMparam, LIMlogl, LIMderiv)
  
  logNormal <- 1 #lognormal: indique si on a une distribution lognormale des effets aleatoires (1) ou Gamma (0)
  
  # Parametres d'integration
  nsim_node <- rep(NA,11)
  nsim_node[1] <- nb.mc # nombre de simulation pour l'integration par Monte carlo, vaut 0 si on ne veut pas faire du MC
  nsim_node[2] <- nb.gh # nombre de points de quadrature a utiliser (preference 5 points pour l'adaptatice et 32 poits pour la non adaptatice)
  nsim_node[3] <- adaptatif # doit-on faire de l'adaptative(1) ou de la non-adaptative(0)
  nsim_node[4] <- int.method# indique la methode d'integration 0=Monte carlo,1= quadrature, 2=quadrature individuel+MC essai, 3=Laplace, 4= monte carlo individuel + quadrature essai
  nsim_node[5] <- nparamfrail
  nsim_node[6] <- 1 # indique si lon fait de la vectorisation dans le calcul integral (1) ou non (0). rmq: la vectorisation permet de reduire le temps de calcul
  nsim_node[7] <- nb.frailty # indique le nombre d'effet aleatoire cas quadrature adaptative
  type.joint <- 1 # type de modele a estimer: 0=joint classique avec un effet aleatoire partage au niveau individuel,1=joint surrogate avec 1 frailty partage indiv et 2 frailties correles essai
  # 2=joint surrogate sans effet aleatoire partage donc deux effets aleatoires a chaque fois"
  nsim_node[8] <- type.joint 
  nsim_node[9] <- nb.gh2 # nombre de point de quadrature a utiliser en cas de non convergence de prefenrence 7 ou 9 pour la pseudo adaptative et 32 pour la non adaptative
  nsim_node[10] <- nb.iterPGH # nombre d'itteration aubout desquelles reestimer les effects aleatoires a posteriori pour la pseude adaptative. si 0 pas de resestimation
  nsim_node[11] <- 1 # model a utiliser pour la generation des donnee en cas de simulation: 1=joint surrogate avec 1 frailty partage indiv, 3=joint frailty copula model
  nsim_node[12] <- 0 # not used: the copula function: 1 = clayton, 2=Gumbel
 # on adapte le nombre de colonne des paramteres estimes au type de modele
  ncol_param_estim <- 24
  nsim_node[13] <- ncol_param_estim # nobre de colenne de la matrice des parametres estimes: depend du type de modele
  
  # Parametres associes au taux de kendall et au bootstrap
  meth.int.kendal <- 4
  method_int_kendal <- meth.int.kendal# methode d'integration pour le taux de kendall: 0= montecarle, 1= quadrature quaussienne classique, 2= approximation de Laplace
  N_MC_kendall <- nb.MC.kendall # nombre de boucle MC pour le calcul du taux de kendal en approximant l'integrale par montye carlo
  nboot_kendal <- nboot.kendall# nombre d'echantillon bootstrap pour le calcul de l'IC du taux de ke,ndall
  Param_kendall_boot <- c(method_int_kendal,N_MC_kendall,nboot_kendal)
  
  # vecteur contenant les taux de kendall issus du bootstrap
  fichier_kendall <- matrix (0,nrow = n_sim1, ncol = 3)
  fichier_R2 <- matrix (0,nrow = n_sim1, ncol = 3)
  
  # nom des fichiers de sortie
  param_estime <- "Parametre_estime.txt"
  param_empirique <- "Parametre_empirique.txt"
  param_empirique_NC <- "Parametre_empirique_NC.txt"
  tableau_rejet <- "tab_rejet.txt"
  NomFichier <- c(kapa,param_estime,param_empirique,param_empirique_NC,tableau_rejet)
  
  #true.init.val # dit si on initialise les parametres avec les vraies(1) valeurs en cas de 
  #simultation, ou des valeurs donnes par default(0). si 0, on estime 4 model de cox a fragilites partages pour les parametres par default
  
  # Parametres initiaux
  #theta.init		 # valeur initiale de theta2
  #sigma.ss.init	 # valeur initiale de sigma ss
  #sigma.tt.init	 # valeur initiale de sigma tt
  #sigma.st.init  # valeur initiale de sigma st
  #gamma.init  # valeur initiale de gamma.ui
  #alpha.init  # valeur initiale de alpha.ui
  #zeta.init  # valeur initiale de zeta_wij
  #betas.init  # valeur initiale de betas
  #betat.init  # valeur initiale de betat
  vbetast = matrix(c(1,1),nrow = 1, ncol = 2) # juste pour besoin de declaration, n'est pas utilise dans cette fonction
  vbetastinit = matrix(c(1,1),nrow = 1, ncol = 2) # juste pour besoin de declaration, n'est pas utilise dans cette fonction
  if(nb.dataset == 1){
    # jeux de donnees (6 colonnes): donnees pour surrogate et death pour true
    if(true.init.val == 2){ # recherche des parametres initiaux
      if(print.iter) cat("+++++++++++initialization of the parameters using reduced models +++++++++++")
      # # estimation of sigma.s and gamma, using an additive gaussian random effects cox model (Rondeau et al. 2008)
      # cox_surr_sigmaS=try(additivePenal(Surv(initTime, timeS, statusS) ~ cluster(trialID) + trt 
      #                                   + slope (trt), correlation = FALSE, data = donnees, n.knots = nz,
      #                                   kappa=kappa0[1]))
      # estimation of gamma, using a shared frailty model (Rondeau et al. 2003)
      cox_surr_sigmaS=try(frailtyPenal(Surv(timeS, statusS) ~ cluster(trialID) + trt 
                                        , data = donnees, n.knots = nz, kappa=kappa0[1], print.times = F))
      
      # # estimation of sigma.tT and alpha, using an additive gaussian random effects cox model (Rondeau et al. 2008)
      # cox_true_sigmaT=try(additivePenal(Surv(initTime, timeT, statusT) ~ cluster(trialID) + trt
      #                                   + slope (trt), correlation = FALSE, data = death, n.knots = nz,
      #                                   kappa=kappa0[2]))
      # estimation of alpha, using a shared fraily model (Rondeau et al. 2003)
      cox_true_sigmaT=try(frailtyPenal(Surv(timeT, statusT) ~ cluster(trialID) + trt
                                        , data = death, n.knots = nz, kappa=kappa0[2], print.times = F))
      
      donnees_death <- merge(donnees,death[,c("patientID","timeT","statusT")])
      # estimation of eta, theta, beta_S and beta_T using a joint frailty model (Rondeau et al. 2007)
      joint_w=try(frailtyPenal(Surv(timeS,statusS) ~ cluster(patientID) + trt + terminal(statusT),
                                  formula.terminalEvent = ~ trt, RandDist = "LogN", 
                                  data = dataUse, n.knots = nz, kappa = kappa0, print.times = F), silent = TRUE)
      
      if(class(cox_surr_sigmaS)=="try-error"){
        if(print.iter) cat("Estimation problem with the shared frailty model: 
              initialization of gamma using the given default value",fill=T)
      }else{
        gamma.init <- cox_surr_sigmaS$sigma2
#        betas.init <- cox_surr_sigmaS$b[length(cox_surr_sigmaS$b)]
      }
      
      if(class(cox_true_sigmaT)=="try-error"){
#        cat("Estimation problem with the additive gaussian random effects cox model: 
#              initialization of sigma.t using the given default value",fill=T)
      }else{
#        betat.init <- cox_true_sigmaT$b[length(cox_true_sigmaT$b)]
      }
      
      if((class(cox_surr_sigmaS)=="try-error") | (class(cox_true_sigmaT)=="try-error")){
        if(print.iter) cat("initialization of alpha using the given default value",fill=T)
      }else{
        alpha.init <- cox_surr_sigmaS$sigma2/cox_true_sigmaT$sigma2
      }
      
      if(class(joint_w)=="try-error"){
        if((class(cox_surr_sigmaS)=="try-error") & (class(cox_true_sigmaT)=="try-error")){
          if(print.iter) cat("Estimation problem with the joint frailty model: 
              initialization of eta, theta, beta_S and beta_T using the given default values",fill=T)
        }else{
            if((class(cox_surr_sigmaS)!="try-error") & (class(cox_true_sigmaT)!="try-error")){
              if(print.iter) cat("Estimation problem with the joint frailty model: 
                initialization of eta and theta using the given default values",fill=T)
              betas.init <- cox_surr_sigmaS$b[length(cox_surr_sigmaS$b)]
              betat.init <- cox_true_sigmaT$b[length(cox_true_sigmaT$b)]
            }else{
              if(class(cox_surr_sigmaS)=="try-error"){
                if(print.iter) cat("Estimation problem with the joint frailty model: 
                initialization of eta, theta and beta_S using the given default values",fill=T)
                betat.init <- cox_true_sigmaT$b[length(cox_true_sigmaT$b)]
              }else{
                if(print.iter) cat("Estimation problem with the joint frailty model: 
                initialization of eta, theta and beta_T using the given default values",fill=T)
                betas.init <- cox_surr_sigmaS$b[length(cox_surr_sigmaS$b)]
              }
            }
        }
      }else{
        zeta.init  <- joint_w$alpha
        theta.init <- joint_w$sigma2
        betas.init <- joint_w$b[length(joint_w$b)-1]
        betat.init <- joint_w$b[length(joint_w$b)]
      }
      true.init.val <- 0 
      # je remets cette variable a 0 car le programme principal ne connait que 0 et 1. le valeur
      # 2 etait juste pour gener l'initialisation avec les models reduits
      if(print.iter) cat("+++++++++++++++++++End initialization+++++++++++++++++++",fill=T)
    }
  }
    
  param_init <- c(theta.init,sigma.ss.init,sigma.tt.init,sigma.st.init,gamma.init,alpha.init,
                  zeta.init,betas.init,betat.init)
  
  revision_echelle <- scale # coefficient pour la division des temps de suivi. permet de reduire l'echelle des temps pour eviter les problemes numeriques en cas d'un nombre eleve de sujet pour certains cluster
  # random.generator <- # generateur des nombre aleatoire, (1) si Random_number() et (2) si uniran(). Random_number() me permet de gerer le seed
  sujet_equi <- equi.subject # dit si on considere la meme proportion de sujet par essai au moment de la generation des donnee (1) ou non (0). dans ce dernier cas remplir les proportion des sujets par essai dansle fichier dedie
  prop_trait <- prop.rand # proportion des patients traites par esssai. si 0 alors on a des proportions variables suivant les essai, remplir le fichier dedie. si non on aura le meme proportion des traites par essai
  
  # parametres de simultation, en cas de simultation
  # gamma1 <- # parametre de la loi gamma
  # gamma2 <- # parametre de la loi gamma
  # theta2 <- # variance des frailties  w_ij S
  eta <- zeta# Zeta associe a w_ij chez les deces
  # gamma.ui <- # variance des frailties associees au risque de base
  # alpha.ui <- # parametre de puissance associe a u_i chez les deces
  theta2_t <- 0.8# variance des frailties  w_ij chez les T, en cas de fragilite correles
  rsqrt_theta <- 0.8# niveau de correlation entre wij_s et wij_t
  gamma.uit <- 0.8# variance des frailties associees au risque de base sur le true
  rsqrt_gamma.ui <- 0.7# niveau de correlation entre us_i et ut_i
  # betas <- # effet fixe du traitement associe au surrogate
  # betat <- # effet fixe du traitement associe au true
  # parametres d'echelle et de forme pour la weibull
  # lambdas <- # lambdas
  # nus <- # nus
  # lambdat <- #   lambdat
  # nut <- # nut
  mode_cens <- 1 # 1= quantille et 2= date fixe
  temps_cens <- time.cens# censure fixe: temps a preciser
  cens0 <- 0.25 # si quentile, proportion des patients censures
  rsqrt <- sqrt(R2)# niveau de correlation souhaite pour les frailties niveau essai
  # sigma.s <- # variance des effest aleatoires au niveau essai en interaction avec le traitement, associee au surrogate
  # sigma.t <- # variance des effest aleatoires au niveau essai en interaction avec le traitement, associee au true
  paramSimul <- c(gamma1, gamma2, theta2, eta, gamma.ui, alpha.ui, theta2_t, rsqrt_theta, gamma.uit,
                  rsqrt_gamma.ui, betas, betat, lambdas, nus, lambdat,nut, mode_cens, temps_cens,
                  cens0, rsqrt, sigma.s, sigma.t, theta.copule)
  
  # Autres parametres de simulation
  weib <- 1# 0= on simule les temps par une loi exponentielle, 1= on simule par une weibull
  # param.weibull <- # parametrisation de la weibull utilisee: 0= parametrisation par default dans le programme de Virginie, 1= parametrisation a l'aide de la fonction de weibull donnee dans le cous de Pierre
  frailty_cor <- 1 # indique si l'on considere pour le modele de simulation deux effets aleatoire correles au niveau essai(=1) ou un effet aleatoire partage(=0) ou encore on simule sans effet aleatoire au niveau essai(=2, model conjoint classique)
  affiche_stat <- 0 # dit si l'on affiche les statistiques des donnees simulees(1) ou non (0)
  seed_ <- 0 # jeux de donnees a retenir pour la validation croisee
  une_donnee <- one.dataset # pour dire si on simule avec un seul jeu de donnees(1) ou pas (0). ceci pour tester le programme d'estimation
  donne_reel <- real.data # dit si 1 a la question precedente dit s'il sagit du jeux de donnees reel (1) ou non (0)
  #gener.only <- # dit si on voudrait seulement generer les donnees(1) ou generer et faire des simulation(0)
  #kappa.use <- # dit si on utilise un kappa a chaque generation de donnee (1) ou le premier kappa pour tous les jeux de donnees(0)
  decoup_simul <- 0# dans le cas ou l'on a decoupe les simulations en plusieurs paquets, donne le nombre de generation de donnees a ne pas considerer avant d'engager les simulations. ceci empeche de reproduire les meme jeux de donnees pour tous les paquets de simulation. vaut 0 si pas de decoupage pevu sinon pour chaque jeux de simulation mettre cette valeur a jour. Exp si 10 paquets de simul pour un total de 100, on affecte 0 pour le premier paquet, 10 pour le second, 20 pour le 3 ieme, ... 90 pour le 10ieme
  aleatoire <- random# dit si on reinitialise la generation des nombre aleatoire avec un environnement different a chaque appel (1) ou non(O).En cas de generation differente, on utilise l'horloge (heure) de l'ordinateur comme graine. Dans ce cas, il n'est pas possible de reproduire les donnees simulees
  nbre_sim <- random.nb.sim# dans le cas ou aleatoire=1, cette variable indique le nombre de generation qui vont etre faites
  graine <- seed # dans le cas ou l'on voudrait avoir la possibilite de reproduire les donnees generees alors on met la variable aleatoire=0 et on donne dans cette variable la graine a utiliser pour la generation
  autreParamSim <- c(weib,param.weibull,frailty_cor,affiche_stat,seed_,une_donnee,donne_reel,gener.only,
                     kappa.use,decoup_simul,aleatoire,nbre_sim,graine,ckappa[1],ckappa[2],pfs)
  
  # autres dichiers de sortie
  # vecteur des pametres
  nva=2 # deux parametres lies aux effets fixes du traitement
  effet <- 0
  if(typeof==0) np <- 2*(nz+2) + nva + nparamfrail
  if(typeof==1) np <- nbintervDC + nbintervR + nva + nparamfrail
  if(typeof==2) np <- 2*nst + nva + effet + nparamfrail
  
  # matrice hessienne
  H_hessOut <- matrix(0,np,np)
  
  # matrice hessienne corigee
  HIHOut <- matrix(0,np,np)
  resOut <- 0 # log-vraisamblance penalisee
  LCV <- c(0,0) # value de LCV(1) et AIC (2)
  
  # Parametre associe au fonctions de risques et survies de base et splines
  # extrait du Joint fortran
  mt11 <- 100
  mt12 <- 100
  
  if(typeof == 1){
    mt1 <- 3*nbintervR
    mt2 <- 3*nbintervDC
  }
  else{
    mt1 <- 100
    mt2 <- 100
  }
  
  # vecteur donnant la taille des tableau pour fortran
  sizeVect=c(np,mt1,mt2,mt11,mt12)
  
  x1Out <- rep(0, mt1) # temps pour la representation des fonction de risque et de la survies surrogate
  lamOut <- matrix(0, nrow = mt1, ncol= 3) # estimates du risque de base surrogate avec IC
  xSu1 <- matrix(0, mt11)
  suOut <- matrix(0, nrow = mt11 , ncol=3) # matrice des estimates de la survie a baseline(avec IC) surrogate et deces
  x2Out <- rep(0, mt2) # temps pour la representation des fonction de risque et de la survies true endpoint
  lam2Out <- matrix(0, nrow = mt2, ncol= 3) # estimates du risque de base deces avec IC
  xSu2 <- matrix(0, mt11)
  su2Out <- matrix(0, nrow = mt12 , ncol = 3) # matrice des estimates de la survie a baseline(avec IC) surrogate et deces
  
  ni <- 0 # nombre d'itteration pour la convergence
  ier <- 0 # informe sur le comportement du modele(-1 = erreur, k = perte de significativite le modele continu, 0 = pas d'erreur)
  istop <- 0 # critere d'arret: 1= le modele a converge, 2= on a attent le nombre max d'itteration, 3= echec inversion de la hessienne, 4= erreur dans les calculs 
  ziOut <- rep(0,nz+6)  # knots for baseline hazard estimated with splines
  Varcov = matrix(0, nrow = 3, ncol = 3) # matrice de variance-covariance de (sigma_S,sigma_ST,sigmaT) obtenue par delta methode a partir de la hesienne, en raison du changement de variable au moment de l'estimation
  dataHessian <- matrix(0, nrow = np, ncol = np) # sauvegarde des matrices hessiennes des differentes simulations 
  dataHessianIH <- matrix(0, nrow = np*n_sim1, ncol = np)
  datab <- matrix(0, nrow = 1, ncol = np) # sauvegarde des vecteurs de parametres des simulation 
  
  
  # print("quelques parametre")
  # print(dim(H_hessOut))
  # print(sizeVect)
  # stop("pour test")
  # proportion de sujet par essai
  prop_i <- rep(1/ntrials1, ntrials1)
  #proportion de sujet traites par essai
  p <- rep(0.5, ntrials1)
  
  vect_kappa <- matrix(0,nrow = n_sim1, ncol = 2)
  
  if(print.iter) 
    affiche.itter <- 1
  else
    affiche.itter <- 0
  
  param_estimes <- NULL
  
  ans <- .Fortran(C_jointsurrogate,
                  as.integer(nsujet1),
                  as.integer(ng),
                  as.integer(ntrials1),
                  as.integer(maxiter),
                  as.integer(nst),
                  as.integer(nparamfrail),
                  as.integer(indice_a_estime),
                  as.integer(param_risque_base),
                  as.integer(nbrevar),
                  as.integer(filtre0),
                  as.matrix(donnees),
                  as.matrix(death),
                  as.double(p),
                  as.double(prop_i),
                  as.integer(n_sim1),
                  EPS2 = as.double(c(LIMparam, LIMlogl, LIMderiv)),
                  as.double(kappa0),
                  as.double(vect_kappa),
                  as.integer(logNormal),
                  nsim_node = as.integer(nsim_node),
                  as.integer(Param_kendall_boot),
                  as.integer(true.init.val),
                  as.double(param_init),
                  as.double(revision_echelle),
                  as.integer(random.generator),
                  as.integer(sujet_equi),
                  as.double(prop_trait),
                  as.double(paramSimul),
                  as.integer(autreParamSim),
                  fichier_kendall = matrix (0,nrow = 1, ncol = 3), # debut section des parametres de sortie
                  fichier_R2 = matrix (0,nrow = 1, ncol = 3),
                  param_estimes = matrix (0,nrow = 1, ncol = ncol_param_estim),
                  as.integer(sizeVect),
                  b = rep(0,np),
                  H_hessOut = matrix(0,np,np),
                  HIHOut = matrix(0,np,np),
                  resOut = 0,
                  LCV = c(0,0),
                  x1Out = rep(0, mt1),
                  lamOut = matrix(0, nrow = mt1, ncol= 3),
                  xSu1 = matrix(0, mt11),
                  suOut = matrix(0, nrow = mt11 , ncol=3),
                  x2Out = rep(0, mt2),
                  lam2Out = matrix(0, nrow = mt2, ncol= 3),
                  xSu2 = matrix(0, mt11),
                  su2Out = matrix(0, nrow = mt12 , ncol = 3),
                  ni= as.integer(0),
                  ier = 0,
                  istop = 0,
                  ziOut = rep(0,nz+6),
                  as.integer(affiche.itter),
                  Varcov = matrix(0, nrow = 3, ncol = 3),
                  dataHessian = matrix(0, nrow = np, ncol = np),
                  dataHessianIH = matrix(0, nrow = np*n_sim1, ncol = np),
                  datab = matrix(0, nrow = 1, ncol = np),
                  as.double(vbetast),
                  as.double(vbetastinit),
                  PACKAGE="frailtypack"
  )
  
  # cat("nombre d'itteration")
  # cat(ans$ni)
  
  # calcul des critere d'evaluation de la surrogacy
  ans$ktau <- data.frame(ans$fichier_kendall)
  names(ans$ktau) <- c("Ktau","inf.95%CI","sup.95%CI")
  
  ans$R2.boot <- data.frame(ans$fichier_R2)
  names(ans$R2.boot) <- c("R2.boot","inf.95%CI","sup.95%CI")
  
  ans$param_estim <- data.frame(ans$param_estim)[,-c(21:23)]
  
  names(ans$param_estim) <- c("theta","SE_theta","zeta","SE_zeta","beta_S","SE_beta_S","beta_T","SE_beta_T","sigma_s",
                              "SE_sigma_s","sigma_t","SE_sigma.t","sigma_sT","SE_sigma_sT","gamma","SE_gamma",
                              "alpha","SE_alpha","R2trial","SE_R2trial","Ktau")
  # reorganisation des donnees
  param.estim <- ans$param_estim[,c("theta","SE_theta","zeta","SE_zeta","gamma","SE_gamma","alpha","SE_alpha","sigma_s",
                                    "SE_sigma_s","sigma_t","SE_sigma.t","sigma_sT","SE_sigma_sT","beta_S","SE_beta_S",
                                    "beta_T","SE_beta_T", "R2trial","SE_R2trial","Ktau")]
  
  param.estim$SE_Ktau <- NA
  
  if(n_sim1==1){
    # calcul des intervalles de confiance des parametres. on ne le fait que si on a un seul jeu de donnees.
    # en cas de simulation, voir les fonction dediees a la synthese des donnees de simulation
    param.estim2 <- data.frame(param.estim[1,c(1,2)])
    names(param.estim2) <- c("Estimate", "Std. Error")
    row.names(param.estim2) <- names(param.estim)[1]
    
    i <- 3
    k <- 2
    while(i <= ncol(param.estim)){
      a=data.frame(param.estim[1,c(i,i+1)])
      names(a)=c("Estimate", "Std. Error")
      param.estim2 <- rbind(param.estim2, a)
      row.names(param.estim2)[k] <- names(param.estim)[i]
      i <- i + 2
      k <- k + 1
    }
    
    param.estim2$"Inf.95%CI"=param.estim2[,1]-1.96*param.estim2[,2]
    param.estim2$"Sup.95%CI"=param.estim2[,1]+1.96*param.estim2[,2]
    
    param.estim2["Ktau",c("Inf.95%CI","Sup.95%CI")] <- ans$ktau[,c(2,3)]
    
    if(indicator.zeta == 0) param.estim2 <- param.estim2[!(row.names(param.estim2) =="zeta"),]
    if(indicator.alpha == 0) param.estim2 <- param.estim2[!(row.names(param.estim2) == "alpha"),]
    if(frail.base == 0) param.estim2 <- param.estim2[!(row.names(param.estim2) == "gamma"),]
    
    param.estim2 <- round(param.estim2, nb.decimal)
    ans$Coefficients <- param.estim2
  }
  
  # resultats a retourner:
  result <- NULL
  result$EPS <- ans$EPS2 
  result$b <- ans$b
  result$varH <- data.frame(ans$H_hessOut)
  result$varHIH <- ans$HIHOut
  result$loglikPenal <- ans$resOut
  result$LCV <- ans$LCV[1]
  result$xS <- ans$x1Out*scale
  result$lamS <- ans$lamOut
  result$survS <- ans$suOut
  result$xT <- ans$x2Out*scale
  result$lamT <- ans$lam2Out
  result$survT <- ans$su2Out
  result$n.iter <- ans$ni
  result$theta <- param.estim$theta
  result$gamma <- param.estim$gamma
  result$alpha <- param.estim$alpha
  result$zeta <- param.estim$zeta
  result$sigma.s <- param.estim$sigma.s
  result$sigma.t <- param.estim$sigma.t
  result$sigma.st <- param.estim$sigma.sT
  result$beta.s <- param.estim$beta_S
  result$beta.t <- param.estim$beta_T
  result$ui <- frail.base
  result$ktau <- ans$ktau
  result$R2.boot <- ans$R2.boot
  result$Coefficients <- ans$Coefficients
  result$kappa  <- kappa0
  result$scale <- scale
  result$data <- dataUse
  result$varcov.Sigma <- ans$Varcov
  result$parameter <- parameter
  result$type.joint <- type.joint
  #result$dataTkendall <- ans$fichier_kendall
  #result$dataR2boot <- ans$fichier_R2
  
  if(is.na(result$n.iter)) {
    result <- NULL # model did not converge 
    print("===Model did not converge===")
  } 
  
  # =====================++++++++++++++++++++++++++++++++++++++++++++++++++++
  # SI L'ON N'EST PAS EN SIMULATION, PENSER A SUPPRIMER LES FICHIERS CREES+++
  # =====================++++++++++++++++++++++++++++++++++++++++++++++++++++

  # file.remove(c("true.txt", "surrogate.txt", "outjoint"))
  # try(file.remove("OutJoint_Result_surrogate.txt"))
  # try(file.remove("kappa_valid_crois.txt"))
  
  if(!is.null(result)) class(result) <- "jointSurroPenal"
  
  # impression du temps de calcul
  if (print.times){
    cost<-(proc.time()-ptm)/60
    cat("The program took", round(cost[3],2), "minutes \n")
  }
  
  return(result)
}




