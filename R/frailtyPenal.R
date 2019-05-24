#' Fit a Shared, Joint or Nested Frailty model
#' 
#' @description{
#' 
#' \if{html}{\bold{Shared Frailty model}
#' 
#' Fit a shared gamma or log-normal frailty model using a semiparametric
#' Penalized Likelihood estimation or parametric estimation on the hazard
#' function. Left-truncated, right-censored data, interval-censored data and
#' strata (up to 6 levels) are allowed. It allows to obtain a non-parametric
#' smooth hazard of survival function. This approach is different from the
#' partial penalized likelihood approach of Therneau et al.
#' 
#' The hazard function, conditional on the frailty term \eqn{\omega}\out{<sub>i</sub>}, of a
#' shared gamma frailty model for the j\out{<sup>th</sup>} subject in the i\out{<sup>th</sup>}
#' group: 
#'
#' {\figure{frailtymodel1.1.png}{options: width="70\%"}}
#' {\figure{frailtymodel1.2.png}{options: width="70\%"}}
#' 
#' where \eqn{\lambda}\out{<sub>0</sub>}(t) is the baseline hazard function, \bold{\eqn{\beta}}
#' the vector of the regression coefficient associated to the covariate vector
#' \bold{\eqn{Z}\out{<sub>ij</sub>}} for the j\out{<sup>th</sup>} individual in the i\out{<sup>th</sup>}
#' group.
#' 
#' Otherwise, in case of a shared log-normal frailty model, we have for the
#' j\out{<sup>th</sup>} subject in the i\out{<sup>th</sup>} group:
#' 
#' {\figure{frailtymodel2.1.png}{options: width="70\%"}}
#' {\figure{frailtymodel2.2.png}{options: width="70\%"}}
#' 
#' From now on, you can also consider time-varying effects covariates in your
#' model, see \code{timedep} function for more details.
#' 
#' \bold{Joint Frailty model}
#' 
#' Fit a joint either with gamma or log-normal frailty model for recurrent and
#' terminal events using a penalized likelihood estimation on the hazard
#' function or a parametric estimation. Right-censored data and strata (up to 6
#' levels) for the recurrent event part are allowed. Left-truncated data is not
#' possible. Joint frailty models allow studying, jointly, survival processes
#' of recurrent and terminal events, by considering the terminal event as an
#' informative censoring.
#' 
#' There is two kinds of joint frailty models that can be fitted with
#' \code{frailtyPenal} :
#' 
#' - The first one (Rondeau et al. 2007) includes a common frailty term to the
#' individuals \eqn{\omega}\out{<sub>i</sub>} for the two rates which will take into account
#' the heterogeneity in the data, associated with unobserved covariates. The
#' frailty term acts differently for the two rates (\eqn{\omega}\out{<sub>i</sub>} for the
#' recurrent rate and \eqn{\omega}\out{<sub>i</sub>}\out{<sup>\eqn{\alpha}</sup>} for the death rate). The
#' covariates could be different for the recurrent rate and death rate.
#' 
#' For the j\out{<sup>th</sup>} recurrence (j=1,...,n\out{<sub>i</sub>}) and the
#' i\out{<sup>th</sup>} subject (i=1,...,G), the joint gamma frailty model
#' for recurrent event hazard function \eqn{r}\out{<sub>ij</sub>}(.) and death rate
#' \eqn{\lambda}\out{<sub>i</sub>} is :
#' 
#' {\figure{frailtymodel3.png}{options: width="70\%"}}
#' 
#' where \eqn{r}\out{<sub>0</sub>}(t) (resp. \eqn{\lambda}\out{<sub>0</sub>}(t)) is the recurrent (resp.
#' terminal) event baseline hazard function, \bold{\eqn{\beta}\out{<sub>1</sub>}} (resp.
#' \bold{\eqn{\beta}\out{<sub>2</sub>}}) the regression coefficient vector, \bold{\eqn{Z}\out{<sub>i</sub>}(t)}
#' the covariate vector. The random effects of frailties \eqn{\omega}\out{<sub>i</sub>} \out{<span>&#126;</span>}  \bold{\eqn{\Gamma}(1/\eqn{\theta},1/\eqn{\theta})} and are
#' iid.
#' 
#' The joint log-normal frailty model will be :
#' 
#' {\figure{frailtymodel4.png}{options: width="70\%"}}
#' {\figure{frailtymodel5.png}{options: width="70\%"}}
#' 
#' - The second one (Rondeau et al. 2011) is quite similar but the frailty term
#' is common to the individuals from a same group. This model is useful for the
#' joint modelling two clustered survival outcomes. This joint models have been
#' developed for clustered semi-competing events. The follow-up of each of the
#' two competing outcomes stops when the event occurs. In this case, j is for
#' the subject and i for the cluster.
#' 
#' {\figure{frailtymodel6.png}{options: width="80\%"}}
#' 
#' It should be noted that in these models it is not recommended to include
#' \eqn{\alpha} parameter as there is not enough information to estimate it and
#' thus there might be convergence problems.
#' 
#' In case of a log-normal distribution of the frailties, we will have :
#' 
#' {\figure{frailtymodel7.png}{options: width="80\%"}}
#' {\figure{frailtymodel8.png}{options: width="80\%"}}
#' 
#' This joint frailty model can also be applied to clustered recurrent events
#' and a terminal event (example on "readmission" data below).
#' 
#' From now on, you can also consider time-varying effects covariates in your
#' model, see \code{timedep} function for more details.
#' 
#' There is a possibility to use a weighted penalized maximum likelihood
#' approach for nested case-control design, in which risk set sampling is
#' performed based on a single outcome (Jazic et al., \emph{Submitted}).
#' 
#' General Joint Frailty model Fit a general joint frailty model for recurrent
#' and terminal events considering two independent frailty terms. The frailty
#' term \eqn{u}\out{<sub>i</sub>} represents the unobserved association between recurrences and
#' death. The frailty term \eqn{v}\out{<sub>i</sub>} is specific to the recurrent event rate.
#' Thus, the general joint frailty model is:
#' 
#' {\figure{frailtymodel9.png}{options: width="90\%"}}
#' 
#' where the \eqn{iid} random effects
#' \bold{\eqn{u}\out{<sub>i</sub>}} \out{&#126;} \bold{\eqn{\Gamma}(1/\eqn{\theta},1/\eqn{\theta})} and the
#' \eqn{iid} random effects
#' \bold{\eqn{v}\out{<sub>i</sub>}} \out{&#126;} \bold{\eqn{\Gamma}(1/\eqn{\eta},1/\eqn{\eta})} are independent
#' from each other.  The joint model is fitted using a penalized likelihood
#' estimation on the hazard. Right-censored data and time-varying covariates
#' \bold{\eqn{Z}\out{<sub>i</sub>}(t)} are allowed.
#' 
#' \bold{Nested Frailty model}
#' 
#' \bold{\emph{Data should be ordered according to cluster and subcluster}}
#' 
#' Fit a nested frailty model using a Penalized Likelihood on the hazard
#' function or using a parametric estimation. Nested frailty models allow
#' survival studies for hierarchically clustered data by including two iid
#' gamma random effects. Left-truncated and right-censored data are allowed.
#' Stratification analysis is allowed (maximum of strata = 2).

#' The hazard function conditional on the two frailties \eqn{v}\out{<sub>i</sub>} and
#' \eqn{w}\out{<sub>ij</sub>} for the k\out{<sup>th</sup>} individual of the j\out{<sup>th</sup>} subgroup of
#' the i\out{<sup>th</sup>} group is :
#' 
#' {\figure{frailtymodel10.png}{options: width="80\%"}}
#' 
#' where \eqn{\lambda}\out{<sub>0</sub>}(t) is the baseline hazard function, \eqn{X}\out{<sub>ijk</sub>}
#' denotes the covariate vector and \eqn{\beta} the corresponding vector of
#' regression parameters.
#' 
#' \bold{Joint Nested Frailty Model}
#' 
#' Fit a joint model for recurrent and terminal events using a penalized
#' likelihood on the hazard functions or a parametric estimation.
#' Right-censored data are allowed but left-truncated data and stratified
#' analysis are not allowed.
#' 
#' Joint nested frailty models allow studying, jointly, survival processes of
#' recurrent and terminal events for hierarchically clustered data, by
#' considering the terminal event as an informative censoring and by including
#' two iid gamma random effects.
#' 
#' The joint nested frailty model includes two shared frailty terms, one for
#' the subgroup (\eqn{u}\out{<sub>fi</sub>}) and one for the group (\eqn{w}\out{<sub>f</sub>}) into the
#' hazard functions. This random effects account the heterogeneity in the data,
#' associated with unobserved covariates. The frailty terms act differently for
#' the two rates (\eqn{u}\out{<sub>fi</sub>}, \eqn{w}\out{<sub>f</sub>}\out{<sup>\eqn{\xi}</sup>} for the recurrent rate and
#' \eqn{u}\out{<sub>fi</sub>}\out{<sup>\eqn{\alpha}</sup>}, \eqn{w}\out{<sub>i</sub>} for the terminal event rate). The covariates
#' could be different for the recurrent rate and death rate.
#' 
#' For the j\out{<sup>th</sup>} recurrence (j = 1, ..., \eqn{n}\out{<sub>i</sub>}) of the i\out{<sup>th</sup>}
#' individual (i = 1, ..., \eqn{m}\out{<sub>f</sub>}) of the \eqn{f}\out{<sup>th</sup>} group (f = 1,...,
#' n), the joint nested gamma frailty model for recurrent event hazard function
#' \eqn{r}\out{<sub>fij</sub>}(.) and for terminal event hazard function \eqn{\lambda}\out{<sub>fi</sub>}
#' is :
#' 
#' {\figure{frailtymodel11.png}{options: width="90\%"}}
#' 
#' where \eqn{r}\out{<sub>0</sub>}(resp. \eqn{\lambda}\out{<sub>0</sub>}) is the recurrent (resp.
#' terminal) event baseline hazard function, \eqn{\beta} (resp. \eqn{\gamma})
#' the regression coefficient vector, \bold{\eqn{X}\out{<sub>fij</sub>}}(t) the covariates
#' vector. The random effects are \eqn{\omega}\out{<sub>f</sub>} \out{<span>&#126;</span>}  \bold{\eqn{\Gamma}(1/\eqn{\eta},1/\eqn{\eta})}
#' and \eqn{u}\out{<sub>fi</sub>} \out{<span>&#126;</span>}  \bold{\eqn{\Gamma}(1/\eqn{\theta},1/\eqn{\theta})}.
#' }
#' \if{latex}{\bold{Shared Frailty model}
#' 
#' Fit a shared gamma or log-normal frailty model using a semiparametric
#' Penalized Likelihood estimation or parametric estimation on the hazard
#' function. Left-truncated, right-censored data, interval-censored data and
#' strata (up to 6 levels) are allowed. It allows to obtain a non-parametric
#' smooth hazard of survival function. This approach is different from the
#' partial penalized likelihood approach of Therneau et al.
#' 
#' The hazard function, conditional on the frailty term \eqn{\omega_i}, of a
#' shared gamma frailty model for the \eqn{j^{th}} subject in the \eqn{i^{th}}
#' group:
#' 
#' \deqn{\lambda_{ij}(t|\omega_i)=\lambda_0(t)\omega_i\exp(\bold{\beta^{'}Z_{ij}})}
#' 
#' \deqn{\omega_i\sim\Gamma\left(\frac{1}{\theta},\frac{1}{\theta}\right)
#' \hspace{0.5cm} \bold{E}(\omega_i)=1
#' \hspace{0.5cm}\bold{Var}(\omega_i)=\theta}
#' 
#' where \eqn{\lambda_0(t)} is the baseline hazard function, \eqn{\bold{\beta}}
#' the vector of the regression coefficient associated to the covariate vector
#' \eqn{\bold{Z_{ij}}} for the \eqn{j^{th}} individual in the \eqn{i^{th}}
#' group.
#' 
#' Otherwise, in case of a shared log-normal frailty model, we have for the
#' \eqn{j^{th}} subject in the \eqn{i^{th}} group:
#' 
#' \deqn{\lambda_{ij}(t|\eta_i)=\lambda_0(t)\exp(\eta_i+\bold{\beta^{'}Z_{ij}})}
#' 
#' \deqn{\eta_i\sim N(0,\sigma^2)}
#' 
#' From now on, you can also consider time-varying effects covariates in your
#' model, see \code{timedep} function for more details.
#' 
#' \bold{Joint Frailty model}
#' 
#' Fit a joint either with gamma or log-normal frailty model for recurrent and
#' terminal events using a penalized likelihood estimation on the hazard
#' function or a parametric estimation. Right-censored data and strata (up to 6
#' levels) for the recurrent event part are allowed. Left-truncated data is not
#' possible. Joint frailty models allow studying, jointly, survival processes
#' of recurrent and terminal events, by considering the terminal event as an
#' informative censoring.
#' 
#' There is two kinds of joint frailty models that can be fitted with
#' \code{frailtyPenal} :
#' 
#' - The first one (Rondeau et al. 2007) includes a common frailty term to the
#' individuals \eqn{(\omega_i)} for the two rates which will take into account
#' the heterogeneity in the data, associated with unobserved covariates. The
#' frailty term acts differently for the two rates ( \eqn{\omega_i} for the
#' recurrent rate and \eqn{\omega_i^{\alpha}} for the death rate). The
#' covariates could be different for the recurrent rate and death rate.
#' 
#' For the \eqn{j^{th}}{j^th} recurrence \eqn{(j=1,...,n_i)} and the
#' \eqn{i^{th}}{i^th} subject \eqn{(i=1,...,G)}, the joint gamma frailty model
#' for recurrent event hazard function \eqn{r_{ij}(.)} and death rate
#' \eqn{\lambda_i(.)} is :
#' 
#' \deqn{\left\{ \begin{array}{ll}
#' r_{ij}(t|\omega_i)=\omega_ir_0(t)\exp(\bold{\beta_1^{'}Z_i(t)}) &
#' \mbox{(Recurrent)} \\
#' \lambda_i(t|\omega_i)=\omega_i^{\alpha}\lambda_0(t)\exp(\bold{\beta_2^{'}Z_i(t)})
#' & \mbox{(Death)} \\ \end{array} \right. }
#' 
#' where \eqn{r_0(t)} (resp. \eqn{\lambda_0(t)}) is the recurrent (resp.
#' terminal) event baseline hazard function, \eqn{\bold{\beta_1}} (resp.
#' \eqn{\bold{\beta_2}}) the regression coefficient vector, \eqn{\bold{Z_i(t)}}
#' the covariate vector. The random effects of frailties
#' \eqn{\omega_i\sim\bold{\Gamma}(\frac{1}{\theta},\frac{1}{\theta})} and are
#' iid.
#' 
#' The joint log-normal frailty model will be :
#' 
#' \deqn{\left\{ \begin{array}{ll}
#' r_{ij}(t|\eta_i)=r_0(t)\exp(\eta_i+\bold{\beta_1^{'}Z_i(t)}) &
#' \mbox{(Recurrent)} \\ \lambda_i(t|\eta_i)=\lambda_0(t)\exp(\alpha
#' \eta_i+\bold{\beta_2^{'}Z_i(t)}) & \mbox{(Death)} \\ \end{array} \right. }
#' 
#' where \deqn{\eta_i\sim N(0,\sigma^2)}
#' 
#' - The second one (Rondeau et al. 2011) is quite similar but the frailty term
#' is common to the individuals from a same group. This model is useful for the
#' joint modelling two clustered survival outcomes. This joint models have been
#' developed for clustered semi-competing events. The follow-up of each of the
#' two competing outcomes stops when the event occurs. In this case, j is for
#' the subject and i for the cluster.
#' 
#' \deqn{\left\{ \begin{array}{ll}
#' r_{ij}(t|u_i)=u_ir_0(t)\exp(\bold{\beta_1^{'}Z_{ij}(t)}) & \mbox{(Time to
#' event)} \\
#' \lambda_{ij}(t|u_i)=u_i^{\alpha}\lambda_0(t)\exp(\bold{\beta_2^{'}Z_{ij}(t)})
#' & \mbox{(Death)} \\ \end{array} \right. }
#' 
#' It should be noted that in these models it is not recommended to include
#' \eqn{\alpha} parameter as there is not enough information to estimate it and
#' thus there might be convergence problems.
#' 
#' In case of a log-normal distribution of the frailties, we will have :
#' 
#' \deqn{\left\{ \begin{array}{ll}
#' r_{ij}(t|v_i)=r_0(t)\exp(v_i+\bold{\beta_1^{'}Z_{ij}(t)}) & \mbox{(Time to
#' event)} \\ \lambda_{ij}(t|v_i)=\lambda_0(t)\exp(\alpha
#' v_i+\bold{\beta_2^{'}Z_{ij}(t)}) & \mbox{(Death)} \\ \end{array} \right. }
#' 
#' where \deqn{v_i\sim N(0,\sigma^2)}
#' 
#' This joint frailty model can also be applied to clustered recurrent events
#' and a terminal event (example on "readmission" data below).
#' 
#' From now on, you can also consider time-varying effects covariates in your
#' model, see \code{timedep} function for more details.
#' 
#' There is a possibility to use a weighted penalized maximum likelihood
#' approach for nested case-control design, in which risk set sampling is
#' performed based on a single outcome (Jazic et al., \emph{Submitted}).
#' 
#' General Joint Frailty model Fit a general joint frailty model for recurrent
#' and terminal events considering two independent frailty terms. The frailty
#' term \eqn{u_i} represents the unobserved association between recurrences and
#' death. The frailty term \eqn{v_i} is specific to the recurrent event rate.
#' Thus, the general joint frailty model is:
#' 
#' \eqn{\left\{ \begin{array}{ll}
#' r_{ij}(t|u_i,v_i)=u_iv_ir_0(t)\exp(\bold{\beta_1^{'}Z_{ij}(t)})
#' =u_iv_ir_{ij}(t) & \mbox{(Recurrent)} \\
#' \lambda_{i}(t|u_i)=u_i\lambda_0(t)\exp(\bold{\beta_1^{'}Z_{i}(t)}) = u_i
#' \lambda_{i}(t) & \mbox{(Death)} \\ \end{array} \right. }
#' 
#' where the \eqn{iid} random effects
#' \eqn{\bold{u_i}\sim\Gamma(\frac{1}{\theta},\frac{1}{\theta})} and the
#' \eqn{iid} random effects
#' \eqn{\bold{v_i}\sim\Gamma(\frac{1}{\eta},\frac{1}{\eta})} are independent
#' from each other.  The joint model is fitted using a penalized likelihood
#' estimation on the hazard. Right-censored data and time-varying covariates
#' \eqn{\bold{Z}_i(t)} are allowed.
#' 
#' \bold{Nested Frailty model}
#' 
#' \bold{\emph{Data should be ordered according to cluster and subcluster}}
#' 
#' Fit a nested frailty model using a Penalized Likelihood on the hazard
#' function or using a parametric estimation. Nested frailty models allow
#' survival studies for hierarchically clustered data by including two iid
#' gamma random effects. Left-truncated and right-censored data are allowed.
#' Stratification analysis is allowed (maximum of strata = 2).
#' 
#' The hazard function conditional on the two frailties \eqn{v_i} and
#' \eqn{w_{ij}} for the \eqn{k^{th}} individual of the \eqn{j^{th}} subgroup of
#' the \eqn{i^{th}} group is :
#' 
#' \deqn{\left\{ \begin{array}{ll}
#' \lambda_{ijk}(t|v_i,w_{ij})=v_iw_{ij}\lambda_0(t)exp(\bold{\beta^{'}X_{ijk}})
#' \\ v_i\sim\Gamma\left(\frac{1}{\alpha},\frac{1}{\alpha}\right)
#' \hspace{0.05cm}i.i.d. \hspace{0.2cm} \bold{E}(v_i)=1
#' \hspace{0.2cm}\bold{Var}(v_i)=\alpha \\
#' w_{ij}\sim\Gamma\left(\frac{1}{\eta},\frac{1}{\eta}\right)\hspace{0.05cm}i.i.d.
#' \hspace{0.2cm} \bold{E}(w_{ij})=1 \hspace{0.2cm} \bold{Var}(w_{ij})=\eta
#' \end{array} \right. }
#' 
#' where \eqn{\lambda_0(t)} is the baseline hazard function, \eqn{X_{ijk}}
#' denotes the covariate vector and \eqn{\beta} the corresponding vector of
#' regression parameters.
#' 
#' \bold{Joint Nested Frailty Model}
#' 
#' Fit a joint model for recurrent and terminal events using a penalized
#' likelihood on the hazard functions or a parametric estimation.
#' Right-censored data are allowed but left-truncated data and stratified
#' analysis are not allowed.
#' 
#' Joint nested frailty models allow studying, jointly, survival processes of
#' recurrent and terminal events for hierarchically clustered data, by
#' considering the terminal event as an informative censoring and by including
#' two iid gamma random effects.
#' 
#' The joint nested frailty model includes two shared frailty terms, one for
#' the subgroup (\eqn{u_{fi}}) and one for the group (\eqn{w_f}) into the
#' hazard functions. This random effects account the heterogeneity in the data,
#' associated with unobserved covariates. The frailty terms act differently for
#' the two rates (\eqn{u_{fi}}, \eqn{w_f^\xi} for the recurrent rate and
#' \eqn{u_{fi}^\alpha, {w_i}} for the terminal event rate). The covariates
#' could be different for the recurrent rate and death rate.
#' 
#' For the \eqn{j^{th}} recurrence (j = 1, ..., \eqn{n_i}) of the \eqn{i^{th}}
#' individual (i = 1, ..., \eqn{m_f}) of the \eqn{f^{th}} group (f = 1, ...,
#' n), the joint nested gamma frailty model for recurrent event hazard function
#' \eqn{r_{fij}}(.) and for terminal event hazard function \eqn{\lambda_{fi}}
#' is :
#' 
#' \deqn{\left\{ \begin{array}{ll} r_{fij}(t|\omega_f, u_{fi}, \bold{X_{fij}})=
#' r_0(t) u_{fi} \omega_f^\xi \exp(\bold{\beta'} \bold{X_{fij}}) &
#' \mbox{(Recurrent)} \\ \lambda_{fi}(t|\omega_f, u_{fi},
#' \bold{X_{fij}})=\lambda_0(t)u_{fi}^\alpha \omega_f \exp(\bold{\gamma'}
#' \bold{X_{fi}}) & \mbox{(Death)} \\ \end{array} \right. }
#' 
#' where \eqn{r_0(t)}(resp. \eqn{\lambda_0(t)}) is the recurrent (resp.
#' terminal) event baseline hazard function, \eqn{\beta} (resp. \eqn{\gamma})
#' the regression coefficient vector, \eqn{\bold{X_{fij}}(t)} the covariates
#' vector. The random effects are \deqn{\omega_f \sim \Gamma \left(
#' \frac{1}{\eta}, \frac{1}{\eta}\right)} and \deqn{ u_{fi} \sim \Gamma \left(
#' \frac{1}{\theta}, \frac{1}{\theta} \right)}
#' }
#' }
#' 
#' @details{
#' Typical usages are for a Cox model
#' \preformatted{frailtyPenal(Surv(time,event)~var1+var2, data, \dots)}
#' 
#' for a shared model
#' \preformatted{frailtyPenal(Surv(time,event)~cluster(group)+var1+var2, data,
#' \dots)}
#' 
#' for a joint model
#' \preformatted{frailtyPenal(Surv(time,event)~cluster(group)+var1+var2+
#' var3+terminal(death), formula.terminalEvent=~ var1+var4, data, \dots)}
#' 
#' for a joint model for clustered data
#' \preformatted{frailtyPenal(Surv(time,event)~cluster(group)+num.id(group2)+
#' var1+var2+var3+terminal(death), formula.terminalEvent=~var1+var4, data,
#' \dots)}
#' 
#' for a joint model for data from nested case-control studies
#' \preformatted{frailtyPenal(Surv(time,event)~cluster(group)+num.id(group2)+
#' var1+var2+var3+terminal(death)+wts(wts.ncc),
#' formula.terminalEvent=~var1+var4, data, \dots)}
#' 
#' for a nested model
#' \preformatted{frailtyPenal(Surv(time,event)~cluster(group)+subcluster(sbgroup)+
#' var1+var2, data, \dots)}
#' 
#' for a joint nested frailty model
#' \preformatted{frailtyPenal(Surv(time,event)~cluster(group)+subcluster(sbgroup)+
#' var1+var2++terminal(death), formula.terminalEvent=~var1+var4, data, \dots)}
#' 
#' The estimated parameter are obtained using the robust Marquardt algorithm
#' (Marquardt, 1963) which is a combination between a Newton-Raphson algorithm
#' and a steepest descent algorithm. The iterations are stopped when the
#' difference between two consecutive log-likelihoods was small
#' \eqn{(<10^{-3})}, the estimated coefficients were stable (consecutive values
#' \eqn{(<10^{-3})}, and the gradient small enough \eqn{(<10^{-3})}. When
#' frailty parameter is small, numerical problems may arise. To solve this
#' problem, an alternative formula of the penalized log-likelihood is used (see
#' Rondeau, 2003 for further details). Cubic M-splines of order 4 are used for
#' the hazard function, and I-splines (integrated M-splines) are used for the
#' cumulative hazard function.
#' 
#' The inverse of the Hessian matrix is the variance estimator and to deal with
#' the positivity constraint of the variance component and the spline
#' coefficients, a squared transformation is used and the standard errors are
#' computed by the \eqn{\Delta}-method (Knight & Xekalaki, 2000). The smooth
#' parameter can be chosen by maximizing a likelihood cross validation
#' criterion (Joly and other, 1998). The integrations in the full log
#' likelihood were evaluated using Gaussian quadrature. Laguerre polynomials
#' with 20 points were used to treat the integrations on \eqn{[0,\infty[}
#' 
#' \bold{INITIAL VALUES}
#' 
#' The splines and the regression coefficients are initialized to 0.1. In case
#' of shared model, the program fits, firstly, an adjusted Cox model to give
#' new initial values for the splines and the regression coefficients. The
#' variance of the frailty term \eqn{\theta} is initialized to 0.1. Then, a
#' shared frailty model is fitted.
#' 
#' In case of a joint frailty model, the splines and the regression
#' coefficients are initialized to 0.5. The program fits an adjusted Cox model
#' to have new initial values for the regression and the splines coefficients.
#' The variance of the frailty term \eqn{\theta} and the coefficient
#' \eqn{\alpha} associated in the death hazard function are initialized to 1.
#' Then, it fits a joint frailty model.
#' 
#' In case of a general joint frailty model we need to initialize the
#' \code{jointGeneral} logical value to \code{TRUE}.
#' 
#' In case of a nested model, the program fits an adjusted Cox model to provide
#' new initial values for the regression and the splines coefficients. The
#' variances of the frailties are initialized to 0.1. Then, a shared frailty
#' model with covariates with only subgroup frailty is fitted to give a new
#' initial value for the variance of the subgroup frailty term. Then, a shared
#' frailty model with covariates and only group frailty terms is fitted to give
#' a new initial value for the variance of the group frailties. In a last step,
#' a nested frailty model is fitted.
#' 
#' In case of a joint nested model, the splines and the regression coefficients
#' are initialized to 0.5 and the variances of the frailty terms \eqn{\eta} and
#' \eqn{\xi} are initialized to 1.  If the option \code{'initialize'} is
#' \code{TRUE}, the program fits a joint frailty model to provide initial
#' values for the splines, covariates coefficients, variance \eqn{\theta} of
#' the frailty terms and \eqn{\alpha}. The variances of the second frailty term
#' (\eqn{\eta}) and the second coefficient \eqn{\xi} are initialized to 1.
#' Then, a joint nested frailty model is fitted.
#' 
#' \bold{NCC DESIGN}
#' 
#' It is possible to fit a joint frailty model for data from nested
#' case-control studies using the approach of weighted penalized maximum
#' likelihood. For this model, only splines can be used for baseline hazards
#' and no time-varying effects of covariates can be included. To accommodate
#' the nested case-control design, the formula for the recurrent events should
#' simply include the special term wts(wts.ncc), where wts.ncc refers to a
#' column of prespecified weights in the data set for every observation.  For
#' details, see Jazic et al., \emph{Submitted} (available on request from the
#' package authors).
#' }
#' 
#' @aliases frailtyPenal waldtest factor.names timedep.names
#' @usage
#' 
#' frailtyPenal(formula, formula.terminalEvent, data, recurrentAG = FALSE,
#' cross.validation = FALSE, jointGeneral,n.knots, kappa, maxit = 300, hazard =
#' "Splines", nb.int, RandDist = "Gamma", betaknots = 1, betaorder = 3,
#' initialize = TRUE, init.B, init.Theta, init.Alpha, Alpha, init.Ksi, Ksi,
#' init.Eta, LIMparam = 1e-3, LIMlogl = 1e-3, LIMderiv = 1e-3, print.times =
#' TRUE)
#' @param formula a formula object, with the response on the left of a
#' \eqn{\sim} operator, and the terms on the right. The response must be a
#' survival object as returned by the 'Surv' function like in survival package.
#' In case of interval-censored data, the response must be an object as
#' returned by the 'SurvIC' function from this package.  Interactions are
#' possible using * or :.
#' @param formula.terminalEvent only for joint and joint nested frailty models
#' : a formula object, only requires terms on the right to indicate which
#' variables are modelling the terminal event.  Interactions are possible using
#' * or :.
#' @param data a 'data.frame' with the variables used in 'formula'.
#' @param recurrentAG Logical value. Is Andersen-Gill model fitted? If so
#' indicates that recurrent event times with the counting process approach of
#' Andersen and Gill is used. This formulation can be used for dealing with
#' time-dependent covariates. The default is FALSE.
#' @param cross.validation Logical value. Is cross validation procedure used
#' for estimating smoothing parameter in the penalized likelihood estimation?
#' If so a search of the smoothing parameter using cross validation is done,
#' with kappa as the seed.  The cross validation is not implemented for several
#' strata, neither for interval-censored data. The cross validation has been
#' implemented for a Cox proportional hazard model, with no covariates. The
#' default is FALSE.
#' @param jointGeneral Logical value. Does the model include two independent
#' random effects? If so, this will fit a general joint frailty model with an
#' association between the recurrent events and a terminal event (explained by
#' the variance \eqn{\theta}) and an association amongst the recurrent events
#' (explained by the variance \eqn{\eta}).
#' @param n.knots integer giving the number of knots to use. Value required in
#' the penalized likelihood estimation.  It corresponds to the (n.knots+2)
#' splines functions for the approximation of the hazard or the survival
#' functions.  We estimate I or M-splines of order 4. When the user set a
#' number of knots equals to k (n.knots=k) then the number of interior knots is
#' (k-2) and the number of splines is (k-2)+order.  Number of knots must be
#' between 4 and 20. (See Note)
#' @param kappa positive smoothing parameter in the penalized likelihood
#' estimation. In a stratified shared model, this argument must be a vector
#' with kappas for both strata.  In a stratified joint model, this argument
#' must be a vector with kappas for both strata for recurrent events plus one
#' kappa for terminal event.  The coefficient kappa of the integral of the
#' squared second derivative of hazard function in the fit (penalized log
#' likelihood). To obtain an initial value for \code{kappa}, a solution is to
#' fit the corresponding shared frailty model using cross validation (See
#' cross.validation).  We advise the user to identify several possible tuning
#' parameters, note their defaults and look at the sensitivity of the results
#' to varying them. Value required. (See Note).
#' @param maxit maximum number of iterations for the Marquardt algorithm.
#' Default is 300
#' @param hazard Type of hazard functions: "Splines" for semiparametric hazard
#' functions using equidistant intervals or "Splines-per" using percentile with
#' the penalized likelihood estimation, "Piecewise-per" for piecewise constant
#' hazard function using percentile (not available for interval-censored data),
#' "Piecewise-equi" for piecewise constant hazard function using equidistant
#' intervals, "Weibull" for parametric Weibull functions. Default is "Splines".
#' In case of \code{jointGeneral = TRUE} or if a joint nested frailty model is
#' fitted, only \code{hazard = "Splines"} can be chosen.
#' @param nb.int Number of time intervals (between 1 and 20) for the parametric
#' hazard functions ("Piecewise-per", "Piecewise-equi"). In a joint model, you
#' need to specify a number of time interval for both recurrent hazard function
#' and the death hazard function (vector of length 2).
#' @param RandDist Type of random effect distribution: "Gamma" for a gamma
#' distribution, "LogN" for a log-normal distribution. Default is "Gamma". Not
#' implemented for nested model. If \code{jointGeneral = TRUE} or if a joint
#' nested frailty model is fitted, the log-normal distribution cannot be
#' chosen.
#' @param betaknots Number of inner knots used for the estimation of B-splines.
#' Default is 1. See 'timedep' function for more details. Not implemented for
#' nested and joint nested frailty models.
#' @param betaorder Order of the B-splines. Default is cubic B-splines (order =
#' 3). See 'timedep' function for more details. Not implemented for nested and
#' joint nested frailty models.
#' @param initialize Logical value, only for joint nested frailty models.
#' Option \code{TRUE} indicates fitting an appropriate standard joint frailty
#' model (without group effect, only the subgroup effect) to provide initial
#' values for the joint nested model. Default is \code{TRUE}.
#' @param init.B A vector of initial values for regression coefficients. This
#' vector should be of the same size as the whole vector of covariates with the
#' first elements for the covariates related to the recurrent events and then
#' to the terminal event (interactions in the end of each component). Default
#' is 0.1 for each (for Cox and shared model) or 0.5 (for joint and joint
#' nested frailty models).
#' @param init.Theta Initial value for variance of the frailties.
#' @param init.Alpha Only for joint and joint nested frailty models : initial
#' value for parameter alpha.
#' @param init.Ksi Only for joint nested frailty model : initial value for
#' parameter \eqn{\xi}.
#' @param init.Eta Only for general joint and joint nested frailty models :
#' initial value for the variance \eqn{\eta} of the frailty \eqn{v_i} (general
#' joint model) and of the frailty \eqn{\omega_i} (joint nested frailty model).
#' @param Alpha Only for joint and joint nested frailty model : input "None" so
#' as to fit a joint model without the parameter alpha.
#' @param Ksi Only for joint nested frailty model : input \code{"None"}
#' indicates a joint nested frailty model without the parameter \eqn{\xi}.
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
#' The following components are included in a 'frailtyPenal' object for each
#' model.
#' 
#' \item{b}{sequence of the corresponding estimation of the coefficients for
#' the hazard functions (parametric or semiparametric), the random effects
#' variances and the regression coefficients.} \item{call}{The code used for
#' the model.} \item{formula}{the formula part of the code used for the model.}
#' \item{coef}{the regression coefficients.} \item{cross.Val }{Logical value.
#' Is cross validation procedure used for estimating the smoothing parameters
#' in the penalized likelihood estimation?} \item{DoF}{Degrees of freedom
#' associated with the "kappa".} \item{groups}{the maximum number of groups
#' used in the fit.} \item{kappa}{ A vector with the smoothing parameters in
#' the penalized likelihood estimation corresponding to each baseline function
#' as components.} \item{loglikPenal}{the complete marginal penalized
#' log-likelihood in the semiparametric case.} \item{loglik}{the marginal
#' log-likelihood in the parametric case.} \item{n}{the number of observations
#' used in the fit.} \item{n.events}{the number of events observed in the fit.}
#' \item{n.iter}{number of iterations needed to converge.}
#' \item{n.knots}{number of knots for estimating the baseline functions in the
#' penalized likelihood estimation.} \item{n.strat}{ number of stratum.}
#' \item{varH}{the variance matrix of all parameters before positivity
#' constraint transformation. Then, the delta method is needed to obtain the
#' estimated variance parameters. That is why some variances don't match with
#' the printed values at the end of the model.} \item{varHIH}{the robust
#' estimation of the variance matrix of all parameters.} \item{x}{matrix of
#' times where both survival and hazard function are estimated. By default
#' seq(0,max(time),length=99), where time is the vector of survival times.}
#' \item{lam}{array (dim=3) of hazard estimates and confidence bands.}
#' \item{surv}{array (dim=3) of baseline survival estimates and confidence
#' bands.} \item{median}{The value of the median survival and its confidence bands. If there are
#' two stratas or more, the first value corresponds to the value for the 
#' first strata, etc.} \item{nbintervR}{Number of intervals (between 1 and 20) for the
#' parametric hazard functions ("Piecewise-per", "Piecewise-equi").}
#' \item{npar}{number of parameters.} \item{nvar}{number of explanatory
#' variables.} \item{LCV}{the approximated likelihood cross-validation
#' criterion in the semiparametric case (with H minus the converged Hessian
#' matrix, and l(.) the full
#' log-likelihood).\deqn{LCV=\frac{1}{n}(trace(H^{-1}_{pl}H) - l(.))}}
#' \item{AIC}{the Akaike information Criterion for the parametric
#' case.\deqn{AIC=\frac{1}{n}(np - l(.))}} \item{n.knots.temp}{initial value
#' for the number of knots.} \item{shape.weib}{shape parameter for the Weibull
#' hazard function.} \item{scale.weib}{scale parameter for the Weibull hazard
#' function.} \item{martingale.res}{martingale residuals for each cluster.}
#' \item{martingaleCox}{martingale residuals for observation in the Cox model.}
#' \item{Frailty}{Logical value. Was model with frailties fitted ?}
#' \item{frailty.pred}{empirical Bayes prediction of the frailty term (ie,
#' using conditional posterior distributions).} \item{frailty.var}{variance of
#' the empirical Bayes prediction of the frailty term (only for gamma frailty
#' models).} \item{frailty.sd}{standard error of the frailty empirical Bayes
#' prediction (only for gamma frailty models).} \item{global_chisq}{a vector
#' with the values of each multivariate Wald test.} \item{dof_chisq}{a vector
#' with the degree of freedom for each multivariate Wald test.}
#' \item{global_chisq.test}{a binary variable equals to 0 when no multivariate
#' Wald is given, 1 otherwise.} \item{p.global_chisq}{a vector with the
#' p_values for each global multivariate Wald test.} \item{names.factor}{Names
#' of the "as.factor" variables.} \item{Xlevels}{vector of the values that
#' factor might have taken.} \item{contrasts}{type of contrast for factor
#' variable.} \item{beta_p.value}{p-values of the Wald test for the estimated
#' regression coefficients.}
#' 
#' The following components are specific to \bold{shared} models.
#' 
#' \item{equidistant}{Indicator for the intervals used the estimation of
#' baseline hazard functions (for splines or pieceiwse-constaant functions) : 1
#' for equidistant intervals ; 0 for intervals using percentile (note:
#' \code{equidistant} = 2 in case of parametric estimation using Weibull
#' distribution).} \item{intcens}{Logical value. Indicator if a joint frailty
#' model with interval-censored data was fitted)} \item{theta}{variance of the
#' gamma frailty parameter \eqn{(\bold{Var}(\omega_i))}} \item{sigma2}{variance
#' of the log-normal frailty parameter \eqn{(\bold{Var}(\eta_i))}}
#' \item{linear.pred}{linear predictor: uses simply "Beta'X" in the cox
#' proportional hazard model or "Beta'X + log w_i" in the shared gamma frailty
#' models, otherwise uses "Beta'X + w_i" for log-normal frailty distribution.}
#' \item{BetaTpsMat}{matrix of time varying-effects and confidence bands (the
#' first column used for abscissa of times)} \item{theta_p.value}{p-value of
#' the Wald test for the estimated variance of the gamma frailty.}
#' \item{sigma2_p.value}{p-value of the Wald test for the estimated variance of
#' the log-normal frailty.}
#' 
#' The following components are specific to \bold{joint} models.
#' \item{intcens}{Logical value. Indicator if a joint frailty model with
#' interval-censored data was fitted)} \item{theta}{variance of the gamma
#' frailty parameter \eqn{(\bold{Var}(\omega_i))} or \eqn{(\bold{Var}(u_i))}}
#' \item{sigma2}{variance of the log-normal frailty parameter
#' \eqn{(\bold{Var}(\eta_i))} or \eqn{(\bold{Var}(v_i))}} \item{eta}{variance
#' of the second gamma frailty parameter in general joint frailty models
#' \eqn{(\bold{Var}(v_i))}} \item{indic_alpha}{indicator if a joint frailty
#' model with \eqn{\alpha} parameter was fitted} \item{alpha}{the coefficient
#' \eqn{\alpha} associated with the frailty parameter in the terminal hazard
#' function.} \item{nbintervR}{Number of intervals (between 1 and 20) for the
#' recurrent parametric hazard functions ("Piecewise-per", "Piecewise-equi").}
#' \item{nbintervDC}{Number of intervals (between 1 and 20) for the death
#' parametric hazard functions ("Piecewise-per", "Piecewise-equi").}
#' \item{nvar}{A vector with the number of covariates of each type of hazard
#' function as components.} \item{nvarRec}{number of recurrent explanatory
#' variables.} \item{nvarEnd}{number of death explanatory variables.}
#' \item{noVar1}{indicator of recurrent explanatory variables.}
#' \item{noVar2}{indicator of death explanatory variables.} \item{xR}{matrix of
#' times where both survival and hazard function are estimated for the
#' recurrent event. By default seq(0,max(time),length=99), where time is the
#' vector of survival times.} \item{xD}{matrix of times for the terminal
#' event.} \item{lamR}{array (dim=3) of hazard estimates and confidence bands
#' for recurrent event.} \item{lamD}{the same value as lamR for the terminal
#' event.} \item{survR}{array (dim=3) of baseline survival estimates and
#' confidence bands for recurrent event.} \item{survD}{the same value as survR
#' for the terminal event.} \item{martingale.res}{martingale residuals for each
#' cluster (recurrent).} \item{martingaledeath.res}{martingale residuals for
#' each cluster (death).} \item{linear.pred}{linear predictor: uses "Beta'X +
#' log w_i" in the gamma frailty model, otherwise uses "Beta'X + eta_i" for
#' log-normal frailty distribution} \item{lineardeath.pred}{linear predictor
#' for the terminal part : "Beta'X + alpha.log w_i" for gamma, "Beta'X +
#' alpha.eta_i" for log-normal frailty distribution} \item{Xlevels}{vector of
#' the values that factor might have taken for the recurrent part.}
#' \item{contrasts}{type of contrast for factor variable for the recurrent
#' part.} \item{Xlevels2}{vector of the values that factor might have taken for
#' the death part.} \item{contrasts2}{type of contrast for factor variable for
#' the death part.} \item{BetaTpsMat}{matrix of time varying-effects and
#' confidence bands for recurrent event (the first column used for abscissa of
#' times of recurrence)} \item{BetaTpsMatDc}{matrix of time varying-effects and
#' confidence bands for terminal event (the first column used for abscissa of
#' times of death)} \item{alpha_p.value}{p-value of the Wald test for the
#' estimated \eqn{\alpha}.} \item{ncc}{Logical value whether nested
#' case-control design with weights was used for the joint model.}
#' 
#' The following components are specific to \bold{nested} models.
#' 
#' \item{alpha}{variance of the cluster effect \eqn{(\bold{Var}(v_{i}))}}
#' \item{eta}{variance of the subcluster effect \eqn{(\bold{Var}(w_{ij}))}}
#' \item{subgroups}{the maximum number of subgroups used in the fit.}
#' \item{frailty.pred.group}{empirical Bayes prediction of the frailty term by
#' group.} \item{frailty.pred.subgroup}{empirical Bayes prediction of the
#' frailty term by subgroup.} \item{linear.pred}{linear predictor: uses "Beta'X
#' + log v_i.w_ij".} \item{subgbyg}{subgroup by group.} \item{n.strat}{A vector
#' with the number of covariates of each type of hazard function as
#' components.} \item{alpha_p.value}{p-value of the Wald test for the estimated
#' variance of the cluster effect.} \item{eta_p.value}{p-value of the Wald test
#' for the estimated variance of the subcluster effect.}
#' 
#' The following components are specific to \bold{joint nested frailty} models.
#' \item{theta}{variance of the subcluster effect \eqn{(\bold{Var}(u_{fi}))}}
#' \item{eta}{variance of the cluster effect \eqn{(\bold{Var}(\omega_f))}}
#' \item{alpha}{the power coefficient \eqn{\alpha} associated with the frailty
#' parameter (\eqn{u_{fi}}) in the terminal event hazard function.}
#' \item{ksi}{the power coefficient \eqn{\xi} associated with the frailty
#' parameter (\eqn{\omega_f}) in the recurrent event hazard function.}
#' \item{indic_alpha}{indicator if a joint frailty model with \eqn{\alpha}
#' parameter was fitted or not.} \item{indic_ksi}{indicator if a joint frailty
#' model with \eqn{\xi} parameter was fitted or not.}
#' \item{frailty.fam.pred}{empirical Bayes prediction of the frailty term by
#' family.} \item{eta_p.value}{p-value of the Wald test for the estimated
#' variance of the cluster effect.} \item{alpha_p.value}{p-value of the Wald
#' test for the estimated power coefficient \eqn{\alpha}.}
#' \item{ksi_p.value}{p-value of the Wald test for the estimated power
#' coefficient \eqn{\xi}.}
#' @note
#' 
#' From a prediction aim, we recommend you to input a data sorted by the group
#' variable with numerical numbers from 1 to n (number of groups). In case of a
#' nested model, we recommend you to input a data sorted by the group variable
#' then sorted by the subgroup variable both with numerical numbers from 1 to n
#' (number of groups) and from 1 to m (number or subgroups). "kappa" and
#' "n.knots" are the arguments that the user have to change if the fitted model
#' does not converge. "n.knots" takes integer values between 4 and 20. But with
#' n.knots=20, the model would take a long time to converge. So, usually, begin
#' first with n.knots=7, and increase it step by step until it converges.
#' "kappa" only takes positive values. So, choose a value for kappa (for
#' instance 10000), and if it does not converge, multiply or divide this value
#' by 10 or 5 until it converges.
#' @seealso \code{\link{SurvIC}}, \code{\link{cluster}},
#' \code{\link{subcluster}}, \code{\link{terminal}}, \code{\link{num.id}},
#' \code{\link{timedep}}
#' @references I. Jazic, S. Haneuse, B. French, G. MacGrogan, and V. Rondeau.
#' Design and analysis of nested case-control studies for recurrent events
#' subject to a terminal event. \emph{Submitted}.
#' 
#' A. Krol, A. Mauguen, Y. Mazroui, A. Laurent, S. Michiels and V. Rondeau
#' (2017). Tutorial in Joint Modeling and Prediction: A Statistical Software
#' for Correlated Longitudinal Outcomes, Recurrent Events and a Terminal Event.
#' \emph{Journal of Statistical Software} \bold{81}(3), 1-52.
#' 
#' V. Rondeau, Y. Mazroui and J. R. Gonzalez (2012). Frailtypack: An R package
#' for the analysis of correlated survival data with frailty models using
#' penalized likelihood estimation or parametric estimation. \emph{Journal of
#' Statistical Software} \bold{47}, 1-28.
#' 
#' Y. Mazroui, S. Mathoulin-Pelissier, P. Soubeyranb and V. Rondeau (2012)
#' General joint frailty model for recurrent event data with a dependent
#' terminalevent: Application to follicular lymphoma data. \emph{Statistics in
#' Medecine}, \bold{31}, 11-12, 1162-1176.
#' 
#' V. Rondeau, J.P. Pignon, S. Michiels (2011). A joint model for the
#' dependance between clustered times to tumour progression and deaths: A
#' meta-analysis of chemotherapy in head and neck cancer. \emph{Statistical
#' methods in medical research} \bold{897}, 1-19.
#' 
#' V. Rondeau, S. Mathoulin-Pellissier, H. Jacqmin-Gadda, V. Brouste, P.
#' Soubeyran (2007). Joint frailty models for recurring events and death using
#' maximum penalized likelihood estimation:application on cancer events.
#' \emph{Biostatistics} \bold{8},4, 708-721.
#' 
#' V. Rondeau, L. Filleul, P. Joly (2006). Nested frailty models using maximum
#' penalized likelihood estimation. \emph{Statistics in Medicine}, \bold{25},
#' 4036-4052.
#' 
#' V. Rondeau, D. Commenges, and P. Joly (2003). Maximum penalized likelihood
#' estimation in a gamma-frailty model. \emph{Lifetime Data Analysis} \bold{9},
#' 139-153.
#' 
#' C.A. McGilchrist, and C.W. Aisbett (1991). Regression with frailty in
#' survival analysis. \emph{Biometrics} \bold{47}, 461-466.
#' 
#' D. Marquardt (1963). An algorithm for least-squares estimation of nonlinear
#' parameters. \emph{SIAM Journal of Applied Mathematics}, 431-441.
#' @keywords models
#' @export
#' @examples
#' 
#' 
#' \dontrun{
#' 
#' ###---  COX proportional hazard model (SHARED without frailties) ---###
#' ###---  estimated with penalized likelihood ---###
#' 
#' data(kidney)
#' frailtyPenal(Surv(time,status)~sex+age,
#' n.knots=12,kappa=10000,data=kidney)
#' 
#' ###---  Shared Frailty model  ---###
#' 
#' frailtyPenal(Surv(time,status)~cluster(id)+sex+age,
#' n.knots=12,kappa=10000,data=kidney)
#' 
#' #-- with an initialisation of regression coefficients
#' 
#' frailtyPenal(Surv(time,status)~cluster(id)+sex+age,
#' n.knots=12,kappa=10000,data=kidney,init.B=c(-1.44,0))
#' 
#' #-- with truncated data
#' 
#' data(dataNested)
#' 
#' frailtyPenal(Surv(t1,t2,event) ~ cluster(group),
#' data=dataNested,n.knots=10,kappa=10000,
#' cross.validation=TRUE,recurrentAG=FALSE)
#' 
#' #-- stratified analysis
#' 
#' data(readmission)
#' frailtyPenal(Surv(time,event)~cluster(id)+dukes+strata(sex),
#' n.knots=10,kappa=c(10000,10000),data=readmission)
#' 
#' #-- recurrentAG=TRUE
#' 
#' frailtyPenal(Surv(t.start,t.stop,event)~cluster(id)+sex+dukes+
#' charlson,data=readmission,n.knots=6,kappa=1e5,recurrentAG=TRUE)
#' 
#' #-- cross.validation=TRUE
#' 
#' frailtyPenal(Surv(t.start,t.stop,event)~cluster(id)+sex+dukes+
#' charlson,data=readmission,n.knots=6,kappa=5000,recurrentAG=TRUE,
#' cross.validation=TRUE)
#' 
#' #-- log-normal distribution
#' 
#' frailtyPenal(Surv(t.start,t.stop,event)~cluster(id)+sex+dukes+
#' charlson,data=readmission,n.knots=6,kappa=5000,recurrentAG=TRUE,
#' RandDist="LogN")
#' 
#' ###--- Joint Frailty model (recurrent and terminal events) ---###
#' 
#' data(readmission)
#' #-- Gap-time
#' modJoint.gap <- frailtyPenal(Surv(time,event)~cluster(id)+sex+dukes+charlson+
#' terminal(death),formula.terminalEvent=~sex+dukes+charlson,
#' data=readmission,n.knots=14,kappa=c(9.55e+9,1.41e+12),
#' recurrentAG=FALSE)
#' 
#' #-- Calendar time
#' modJoint.calendar <- frailtyPenal(Surv(t.start,t.stop,event)~cluster(id)+
#' sex+dukes+charlson+terminal(death),formula.terminalEvent=~sex
#' +dukes+charlson,data=readmission,n.knots=10,kappa=c(9.55e9,1.41e12),
#' recurrentAG=TRUE)
#' 
#' #-- without alpha parameter
#' modJoint.gap <- frailtyPenal(Surv(time,event)~cluster(id)+sex+dukes+charlson+
#' terminal(death),formula.terminalEvent=~sex+dukes+charlson,
#' data=readmission,n.knots=10,kappa=c(9.55e9,1.41e12),
#' recurrentAG=FALSE,Alpha="None")
#' 
#' #-- log-normal distribution
#' 
#' modJoint.log <- frailtyPenal(Surv(t.start,t.stop,event)~cluster(id)+sex
#' +dukes+charlson+terminal(death),formula.terminalEvent=~sex
#' +dukes+charlson,data=readmission,n.knots=10,kappa=c(9.55e9,1.41e12),
#' recurrentAG=TRUE,RandDist="LogN")
#' 
#' ###--- Joint frailty model for NCC data ---###
#' data(dataNCC)
#' modJoint.ncc <- frailtyPenal(Surv(t.start,t.stop,event)~cluster(id)+cov1
#' +cov2+terminal(death)+wts(ncc.wts), formula.terminalEvent=~cov1+cov2,
#' data=dataNCC,n.knots=8,kappa=c(1.6e+10, 5.0e+03),recurrentAG=TRUE, RandDist="LogN") 
#' 
#' 
#' ###--- Joint Frailty model for clustered data ---###
#' 
#' #-- here is generated cluster (5 clusters)
#' readmission <- transform(readmission,group=id%%5+1)
#' 
#' #-- exclusion all recurrent events --#
#' #--  to obtain framework of semi-competing risks --#
#' readmission2 <- subset(readmission, (t.start == 0 & event == 1) | event == 0)
#' 
#' joi.clus.gap <- frailtyPenal(Surv(time,event)~cluster(group)+
#' num.id(id)+dukes+charlson+sex+chemo+terminal(death),
#' formula.terminalEvent=~dukes+charlson+sex+chemo,
#' data=readmission2,recurrentAG=FALSE, n.knots=8,
#' kappa=c(1.e+10,1.e+10) ,Alpha="None")
#' 
#' 
#' ###--- General Joint model (recurrent and terminal events) 
#' with 2 covariates ---###
#' 
#' data(readmission)
#' modJoint.general <- frailtyPenal(Surv(time,event) ~ cluster(id) + dukes +
#' charlson + sex  + chemo + terminal(death), 
#' formula.terminalEvent = ~ dukes + charlson + sex + chemo,
#' data = readmission, jointGeneral = TRUE,  n.knots = 8, 
#' kappa = c(2.11e+08, 9.53e+11))
#' 
#' 
#' ###--- Nested Frailty model ---###
#' 
#' ##***** WARNING *****##
#' # Data should be ordered according to cluster and subcluster
#' 
#' data(dataNested)
#' modClu <- frailtyPenal(Surv(t1,t2,event)~cluster(group)+
#' subcluster(subgroup)+cov1+cov2,data=dataNested,
#' n.knots=8,kappa=50000)
#' 
#' modClu.str <- frailtyPenal(Surv(t1,t2,event)~cluster(group)+
#' subcluster(subgroup)+cov1+strata(cov2),data=dataNested,
#' n.knots=8,kappa=c(50000,50000))
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
#' JoiNes_GapSpline <- frailtyPenal(formula = Surv(time, event) 
#' ~ subcluster(id) + cluster(group) + dukes + terminal(death), 
#' formula.terminalEvent = ~dukes, data = readmissionNested, 
#' recurrentAG = FALSE, n.knots = 8, kappa = c(9.55e+9, 1.41e+12), 
#' initialize = TRUE, init.Alpha = 1.091, Ksi = "None")
#' 
#' }
#' 
#' 
"frailtyPenal" <-
  function (formula, formula.terminalEvent, data, recurrentAG=FALSE, cross.validation=FALSE, jointGeneral, n.knots, kappa,maxit=300, 
            hazard="Splines", nb.int, RandDist="Gamma", betaknots=1,betaorder=3, initialize=TRUE, init.B, init.Theta, init.Alpha, Alpha, init.Ksi, Ksi, init.Eta,
            LIMparam=1e-3, LIMlogl=1e-3, LIMderiv=1e-3, print.times=TRUE){

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
    
    if (missing(jointGeneral)) jointGeneral<-FALSE
    if (!missing(init.Eta) & jointGeneral)  init.Alpha <- init.Eta
    
    # al suppression de l'argument joint
    if (!missing(formula.terminalEvent)) joint <- TRUE
    else joint <- FALSE
    if ((!missing(Alpha) | !missing(init.Alpha)) & !joint) stop("init.Alpha and Alpha parameters belong to joint frailty model")
    
    
    #ad 15/02/12 :add Audrey
    m2 <- match.call()
    m2$formula <- m2$formula.terminalEvent <- m2$recurrentAG <- m2$cross.validation <- m2$n.knots <- m2$kappa <- m2$maxit <- m2$hazard <- m2$nb.int <- m2$RandDist <- m2$betaorder <- m2$betaknots <- m2$init.B <- m2$LIMparam <- m2$LIMlogl <- m2$LIMderiv <- m2$print.times <- m2$init.Theta <- m2$init.Alpha <- m2$Alpha <- m2$init.Ksi <- m2$Ksi <- m2$init.Eta <- m2$Eta <- m2$initialize <- m2$... <- NULL
    Names.data <- m2$data
    
    #### Betaknots et betaorder ####
    if (betaknots > 10) stop("Number of knots for beta(t) greater than 10 is useless, please choose a number between 0 and 10 (3 is optimal)")
    if ((betaorder == 0)|(betaorder > 4)) stop("B-splines order for beta(t) must be a number between 1 and 4 (3 is optimal)")
    
    #### Frailty distribution specification ####
    if (!(RandDist %in% c("Gamma","LogN"))) { stop("Only 'Gamma' and 'LogN' distributions for frailties are allowed") }
    logNormal <- switch(RandDist,"Gamma"=0,"LogN"=1)
    
    if (RandDist=="LogN" & jointGeneral==TRUE)        stop("Log normal distribution is not available for the Joint General Model !")
    if ((hazard!="Splines") & jointGeneral== TRUE)    stop("No general joint frailty model allowed here! Only 'Splines' is allowed")
    
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
          if (!missing(nb.int)){
            stop("When the hazard function equals 'Splines' or 'Weibull', 'nb.int' argument must be deleted.")
          }
          size1 <- 100
          size2 <- 100
          equidistant <- 1
          nbintervR <- 0
          nbintervDC <- 0
        }
        ### Weibull
        if (typeof == 2){
          if (!missing(nb.int)){
            stop("When the hazard function equals 'Splines' or 'Weibull', 'nb.int' argument must be deleted.")
          }
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
          if (haztemp == "Splines-per") equidistant <- 0
          else equidistant <- 1
          size1 <- 100
          size2 <- 100
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
    
    #AD:
    if (missing(formula))stop("The argument formula must be specified in any model")
    if(class(formula)!="formula")stop("The argument formula must be a formula")
    
    if(typeof == 0){
      #AD:
      if (missing(n.knots))stop("number of knots are required")   
      #AD:	 
      n.knots.temp <- n.knots	
      #AD
      if (n.knots<4) n.knots<-4
      if (n.knots>20) n.knots<-20
      
      if (missing(kappa))stop("smoothing parameter (kappa) is required")
      #AD:
      
      if (length(kappa)>1 & cross.validation){
        stop("The cross validation is not implemented for two strata or more")
      }
      
      if (joint & cross.validation){
        stop("The cross validation is not implemented for the joint model")  
      }  
    }else{
      if (!(missing(n.knots)) || !(missing(kappa)) || !(missing(cross.validation))){
        stop("When parametric hazard function is specified, 'kappa', 'n.knots' and 'cross.validation' arguments must be deleted.")
      }
      n.knots <- 0
      kappa <- 0
      crossVal <- 0      
    }
    call <- match.call()	
    
    m <- match.call(expand.dots = FALSE) # recupere l'instruction de l'utilisateur	
    
    m$formula.terminalEvent <- m$n.knots <- m$recurrentAG <- m$cross.validation <- m$jointGeneral <- m$kappa <- m$maxit <- m$hazard <- m$nb.int <- m$RandDist <- m$betaorder <- m$betaknots <- m$init.B <- m$LIMparam <- m$LIMlogl <- m$LIMderiv <-  m$print.times <- m$init.Theta <- m$init.Alpha <- m$Alpha <- m$init.Ksi <- m$Ksi <- m$init.Eta <- m$Eta <- m$initialize <- m$... <- NULL    
    special <- c("strata", "cluster", "subcluster", "terminal","num.id","timedep", "wts") #wts for weights (ncc design) ncc - nested case-control
    
    Terms <- if (missing(data)){ 
      terms(formula, special)
    }else{
      terms(formula, special, data = data)
    }
    
    ord <- attr(Terms, "order") # longueur de ord=nbre de var.expli
    #if (length(ord) & any(ord != 1))stop("Interaction terms are not valid for this function")
    #si pas vide tous si il ya au moins un qui vaut 1 on arrete
    
    m$formula <- Terms    
    
    m[[1]] <- as.name("model.frame") # m[[1]]=frailtypenal, il le remplace par model.frame en fait    
    
    m <- eval(m, sys.parent()) #ici la classe de m est un data.frame donc il recupere ce qu'on lui donne en argument
    
    cluster <- attr(Terms, "specials")$cluster # (indice) nbre de var qui sont en fonction de cluster()
    # al 13/02/14 : suppression de l'argument Frailty
    if (length(cluster)) Frailty <- TRUE
    else                 Frailty <- FALSE
    
    subcluster <- attr(Terms, "specials")$subcluster #nbre de var qui sont en fonction de subcluster()	
    
    wts <- attr(Terms, "specials")$wts #nbre de var qui sont en fonction de wts()	
    
    # booleen pour voir si l'objet Y est reponse avant tri des donnees Surv ou SurvIC
    classofY <- attr(model.extract(m, "response"),"class")
    # attention le package pec rajoute un element dans l'attribut "class" des objets de survie
    if (length(classofY)>1) classofY <- classofY[2]
    
    typeofY <- attr(model.extract(m, "response"),"type") # type de reponse : interval etc..
    
    # # # We retrieve here the t0 name give in the 'data' dataframe (needed for checking if left truncated data)
    # # # vart0 <- dimnames(m)[[2]][1]	
    # # # vart0 <- gsub('(.*\\()(.*)', '\\2', strsplit(vart0,",")[[1]][1])
    
    
    #Al : tri du jeu de donnees par cluster croissant
    if (length(cluster)){
      tempc <- untangle.specials(Terms, "cluster", 1:10)
      ord <- attr(Terms, "order")[tempc$terms] # le type ordinal de la variable (0,1 ou 2)
      if (any(ord > 1))stop("Cluster can not be used in an interaction")
      m <- m[order(m[,tempc$vars]),] # soit que des nombres, soit des caracteres
      ordre <- as.integer(row.names(m)) # recupere l'ordre du data set
      cluster <- strata(m[, tempc$vars], shortlabel = TRUE)
      uni.cluster <- unique(cluster)
    }
    
    
    #Verification de l'utilisation du parametre Ksi pour le nested joint modele
    if ((!missing(Ksi) | !missing (init.Ksi)) & (!joint & !length(subcluster))) stop("init.Ksi and Ksi parameters belong only to nested joint frailty model")
    
    # verification de la structure nested si besoin
    if (length(subcluster) && Frailty == TRUE){
      tempsub <- untangle.specials(Terms, "subcluster", 1:10)
      ordsub <- attr(Terms, "order")[tempsub$terms] # type ordinal de la variable
      if (any(ordsub > 1))stop("subcluster can not be used in an interaction")
      if (any(ifelse(apply(ifelse(table(m[,tempsub$vars],m[,tempc$vars])>0,1,0),1,sum)==1,FALSE,TRUE))){
        stop("nested structure is necessary to fit a nested model")
      }
      
      # tri par ordre croissant de subcluster a l'interieur des clusters
      m <- m[order(m[,tempc$vars],m[,tempsub$vars]),]
      ordretmp <- order(m[,tempsub$vars])
      
      subcluster <- strata(m[, tempsub$vars], shortlabel = TRUE)
      ordre <- as.integer(row.names(m))
      
      subcluster <- as.integer(subcluster) # a determiner si il y en a besoin
      curr <- subcluster[1]
      subcluster[1] <- 1
      for (i in 2:length(subcluster)) {
        if (subcluster[i] == curr) { subcluster[i] <- subcluster[i-1] }
        else {
          curr <- subcluster[i]
          subcluster[i] <- subcluster[i-1] + 1
        }
      }
    }
    
    #Al
    if (NROW(m) == 0)stop("No (non-missing) observations") #nombre ligne different de 0
    
    Y <- model.extract(m, "response") # objet de type Surv =Time
    
    if (classofY == "SurvIC") intcens <- TRUE # booleen censure par intervalle
    else intcens <- FALSE
    
    if (intcens == TRUE) {
      if (classofY != "SurvIC") stop("When interval censoring, must use the SurvIC fonction")
    } else {
      #if (!inherits(Y, "Surv")) stop("Response must be a survival object") #test si c bien un objet de type "Surv"
      if (classofY != "Surv") stop("Response must be a survival object")
    }
    
    ll <- attr(Terms, "term.labels")#liste des variables explicatives
    
    ##########################################################################################################################
    # add Myriam 11/05/2016 Checking left-truncating data 
    # # # troncat <- function(cluster, time0){
    # # # mini = time0[1]
    # # # indiv <- cluster[1]
    # # # for (i in 1:length(time0)){
    # # # if (cluster[i] != indiv){
    # # # if (mini != 0) stop('Sorry but left-troncature is not allowed for joint or joint nested frailty models')
    # # # else {
    # # # indiv <- cluster[i]
    # # # mini <- time0[i]
    # # # }
    # # # }
    # # # else{
    # # # if (mini > time0[i]) mini <- time0[i]
    # # # }	
    # # # }
    # # # }
    # # # Timet0 <- data[vart0]
    # # # if(length(cluster)){		
    # # # varclust <- gsub('(.*\\()(.*)\\)', '\\2', tempc$vars) # pour recuperer le nom de la variable definie par cluster()
    # # # }
    #verification pour le modele joint nested
    # # # if (length(subcluster) & joint){
    # # # varsubclust <- gsub('(.*\\()(.*)\\)', '\\2', tempsub$vars)
    # # # Timet0 <- Timet0[order(data[,varclust], data[,varsubclust]),]
    # # # if (typeofY == "counting") troncat(as.numeric(as.vector(data[order(data[,varclust], data[,varsubclust]),varsubclust])), Timet0)
    # # # }
    # verification pour le modele joint
    # # # if (!length(subcluster) & joint){
    # # # Timet0 <- list(Timet0)[[1]][order(data[,varclust]),]
    # # # #data[,varclust] = data[order(data[,varclust]), varclust]
    # # # if (typeofY == "counting") troncat(as.numeric(as.vector(data[order(data[,varclust]), varclust])), Timet0)
    # # # }
    
    # # # NOTE : 07-03-17 : this operation makes errors in frailtyPenal execution when the user put in surv() formula 
    # # # covariates with a such name : 'data$var' or 'data[,2]' or 'data[,"var"]'
    # # # So, for the moment we disable this process, any left-truncated data will not be  detected by frailtypack 
    # # # and calculation will done with  bias because left truncating is not took into account joint and joint nested modelling
    ##############################################################################################################################    
    #=========================================================>
    
    mt <- attr(m, "terms") #m devient de class "formula" et "terms"
    
    X <- if (!is.empty.model(mt)) model.matrix(mt, m) #, contrasts) #idem que mt sauf que ici les factor sont divise en plusieurs variables
    
    ind.place <- unique(attr(X,"assign")[duplicated(attr(X,"assign"))]) ### unique : changement au 25/09/2014
    
    vec.factor <- NULL
    vec.factor <- c(vec.factor,ll[ind.place])
    
    #=========================================================>
    # On determine le nombre de categorie pour chaque var categorielle
    strats <- attr(Terms, "specials")$strata #nbre de var qui sont en fonction de strata()
    cluster <- attr(Terms, "specials")$cluster #nbre de var qui sont en fonction de cluster()
    num.id <- attr(Terms, "specials")$num.id #nbre de var qui sont en fonction de patkey()
    vartimedep <- attr(Terms, "specials")$timedep #nbre de var en fonction de timedep()
    wts <- attr(Terms, "specials")$wts #nbre de var en fonction de wts()....redundant but keeping with cluster()
    
    #booleen pour savoir si au moins une var depend du tps
    if (is.null(vartimedep)) timedep <- 0
    else timedep <- 1
    
    if (intcens & (equidistant == 0)) stop("You can not fit a model with a baseline hazard function estimated using percentiles and interval censoring")
    if (intcens & timedep) stop("You can not use time-varying effect covariates with interval censoring")
    if (intcens & cross.validation) stop("You can not do cross validation with interval censoring")
    if (intcens & logNormal) stop("It is currently impossible to fit a model with interval censoring and a log normal distribution of the frailties")
    if (timedep & logNormal) stop("You can not use time-varying effect covariates with a log normal distribution of the frailties")
    
    if (intcens & joint & length(subcluster)) stop("Cox model doesn't support intervall-censored survival data")
    if (intcens & !joint & length(subcluster)) stop ("Nested model doesn't support intervall-censored survival data")
    
    if(is.null(num.id)){
      joint.clust <- 1
    }else{
      joint.clust <- 0
      if (!joint) stop("num.id function can only be used with joint models")
    }
    if (jointGeneral==TRUE) joint.clust <- 2
    if (joint & length(subcluster)) joint.clust <- 3
    
    subcluster <- attr(Terms, "specials")$subcluster #nbre de var qui sont en fonction de subcluster()
    
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
    
    if (length(wts)){
      ll <- ll[-grep("wts",ll)]
    }
    
    #   plus besoin de as.factor() pour afficher le test de Wald global
    if (length(grep("strata",vec.factor))) vec.factor <- vec.factor[-grep("strata",vec.factor)]
    if (length(grep("cluster",vec.factor))) vec.factor <- vec.factor[-grep("cluster",vec.factor)]
    if (length(grep("subcluster",vec.factor))) vec.factor <- vec.factor[-grep("subcluster",vec.factor)]
    if (length(grep("num.id",vec.factor))) vec.factor <- vec.factor[-grep("num.id",vec.factor)]
    if (length(grep("wts",vec.factor))) vec.factor <- vec.factor[-grep("wts",vec.factor)]
    
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
      }
    } ) 
    
    if(length(grep("terminal",ll))>0){
      ind.place <- grep(paste(vec.factor,collapse="|"),ll[-grep("terminal",ll)])
    }
    else{
      ind.place <- grep(paste(vec.factor,collapse="|"),ll)
    }  
    if(length(vec.factor) > 0){
      vect.fact <- attr(X,"dimnames")[[2]]
      
      #vect.fact <- vect.fact[grep("factor",vect.fact)]
      vect.fact <- vect.fact[grep(paste(vec.factor,collapse="|"),vect.fact)]
      
      # 		vect.fact <-apply(matrix(vect.fact,ncol=1,nrow=length(vect.fact)),MARGIN=1,FUN=function(x){
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
                which <- i
              }
            }
          }
          occur[i] <- length.grep
          
        }else{         
          if(length(vect.fact[-which.interaction])>0){occur[i] <- length(grep(vec.factor[i],vect.fact[-which.interaction]))
          }else{occur[i] <- length(grep(vec.factor[i],vect.fact))}
        }
      }
    }
    
    if (joint & length(subcluster)){
      if(initialize) initialize <- 1
      else initialize <- 0
    }
    else{
      if (!missing (initialize)) stop ('Sorry but \'initialize\' option is only available for joint nested frailty model')
    }
    
    
    #=========================================================>
    
    terminalEvent <- attr(Terms, "specials")$terminal #nbre de var qui sont en fonction de terminal()
    
    dropx <- NULL
    
    if (length(cluster) & Frailty == TRUE){
      tempc <- untangle.specials(Terms, "cluster", 1:10)
      ord <- attr(Terms, "order")[tempc$terms]
      if (any(ord > 1))stop("Cluster can not be used in an interaction")      
      cluster <- strata(m[, tempc$vars], shortlabel = TRUE)    
      dropx <- tempc$terms
      uni.cluster<-unique(cluster)
    }
    else if (!length(cluster) & Frailty == TRUE){      
      stop("grouping variable is needed")      
    }
    else if (length(cluster) & Frailty == FALSE){
      stop("cluster not necessary for proportional hazard model")
    }
    else if (!length(cluster) & Frailty == FALSE){
      cluster <- 1:nrow(m) #nrow(data) # valeurs inutiles pour un modele de Cox
      uni.cluster <- 1:nrow(m) #nrow(data)
    }
    
    if (!missing(RandDist) & (Frailty == FALSE)){
      stop("RandDist not necessary for proportional hazard model")
    }
    
    if (length(num.id)){
      temppat <- untangle.specials(Terms, "num.id", 1:10)
      num.id <- m[,temppat$vars]
      dropx <- c(dropx,temppat$terms)
    }
    
    if(length(uni.cluster)==1){ 
      stop("grouping variable must have more than 1 level")
    }  
    
    if (length(subcluster)){
      tempsub <- untangle.specials(Terms, "subcluster", 1:10)
      ordsub <- attr(Terms, "order")[tempsub$terms]
      if (any(ordsub > 1))stop("subcluster can not be used in an interaction")
      subcluster <- strata(m[, tempsub$vars], shortlabel = TRUE)
      dropx <- c(dropx,tempsub$terms)
      uni.subcluster<-unique(subcluster)
      #if (joint)stop("joint model is not implemented for nested model")      
      if(length(uni.subcluster)==1){
        stop("subcluster variable must have more than 1 level")
      }
      
    }
    #AD:	
    if (length(cluster) == length(subcluster)){
      if (all(all.equal(cluster,subcluster)==T)){
        stop("'Subgroup' variable and 'group' variable need to be different")
      }
    }
    #AD:	
    if (length(strats)){      
      temp <- untangle.specials(Terms, "strata", 1)
      dropx <- c(dropx, temp$terms)
      if (length(temp$vars) == 1)strata.keep <- m[[temp$vars]]
      else strata.keep <- strata(m[, temp$vars], shortlabel = TRUE)
      strats <- as.numeric(strata.keep)
      uni.strat<-length(unique(strats))      
      if (missing(kappa)) stop("smoothing parameter (kappa) is required")      
      if (length(subcluster)){
        if (uni.strat > 2) {
          if (joint) stop("maximum number of strata for nested frailty joint model is 2")
          else stop("maximum number of strata for nested frailty joint model is 2")
        }
      }
      else{
        if (uni.strat > 6) stop("maximum number of strata is 6")
      }
      if ((uni.strat > 2) & (intcens)) stop("maximum number of strata for interval censored data is 2")
      if ((uni.strat > 2) & (timedep)) stop("maximum number of strata for time-varying effect of covariates is 2")      
    }else{
      uni.strat<-1
      strats <- rep(1,nrow(data))
    }    
    if (!joint & (typeof==0) & (length(kappa)!=uni.strat)) stop("wrong length of argument 'kappa' for the current stratification")
    
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
    
    
    #IJ 2018:
    if (length(wts)){      
      tempwts <- untangle.specials(Terms, "wts", 1:10)    
      ord <- attr(Terms, "order")[tempwts$terms] # ord[6]=1 ici dans notre exemple      
      if (any(ord > 1))stop("Weights can not be used in an interaction")
      dropx <- c(dropx,tempwts$terms) # vecteur de position
      weights.vec <- strata(m[, tempwts$vars], shortlabel = TRUE)
      weights.vec <- as.numeric(as.character(weights.vec))      
    } else{
      weights.vec <- rep(1, nrow(data))
    }
    weights.agg <- aggregate(weights.vec, by=list(cluster), FUN=function(x) x[length(x)])[,2]	
    
    
    #type <- attr(Y, "type")
    type <- typeofY 
    if (type != "right" && type != "counting" && type != "interval" && type != "intervaltronc") { # Cox supporte desormais la censure par intervalle
      stop(paste("Cox model doesn't support \"", type, "\" survival data", sep = ""))
    }    
    #	if ((type == "interval" || type == "interval2" || type == "intervaltronc") && intcens == FALSE) { # rajout
    #		stop("You are trying to do interval censoring without intcens = TRUE")
    #	}    
    if (type != "counting" && recurrentAG) {
      stop("recurrentAG needs counting process formulation")
    }    
    if (intcens == TRUE & recurrentAG == TRUE) {
      stop("recurrentAG invalid for interval censored data")
    }
    
    #drop contient les position liees au fonction() ic ex:cluster(id) et terminal(death)    
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
      noVar1 <- 1
    }else{
      X <- X[, -1, drop = FALSE]
      noVar1 <- 0
    }
    # on enleve ensuite la premiere colonne correspondant a id 
    nvar<-ncol(X) #nvar==1 correspond a 2 situations:    
    # au cas ou on a aucune var explicative dans la partie rec, mais X=0
    # cas ou on a 1seul var explicative, ici X est en general different de 0    
    varnotdep <- colnames(X)[-grep("timedep",colnames(X))]
    vardep <- colnames(X)[grep("timedep",colnames(X))]
    vardep <- apply(matrix(vardep,ncol=1,nrow=length(vardep)),1,timedep.names)    
    if (length(intersect(varnotdep,vardep)) != 0) {
      stop("A variable is both used as a constant and time-varying effect covariate")
    }
    nvartimedep <- length(vardep)    
    filtretps <- rep(0,nvar)
    filtretps[grep("timedep",colnames(X))] <- 1    
    var<-matrix(c(X),nrow=nrow(X),ncol=nvar)
    
    n<-nrow(X)
    
    #add Alexandre 04/06/2012
    #lire les donnees differemment si censure par intervalle
    if (intcens==TRUE) {
      if (type=="intervaltronc") {
        tt0 <- Y[,1]
        tt1 <- Y[,2]
        ttU <- Y[,3]
        cens <- Y[,4]
      } 
      else {
        tt0 <- rep(0,n)
        tt1 <- Y[,1]
        ttU <- Y[,2]
        cens <- Y[,3]
        tt1[tt1==0] <- 0.1
      }
    } 
    else {
      if (type=="right"){
        tt0 <- rep(0,n)
        tt1 <- Y[,1]
        cens <- Y[,2]
        ttU <- Y[,1] # rajouter quand meme dans frailPenal mais ne sera pas utilise
      } 
      else {
        tt0 <- Y[,1]
        tt1 <- Y[,2]
        cens <- Y[,3]
        ttU <- Y[,2] # ne sera pas pris en compte dans le tri des temps de survie dans frailtypack.f90
      }                   # attention ne pas mettre de 0 sinon en cas de left trunc probleme dans la logV
    }    
    if (min(cens)==0) cens.data<-1
    if (min(cens)==1 && max(cens)==1) cens.data<-0    
    AG<-ifelse(recurrentAG,1,0)
    if (typeof == 0){
      crossVal<-ifelse(cross.validation,0,1)
    }  
    #=======================================>
    #======= Construction du vecteur des indicatrice
    if(length(vec.factor) > 0){
      #		ind.place <- ind.place -1
      k <- 0
      for(i in 1:length(vec.factor)){
        ind.place[i] <- ind.place[i]+k
        k <- k + occur[i]-1
      }
    }
    
    #==================================
    # Begin SHARED MODEL
    #
    
    if (!joint & !length(subcluster)){
      if((equidistant %in% c(0,1)) & (typeof == 1)){
        if (missing(nb.int)) stop("Number of time interval 'nb.int' is required")
        if (length(nb.int) != 1) stop("Wrong length of number of time interval argument 'nb.int'")
        if (class(nb.int) != "numeric") stop("The argument 'nb.int' must be a numeric")
        if ((nb.int < 1)) stop("Number of time interval 'nb.int' must be between 1 and 20")
        if (nb.int > 20){
          nb.int <- 20
          indic.nb.int <- 1 # equals 1 for nb.int > 20
        }
        else{
          indic.nb.int <- 0 # equals 0 for nb.int < 20
        }
        nbintervR <- nb.int
        size1 <- 3*nbintervR
      }
      if ((typeof == 0) | (typeof == 2)) indic.nb.int <- 0      
      if (sum(as.double(var))==0) nvar <- 0      
      if (timedep==0){
        npbetatps <- 0
      }else{
        npbetatps <- (betaknots+betaorder-1)*nvartimedep
      }      
      np <- switch(as.character(typeof),
                   "0"=(as.integer(uni.strat) * (as.integer(n.knots) + 2) + as.integer(nvar) + as.integer(Frailty)) + npbetatps,
                   "1"=(as.integer(uni.strat) * nbintervR + nvar + as.integer(Frailty)) + npbetatps,
                   "2"=(as.integer(uni.strat) * 2 + nvar + as.integer(Frailty)) + npbetatps)
      
      # traitement de l'initialisation du Beta rentre par l'utilisateur
      Beta <- rep(0,np)
      if (!missing(init.B)) {
        if (length(init.B) != nvar) stop("Wrong number of regression coefficients in init.B")
        if (timedep) stop("You can hardly know initial regression coefficient while time-varying effect")
        Beta <- c(rep(0,np-nvar),init.B)        
      }
      if (!missing(init.Theta)) {
        if (!is.numeric(init.Theta)) stop("init.Theta must be numeric")
        if (Frailty==FALSE) stop("init.Theta does not exist in a Cox proportional hazard model")
        Beta[np-nvar] <- init.Theta
      }
      
      xSuT <- matrix(0,nrow=100,ncol=uni.strat)
      if (typeof==0){
        mt1 <- size1
      }else{
        mt1 <- 100
      }
      size2 <- mt1
      
      flush.console()
      if (print.times){
        ptm<-proc.time()
        cat("\n")
        cat("Be patient. The program is computing ... \n")
      }
      
      ans <- .Fortran(C_frailpenal,                      
                      as.integer(n),
                      as.integer(length(uni.cluster)),
                      as.integer(cens.data),
                      as.integer(uni.strat),
                      as.integer(Frailty),
                      as.integer(n.knots),
                      as.double(kappa),
                      as.double(tt0),
                      as.double(tt1),
                      as.integer(cens),
                      
                      as.integer(cluster),
                      as.integer(nvar),
                      as.double(strats),
                      as.double(var),
                      as.integer(AG),
                      as.integer(noVar1),
                      as.integer(maxit),
                      as.integer(crossVal),
                      np=as.integer(np),
                      b=as.double(Beta),
                      
                      as.double(matrix(0,nrow=np,ncol=np)),
                      as.double(matrix(0,nrow=np,ncol=np)),
                      as.double(0),
                      LCV=as.double(rep(0,2)),
                      as.double(matrix(0,nrow=size1,ncol=uni.strat)),
                      as.double(array(0,dim=c(size1,3,uni.strat))), #lam
                      xSuT=as.double(xSuT),
                      as.double(array(0,dim=c(size2,3,uni.strat))),
                      as.integer(typeof),
                      as.integer(equidistant),
                      
                      as.integer(nbintervR),
                      as.integer(size1),
                      as.integer(0),
                      as.integer(0),
                      as.integer(0),
                      as.double(c(0,0)),
                      as.double(0),
                      istop=as.integer(0),
                      shape.weib=as.double(rep(0,2)),
                      scale.weib=as.double(rep(0,2)),
                      
                      as.integer(mt1),
                      zi=as.double(rep(0,(n.knots+6))),
                      martingale.res=as.double(rep(0,as.integer(length(uni.cluster)))),
                      martingaleCox=as.double(rep(0,n)),
                      frailty.pred=as.double(rep(0,as.integer(length(uni.cluster)))),
                      frailty.var=as.double(rep(0,as.integer(length(uni.cluster)))),
                      frailty.sd=as.double(rep(0,as.integer(length(uni.cluster)))),
                      linear.pred=as.double(rep(0,n)),
                      time=as.double(rep(0,(nbintervR+1))),
                      as.integer(intcens), # rajout
                      
                      as.double(ttU), # rajout
                      logNormal=as.integer(logNormal),
                      timedep=as.integer(timedep),
                      as.integer(betaknots),
                      as.integer(betaorder),
                      as.integer(filtretps),
                      BetaTpsMat=as.double(matrix(0,nrow=101,ncol=1+4*nvartimedep)),
                      EPS=as.double(c(LIMparam,LIMlogl,LIMderiv))
      )#,
      #PACKAGE = "frailtypack") # 58 arguments
      #AD:      
      if (ans$istop == 4){
        warning("Problem in the loglikelihood computation. The program stopped abnormally. Please verify your dataset. \n")
      }      
      if (ans$istop == 2){
        warning("Model did not converge.")
      }
      if (ans$istop == 3){
        warning("Matrix non-positive definite.")
      }
      
      #AD:      
      if (noVar1 == 1) nvar<-0      
      np <- ans[[19]]
      fit <- NULL
      fit$b <- ans$b
      fit$na.action <- attr(m, "na.action")
      fit$call <- call
      fit$n <- n
      fit$groups <- length(uni.cluster)
      fit$n.events <- ans[[34]]
      #Al:
      if(dim(table(cens,cluster))[1]==2)	fit$n.eventsbygrp <- table(cens,cluster)[2,] 
      
      if(as.character(typeof)=="0"){
        fit$logLikPenal <- ans[[23]]
      }else{
        fit$logLik <- ans[[23]]
      }
      
      if (Frailty) {
        fit$coef <- ans[[20]][(np - nvar - npbetatps + 1):np]
        if (logNormal == 0){fit$theta <- (ans[[20]][np - nvar - npbetatps])^2
        }else{fit$sigma2 <- (ans[[20]][np - nvar - npbetatps])^2
        }
      }
      if (!Frailty) {
        if (logNormal == 0) fit$theta <- NULL
        else fit$sigma2 <- NULL
      }
      if (noVar1 == 1) {
        fit$coef <- NULL
      } 
      else{
        fit$coef <- ans[[20]][(np - nvar - npbetatps + 1):np]
        noms <- factor.names(colnames(X))
        if(length(grep(":",noms))>0)noms <- factor.names(noms)   
        if (timedep == 1){ # on enleve les parametres des B-splines qui ne serviront pas a l'utilisateur
          while (length(grep("timedep",noms))!=0){
            pos <- grep("timedep",noms)[1]
            noms <- noms[-pos]
            fit$coef <- fit$coef[-(pos:(pos+betaorder+betaknots-1))]
          }
        }
        names(fit$coef) <- noms
      }      
      temp1 <- matrix(ans[[21]], nrow = np, ncol = np)
      temp2 <- matrix(ans[[22]], nrow = np, ncol = np)
      if (Frailty) {
        fit$varTheta <- c(temp1[(np - nvar - npbetatps),(np - nvar - npbetatps)],temp2[(np - nvar - npbetatps),(np - nvar - npbetatps)])
        fit$varTheta <- ((2*ans[[20]][np - nvar - npbetatps])^2)*fit$varTheta # delta-method
        if (logNormal == 0)	fit$theta_p.value <- 1 - pnorm(fit$theta/sqrt(fit$varTheta[1]))
        else fit$sigma2_p.value <- 1 - pnorm(fit$sigma2/sqrt(fit$varTheta[1]))
        
      }
      
      #AD:modification des dimensions des tableaux
      if(nvar > 0){        
        fit$varH <- temp1[(np - nvar - npbetatps + 1):np, (np - nvar - npbetatps + 1):np]
        fit$varHIH <- temp2[(np - nvar - npbetatps + 1):np, (np - nvar - npbetatps + 1):np]
        noms <- factor.names(colnames(X))
        if(length(grep(":",noms))>0)noms <- factor.names(noms)
        if (timedep == 1){ # on enleve les variances des parametres des B-splines
          while (length(grep("timedep",noms))!=0){
            pos <- grep("timedep",noms)[1]
            noms <- noms[-pos]
            fit$varH <- fit$varH[-(pos:(pos+betaorder+betaknots-1)),-(pos:(pos+betaorder+betaknots-1))]
            fit$varHIH <- fit$varHIH[-(pos:(pos+betaorder+betaknots-1)),-(pos:(pos+betaorder+betaknots-1))]
          }
        }
      }      
      fit$varHtotal <- temp1 # new Al: 20/06/13
      fit$varHIHtotal <- temp2      
      fit$formula <- formula(Terms)
      
      fit$x <- matrix(ans[[25]], nrow = size1, ncol = uni.strat)
      fit$lam <- if(typeof == 1){array(ans[[26]][seq(1,length(ans[[26]]),3)], dim = c(nb.int,3,uni.strat))} else{array(ans[[26]], dim = c(size1,3,uni.strat))} # Le lam s'ecrit differemment selon la fonction de hasard (Piecewise selon le nombre d'intervalles specifies, Weibull et Splines selon size1 = 100)
      fit$xSu <- matrix(ans$xSuT, nrow = 100, ncol = uni.strat)
      fit$surv <- array(ans[[28]], dim = c(size2,3,uni.strat))
      
      fit$type <- type
      fit$n.strat <- uni.strat
      fit$n.iter <- ans[[33]]
      
      median <- NULL
      for (i in (1:fit$n.strat)) median[i] <- ifelse(typeof==0, minmin(fit$surv[,1,i],fit$x), minmin(fit$surv[,1,i],fit$xSu))
      lower <- NULL
      for (i in (1:fit$n.strat)) lower[i] <- ifelse(typeof==0, minmin(fit$surv[,2,i],fit$x), minmin(fit$surv[,2,i],fit$xSu))
      upper <- NULL
      for (i in (1:fit$n.strat)) upper[i] <- ifelse(typeof==0, minmin(fit$surv[,3,i],fit$x), minmin(fit$surv[,3,i],fit$xSu))
      fit$median <- cbind(lower,median,upper)
      
      if (typeof == 0){
        fit$n.knots<-n.knots
        if (uni.strat > 1) fit$kappa <- ans[[36]]
        else fit$kappa <- ans[[36]][1]
        fit$DoF <- ans[[37]]
        fit$cross.Val<-cross.validation
        fit$n.knots.temp <- n.knots.temp
        fit$zi <- ans$zi
      }
      if(typeof == 1) fit$time <- ans$time
      #AD:      
      fit$LCV <- ans$LCV[1]
      fit$AIC <- ans$LCV[2]
      fit$npar <- np
      fit$nvar <- nvar
      fit$noVar1 <- noVar1
      fit$indic.nb.int <- indic.nb.int
      #AD:      
      if(ans[[35]]==2000) stop("The cross validation procedure cannot be finished. Try to change
                               either the number of knots or the seed for kappa parameter")
      
      fit$typeof <- typeof
      fit$equidistant <- equidistant
      fit$nbintervR <- nbintervR
      fit$istop <- ans$istop
      
      fit$AG <- recurrentAG
      fit$intcens <- intcens # rajout
      fit$logNormal <- ans$logNormal
      
      fit$shape.weib <- ans$shape.weib
      fit$scale.weib <- ans$scale.weib
      fit$Names.data <- Names.data
      if (Frailty) fit$Names.cluster <- Names.cluster
      fit$Frailty <- Frailty
      if (Frailty){
        fit$martingale.res <- ans$martingale.res
        fit$frailty.pred <- ans$frailty.pred
        if (logNormal==0){
          fit$frailty.var <- ans$frailty.var
          fit$frailty.sd <- ans$frailty.sd
        }
      }
      else{
        fit$martingaleCox <- ans$martingaleCox
      }
      if (Frailty) fit$linear.pred <- ans$linear.pred[order(ordre)] # pour remettre dans le bon ordre
      else fit$linear.pred <- ans$linear.pred      
      fit$BetaTpsMat <- matrix(ans$BetaTpsMat,nrow=101,ncol=1+4*nvartimedep)
      fit$nvartimedep <- nvartimedep      
      fit$Names.vardep <- vardep      
      fit$EPS <- ans$EPS      
      #
      #========================= Test de Wald pour shared
      
      if(ans$istop==1){
        if ((length(vec.factor) > 0) & (timedep == 0)){
          Beta <- ans[[20]][(np - nvar + 1):np]
          VarBeta <- fit$varH
          nfactor <- length(vec.factor)
          p.wald <- rep(0,nfactor)
          
          if(fit$istop == 1) fit$global_chisq <- waldtest(N=nvar,nfact=nfactor,place=ind.place,modality=occur,b=Beta,Varb=VarBeta)
          else fit$global_chisq <- 0 
          
          fit$dof_chisq <- occur
          fit$global_chisq.test <- 1
          # Calcul de pvalue globale
          for(i in 1:length(vec.factor)){
            p.wald[i] <- signif(1 - pchisq(fit$global_chisq[i], occur[i]), 3)
          }
          fit$p.global_chisq <- p.wald
          fit$names.factor <- vec.factor
        }
        else{
          fit$global_chisq.test <- 0
        }
      }
      
      if(!is.null(fit$coef)){
        if(nvar != 1){
          seH <- sqrt(diag(fit$varH))
        }else{
          seH <- sqrt(fit$varH)
        }
        fit$beta_p.value <- 1 - pchisq((fit$coef/seH)^2, 1)
      }
      
      #===============================================
      if (length(Xlevels) >0)fit$Xlevels <- Xlevels
      fit$contrasts <- contr.save
      attr(fit,"joint")<-joint
      attr(fit,"subcluster")<-FALSE
      class(fit) <- "frailtyPenal"
      
    }  # End SHARED MODEL
    
    
    #
    # Begin JOINT MODEL
    #
    
    if (joint & !length(subcluster)){
      # Preparing data ...
      #AD:
      if(Frailty =="FALSE"){
        stop("For joint frailty models, 'Frailty' must be equal to 'TRUE' ")
      }
      #AD
      if (classofY == "Surv"){
        if (!recurrentAG){
          if(joint.clust==0){    
            tempdc <- aggregate(tt1,by=list(num.id),FUN=sum)[,2]
            lignedc0 <- length(tempdc)
            tempdc <- cbind(rep(0,lignedc0),tempdc)
            clusterdc <- aggregate(cluster,by=list(num.id),FUN=function(x) x[length(x)])[,2]
            tt1.death <- 0
            tt0.death <- 0
            
            tt1 <- aggregate(tt1,by=list(num.id), FUN=function(x) x[1])[,2]
            tt0 <- aggregate(tt0,by=list(num.id), FUN=function(x) x[1])[,2]
            cluster <- aggregate(cluster,by=list(num.id), FUN=function(x) x[1])[,2]
            table <- as.data.frame(cbind(tt0,tt1,cluster))
            table <- table[order(table$cluster),]
            
            cluster <- table$cluster
            tt1 <- table$tt1
            tt0 <- table$tt0            
            
            n <- length(tt0)
            uni.cluster<-unique(num.id)#unique(cluster)          
          }
          else{
            tt1.death<-aggregate(tt1,by=list(cluster),FUN=sum)[,2]
            tt0.death<-rep(0,length(tt1.death))
            clusterdc <- 0
            lignedc0 <- 0
            tempdc <- 0
          }
        }
        else{
          if(joint.clust==0){
            #tempdc <- aggregate(tt1,by=list(num.id,cluster),FUN=sum)[,2]
            tempdc<-aggregate(tt1,by=list(num.id),FUN=function(x) x[length(x)])[,2]
            lignedc0 <- length(tempdc)
            tempdc <- cbind(rep(0,lignedc0),tempdc)
            clusterdc <- aggregate(cluster,by=list(num.id),FUN=function(x) x[length(x)])[,2]
            tt1.death <- 0
            tt0.death <- 0
            
            tt1 <- aggregate(tt1,by=list(num.id), FUN=function(x) x[1])[,2]
            tt0 <- aggregate(tt0,by=list(num.id), FUN=function(x) x[1])[,2]
            cluster <- aggregate(cluster,by=list(num.id), FUN=function(x) x[1])[,2]
            table <- as.data.frame(cbind(tt0,tt1,cluster))
            table <- table[order(table$cluster),]
            
            cluster <- table$cluster
            tt1 <- table$tt1
            tt0 <- table$tt0            
            
            n <- length(tt0)
            uni.cluster<-unique(num.id)#unique(cluster)
            
          }else{
            tt1.death<-aggregate(tt1,by=list(cluster),FUN=function(x) x[length(x)])[,2]
            tt0.death<-rep(0,length(tt1.death))
            clusterdc <- 0
            lignedc0 <- 0
            tempdc <- 0            
          }
        }
      }
      else{ # Interval censoring
        if (recurrentAG == TRUE) stop("You can't fit joint models on interval-censored data with recurrentAG = TRUE")
        if (joint.clust==0){
          tempdc0 <- aggregate(tt0,by=list(num.id),FUN=function(x) x[length(x)])[,2]
          tempdc <- aggregate(tt1,by=list(num.id),FUN=function(x) x[length(x)])[,2]
          
          lignedc0 <- length(tempdc)
          tempdc <- cbind(tempdc0,tempdc)
          clusterdc <- aggregate(cluster,by=list(num.id),FUN=function(x) x[length(x)])[,2]
          tt1.death <- 0
          tt0.death <- 0
          
          # prendre en compte seulement un evenement pour le joint cluster
          tt0 <- aggregate(tt0,by=list(num.id),FUN=function(x) x[1])[,2]
          tt1 <- aggregate(tt1,by=list(num.id),FUN=function(x) x[1])[,2]
          ttU <- aggregate(ttU,by=list(num.id),FUN=function(x) x[1])[,2]
          cens <- aggregate(cens,by=list(num.id),FUN=function(x) x[1])[,2]
          cluster <- aggregate(cluster,by=list(num.id),FUN=function(x) x[1])[,2]
          
          if (!is.null(ncol(var))){ # if more than one covariate 
            varAG<-aggregate(var[,1],by=list(num.id), FUN=function(x) x[1])[,2]
            if (ncol(var)>1){
              for (i in 2:ncol(var)){
                varAG.i<-aggregate(var[,i],by=list(num.id), FUN=function(x) x[1])[,2]
                varAG<-cbind(varAG,varAG.i)
              }
            }
            var<-varAG
          }
          else{
            var<-aggregate(var,by=list(num.id), FUN=function(x) x[1])[,2]
          }
          nobs <- n
          n <- length(tt0)			  
        }
        else{
          tt1.death<-aggregate(tt1,by=list(cluster),FUN=function(x) x[length(x)])[,2]
          tt0.death<-rep(0,length(tt1.death))
          clusterdc <- 0
          lignedc0 <- 0
          tempdc <- 0
          
          # prendre en compte seulement un evenement pour le joint
          tt0 <- aggregate(tt0,by=list(cluster),FUN=function(x) x[1])[,2]
          tt1 <- aggregate(tt1,by=list(cluster),FUN=function(x) x[1])[,2]
          ttU <- aggregate(ttU,by=list(cluster),FUN=function(x) x[1])[,2]
          cens <- aggregate(cens,by=list(cluster),FUN=function(x) x[1])[,2]
          if (!is.null(ncol(var))){ # if more than one covariate 
            varAG<-aggregate(var[,1],by=list(cluster), FUN=function(x) x[1])[,2]
            if (ncol(var)>1){
              for (i in 2:ncol(var)){
                varAG.i<-aggregate(var[,i],by=list(cluster), FUN=function(x) x[1])[,2]
                varAG<-cbind(varAG,varAG.i)
              }
            }
            var<-varAG
          }
          else{
            var<-aggregate(var,by=list(cluster), FUN=function(x) x[1])[,2]
          }
          nobs <- n
          n <- length(tt0)
        }
      }      
      Terms2 <- if (missing(data)){
        if (!missing(formula.terminalEvent))terms(formula.terminalEvent, special)
      }
      else{
        if (!missing(formula.terminalEvent))terms(formula.terminalEvent, special, data = data)
      }
      #AD:
      if (!missing(formula.terminalEvent)){
        ord2 <- attr(Terms2, "order")        
        if (length(ord2) & any(ord2 != 1)){
          #       		stop("Interaction terms are not valid for terminal event formula")
        }
      }
      #AD:
      
      #AD:Joint model needs "terminal()"
      if (ind.terminal){
        if(joint.clust==0 ){
          icdc00 <- aggregate(terminal,by=list(num.id),FUN=function(x) x[length(x)])[,2] #+aggregate(cens,by=list(num.id),FUN=function(x) x[length(x)])[,2]
          terminalEvent <- 0
        }
        else{
          terminalEvent<-aggregate(terminal,by=list(cluster),FUN=function(x) x[length(x)])[,2]
          icdc00 <- 0
        }
      }
      else{
        stop(" Joint frailty model miss specified ")
      }
      #AD:     
      
      # terminalEvent might be 0-1
      if(joint.clust==0){
        if (!all(icdc00%in%c(2,1,0))){
          stop("terminal must contain a variable coded 0-1 and a non-factor variable")
        }
      }
      else{
        if (!all(terminalEvent%in%c(2,1,0))){
          stop("terminal must contain a variable coded 0-1 and a non-factor variable")
        }
      }      
      m2 <- match.call(expand.dots = FALSE)
      ## AD: modified 20 06 2011, for no covariates on terminal event part
      if (missing(formula.terminalEvent)){
        m2$n.knots <- m2$recurrentAG <- m2$cross.validation <- m2$kappa <- m2$maxit <- m2$hazard <- m2$nb.int <- m2$RandDist <- m2$betaorder <- m2$betaknots <- m2$init.B <- m2$LIMparam <- m2$LIMlogl <- m2$LIMderiv <- m2$print.times <- m2$init.Theta <- m2$init.Alpha <- m2$Alpha <- m2$init.Ksi <- m2$Ksi <- m2$init.Eta <- m2$Eta <- m2$initialize <- m2$... <- NULL
      }else{
        m2$formula.terminalEvent <- m2$n.knots <- m2$recurrentAG <- m2$cross.validation <- m2$jointGeneral<- m2$kappa <- m2$maxit <- m2$hazard <- m2$nb.int <- m2$RandDist <- m2$betaorder <- m2$betaknots <- m2$init.B <- m2$LIMparam <- m2$LIMlogl <- m2$LIMderiv <- m2$print.times <- m2$init.Theta <- m2$init.Alpha <- m2$Alpha <- m2$init.Ksi <- m2$Ksi <- m2$init.Eta <- m2$Eta <- m2$initialize <- m2$... <- NULL
      }     
      
      m2$formula <- Terms2
      m2[[1]] <- as.name("model.frame")
      m2 <- eval(m2, sys.parent()) #ici il prend les colonne associe au var.exp, si rien il prend tout     
      
      match.noNA<-dimnames(m2)[[1]]%in%dimnames(m)[[1]]#masque logique pour ne pas avoir de NA
      
      m2<-m2[match.noNA, ,drop=FALSE]#m2 inchanger si pas de NA
      
      if (!missing(formula.terminalEvent))newTerms2<-Terms2
      
      #=========================================================>
      if (!missing(formula.terminalEvent)){
        X2 <- model.matrix(newTerms2, m2)
        lldc <- attr(newTerms2,"term.labels")
        #ind.placedc <- grep("factor",lldc)
        ind.placedc <- unique(attr(X2,"assign")[duplicated(attr(X2,"assign"))])#changement unique le 26/09/2014
        vec.factordc <- NULL
        vec.factordc <- c(vec.factordc,lldc[ind.placedc])        
        mat.factordc <- matrix(vec.factordc,ncol=1,nrow=length(vec.factordc))
        # Fonction servant a prendre les termes entre "as.factor"
        vec.factordc <-apply(mat.factordc,MARGIN=1,FUN=function(x){
          if (length(grep("factor",x))>0){
            if(length(grep(":",x))>0){
              if(grep('\\(',unlist(strsplit(x,split="")))[1]<grep(":",unlist(strsplit(x,split="")))[1]  && length(grep('\\(',unlist(strsplit(x,split=""))))==1){						
                pos1 <- grep("r",unlist(strsplit(x,split="")))[1]+2
                pos2 <- grep(":",unlist(strsplit(x,split="")))[1]-2
                pos3 <- grep(":",unlist(strsplit(x,split="")))[1]
                pos4 <- length(unlist(strsplit(x,split="")))
                return(paste(substr(x,start=pos1,stop=pos2),substr(x,start=pos3,stop=pos4),sep=""))
              }
              else if(grep("\\(",unlist(strsplit(x,split="")))[1]>grep(":",unlist(strsplit(x,split="")))[1] && length(grep('\\(',unlist(strsplit(x,split=""))))==1){
                pos2 <- grep(":",unlist(strsplit(x,split="")))[1]
                pos3 <- grep("\\(",unlist(strsplit(x,split="")))[1]+1
                pos4 <- length(unlist(strsplit(x,split="")))-1
                return(paste(substr(x,start=1,stop=pos2),substr(x,start=pos3,stop=pos4),sep=""))
              }
              else{#both factors
                pos1 <- grep("r",unlist(strsplit(x,split="")))[1]+2
                pos2 <- grep(":",unlist(strsplit(x,split="")))[1]-2
                pos3 <- grep("\\(",unlist(strsplit(x,split="")))[2]+1
                pos4 <- length(unlist(strsplit(x,split="")))-1
                return(paste(substr(x,start=pos1,stop=pos2),":",substr(x,start=pos3,stop=pos4),sep=""))
              }
            }
            else{
              pos1 <- grep("r",unlist(strsplit(x,split="")))[1]+2
              pos2 <- length(unlist(strsplit(x,split="")))-1
              return(substr(x,start=pos1,stop=pos2))
            }
          }
          else return(x)
        }  )        
        # On determine le nombre de categorie pour chaque var categorielle
        if(length(vec.factordc) > 0){
          vect.factdc <- attr(X2,"dimnames")[[2]]
          vect.factdc <- vect.factdc[grep(paste(vec.factordc,collapse="|"),vect.factdc)]			  
          occurdc <- rep(0,length(vec.factordc))
          #    for(i in 1:length(vec.factordc)){
          #      occurdc[i] <- length(grep(vec.factordc[i],vect.factdc))
          #    }
          #  }
          interactiondc<-as.vector(apply(matrix(vect.factdc,nrow=length(vect.factdc)),MARGIN=1,FUN=function(x){length(grep(":",unlist(strsplit(x,split=""))))}))
          which.interactiondc <- which(interactiondc==1)				
          for(i in 1:length(vec.factordc)){				  
            if(length(grep(":",unlist(strsplit(vec.factordc[i],split=""))))>0){					
              pos <- grep(":",unlist(strsplit(vec.factordc[i],split="")))
              length.grep <- 0
              for(j in 1:length(vect.factdc)){
                if(j%in%which.interactiondc){
                  if(length(grep(substr(vec.factordc[i],start=1,stop=pos-1),vect.factdc[j]))>0 && length(grep(substr(vec.factordc[i],start=pos+1,stop=length(unlist(strsplit(vec.factordc[i],split="")))),vect.factdc[j]))>0){
                    length.grep <- length.grep + 1
                    which <- i
                  }
                }
              }
              occurdc[i] <- length.grep					
            }
            else{				
              if(length(vect.factdc[-which.interactiondc])>0) occurdc[i] <- length(grep(vec.factordc[i],vect.factdc[-which.interactiondc]))
              else occurdc[i] <- length(grep(vec.factordc[i],vect.factdc))
            }
          }
        }   
        #=========================================================>
        assign <- lapply(attrassign(X2, newTerms2)[-1], function(x) x - 1)
        Xlevels2 <- .getXlevels(newTerms2, m2)
        contr.save2 <- attr(X2, 'contrasts')
        #========================================>
        if(length(vec.factordc) > 0){
          positiondc <- unlist(assign,use.names=F)
        }
        #========================================>
        if (ncol(X2) == 1){
          X2<-X2-1
          noVar2 <- 1
        }
        else{
          X2 <- X2[, -1, drop = FALSE]
          noVar2 <- 0
        }			
        nvar2 <- ncol(X2)
        #       if(sum(ord)>length(ord)){
        #          for(i in 1:length(ord)){
        #            if(ord[i]>1){
        #              name_v1 <- strsplit(as.character(lldc[i]),":")[[1]][1]
        #              name_v2 <- strsplit(as.character(lldc[i]),":")[[1]][2]
        #              if(length(grep("as.factor",name_v1))>0){name_v1<-substring(name_v1,11,nchar(name_v1)-1)
        #                                                      v1 <- as.factor(data[,names(data)==name_v1])}
        #              else{v1 <- data[,names(data)==name_v1]}
        #              if(length(grep("as.factor",name_v2))>0){name_v2<-substring(name_v2,11,nchar(name_v2)-1)
        #                                                      v2 <- as.factor(data[,names(data)==name_v2])}
        #              else{v2 <- data[,names(data)==name_v2]}
        #              
        #              if(is.factor(v1) && length(levels(v1))>2)stop("Interactions not allowed for factors with 3 or more levels (yet)")
        #              if(is.factor(v2) && length(levels(v2))>2)stop("Interactions not allowed for factors with 3 or more levels (yet)")
        #              
        #              
        #            }
        #          }
        #        }
        vartimedep2 <- attr(Terms2, "specials")$timedep #nbre de var en fonction de timedep()
        
        # verifier qu'il y ait du timedep dans la deuxieme formule
        if (!is.null(vartimedep2)) timedep <- 1			
        varnotdep2 <- colnames(X2)[-grep("timedep",colnames(X2))]
        vardep2 <- colnames(X2)[grep("timedep",colnames(X2))]
        vardep2 <- apply(matrix(vardep2,ncol=1,nrow=length(vardep2)),1,timedep.names)			
        if (length(intersect(varnotdep2,vardep2)) != 0) {
          stop("A variable is both used as a constant and time-varying effect covariate in the formula of terminal event")
        }			
        nvartimedep2 <- length(vardep2)			
        filtretps2 <- rep(0,nvar2)
        filtretps2[grep("timedep",colnames(X2))] <- 1			
        vardc.temp<-matrix(c(X2),nrow=nrow(X2),ncol=nvar2)	
        
        if(is.null(nrow(m2))){
          if (length(m2) != nrow(m)){
            stop(" There are missing values in the covariates modelling the terminal event. \n Prepare data only with complete cases")
          }
        }
        else{		  
          if (nrow(m2) != nrow(m)){
            stop(" There are missing values in the covariates modelling the terminal event. \n Prepare data only with complete cases")
          }          
        } 
        
        if(joint.clust==0 ){ 
          if (!is.null(ncol(vardc.temp))){
            vaxdc00<-aggregate(vardc.temp[,1],by=list(num.id), FUN=function(x) x[length(x)])[,2]  # num.id au lieu de cluster           
            if (ncol(vardc.temp)>1){
              for (i in 2:ncol(vardc.temp)){
                vaxdc00.i<-aggregate(vardc.temp[,i],by=list(num.id), FUN=function(x) x[length(x)])[,2]
                vaxdc00<-cbind(vaxdc00,vaxdc00.i)
              }
            }
          }
          else{
            vaxdc00<-aggregate(vardc.temp,by=list(num.id), FUN=function(x) x[length(x)])[,2]
          }
          vardc <- 0
        }
        else{
          if (!is.null(ncol(vardc.temp))){
            vardc<-aggregate(vardc.temp[,1],by=list(cluster), FUN=function(x) x[length(x)])[,2]            
            if (ncol(vardc.temp)>1){              
              for (i in 2:ncol(vardc.temp)){
                vardc.i<-aggregate(vardc.temp[,i],by=list(cluster), FUN=function(x) x[length(x)])[,2]
                vardc<-cbind(vardc,vardc.i)
              }
            }
          }
          else{
            vardc<-aggregate(vardc.temp,by=list(cluster), FUN=function(x) x[length(x)])[,2]
          }
          vaxdc00 <- 0
        }
      }
      else{
        noVar2 <- 1
        vardc<-0
      }      
      if ((classofY == "SurvIC") & (joint.clust==1 || joint.clust== 2) & (recurrentAG=FALSE)) cluster <- aggregate(cluster,by=list(cluster),FUN=function(x) x[1])[,2]
      # 	cluster <- aggregate(cluster,by=list(cluster),FUN=function(x) x[1])[,2]
      nvarRec<-nvar      
      if (!missing(formula.terminalEvent)){
        #=======================================>
        #======= Construction du vecteur des indicatrice
        
        if(length(vec.factordc) > 0){
          k <- 0
          for(i in 1:length(vec.factordc)){
            ind.placedc[i] <- ind.placedc[i]+k
            k <- k + occurdc[i]-1
          }
        }        
        #==================================
        if(joint.clust==1 || joint.clust==2){
          if(is.null(nrow(vardc))){
            nvarEnd<-1
          }
          else{
            nvarEnd<-ncol(vardc)
          }
        }
        else{
          if(is.null(nrow(vaxdc00))){
            nvarEnd<-1
          }
          else{
            nvarEnd<-ncol(vaxdc00)
          }
        }
      }
      else nvarEnd<-0      
      if (sum(as.double(var))==0) nvarRec <- 0
      if ((joint.clust==0) & sum(as.double(vaxdc00))==0) nvarEnd <- 0
      if ((joint.clust==1 || joint.clust==2 ) & sum(as.double(vardc))==0) nvarEnd <- 0
      
      nvar<-nvarRec+nvarEnd
      
      # ... end preparing data
      #AD:
      effet <- 1
      indic_alpha <- 1
      indic_xi <- 3
      
      if (!missing(Alpha)){ # new : joint more flexible alpha = 0
        if (Alpha=="None") indic_alpha <- 0
        else stop("Alpha can only take 'None' as a value in this version of frailtypack package")
      }      
      nst <- uni.strat #2
      
      indices <- c(indic_alpha, indic_xi)
      
      if((equidistant %in% c(0,1)) & (typeof == 1)){
        if (missing(nb.int)) stop("Number of time interval for recurrences and terminal event 'nb.int' is required")
        if (length(nb.int) != 2) stop("The length of argument 'nb.int' should be 2. Must indicate for both recurrent events and terminal event.")
        if (class(nb.int) != "numeric") stop("The argument 'nb.int' must be a numeric")
        if (nb.int[1] < 1) stop("Number of time interval 'nb.int' must be between 1 and 20")
        if (nb.int[2] < 1) stop("Number of time interval 'nb.int' must be between 1 and 20")
        
        if (nb.int[1] > 20){
          nb.int[1] <- 20
          indic.nb.int1 <- 1 # equals 1 for nb.int1 > 20
        }
        else{
          indic.nb.int1 <- 0 # equals 1 for nb.int1 < 20
        }			
        if (nb.int[2] > 20){
          nb.int[2] <-20
          indic.nb.int2 <- 1 # equals 1 for nb.int1 > 20
        }
        else{
          indic.nb.int2 <- 0 # equals 1 for nb.int1 < 20
        }			
        nbintervR <- nb.int[1]
        size1 <- 3*nbintervR
        nbintervDC <- nb.int[2]
        size2 <- 3*nbintervDC
      }
      if ((typeof == 0) | (typeof == 2)){
        indic.nb.int1 <- 0
        indic.nb.int2 <- 0
      }      
      if (timedep==0){
        npbetatps1 <- 0
        npbetatps2 <- 0
      }
      else{
        npbetatps1 <- (betaknots+betaorder-1)*nvartimedep
        npbetatps2 <- (betaknots+betaorder-1)*nvartimedep2
      }
      npbetatps <- npbetatps1 + npbetatps2      
      np <- switch(as.character(typeof),
                   "0"=((nst+1) * (n.knots + 2) + nvarRec + nvarEnd + effet + indic_alpha + npbetatps),
                   "1"=(nst*nbintervR + nbintervDC + nvarRec + nvarEnd + effet + indic_alpha + npbetatps),
                   "2"=(2*(nst+1) + nvarRec + nvarEnd + effet + indic_alpha + npbetatps))
      
      if (all(all.equal(as.numeric(cens),terminal)==T)){
        stop("'Recurrent event' variable and 'Terminal event' variable need to be different")
      }
      
      # traitement de l'initialisation du Beta rentre par l'utilisateur
      Beta <- rep(0,np)
      if (!missing(init.B)){
        if (length(init.B) != nvar) stop("Wrong number of regression coefficients in init.B")
        if (timedep) stop("You can hardly know initial regression coefficient while time-varying effect")
        Beta <- c(rep(0,np-nvar),init.B)
      }
      if (!missing(init.Theta)){
        if (!is.numeric(init.Theta)) stop("init.Theta must be numeric")
        Beta[np-nvar-indic_alpha] <- init.Theta
      }
      if (!missing(init.Alpha)){
        if (!missing(Alpha)) stop("You can not both initialize alpha parameter and fit a joint model without it")
        if (!is.numeric(init.Alpha)) stop("init.Alpha must be numeric")
        Beta[np-nvar] <- init.Alpha
      }      
      xSu1 <- rep(0,100)
      xSu2 <- rep(0,100)
      if (typeof==0){
        mt11 <- size1
        mt12 <- size2
      }
      else{
        mt11 <- 100
        mt12 <- 100
      }      
      initialize <- 1      
      #npinit <- switch(as.character(typeof),
      # "0"=((n.knots + 2) + nvarRec + effet),
      # "1"=(nbintervR + nvarRec + nvarEnd + effet),
      # "2"=(2 + nvarRec + effet))      
      if ((uni.strat > 1 || joint.clust==2) & (joint.clust==0)) stop("stratification for clustered joint model is not yet allowed")
      if ((uni.strat > 1 || joint.clust==2) & (intcens)) stop("stratification for joint model with interval censored data is not yet allowed")
      if ((uni.strat > 1 || joint.clust==2) & (timedep)) stop("stratification for joint model and time-varying effect of covariates are not yet allowed")
      
      if ((typeof==0) & (length(kappa)!=(uni.strat+1))) stop("wrong length of argument kappa")
      
      ordretmp <- rep(0, n)
      
      flush.console()
      if (print.times){
        ptm<-proc.time()
        cat("\n")
        cat("Be patient. The program is computing ... \n")
      }  
      
      ans <- .Fortran(C_joint,                      
                      as.integer(n),
                      as.integer(c(length(uni.cluster),0, uni.strat)), #IJ 
                      #as.integer(uni.strat),
                      as.integer(strats),
                      as.integer(lignedc0),###
                      as.integer(n.knots),
                      #k0=as.double(kappa), # joint intcens,tps,cluster
                      axT=as.double(kappa), # joint avec generalisation de strate
                      as.double(tt0),
                      as.double(tt1),
                      
                      as.integer(cens),
                      as.integer(cluster),
                      as.integer(clusterdc),
                      as.integer(0),###
                      as.double(tt0.death),
                      as.double(tt1.death),
                      as.integer(terminalEvent),
                      as.double(tempdc),###
                      as.integer(icdc00),###
                      as.integer(nvarRec),
                      as.double(var),
                      
                      as.integer(nvarEnd),
                      as.double(vardc),
                      as.double(vaxdc00),###
                      as.integer(c(noVar1,noVar2)),
                      #as.integer(noVar2),
                      as.double(weights.agg),
                      as.integer(maxit),
                      np=as.integer(np),
                      b=as.double(Beta),
                      H=as.double(matrix(0,nrow=np,ncol=np)),
                      HIH=as.double(matrix(0,nrow=np,ncol=np)),
                      
                      loglik=as.double(0),
                      LCV=as.double(rep(0,2)),
                      xR=as.double(matrix(0,nrow=size1,ncol=uni.strat)),
                      lamR=as.double(array(0,dim=c(size1,3,uni.strat))),
                      xSuR=as.double(xSu1),
                      survR=as.double(array(0,dim=c(mt11,3,uni.strat))),
                      xD=as.double(rep(0,size2)),
                      lamD=as.double(matrix(0,nrow=size2,ncol=3)),
                      xSuD=as.double(xSu2),
                      survD=as.double(matrix(0,nrow=mt12,ncol=3)),
                      
                      as.integer(c(typeof, equidistant)),
                      #as.integer(equidistant),
                      as.integer(c(nbintervR, nbintervDC)),
                      #as.integer(nbintervDC),
                      as.integer(c(size1,size2,mt11,mt12)),###
                      counts=as.integer(c(0,0,0,0)),					  
                      IerIstop = as.integer(c(0,0)),
                      #ier=as.integer(0),
                      #istop=as.integer(0),
                      paraweib=as.double(rep(0,4)),
                      #			shape.weib=as.double(rep(0,2)),
                      #			scale.weib=as.double(rep(0,2)),
                      MartinGale=as.double(matrix(0,nrow=length(uni.cluster),ncol=5)),###
                      
                      linear.pred=as.double(rep(0,n)),
                      lineardc.pred=as.double(rep(0,as.integer(length(uni.cluster)))),
                      zi=as.double(rep(0,(n.knots+6))),
                      time=as.double(rep(0,(nbintervR+1))),
                      timedc=as.double(rep(0,(nbintervDC+1))),
                      # 			kendall=as.double(matrix(0,nrow=4,ncol=2)),
                      #			as.integer(initialize),
                      #			as.integer(npinit),
                      #			Bshared=as.double(rep(0,npinit)),
                      linearpredG=as.double(rep(0,lignedc0)),
                      joint.clust=as.integer(joint.clust),
                      as.integer(intcens),
                      as.integer(indices),
                      #as.integer(indic_alpha),
                      #as.integer(indic_xi),
                      # censure par intervalle, indic_alpha
                      as.double(ttU),
                      as.integer(ordretmp),
                      as.integer(initialize),
                      
                      logNormal=as.integer(logNormal),
                      paratps=as.integer(c(timedep,betaknots,betaorder)),
                      as.integer(c(filtretps,filtretps2)),
                      BetaTpsMat=as.double(matrix(0,nrow=101,ncol=1+4*nvartimedep)),
                      BetaTpsMatDc=as.double(matrix(0,nrow=101,ncol=1+4*nvartimedep2)),
                      EPS=as.double(c(LIMparam,LIMlogl,LIMderiv))
      )#,
      #PACKAGE = "frailtypack") # 65 arguments
      
      MartinGale <- matrix(ans$MartinGale,nrow=as.integer(length(uni.cluster)),ncol=5) 
      
      istop <- ans$IerIstop[2]
      
      if (istop == 4){
        warning("Problem in the loglikelihood computation. The program stopped abnormally. Please verify your dataset. \n")
      }      
      if (istop == 2){
        warning("Model did not converge.")
      }
      if (istop == 3){
        warning("Matrix non-positive definite.")
      }       
      #AD:
      if (noVar1==1 & noVar2==1) nvar<-0
      #AD:      
      np <- ans$np
      fit <- NULL
      fit$b <- ans$b
      fit$na.action <- attr(m, "na.action")
      fit$call <- call
      if (classofY == "SurvIC"){
        fit$n <- nobs
        if (typeofY == "intervaltronc") fit$indic.trunc <- 1
        else fit$indic.trunc <- 0
      }else{
        fit$n <- n
      }      
      if (joint.clust == 0) fit$ind <- lignedc0
      fit$groups <- length(uni.cluster)
      fit$n.events <- ans$counts[2]
      fit$n.deaths <- ans$counts[3]
      fit$n.censored<-ans$counts[4]      
      if(as.character(typeof)=="0"){
        fit$logLikPenal <- ans$loglik
      }
      else{
        fit$logLik <- ans$loglik
      }
      #AD:
      fit$LCV <- ans$LCV[1]
      fit$AIC <- ans$LCV[2]
      #AD:      
      if (logNormal == 0) fit$theta <- ans$b[np - nvar - npbetatps - indic_alpha]^2
      else fit$sigma2 <- ans$b[np - nvar - npbetatps - indic_alpha]^2      
      if (indic_alpha == 1) fit$alpha <- ans$b[np - nvar - npbetatps]
      if (joint.clust==2) fit$eta <- ans$b[np - nvar]^2
      fit$npar <- np      
      #AD:
      if ((noVar1==1 & noVar2==1)){
        fit$coef <- NULL
      }else{
        fit$coef <- ans$b[(np - nvar - npbetatps + 1):np]      
        noms <- c(factor.names(colnames(X)),factor.names(colnames(X2)))  
        if(length(grep(":",noms))>0)noms <- factor.names(noms)      
        if (timedep == 1){
          while (length(grep("timedep",noms))!=0){
            pos <- grep("timedep",noms)[1]
            noms <- noms[-pos]
            fit$coef <- fit$coef[-(pos:(pos+betaorder+betaknots-1))]
          }
        }
        names(fit$coef) <- noms
        #	if (missing(formula.terminalEvent)){
        #	   names(fit$coef) <- c(factor.names(colnames(X)))
        #	}else{
        #          names(fit$coef) <- c(factor.names(colnames(X)), factor.names(colnames(X2)))
        #      }
      }
      
      #AD:
      temp1 <- matrix(ans$H, nrow = np, ncol = np)
      temp2 <- matrix(ans$HIH, nrow = np, ncol = np)
      
      #Al:       
      fit$varHtotal <- temp1
      fit$varHIHtotal <- temp2      
      fit$varH <- temp1[(np - nvar - npbetatps - indic_alpha):np, (np - nvar - npbetatps - indic_alpha):np]
      fit$varHIH <- temp2[(np - nvar - npbetatps - indic_alpha):np, (np - nvar - npbetatps - indic_alpha):np]
      
      if (indic_alpha == 1)fit$alpha_p.value <- 1 - pchisq((fit$alpha/sqrt(diag(fit$varH))[2])^2,1)
      
      if (logNormal == 0){seH.frail <- sqrt(((2 * (fit$theta^0.5))^2) * diag(fit$varH)[1])
      fit$theta_p.value <- 1 - pnorm(fit$theta/seH.frail)
      }else{ seH.frail <- sqrt(((2 * (fit$sigma2^0.5))^2) * diag(fit$varH)[1])
      fit$sigma2_p.value <- 1 - pnorm(fit$sigma2/seH.frail)
      }
      
      
      if (indic_alpha == 1) noms <- c("theta","alpha",factor.names(colnames(X)),factor.names(colnames(X2)))
      else noms <- c("theta",factor.names(colnames(X)),factor.names(colnames(X2)))
      if(length(grep(":",noms))>0)noms <- factor.names(noms)
      if (timedep == 1){ # on enleve les variances des parametres des B-splines
        while (length(grep("timedep",noms))!=0){
          pos <- grep("timedep",noms)[1]
          noms <- noms[-pos]
          fit$varH <- fit$varH[-(pos:(pos+betaorder+betaknots-1)),-(pos:(pos+betaorder+betaknots-1))]
          fit$varHIH <- fit$varHIH[-(pos:(pos+betaorder+betaknots-1)),-(pos:(pos+betaorder+betaknots-1))]
        }
      }
      fit$nvar<-c(nvarRec,nvarEnd)
      fit$nvarnotdep<-c(nvarRec-nvartimedep,nvarEnd-nvartimedep2)
      #	fit$formula <- formula(Terms)
      
      fit$xR <- matrix(ans$xR, nrow = size1, ncol = uni.strat)
      fit$lamR <- if(typeof == 1){array(ans$lamR[seq(1,length(ans$lamR),3)], dim = c(nb.int[1],3,uni.strat))} else{array(ans$lamR, dim = c(size1,3,uni.strat))} 
      fit$xSuR <- matrix(ans$xSuR, nrow = 100, ncol = uni.strat)
      fit$survR <- array(ans$survR, dim = c(mt11,3,uni.strat))
      
      fit$xD <- ans$xD
      fit$lamD <- if(typeof == 1) {matrix(ans$lamD[seq(1,length(ans$lamD),3)], nrow = nb.int[2], ncol = 3)} else{matrix(ans$lamD, nrow = size2, ncol = 3)}
      fit$xSuD <- ans$xSuD
      fit$survD <- matrix(ans$survD, nrow = mt12, ncol = 3)
      
      
      fit$type <- type
      fit$n.strat <- uni.strat
      fit$n.iter <- ans$counts[1]
      fit$typeof <- typeof
      if (typeof == 0){
        fit$n.knots<-n.knots
        fit$kappa <- ans$axT
        fit$cross.Val<-cross.validation
        fit$n.knots.temp <- n.knots.temp
        fit$zi <- ans$zi
      }
      if(typeof == 1){
        fit$time <- ans$time
        fit$timedc <- ans$timedc
      }
      
      medianR <- NULL
      for (i in (1:fit$n.strat)) medianR[i] <- ifelse(typeof==0, minmin(fit$survR[,1,i],fit$xR), minmin(fit$survR[,1,i],fit$xSuR))
      lowerR <- NULL
      for (i in (1:fit$n.strat)) lowerR[i] <- ifelse(typeof==0, minmin(fit$survR[,2,i],fit$xR), minmin(fit$survR[,2,i],fit$xSuR))
      upperR <- NULL
      for (i in (1:fit$n.strat)) upperR[i] <- ifelse(typeof==0, minmin(fit$survR[,3,i],fit$xR), minmin(fit$survR[,3,i],fit$xSuR))
      fit$medianR <- cbind(lowerR,medianR,upperR)
      
      medianD <- ifelse(typeof==0, minmin(fit$survD[,1],fit$xD), minmin(fit$survD[,1],fit$xSuD))
      lowerD <- ifelse(typeof==0, minmin(fit$survD[,3],fit$xD), minmin(fit$survD[,3],fit$xSuD))
      upperD <- ifelse(typeof==0, minmin(fit$survD[,2],fit$xD), minmin(fit$survD[,2],fit$xSuD))
      fit$medianD <- cbind(lowerD,medianD,upperD)
      
      #AD:
      fit$noVar1 <- noVar1
      fit$noVar2 <- noVar2
      fit$nbintervR <- nbintervR
      fit$nbintervDC <- nbintervDC
      fit$nvarRec <- nvarRec
      fit$nvarEnd <- nvarEnd
      fit$istop <- istop
      fit$indic.nb.intR <- indic.nb.int1
      fit$indic.nb.intD <- indic.nb.int2
      fit$shape.weib <- ans$paraweib[1:2]#ans$shape.weib
      fit$scale.weib <- ans$paraweib[3:4]#ans$scale.weib
      #AD:
      
      # verif que les martingales ont ete bien calculees
      msg <- "Problem in the estimation of the random effects (perhaps high number of events in some clusters)"
      if (Frailty){
        if (any(MartinGale[,1]==0)){
          fit$martingale.res <- msg
          fit$martingaledeath.res <- msg          
          fit$frailty.pred <- msg
          #		fit$frailty.var <- msg          
          fit$linear.pred <- msg
          fit$lineardeath.pred <- msg
        }
        else{
          fit$martingale.res <- MartinGale[,1]#ans$martingale.res
          fit$martingaledeath.res <- MartinGale[,2]#ans$martingaledc.res          
          fit$frailty.pred <- MartinGale[,3]#ans$frailty.pred
          #		fit$frailty.var <- MartinGale[,4]#ans$frailty.var          
          fit$linear.pred <- ans$linear.pred[order(ordre)]
          if (joint.clust==0) fit$lineardeath.pred <- ans$linearpredG 
          else fit$lineardeath.pred <- ans$lineardc.pred 
        }
      }      
      #    if (joint.clust==0){
      #        fit$kendall <- matrix(ans$kendall,nrow=4,ncol=2)
      #    }
      fit$joint.clust <- ans$joint.clust
      fit$AG <- recurrentAG
      fit$intcens <- intcens # rajout
      
      fit$indic_alpha <- indic_alpha
      if(joint.clust==2)fit$indic_alpha <- 0
      fit$logNormal <- ans$logNormal
      fit$BetaTpsMat <- matrix(ans$BetaTpsMat,nrow=101,ncol=1+4*nvartimedep)
      fit$BetaTpsMatDc <- matrix(ans$BetaTpsMatDc,nrow=101,ncol=1+4*nvartimedep2)
      fit$nvartimedep <- c(nvartimedep,nvartimedep2)
      
      fit$Names.vardep <- vardep
      fit$Names.vardepdc <- vardep2
      
      fit$EPS <- ans$EPS     
      
      #================================> For the reccurrent
      #========================= Test de Wald
      
      if ((length(vec.factor) > 0) & (timedep == 0)){
        Beta <- ans$b[(np - nvar + 1):np]
        if (indic_alpha == 1) VarBeta <- diag(diag(fit$varH)[-c(1,2)])
        else VarBeta <- diag(diag(fit$varH)[-c(1)])
        nfactor <- length(vec.factor)
        p.wald <- rep(0,nfactor)
        ntot <- nvarEnd + nvarRec
        
        if(fit$istop == 1) fit$global_chisq <- waldtest(N=nvarRec,nfact=nfactor,place=ind.place,modality=occur,b=Beta,Varb=VarBeta,Llast=nvarEnd,Ntot=ntot)
        else fit$global_chisq <- 0
        
        fit$dof_chisq <- occur
        fit$global_chisq.test <- 1
        # Calcul de pvalue globale
        for(i in 1:length(vec.factor)){
          p.wald[i] <- signif(1 - pchisq(fit$global_chisq[i], occur[i]), 3)
        }
        fit$p.global_chisq <- p.wald
        fit$names.factor <- vec.factor
      }
      else{
        fit$global_chisq.test <- 0
      }
      
      #================================> For the death
      #========================= Test de Wald
      
      if (!missing(formula.terminalEvent)){
        if ((length(vec.factordc) > 0) & (timedep == 0)){
          Beta <- ans$b[(np - nvar + 1):np]
          #if (indic_alpha == 1) VarBeta <- diag(diag(fit$varH)[-c(1,2)])
          #else VarBeta <- diag(diag(fit$varH)[-c(1)])
          if(indic_alpha == 1) VarBeta <- fit$varH[3:dim(fit$varH)[1],3:dim(fit$varH)[2]]
          else VarBeta <- fit$varH[2:dim(fit$varH)[1],2:dim(fit$varH)[2]]
          nfactor <- length(vec.factordc)
          p.walddc <- rep(0,nfactor)
          ntot <- nvarEnd + nvarRec
          
          if(fit$istop == 1)	fit$global_chisq_d <- waldtest(N=nvarEnd,nfact=nfactor,place=ind.placedc,modality=occurdc,b=Beta,Varb=VarBeta,Lfirts=nvarRec,Ntot=ntot)
          else fit$global_chisq_d <- 0
          
          fit$dof_chisq_d <- occurdc
          fit$global_chisq.test_d <- 1
          # Calcul de pvalue globale
          for(i in 1:length(vec.factordc)){
            p.walddc[i] <- signif(1 - pchisq(fit$global_chisq_d[i], occurdc[i]), 3)
          }
          fit$p.global_chisq_d <- p.walddc
          fit$names.factordc <- vec.factordc
        }
        else{
          fit$global_chisq.test_d <- 0
        }
      }
      else{
        fit$global_chisq.test_d <- 0
      }
      if(!is.null(fit$coef)){
        if (indic_alpha == 1 || fit$joint.clust==2) {
          seH <- sqrt(diag(fit$varH))[-c(1,2)]
        }else{
          seH <- sqrt(diag(fit$varH))[-1]
        }
        fit$beta_p.value <- 1 - pchisq((fit$coef/seH)^2, 1)
      }
      
      if (length(Xlevels) >0)	fit$Xlevels <- Xlevels
      fit$contrasts <- contr.save
      if (length(Xlevels2) >0) fit$Xlevels2 <- Xlevels2
      fit$contrasts2 <- contr.save2
      
      #NCC design
      fit$ncc <- FALSE
      if(length(wts))fit$ncc <- TRUE
      
      
      fit$formula <- formula
      fit$formula.terminalEvent <- formula.terminalEvent
      
      attr(fit,"joint")<-joint
      attr(fit,"subcluster")<-FALSE
      class(fit) <- "jointPenal"
    }  # End JOINT MODEL    
    
    #
    # Begin NESTED MODEL
    #
    
    effet <- 1
    # Modified ML 24/03/2015 for Nested Joint model
    if (length(subcluster) & !joint){
      if (logNormal == 1) stop("Nested model not implemented yet for log normal distribution of frailties")      
      if((equidistant %in% c(0,1)) & (typeof == 1)){
        if (missing(nb.int)) stop("Number of time interval 'nb.int' is required")
        if (length(nb.int) != 1) stop("Wrong length of number of time interval argument 'nb.int'")
        if (class(nb.int) != "numeric") stop("The argument 'nb.int' must be a numeric")
        if (nb.int < 1) stop("Number of time interval 'nb.int' must be between 1 and 20")
        if (nb.int > 20){
          nb.int <-20
          indic.nb.int <- 1 # equals 1 for nb.int > 20
        }
        else{
          indic.nb.int <- 0 # equals 1 for nb.int < 20
        }
        nbintervR <- nb.int
        size1 <- 3*nbintervR
      }
      if ((typeof == 0) | (typeof == 2)) indic.nb.int <- 0
      if (sum(as.double(var))==0) nvar <- 0      
      np <- switch(as.character(typeof),
                   "0"=(as.integer(uni.strat) * (as.integer(n.knots) + 2) + as.integer(nvar) + 2 * as.integer(Frailty)),
                   
                   "1"=(as.integer(uni.strat) * nbintervR + nvar + 2 * as.integer(Frailty)),
                   
                   "2"=(as.integer(uni.strat) * 2 + nvar + 2 * as.integer(Frailty)))      
      xSu1 <- rep(0,100)
      xSu2 <- rep(0,100)
      if (typeof==0){
        mt1 <- size1
      }
      else{
        mt1 <- 100
      }
      size2 <- mt1      
      if (length(kappa)==1) kappa <- c(kappa,0)      
      ########### group and subgroup
      grpe <- function(g){        
        grp <- unique(g)        
        res <- rep(0,length(grp))        
        for(i in 1:length(res)){
          res[i] = sum(grp[i]==g)
        }
        return(res)
      }      
      grp <- grpe(as.integer(cluster))      
      subgrpe <- function(g,sg){        
        j <- 0
        k <- 0
        res <- rep(0,length(g))			
        for(i in 1:length(g)){
          k <- k + g[i]
          j <- j + 1
          temp <- sg[j:k]
          res[i] <- length(grpe(temp))
          j <- k
        }
        return(res)
      }      
      subgbyg <- subgrpe(grp,as.integer(subcluster))      
      maxng <- max(subgbyg)
      ngg <- length(uni.cluster)    
      
      flush.console()
      if (print.times){
        ptm<-proc.time()
        cat("\n")
        cat("Be patient. The program is computing ... \n")
      }	
      
      
      ans <- .Fortran(C_nested,
                      as.integer(n),
                      as.integer(length(uni.cluster)),
                      as.integer(length(uni.subcluster)),
                      as.integer(uni.strat),
                      as.integer(n.knots),
                      as.double(kappa),
                      as.double(tt0),
                      as.double(tt1),
                      as.integer(cens),
                      as.integer(cluster),
                      
                      as.integer(subcluster),
                      as.integer(nvar),
                      as.double(strats),
                      as.double(var),
                      as.integer(AG),
                      as.integer(noVar1),
                      as.integer(maxit),
                      as.integer(crossVal),
                      as.integer(np),
                      as.integer(maxng),
                      
                      b=as.double(rep(0,np)),
                      H=as.double(matrix(0,nrow=np,ncol=np)),
                      HIH=as.double(matrix(0,nrow=np,ncol=np)),
                      loglik=as.double(0),
                      LCV=as.double(rep(0,2)),
                      x1=as.double(rep(0,size1)),
                      lam=as.double(matrix(0,nrow=size1,ncol=3)),
                      xSu1=as.double(xSu1),
                      surv=as.double(matrix(0,nrow=size2,ncol=3)),
                      x2=as.double(rep(0,size1)),
                      
                      lam2=as.double(matrix(0,nrow=size1,ncol=3)),
                      xSu2=as.double(xSu2),
                      surv2=as.double(matrix(0,nrow=size2,ncol=3)),
                      as.integer(typeof),
                      as.integer(equidistant),
                      as.integer(nbintervR),
                      as.integer(size1),
                      ni=as.integer(0),
                      cpt=as.integer(0),
                      ier=as.integer(0),
                      
                      k0=as.double(c(0,0)),
                      ddl=as.double(0),
                      istop=as.integer(0),
                      shape.weib=as.double(rep(0,2)),
                      scale.weib=as.double(rep(0,2)),
                      as.integer(mt1),
                      zi=as.double(rep(0,(n.knots+6))),
                      time=as.double(rep(0,(nbintervR+1))),
                      martingale.res=as.double(rep(0,as.integer(length(uni.subcluster)))),
                      frailty.pred.group=as.double(rep(0,as.integer(length(uni.cluster)))),
                      
                      frailty.pred.subgroup=as.double(matrix(0,nrow=ngg,ncol=maxng)),
                      frailty.var.group=as.double(rep(0,as.integer(length(uni.cluster)))),
                      frailty.var.subgroup=as.double(matrix(0,nrow=ngg,ncol=maxng)),
                      frailty.sd.group=as.double(rep(0,as.integer(length(uni.cluster)))),
                      frailty.sd.subgroup=as.double(matrix(0,nrow=ngg,ncol=maxng)),
                      linear.pred=as.double(rep(0,n)),
                      EPS=as.double(c(LIMparam,LIMlogl,LIMderiv))
      )#,
      #PACKAGE = "frailtypack") # 57 arguments
      
      if (ans$istop == 4){
        warning("Problem in the loglikelihood computation. The program stopped abnormally. Please verify your dataset. \n")
      }      
      if (ans$istop == 2){
        warning("Model did not converge.")
      }
      if (ans$istop == 3){
        warning("Matrix non-positive definite.")
      }      
      nst <- as.integer(uni.strat)      
      if (noVar1 == 1) nvar<-0      
      np <- np
      fit <- NULL
      fit$b <- ans$b
      fit$na.action <- attr(m, "na.action")
      fit$call <- call
      fit$n <- n
      fit$groups <- length(uni.cluster)
      fit$subgroups <- length(uni.subcluster)
      fit$n.events <- ans$cpt
      if(as.character(typeof)=="0"){
        fit$logLikPenal <- ans$loglik
      }
      else{
        fit$logLik <- ans$loglik
      }      
      fit$alpha<-ans$b[np-nvar-1]^2
      fit$eta<-ans$b[np-nvar]^2      
      if (noVar1 == 1) {
        fit$coef <- NULL
      }
      else{
        fit$coef <- ans$b[(np - nvar + 1):np]
        names(fit$coef) <- factor.names(colnames(X))
      }      
      temp1 <- matrix(ans$H, nrow = np, ncol = np)
      temp2 <- matrix(ans$HIH, nrow = np, ncol = np)
      
      fit$varH <- temp1[(np - nvar - 1):np, (np - nvar - 1):np]
      fit$varHIH <- temp2[(np - nvar - 1):np, (np - nvar - 1):np]
      
      seH.alpha <- sqrt(((2 * (fit$alpha^0.5))^2) * diag(fit$varH)[1])
      fit$alpha_p.value <- 1 - pnorm(fit$alpha/seH.alpha)
      
      seH.eta <- sqrt(((2 * (fit$eta^0.5))^2) * diag(fit$varH)[2])
      fit$eta_p.value <- 1 - pnorm(fit$eta/seH.eta)
      
      
      #	fit$formula <- formula(Terms)
      
      fit$x <- cbind(ans$x1,ans$x2)
      fit$lam <- if(typeof == 1) {array(c(ans$lam[seq(1,length(ans$lam),3)],ans$lam2[seq(1,length(ans$lam2),3)]), dim=c(nb.int,3,2))} else{array(c(ans$lam,ans$lam2), dim=c(size1,3,2))}
      fit$xSu <- cbind(ans$xSu1,ans$xSu2)
      fit$surv <- array(c(ans$surv,ans$surv2), dim=c(size2,3,2))
      
      fit$type <- type
      fit$n.strat <- uni.strat
      fit$n.iter <- ans$ni
      fit$typeof <- typeof
      fit$noVar1 <- noVar1
      
      median <- NULL
      for (i in (1:fit$n.strat)) median[i] <- ifelse(typeof==0, minmin(fit$surv[,1,i],fit$x), minmin(fit$surv[,1,i],fit$xSu))
      lower <- NULL
      for (i in (1:fit$n.strat)) lower[i] <- ifelse(typeof==0, minmin(fit$surv[,3,i],fit$x), minmin(fit$surv[,3,i],fit$xSu))
      upper <- NULL
      for (i in (1:fit$n.strat)) upper[i] <- ifelse(typeof==0, minmin(fit$surv[,2,i],fit$x), minmin(fit$surv[,2,i],fit$xSu))
      fit$median <- cbind(lower,median,upper)
      
      if (typeof == 0){
        fit$n.knots<-n.knots
        if (uni.strat > 1) fit$kappa <- ans$k0
        else fit$kappa <- ans$k0[1]
        fit$DoF <- ans$ddl
        fit$cross.Val<-cross.validation
        fit$zi <- ans$zi
      }
      if(typeof == 1)fit$time <- ans$time
      #AD:
      fit$nbintervR <- nbintervR
      fit$nvar <- nvar
      fit$LCV <- ans$LCV[1]
      fit$AIC <- ans$LCV[2]
      fit$npar <- np
      fit$nst <- nst
      if (typeof == 0){
        #     	fit$indic.Kappa2 <- indic.Kappa2
        fit$n.knots.temp <- n.knots.temp
      }
      fit$indic.nb.int <- indic.nb.int
      fit$istop <- ans$istop
      fit$shape.weib <- ans$shape.weib
      fit$scale.weib <- ans$scale.weib
      fit$AG <- recurrentAG
      fit$EPS <- ans$EPS
      
      #   if (Frailty){      
      msg <- "Problem in the estimation of the random effects (perhaps high number of events in some clusters)"
      if (any(ans$martingale.res==0)){        
        fit$martingale.res <- msg
        fit$frailty.pred.group <- msg
        fit$frailty.pred.subgroup <- msg
        fit$linear.pred <- msg        
      }
      else{        
        fit$martingale.res <- ans$martingale.res
        fit$frailty.pred.group <- ans$frailty.pred.group        
        nom1 <- paste("g_",c(1:ngg),sep="")
        nom2 <- paste("sub_g",c(1:maxng))
        
        frailty.pred.subgroup <- as.data.frame(matrix(round(ans$frailty.pred.subgroup,6),ncol=maxng))
        rownames(frailty.pred.subgroup) <- nom1
        colnames(frailty.pred.subgroup) <- nom2
        for (i in 1:ngg) {
          if (subgbyg[i] < max(subgbyg)) {
            frailty.pred.subgroup[i,(subgbyg[i]+1):max(subgbyg)] <- "."
          }
        }
        fit$frailty.pred.subgroup <- frailty.pred.subgroup			
        #	fit$frailty.var.group <- ans$frailty.var.group			
        # 	frailty.var.subgroup <- as.data.frame(matrix(round(ans$frailty.var.subgroup,6),nc=maxng))
        # 	rownames(frailty.var.subgroup) <- nom1
        # 	colnames(frailty.var.subgroup) <- nom2
        # 
        # 	if(sum(which(subgbyg < max(subgbyg)))>0)frailty.var.subgroup[which(subgbyg < max(subgbyg)),(subgbyg[which(subgbyg < max(subgbyg))]+1):max(subgbyg)] <- "."
        
        #	fit$frailty.var.subgroup <- frailty.var.subgroup			
        #	fit$frailty.sd.group <- ans$frailty.sd.group			
        #	frailty.sd.subgroup <- as.data.frame(matrix(round(ans$frailty.sd.subgroup,6),nc=maxng))
        #	rownames(frailty.sd.subgroup) <- nom1
        #	colnames(frailty.sd.subgroup) <- nom2
        #	if(sum(which(subgbyg < max(subgbyg)))>0)frailty.sd.subgroup[which(subgbyg < max(subgbyg)),(subgbyg[which(subgbyg < max(subgbyg))]+1):max(subgbyg)] <- "."
        
        #	fit$frailty.sd.subgroup <- frailty.sd.subgroup
        fit$linear.pred <- ans$linear.pred[order(ordre)]
      }
      fit$subgbyg <- subgbyg
      #   }
      if(ans$ier==2000)
        stop("The cross validation procedure cannot be finished. Try to change
             either the number of knots or the seed for kappa parameter")
      
      
      #========================= Test de Wald pour nested
      
      if(length(vec.factor) > 0){
        Beta <- ans$b[(np - nvar + 1):np]
        #VarBeta <- fit$varH[2:(nvar+1),2:(nvar+1)]	
        VarBeta <- fit$varH[3:dim(fit$varH)[1], 3:dim(fit$varH)[2]]
        nfactor <- length(vec.factor)
        p.wald <- rep(0,nfactor)
        
        if (fit$istop == 1) fit$global_chisq <- waldtest(N=nvar,nfact=nfactor,place=ind.place,modality=occur,b=Beta,Varb=VarBeta)
        else fit$global_chisq <- 0
        
        fit$dof_chisq <- occur
        fit$global_chisq.test <- 1
        # Calcul de pvalue globale
        for(i in 1:length(vec.factor)){
          p.wald[i] <- signif(1 - pchisq(fit$global_chisq[i], occur[i]), 3)
        }
        fit$p.global_chisq <- p.wald
        fit$names.factor <- vec.factor
      }
      else{
        fit$global_chisq.test <- 0
      }
      
      if(!is.null(fit$coef)){
        if(length(fit$coef)>1) seH <- sqrt(diag(fit$varH))[-c(1:2)]
        else seH <- sqrt(fit$varH)[-c(1:2)]
        
        fit$beta_p.value <- 1 - pchisq((fit$coef/seH)^2, 1)
      }
      
      if (length(Xlevels) >0) fit$Xlevels <- Xlevels
      fit$contrasts <- contr.save
      
      fit$formula <- formula
      
      attr(fit,"joint")<-joint
      attr(fit,"subcluster")<-TRUE
      class(fit) <- "nestedPenal"      
    } # End NESTED MODEL
    
    # ===================================================
    # BEGIN JOINT NESTED MODEL 
    # ===================================================	
    if (length(subcluster) & joint){			
      #if (!Frailty) stop("For joint nested frailty models, 'Frailty' must be equal to 'TRUE' !")
      if (missing(formula.terminalEvent)) stop ("For joint nested frailty model, 'formula.terminalEvent' is required !")
      if (joint.clust != 3) stop("Argument is mispecified for joint nested frailty model. Please look at the frailtyPenal documentation.")
      if (logNormal) stop("Sorry but log normal distribution is not available for joint nested frailty models")
      if (uni.strat > 1) stop("Sorry but stratification for joint nested frailty model is not allowed")
      
      if (classofY == "Surv")
      {
        if (!recurrentAG)
        {
          tt1.death <- aggregate(tt1, by=list(subcluster), FUN=sum)[,2]
          tt0.death <- rep(0,length(tt1.death))
          clusterdc <- 0
          lignedc0 <- 0
          tempdc <- 0
        }
        else{
          tt1.death <- aggregate(tt1, by=list(subcluster), FUN=function(x) x[length(x)])[,2] # Myriam ! modifie sort()
          tt0.death <- rep(0,length(tt1.death))
          clusterdc <- 0
          lignedc0 <- 0
          tempdc <- 0
          
        }
        #stop ("The counting process approach of Andersen and Gill with a calendar timescale for recurrent event times is not yet allowed")}
      }
      else stop("Interval-censored data are not allowed for joint nested model.")
      
      Terms2 <- if(missing(data))
      {
        terms(formula.terminalEvent, special)
      }else{
        terms(formula.terminalEvent, special, data = data)
      }
      ord2 <- attr(Terms2, "order")
      if (length(ord2) & any(ord2 != 1)) stop("Interaction terms are not valid for terminalEvent formula")
      
      if (ind.terminal)
      {
        terminalEvent<-aggregate(terminal,by=list(subcluster),FUN=function(x) x[length(x)])[,2]
        icdc00 <- 0
      }else{
        stop("Nested joint frailty model mispecified. Please look at the frailtyPenal documentation.")
      }
      
      
      # TerminalEvent doit etre de la forme 0-1
      if(!all(terminalEvent %in% c(2,1,0))) stop("'terminal' must contain a variable coded 0-1 and a non-factor variable")
      m2 <- match.call(expand.dots = FALSE)
      
      m2$formula.terminalEvent <- m2$n.knots <- m2$recurrentAG <- m2$cross.validation <- m2$jointGeneral<- m2$kappa <- m2$maxit <- m2$hazard <- m2$nb.int <- m2$RandDist <- m2$betaorder <- m2$betaknots <- m2$init.B <- m2$LIMparam <- m2$LIMlogl <- m2$LIMderiv <- m2$print.times <- m2$init.Theta <- m2$init.Alpha <- m2$Alpha <- m2$init.Ksi <- m2$Ksi <- m2$init.Eta <- m2$Eta <- m2$initialize <- m2$... <- NULL
      
      m2$formula <- Terms2
      m2[[1]] <- as.name("model.frame")
      m2 <- eval(m2, sys.parent())
      
      match.noNA <- dimnames(m2)[[1]]%in%dimnames(m)[[1]]
      m2 <- m2[match.noNA, ,drop=FALSE]
      
      newTerms2 <- Terms2
      
      #if (!missing(formula.terminalEvent))		{
      X2 <- model.matrix(newTerms2, m2)
      lldc <- attr(newTerms2,"term.labels")
      ind.placedc <- unique(attr(X2,"assign")[duplicated(attr(X2,"assign"))])#changement unique le 26/09/2014
      
      vec.factordc <- NULL
      vec.factordc <- c(vec.factordc,lldc[ind.placedc])			
      mat.factordc <- matrix(vec.factordc,ncol=1,nrow=length(vec.factordc))
      
      # Fonction servant a recuperer les termes dependants de "as.factor" en ce qui concerne l'evenement terminal
      vec.factordc <-apply(mat.factordc,MARGIN=1,FUN=function(x){
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
            return(substr(x,start=pos1,stop=pos2))
          }
        }else{
          return(x)
        }
      }
      )
      
      # On determine le nombre de categories pour chaque var categorielle
      if(length(vec.factordc) > 0)
      {
        vect.factdc <- attr(X2,"dimnames")[[2]]
        vect.factdc <- vect.factdc[grep(paste(vec.factordc,collapse="|"),vect.factdc)]			  
        occurdc <- rep(0,length(vec.factordc))
        
        interactiondc<-as.vector(apply(matrix(vect.factdc,nrow=length(vect.factdc)),MARGIN=1,FUN=function(x){length(grep(":",unlist(strsplit(x,split=""))))}))
        which.interactiondc <- which(interactiondc==1)
        
        for(i in 1:length(vec.factordc))
        {			  
          if(length(grep(":",unlist(strsplit(vec.factordc[i],split=""))))>0)
          {				
            pos <- grep(":",unlist(strsplit(vec.factordc[i],split="")))
            length.grep <- 0
            for(j in 1:length(vect.factdc))
            {
              if(j%in%which.interactiondc)
              {
                if(length(grep(substr(vec.factordc[i],start=1,stop=pos-1),vect.factdc[j]))>0 && length(grep(substr(vec.factordc[i],start=pos+1,stop=length(unlist(strsplit(vec.factordc[i],split="")))),vect.factdc[j]))>0)
                {
                  length.grep <- length.grep + 1
                  which <- i
                }
              }
            }
            occurdc[i] <- length.grep				
          }else{				
            if(length(vect.factdc[-which.interactiondc])>0)	occurdc[i] <- length(grep(vec.factordc[i],vect.factdc[-which.interactiondc]))
            else occurdc[i] <- length(grep(vec.factordc[i],vect.factdc))
          }
        }
      }	   
      #=========================================================>			
      assign <- lapply(attrassign(X2, newTerms2)[-1], function(x) x - 1)
      Xlevels2 <- .getXlevels(newTerms2, m2)
      contr.save2 <- attr(X2, 'contrasts')
      #========================================>
      if(length(vec.factordc) > 0){
        positiondc <- unlist(assign,use.names=F)
      }
      #========================================>
      if (ncol(X2) == 1)
      {
        X2<-X2-1
        noVar2 <- 1
      }else{
        X2 <- X2[, -1, drop = FALSE]
        noVar2 <- 0
      }			
      nvar2 <- ncol(X2)	
      
      if ((!is.null(attr(Terms2, "specials")$timedep))|timedep) stop("Time-dependant covariates are not allowed for joint nested frailty model.")		
      
      #vartimedep2 <- if (!is.null(attr(Terms2, "specials")$timedep)) stop("Nested frailty joint model not yet able to account time-dependant covariate")
      # if (!is.null(vartimedep2)) timedep <- 1
      #varnotdep2 <- colnames(X2)[-grep("timedep", colnames(X2))]
      vardep2 <- colnames(X2)[grep("timedep", colnames(X2))]
      vardep2 <- apply(matrix(vardep2, ncol=1, nrow=length(vardep2)), 1, timedep.names)
      
      nvartimedep2 <- length(vardep2)
      filtretps2 <- rep(0, nvar2)		
      #filtretps2[grep("timedep", colnames(X2))] <- 1	
      
      vardc.temp<-matrix(c(X2),nrow=nrow(X2),ncol=nvar2)   
      
      #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      #		Traitement des valeurs manquantes
      #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      if(is.null(nrow(m2)))
      {
        if (length(m2) != nrow(m)){
          stop("There are missing values in the covariates modelling the terminal event. \n Prepare data only with complete cases")
        }
      }else{          
        if (nrow(m2) != nrow(m)){
          stop(" There are missing values in the covariates modelling the terminal event. \n Prepare data only with complete cases")
        }          
      }        
      
      if (!is.null(ncol(vardc.temp)))
      {
        vardc<-aggregate(vardc.temp[,1],by=list(sort(subcluster)), FUN=function(x) x[length(x)])[,2]
        if (ncol(vardc.temp)>1)
        {
          for (i in 2:ncol(vardc.temp))
          {
            vardc.i<-aggregate(vardc.temp[,i],by=list(sort(subcluster)), FUN=function(x) x[length(x)])[,2]
            vardc<-cbind(vardc,vardc.i)
          }
        }
      }else vardc<-aggregate(vardc.temp,by=list(sort(subcluster)), FUN=function(x) x[length(x)])[,2]
      
      vaxdc00 <- 0
      
      #} # if ! missing formula.terminalEvent
      
      nvarRec <- nvar
      
      #=======================================>
      #======= Construction du vecteur des indicatrice
      
      if(length(vec.factordc) > 0){
        k <- 0
        for(i in 1:length(vec.factordc)){
          ind.placedc[i] <- ind.placedc[i]+k
          k <- k + occurdc[i]-1
        }
      }
      
      if(is.null(nrow(vardc))) nvarEnd <- 1
      else nvarEnd <- ncol(vardc)
      
      if (sum(as.double(var)) == 0) nvarRec <- 0
      if (sum(as.double(vardc)) == 0) nvarEnd <- 0
      nvar <- nvarRec + nvarEnd
      
      # ... End preparing data 
      
      effet <- 2
      indic_alpha <- 1
      indic_xi <- 1
      
      if (!missing(Alpha))
      {
        if (Alpha == "None") indic_alpha <- 0
        else stop("Alpha can only take 'None' as a value in this version of frailtypack package ")
      }
      if (!missing(Ksi))
      {	
        if (Ksi == "None") indic_xi <- 0
        else stop("Ksi can only take 'None' as a value in this version of frailtypack package")
      }
      
      indices <- c(indic_alpha, indic_xi)
      
      nst <- uni.strat 
      if (typeof == 1) stop ("!Warning! Piecewise baseline hazard function is not yet allowed for joint nested frailty model.")
      else {
        indic.nb.int1 <- 0
        indic.nb.int2 <- 0
      }		
      np <- switch(as.character(typeof),
                   '0' = ((nst+1) * (n.knots + 2) + nvar + effet + indic_alpha + indic_xi), #splines
                   '2' = (2*(nst+1) + nvar + effet + indic_alpha + indic_xi)) # weibull
      if (all(all.equal(as.numeric(cens), terminal) == T)) stop("'Recurrent event' variable and 'Terminal event' variable need to be different")
      
      # Initialisation des parametres beta 
      #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      Beta <- rep(0,np)		
      if(!missing(init.B)){
        if (length(init.B) != nvar) stop ("Wrong number of regression coefficient in init.B")
        Beta <- c(rep(0, np-nvar), init.B) 
      }		
      if (!missing(init.Ksi)){
        if(!missing(Ksi)) stop("You can't both initialize Ksi parameter and fit a joint model without Ksi in joint nested frailty model.")
        if(!is.numeric(init.Ksi)) stop ("'init.Ksi' must be numeric")
        Beta[np-nvar] <- init.Ksi
      }		
      if (!missing(init.Alpha)){
        if(!missing(Alpha)) stop("You can't both initialize alpha parameter and fit a joint model without it")
        if(!is.numeric(init.Alpha)) stop ("'init.Alpha' must be numeric")
        Beta[np-nvar-indic_xi] <- init.Alpha 
      }		
      if (!missing(init.Eta)){
        if (!is.numeric(init.Eta)) stop("'init.Eta' must be numeric")
        Beta[np-nvar-indic_xi-indic_alpha] <- init.Eta
      }		
      if (!missing(init.Theta)){
        if (!is.numeric(init.Theta)) stop("'init.Theta' must be numeric")
        Beta[np-nvar-indic_xi-indic_alpha-1] <- init.Theta
      }
      #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
      # Clusterfamdc pour les deces : doit avoir autant de valeur que d'individus : 
      clusterfamdc <- aggregate(cluster, by=list(subcluster), FUN = function(x) x[length(x)])[,2]
      
      xSu1 <- rep(0,100)
      xSu2 <- rep(0,100)
      
      if(typeof==0){
        mt11 <- size1
        mt12 <- size2
      }else{
        mt11 <- 100
        mt12 <- 100
      }
      
      if ((typeof == 0) & (length(kappa) != (uni.strat+1))) stop ("Wrong length of argument kappa")
      
      weights.vec <- rep(1, nrow(data))
      weights.agg <- aggregate(weights.vec, by=list(cluster), FUN=function(x) x[length(x)])[,2]	
      
      #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      # Fin de preparation des arguments 
      #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%	
      
      flush.console()
      if (print.times){
        ptm<-proc.time()
        cat("\n")
        cat("Be patient. The program is computing ... \n")
      }	
      
      ans <- .Fortran(C_joint, 
                      as.integer(n), 		
                      as.integer(c(length(uni.subcluster),length(uni.cluster),uni.strat)),
                      as.integer(strats), 
                      as.integer(lignedc0), 
                      as.integer(n.knots),
                      #k0=as.double(kappa)
                      axT = as.double(kappa),
                      as.double(tt0), 
                      as.double(tt1),
                      
                      as.integer(cens), 
                      as.integer(subcluster), 
                      as.integer(clusterdc), 
                      as.integer(clusterfamdc),				
                      as.double(tt0.death), 
                      as.double(tt1.death),
                      as.integer(terminalEvent),
                      as.integer(tempdc),
                      as.integer(icdc00),
                      as.integer(nvarRec),
                      as.double(var),
                      
                      as.integer(nvarEnd),
                      as.double(vardc),
                      as.double(vaxdc00),
                      as.integer(c(noVar1, noVar2)),
                      as.double(weights.agg),
                      #as.integer(noVar2),
                      as.integer(maxit),				
                      np = as.integer(np),
                      b = as.double(Beta),
                      H_hessOut = as.double(matrix(0, nrow = np, ncol = np)),
                      HIHOut = as.double(matrix(0, nrow=np, ncol=np)),
                      
                      loglik = as.double(0), 
                      LCV = as.double(rep(0,2)),
                      xR = as.double(matrix(0, nrow=size1, ncol=uni.strat)), 	
                      lamR = as.double(array(0, dim=c(size1, 3, uni.strat))), 
                      xSuR = as.double(xSu1),
                      survR = as.double(array(0, dim=c(mt11, 3, uni.strat))),
                      xD = as.double(rep(0, size2)),
                      lamD = as.double(matrix(0, nrow=size2, ncol=3)),
                      xSuD = as.double(xSu2),
                      survD = as.double(matrix(0, nrow=mt12, ncol=3)),
                      
                      as.integer(typeof,equidistant),
                      as.integer(c(nbintervR,nbintervDC)),
                      #as.integer(nbintervDC),
                      as.integer(c(size1,size2, mt11, mt12)),
                      counts = as.integer(c(0,0,0,0)),				
                      IerIstop = as.integer(c(0,0)),
                      #ier = as.integer(0),
                      #istop = as.integer(0),
                      paraweib = as.double(rep(0,4)),
                      MartinGales = as.double(matrix(0,nrow=length(uni.subcluster), ncol=5)),
                      
                      linear.pred = as.double(rep(0,n)),
                      lineardc.pred = as.double(rep(0,as.integer(length(uni.subcluster)))),
                      zi = as.double(rep(0,(n.knots+6))),
                      time = as.double(rep(0,(nbintervR+1))),
                      timedc = as.double(rep(0,(nbintervDC+1))),
                      linearpredG = as.double(rep(0, lignedc0+1)),	
                      typeJoint0 = as.integer(joint.clust),
                      as.integer(intcens),
                      as.integer(indices),
                      #as.integer(indic_alpha),
                      #as.integer(indic_xi),
                      as.double(ttU),
                      as.integer(ordretmp),
                      as.integer(initialize),
                      
                      logNormal0 = as.integer(logNormal),		
                      paratps = as.integer(c(timedep, betaknots, betaorder)),
                      as.integer(c(filtretps, filtretps2)),
                      BetaTpsMat = as.double(matrix(0,nrow=101, ncol=1+4*nvartimedep)),
                      BetaTpsMatDc = as.double(matrix(0,nrow=101, ncol=1+4*nvartimedep2)),
                      EPS = as.double(c(LIMparam, LIMlogl, LIMderiv))
      )#,
      #PACKAGE = "frailtypack") #65 arguments
      
      ###########################################
      ### Verification de l'execution du code ###
      ###########################################
      ier <- ans$IerIstop[1]
      istop <- ans$IerIstop[2]
      if (istop == 4) warning("Problem in the likelihood computation. The program stopped abnormally. Please verify your dataset.")
      if (istop == 2) warning("Model did not converge.")
      if (istop == 3) warning("Matrix non-positive definite.")
      
      MartinGale <- matrix(ans$MartinGales, nrow=as.integer(length(uni.subcluster)), ncol=5)
      
      fit <- NULL
      
      fit$b <- ans$b
      fit$na.action <- attr(m,"na.action")
      fit$call <- call
      fit$n <- n
      fit$groups <- length(uni.cluster)
      fit$subgroups <- length(uni.subcluster)
      fit$n.events <- ans$counts[2]
      fit$n.death <- ans$counts[3]
      fit$n.censored <- ans$counts[4]
      fit$npar <- np
      fit$nst <- nst
      fit$LCV <- ans$LCV[1]
      fit$AIC <- ans$LCV[2]
      
      fit$nvar<-c(nvarRec,nvarEnd)
      #fit$nvarnotdep<-c(nvarRec-nvartimedep,nvarEnd-nvartimedep2)
      #	fit$formula <- formula(Terms)
      
      fit$xR <- matrix(ans$xR, nrow = size1, ncol = uni.strat)
      fit$lamR <- array(ans$lamR, dim = c(size1,3,uni.strat))
      fit$xSuR <- matrix(ans$xSuR, nrow = 100, ncol = uni.strat)
      fit$survR <- array(ans$survR, dim = c(mt11,3,uni.strat))
      
      fit$xD <- ans$xD
      fit$lamD <- matrix(ans$lamD, nrow = size2, ncol = 3)
      fit$xSuD <- ans$xSuD
      fit$survD <- matrix(ans$survD, nrow = mt12, ncol = 3)      
      
      fit$type <- type
      fit$n.strat <- uni.strat
      fit$n.iter <- ans$counts[1]
      fit$typeof <- typeof
      
      medianR <- NULL
      for (i in (1:fit$n.strat)) medianR[i] <- ifelse(typeof==0, minmin(fit$survR[,1,i],fit$xR), minmin(fit$survR[,1,i],fit$xSuR))
      lowerR <- NULL
      for (i in (1:fit$n.strat)) lowerR[i] <- ifelse(typeof==0, minmin(fit$survR[,2,i],fit$xR), minmin(fit$survR[,2,i],fit$xSuR))
      upperR <- NULL
      for (i in (1:fit$n.strat)) upperR[i] <- ifelse(typeof==0, minmin(fit$survR[,3,i],fit$xR), minmin(fit$survR[,3,i],fit$xSuR))
      fit$medianR <- cbind(lowerR,medianR,upperR)
      
      medianD <- ifelse(typeof==0, minmin(fit$survD[,1],fit$xD), minmin(fit$survD[,1],fit$xSuD))
      lowerD <- ifelse(typeof==0, minmin(fit$survD[,3],fit$xD), minmin(fit$survD[,3],fit$xSuD))
      upperD <- ifelse(typeof==0, minmin(fit$survD[,2],fit$xD), minmin(fit$survD[,2],fit$xSuD))
      fit$medianD <- cbind(lowerD,medianD,upperD)
      
      fit$noVar1 <- noVar1
      fit$noVar2 <- noVar2
      fit$nbintervR <- nbintervR
      fit$nbintervDC <- nbintervDC
      fit$nvarRec <- nvarRec
      fit$nvarEnd <- nvarEnd
      fit$istop <- istop
      fit$indic.nb.intR <- indic.nb.int1
      fit$indic.nb.intD <- indic.nb.int2
      fit$shape.weib <- ans$paraweib[1:2]#ans$shape.weib
      fit$scale.weib <- ans$paraweib[3:4]#ans$scale.weib
      
      fit$joint.clust <- ans$typeJoint0
      fit$AG <- recurrentAG
      fit$intcens <- intcens 
      
      fit$indic_alpha <- indic_alpha
      fit$indic_ksi <- indic_xi
      fit$logNormal <- ans$logNormal0
      #fit$BetaTpsMat <- matrix(ans$BetaTpsMat,nrow=101,ncol=1+4*nvartimedep) !! Pas encore MEP pour time-varying covariables
      #fit$BetaTpsMatDc <- matrix(ans$BetaTpsMatDc,nrow=101,ncol=1+4*nvartimedep2)
      #fit$nvartimedep <- c(nvartimedep,nvartimedep2)
      
      fit$Names.vardep <- vardep
      fit$Names.vardepdc <- vardep2      
      fit$EPS <- ans$EPS
      
      if(typeof == 0) fit$logLikPenal <- ans$loglik
      else fit$logLik <- ans$loglik 	
      
      fit$theta <- ans$b[np-nvar-indic_xi-indic_alpha-1]^2	
      fit$eta <- ans$b[np-nvar-indic_alpha-indic_xi]^2
      
      if (indic_alpha == 1) fit$alpha <- ans$b[np-nvar-indic_xi]		
      if (indic_xi == 1) fit$ksi <- ans$b[np-nvar]
      
      
      if (noVar1==1 & noVar2==1) fit$coef <- NULL
      else {
        fit$coef <- ans$b[(np-nvar+1):np]			
        noms <- c(factor.names(colnames(X)),factor.names(colnames(X2)))			
        if(length(grep(":", noms))>0) noms <- factor.names(noms)			
        #if timedep == 1			
        names(fit$coef) <- noms
      }	
      
      temp1 <- matrix(ans$H_hessOut, nrow = np, ncol = np)
      temp2 <- matrix(ans$HIHOut, nrow = np, ncol = np)	
      
      fit$varHtotal <- temp1
      fit$varHIHtotal <- temp2
      
      fit$varH <- temp1[(np-nvar-indic_xi-indic_alpha-1):np, (np-nvar-indic_xi-indic_alpha-1):np]
      fit$varHIH <- temp2[(np-nvar-indic_xi-indic_alpha-1):np, (np-nvar-indic_xi-indic_alpha-1):np]
      
      seH.theta <- sqrt(((2 * (fit$theta^0.5))^2) * diag(fit$varH)[1])
      fit$theta_p.value <- 1 - pnorm(fit$theta/seH.theta)
      
      seH.eta <- sqrt(((2 * (fit$eta^0.5))^2) * diag(fit$varH)[2])
      fit$eta_p.value <- 1 - pnorm(fit$eta/seH.eta)
      
      if (indic_alpha == 1) fit$alpha_p.value <- 1 - pchisq((fit$alpha/sqrt(diag(fit$varH))[3])^2,1)
      if (indic_xi == 1)fit$ksi_p.value <- 1 - pchisq((fit$ksi/sqrt(diag(fit$varH))[3+indic_alpha])^2,1)
      
      
      
      if (indic_alpha == 1) noms <- c("theta","alpha",factor.names(colnames(X)),factor.names(colnames(X2)))
      else noms <- c("theta",factor.names(colnames(X)),factor.names(colnames(X2)))
      if(length(grep(":",noms))>0)noms <- factor.names(noms)		
      
      if (typeof == 0){
        fit$n.knots<-n.knots
        fit$kappa <- ans$axT
        fit$cross.Val<-cross.validation
        fit$n.knots.temp <- n.knots.temp
        fit$zi <- ans$zi
      }
      # if constante par morceaux      
      # verif que les martingales ont ete bien calculees
      #msg <- "Problem in the estimation of the random effects (perhaps high number of events in some clusters)"
      #if (Frailty){
      #	if (any(MartinGale[,1]==0)){
      #		fit$martingale.res <- msg
      #		fit$martingaledeath.res <- msg          
      #		fit$frailty.pred <- msg   
      #		fit$linear.pred <- msg
      #		fit$lineardeath.pred <- msg
      #	}else{
      #		fit$martingale.res <- MartinGale[,1]
      #		fit$martingaledeath.res <- MartinGale[,2]          
      #		fit$frailty.pred <- MartinGale[,3]          
      #		fit$linear.pred <- ans$linearpred[order(ordre)]
      #		fit$lineardeath.pred <- ans$linearpreddc
      #	}
      #}    
      fit$martingale.res <- MartinGale[,1] # seulement pour les familles
      fit$martingaledeath.res <- MartinGale[,2]		
      fit$frailty.pred <- MartinGale[,3]   
      fit$frailty.fam.pred <- MartinGale[1:length(uni.cluster),5] ### AK 12/12/2016 (family frailties)
      fit$linear.pred <- ans$linear.pred[order(ordre)]# a faire
      fit$lineardeath.pred <- ans$lineardc.pred
      
      #================================> For the Recurrent Event
      #========================= Test de Wald
      
      if (length(vec.factor) > 0){
        Beta <- ans$b[(np - nvar + 1):np]
        
        if (indic_alpha == 1) {
          #if (indic_xi ==1) VarBeta <- diag(diag(fit$varH)[-c(1,2,3,4)])
          #else VarBeta <- diag(diag(fit$varH)[-c(1,2,3)])
          if(indic_xi == 1) VarBeta <- fit$varH[5:dim(fit$varH)[1],5:dim(fit$varH)[2]]
          else VarBeta <- fit$varH[4:dim(fit$varH)[1],4:dim(fit$varH)[2]]
        }else{
          #if (indic_xi == 1) VarBeta <- diag(diag(fit$varH)[-c(1,2,3)])
          #else VarBeta <- diag(diag(fit$varH)[-c(1,2)])
          if(indic_xi == 1) VarBeta <- fit$varH[4:dim(fit$varH)[1],4:dim(fit$varH)[2]]
          else VarBeta <- fit$varH[3:dim(fit$varH)[1],3:dim(fit$varH)[2]]
        }	
        
        nfactor <- length(vec.factor)
        p.wald <- rep(0,nfactor)
        ntot <- nvarEnd + nvarRec
        
        if(fit$istop == 1) fit$global_chisq <- waldtest(N=nvarRec,nfact=nfactor,place=ind.place,modality=occur,b=Beta,Varb=VarBeta,Llast=nvarEnd,Ntot=ntot)		 
        else fit$global_chisq <- 0
        
        fit$dof_chisq <- occur
        fit$global_chisq.test <- 1
        # Calcul de pvalue globale
        for(i in 1:length(vec.factor)){
          p.wald[i] <- signif(1 - pchisq(fit$global_chisq[i], occur[i]), 3)
        }
        fit$p.global_chisq <- p.wald
        fit$names.factor <- vec.factor
      }else{
        fit$global_chisq.test <- 0
      }
      
      #================================> For the death
      #========================= Test de Wald      
      if (length(vec.factordc) > 0){
        Beta <- ans$b[(np-nvar+1):np]
        
        if (indic_alpha == 1) {
          if(indic_xi ==1) VarBeta <- diag(diag(fit$varH)[-c(1,2,3,4)])
          else VarBeta <- diag(diag(fit$varH)[-c(1,2,3)])
        }else{
          if(indic_xi == 1) VarBeta <- diag(diag(fit$varH)[-c(1,2,3)])
          else VarBeta <- diag(diag(fit$varH)[-c(1,2)])
        }			
        nfactor <- length(vec.factordc)
        p.walddc <- rep(0,nfactor)
        ntot <- nvarEnd + nvarRec
        
        if(fit$istop == 1) fit$global_chisq_d <- waldtest(N=nvarEnd,nfact=nfactor,place=ind.placedc,modality=occurdc,b=Beta,Varb=VarBeta,Lfirts=nvarRec,Ntot=ntot)
        else fit$global_chisq_d <- 0 
        
        fit$dof_chisq_d <- occurdc
        fit$global_chisq.test_d <- 1
        # Calcul de pvalue globale
        for(i in 1:length(vec.factordc)){
          p.walddc[i] <- signif(1 - pchisq(fit$global_chisq_d[i], occurdc[i]), 3)
        }
        fit$p.global_chisq_d <- p.walddc
        fit$names.factordc <- vec.factordc
      }else fit$global_chisq.test_d <- 0
      
      
      if(!is.null(fit$coef)){
        
        if (indic_alpha == 1 & indic_xi == 1) seH <- sqrt(diag(fit$varH))[-c(1:4)]
        else if (indic_alpha == 0 | indic_xi == 0) seH <- sqrt(diag(fit$varH))[-c(1:3)]
        else if(indic_alpha == 0 & indic_xi == 0) seH <- sqrt(diag(fit$varH))[-c(1:2)]
        
        
        fit$beta_p.value <- 1 - pchisq((fit$coef/seH)^2, 1)
      }
      
      if (length(Xlevels) >0)	fit$Xlevels <- Xlevels
      fit$contrasts <- contr.save
      if (length(Xlevels2) >0) fit$Xlevels2 <- Xlevels2
      fit$contrasts2 <- contr.save2
      
      fit$formula <- formula
      fit$formula.terminalEvent <- formula.terminalEvent
      
      attr(fit,"joint")<-joint
      attr(fit,"subcluster")<-TRUE
      class(fit) <- "jointNestedPenal"
      
      # ===================================================
      # END NESTED JOINT MODEL 
      # ===================================================		
    }
    
    if (print.times){
      cost<-proc.time()-ptm
      cat("The program took", round(cost[3],2), "seconds \n")
    }
    fit	
  }
