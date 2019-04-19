#' Prediction probabilities for Cox proportional hazard, Shared, Joint frailty
#' models, Joint models for longitudinal data and a terminal event and
#' Trivariate joint model for longitudinal data, recurrent events and a
#' terminal event (linear and non-linear).
#' 
#' @description{ 
#' 
#' \ifelse{html}{\bold{For Cox proportional hazard model}
#' 
#' A predictive probability of event between t and horizon time t+w, with w the
#' window of prediction. 
#' 
#' {\figure{prediction1.png}{options: width="100\%"}}
#' 
#' \bold{For Gamma Shared Frailty model for clustered (not recurrent) events}
#' 
#' Two kinds of predictive probabilities can be calculated:
#' 
#' - a conditional predictive probability of event between t and horizon time
#' t+w, i.e. given a specific group
#' 
#' {\figure{prediction2.png}{options: width="100\%"}}
#' 
#' - a marginal predictive probability of event between t and horizon time t+w,
#' i.e. averaged over the population
#' 
#' {\figure{prediction3.png}{options: width="100\%"}}
#' 
#' \bold{For Gaussian Shared Frailty model for clustered (not recurrent)
#' events}
#' 
#' Two kinds of predictive probabilities can be calculated:
#' 
#' - a conditional predictive probability of event between t and horizon time
#' t+w, i.e. given a specific group and given a specific Gaussian random effect
#' \eqn{\eta}
#' 
#' {\figure{prediction4.png}{options: width="100\%"}}
#' 
#' - a marginal predictive probability of event between t and horizon time t+w,
#' i.e. averaged over the population
#' 
#' {\figure{prediction5.png}{options: width="100\%"}}
#' 
#' \bold{For Gamma Shared Frailty model for recurrent events}
#' 
#' Two kinds of predictive probabilities can be calculated:
#' 
#' - A marginal predictive probability of event between t and horizon time t+w,
#' i.e. averaged over the population.
#' 
#' {\figure{prediction6.png}{options: width="100\%"}}
#' 
#' - a conditional predictive probability of event between t and horizon time
#' t+w, i.e. given a specific individual.
#' 
#' This prediction method is the same as the conditional gamma prediction
#' method applied for clustered events (see formula \eqn{P}\out{<sup>cond</sup>} before).
#' 
#' \bold{For Gaussian Shared Frailty model for recurrent events}
#' 
#' Two kinds of predictive probabilities can be calculated:
#' 
#' - A marginal predictive probability of event between t and horizon time t+w,
#' i.e. averaged over the population.
#' 
#' {\figure{prediction7.png}{options: width="100\%"}}
#' 
#' - a conditional predictive probability of event between t and horizon time
#' t+w, i.e. given a specific individual.
#' 
#' This prediction method is the same as the conditional Gaussian prediction
#' method applied for clustered events (see formula \eqn{P}\out{<sup>cond</sup>} before).
#' 
#' It is possible to compute all these predictions in two ways on a scale of
#' times : - either you want a cumulative probability of developing the event
#' between t and t+w (with t fixed, but with a varying window of prediction w);
#' - either you want at a specific time the probability to develop the event in
#' the next w (ie, for a varying prediction time t, but for a fixed window of
#' prediction). See Details.
#' 
#' \bold{For Joint Frailty model}
#' 
#' Prediction for two types of event can be calculated : for a terminal event
#' or for a new recurrent event, knowing patient's characteristics.
#' 
#' \bold{ - Prediction of death knowing patients' characteristics : }
#' 
#' It is to predict the probability of death in a specific time window given
#' the history of patient i before the time of prediction t. The history
#' H\out{<sub>i</sub>}\out{<sup>J,l</sup>}, (l=1,2) is the information on covariates before time
#' t, but also the number of recurrences and the time of occurences. Three
#' types of marginal probabilities are computed:
#' 
#' - a prediction of death between t and t+w given that the patient had exactly
#' J recurrences (H\out{<sub>i</sub>}\out{<sup>J,1</sup>}) before t
#' 
#' {\figure{prediction8.png}{options: width="100\%"}}
#' 
#' - a prediction of death between t and t+w given that the patient had at
#' least J recurrences (H\out{<sub>i</sub>}\out{<sup>J,2</sup>}) before t
#' 
#' {\figure{prediction9.png}{options: width="100\%"}}
#' 
#' - a prediction of death between t and t+w considering the recurrence history
#' only in the parameters estimation. It corresponds to the average probability
#' of death between t and t+w for a patient with these given characteristics.
#' 
#' {\figure{prediction10.png}{options: width="100\%"}}
#' 
#' \bold{ - Prediction of risk of a new recurrent event knowing patients'
#' characteristics : }
#' 
#' It is to predict the probability of a new recurrent event in a specific time
#' window given the history of patient i before the time of prediction t. The
#' history H\out{<sub>i</sub>}\out{<sup>J</sup>} is the information on covariates before time t, but also
#' the number of recurrences and the time of occurences. The marginal
#' probability computed is a prediction of a new recurrent event between t and
#' t+w given that the patient had exactly J recurrences (\eqn{H}\out{<sub>i</sub>}\out{<sup>J</sup>}) before t:
#' 
#' {\figure{prediction11.1.png}{options: width="100\%"}}
#' {\figure{prediction11.2.png}{options: width="100\%"}}
#' 
#' It is possible to compute all these predictions in two ways : - either you
#' want a cumulative probability of developing the event between t and t+w
#' (with t fixed, but with a varying window of prediction w); - either you want
#' at a specific time the probability to develop the event in the next w (ie,
#' for a varying prediction time t, but for a fixed window of prediction). See
#' Details.
#' 
#' With Gaussian frailties (\eqn{\eta}), the same expressions are used but with
#' u\out{<sub>i</sub>}\out{<sup>J</sup>} replaced by \eqn{exp}(J\\eqn{eta}\out{<sub>i</sub>}) and \eqn{g(\eta)}
#' corresponds to the Gaussian distribution.
#' 
#' \bold{For Joint Nested Frailty models}
#' 
#' Prediction of the probability of developing a terminal event between t and
#' t+w for subject i who survived by time t based on the visiting and disease
#' histories of their own and other family members observed by time t.
#' 
#' Let (Y\out{<sub>fi</sub>}\out{<sup>R</sup>}(t)) be the history of subject i in family f,
#' before time t, which includes all the recurrent events and covariate
#' information. For disease history, let T\out{<sub>fi</sub>}\out{<sup>D</sup>}(t) = min(T\out{<sub>fi</sub>},t) be
#' the observed time to an event before t ; \eqn{\delta}\out{<sub>fi</sub>}\out{<sup>D</sup>}(t) the disease
#' indicator by time t and X\out{<sub>fi</sub>}\out{<sup>D</sup>}(t) the covariate information
#' observed up to time t. We define the family history of subject i in
#' family f by {\figure{prediction12.png}{options: width="100\%"}}
#' 
#' which includes the visiting and disease history of all subjects except for
#' subject i in family f as well as their covariate information by
#' time t.
#' 
#' The prediction probability can be written as :
#' 
#' {\figure{prediction13.1.png}{options: width="100\%"}}
#' {\figure{prediction13.2.png}{options: width="100\%"}}
#' 
#' \bold{For Joint models for longitudinal data and a terminal event}
#' 
#' The predicted probabilities are calculated in a specific time window given
#' the history of biomarker measurements before the time of prediction t
#' (\out{&#x1B4;}\out{<sub>i</sub>}(t)). The probabilities are conditional also on
#' covariates before time t and that the subject was at risk at t.  The
#' marginal predicted probability of the terminal event is
#' 
#' {\figure{prediction14.png}{options: width="100\%"}}
#' 
#' These probabilities can be calculated in several time points with fixed time
#' of prediction t and varying window w or with fixed window w and varying time
#' of prediction t. See Details for an example of how to construct time
#' windows.
#' 
#' \bold{For Trivariate joint models for longitudinal data, recurrent events
#' and a terminal event}
#' 
#' The predicted probabilities are calculated in a specific time window given
#' the history of biomarker measurements \out{&#x1B4;}\out{<sub>i</sub>}(t) and recurrences
#' H\out{<sub>i</sub>}\out{<sup>J,1</sup>} (complete history of recurrences with known J number
#' of observed events) before the time of prediction t. The probabilities are
#' conditional also on covariates before time t and that the subject was at
#' risk at t.  The marginal predicted probability of the terminal event is
#' 
#' {\figure{prediction15.png}{options: width="100\%"}}
#' 
#' The biomarker history can be represented using a linear (\code{trivPenal})
#' or non-linear mixed-effects model (\code{trivPenalNL}).
#' 
#' These probabilities can be calculated in several time points with fixed time
#' of prediction t and varying window w or with fixed window w and varying time
#' of prediction t. See Details for an example of how to construct time
#' windows.
#' 
#' }{\bold{For Cox proportional hazard model}
#' 
#' A predictive probability of event between t and horizon time t+w, with w the
#' window of prediction.
#' 
#' \deqn{P(t,t+w)=\frac{S_i(t)-S_i(t+w)}{S_i(t)}=1-\left(\frac{S_0(t+w)}{S_0(t)}\right)^{\exp(\beta'Z_i)}}
#' 
#' \bold{For Gamma Shared Frailty model for clustered (not recurrent) events}
#' 
#' Two kinds of predictive probabilities can be calculated:
#' 
#' - a conditional predictive probability of event between t and horizon time
#' t+w, i.e. given a specific group
#' 
#' \deqn{P^{cond}(t,t+w)=\frac{S_{ij}(t|u_i)-S_{ij}(t+w|u_i)}{S_{ij}(t|u_i)}=1-\left(\frac{S_0(t+w)}{S_0(t)}\right)^{u_i\exp(\beta'Z_{ij})}}
#' 
#' - a marginal predictive probability of event between t and horizon time t+w,
#' i.e. averaged over the population
#' 
#' \deqn{P^{marg}(t,t+w)=1-\left(\frac{1+\theta
#' H_0(t)\exp(\beta'Z_{ij})}{1+\theta H_0(t+w)\exp(\beta'Z_{ij})}\right)^{1/
#' \theta}}
#' 
#' \bold{For Gaussian Shared Frailty model for clustered (not recurrent)
#' events}
#' 
#' Two kinds of predictive probabilities can be calculated:
#' 
#' - a conditional predictive probability of event between t and horizon time
#' t+w, i.e. given a specific group and given a specific Gaussian random effect
#' \eqn{\eta}
#' 
#' \deqn{P^{cond}(t,t+w)=\frac{S_{ij}(t|\eta_i)-S_{ij}(t+w|\eta_i)}{S_{ij}(t|\eta_i)}=1-\left(\frac{S_0(t+w)}{S_0(t)}\right)^{\exp(\eta_i+\beta'Z_{ij})}}
#' 
#' - a marginal predictive probability of event between t and horizon time t+w,
#' i.e. averaged over the population
#' 
#' \deqn{P^{marg}(t,t+w)=\frac{\int_{-\infty}^{+\infty}(S_{ij}(t|\eta_i)-S_{ij}(t+w|\eta_i))g(\eta)d\eta}{\int_{-\infty}^{+\infty}S_{ij}(t)g(\eta)d\eta}}
#' 
#' \bold{For Gamma Shared Frailty model for recurrent events}
#' 
#' Two kinds of predictive probabilities can be calculated:
#' 
#' - A marginal predictive probability of event between t and horizon time t+w,
#' i.e. averaged over the population.
#' 
#' \deqn{P^{marg}(t,t+w)=\frac{\int_{0}^{+\infty}(S_{i(J+1)}(t|u_i) -
#' S_{ij}(t+w|u_i))\cdot(u_i)^J S_{ij}(X_{iJ}|u_i)
#' g(u)du}{\int_{0}^{+\infty}S_{i(J+1)}(t|u_i) (u_i)^J
#' S_{i(J+1)}(X_{iJ}|u_i))g(u)du}}
#' 
#' - a conditional predictive probability of event between t and horizon time
#' t+w, i.e. given a specific individual.
#' 
#' This prediction method is the same as the conditional gamma prediction
#' method applied for clustered events (see formula \deqn{P^{cond}} before).
#' 
#' \bold{For Gaussian Shared Frailty model for recurrent events}
#' 
#' Two kinds of predictive probabilities can be calculated:
#' 
#' - A marginal predictive probability of event between t and horizon time t+w,
#' i.e. averaged over the population.
#' 
#' \deqn{P^{marg}(t,t+w)=\frac{\int_{0}^{+\infty}(S_{i(J+1)}(t|\eta_i) -
#' S_{ij}(t+w|\eta_i))\cdot \exp(J\eta_i) S_{ij}(X_{iJ}|\eta_i)
#' g(\eta)d\eta}{\int_{0}^{+\infty}S_{i(J+1)}(t|\eta_i) \exp(J\eta_i)
#' S_{i(J+1)}(X_{iJ}|\eta_i))g(\eta)d\eta}}
#' 
#' - a conditional predictive probability of event between t and horizon time
#' t+w, i.e. given a specific individual.
#' 
#' This prediction method is the same as the conditional Gaussian prediction
#' method applied for clustered events (see formula \deqn{P^{cond}} before).
#' 
#' It is possible to compute all these predictions in two ways on a scale of
#' times : - either you want a cumulative probability of developing the event
#' between t and t+w (with t fixed, but with a varying window of prediction w);
#' - either you want at a specific time the probability to develop the event in
#' the next w (ie, for a varying prediction time t, but for a fixed window of
#' prediction). See Details.
#' 
#' \bold{For Joint Frailty model}
#' 
#' Prediction for two types of event can be calculated : for a terminal event
#' or for a new recurrent event, knowing patient's characteristics.
#' 
#' \bold{ - Prediction of death knowing patients' characteristics : }
#' 
#' It is to predict the probability of death in a specific time window given
#' the history of patient i before the time of prediction t. The history
#' \eqn{H_i^{J,l}}, (\eqn{l=1,2}) is the information on covariates before time
#' t, but also the number of recurrences and the time of occurences. Three
#' types of marginal probabilities are computed:
#' 
#' - a prediction of death between t and t+w given that the patient had exactly
#' J recurrences (\eqn{H_i^{J,1}}) before t
#' 
#' \deqn{P^1(t,t+w)=P(D_i \le
#' t+w|D_i>t,H_i^{J,1})=\frac{\int_0^\infty[S_i^D(t)-S_i^D(t+w)](u_i)^JS_{i(J+1)}^R(t)g(u)du_i}{\int_0^\infty
#' S_i^D(t)(u_i)^JS_{i(J+1)}^R(t)g(u)du_i}}
#' 
#' - a prediction of death between t and t+w given that the patient had at
#' least J recurrences (\eqn{H_i^{J,2}}) before t
#' 
#' \deqn{P^2(t,t+w)=P(D_i \le
#' t+w|D_i>t,H_i^{J,2})=\frac{\int_0^\infty[S_i^D(t)-S_i^D(t+w)](u_i)^JS_{iJ}^R(X_{iJ})g(u)du_i}{\int_0^\infty
#' S_i^D(t)(u_i)^JS_{iJ}^R(X_{iJ})g(u)du_i}}
#' 
#' - a prediction of death between t and t+w considering the recurrence history
#' only in the parameters estimation. It corresponds to the average probability
#' of death between t and t+w for a patient with these given characteristics.
#' 
#' \deqn{P^3(t,t+w)=P(D_i \le
#' t+w|D_i>t)=\frac{\int_0^\infty[S_i^D(t)-S_i^D(t+w)]g(u)du_i}{\int_0^\infty
#' S_i^D(t)g(u)du_i}}
#' 
#' \bold{ - Prediction of risk of a new recurrent event knowing patients'
#' characteristics : }
#' 
#' It is to predict the probability of a new recurrent event in a specific time
#' window given the history of patient i before the time of prediction t. The
#' history \eqn{H_i^J} is the information on covariates before time t, but also
#' the number of recurrences and the time of occurences.  The marginal
#' probability computed is a prediction of a new recurrent event between t and
#' t+w given that the patient had exactly J recurrences (\eqn{H_i^J}) before t:
#' 
#' \deqn{P^R(t,t+w)=P(X_{i(j+1)} \le t+w|X_{i(j+1)}>t, D_i>t, H_i^J)=}
#' 
#' \deqn{\frac{\int_0^\infty[S_{i(J+1)}^R(t)-S_{i(J+1)}^R(t+w)]S_i^D(t)(u_i)^JS_{i(J+1)}^R(X_{ij})g(u)du_i}{\int_0^\infty
#' S_{i(J+1)}^R(t)S_i^D(t)(u_i)^JS_{i(J+1)}^R(X_{ij})g(u)du_i}}
#' 
#' It is possible to compute all these predictions in two ways : - either you
#' want a cumulative probability of developing the event between t and t+w
#' (with t fixed, but with a varying window of prediction w); - either you want
#' at a specific time the probability to develop the event in the next w (ie,
#' for a varying prediction time t, but for a fixed window of prediction). See
#' Details.
#' 
#' With Gaussian frailties (\eqn{\eta}), the same expressions are used but with
#' \eqn{{u_i}^{J}} replaced by \eqn{\exp(J\eta_i)} and \eqn{g(\eta)}
#' corresponds to the Gaussian distribution.
#' 
#' \bold{For Joint Nested Frailty models}
#' 
#' Prediction of the probability of developing a terminal event between t and
#' t+w for subject i who survived by time t based on the visiting and disease
#' histories of their own and other family members observed by time t.
#' 
#' Let (\eqn{Y_{fi}^R(t)}) be the history of subject \eqn{i} in family \eqn{f},
#' before time \eqn{t}, which includes all the recurrent events and covariate
#' information. For disease history, let \eqn{T_{fi}^D(t) = min(T_{fi},t)} be
#' the observed time to an event before t ; \eqn{\delta_{fi}^D(t)} the disease
#' indicator by time \eqn{t} and \eqn{X_{fi}^D(t)} the covariate information
#' observed up to time t. We define the family history of subject \eqn{i} in
#' family \eqn{f} by \deqn{H_{f(-i)}(t) = \{Y_{fl}^R(t), T_{fl}^D(t),
#' \delta_{fl}^D(t), X_{fl}^D(t), \forall l \in \{1, ..., i-1, i+1, ..., m_f
#' \}\} }
#' 
#' which includes the visiting and disease history of all subjects except for
#' subject \eqn{i} in family \eqn{f} as well as their covariate information by
#' time \eqn{t}.
#' 
#' The prediction probability can be written as :
#' 
#' \deqn{ P(T_{fi}^D < t+s|T_{fi}^D > t, Y_i(t), H_{(f-i)}(t)) = }
#' 
#' \deqn{ \frac{ \int \int P(t<T_{fi}^D < t+s| X_{fi}^D, \omega_{fi})
#' P(Y_i(t)|X_{fi}^R(t), \omega_{i}) P(H_{f(-i)}(t)|X_{f(-i)}(t), \omega_{fi})
#' g_{ui} g_{\omega f} }{\int \int P(T_{fi}^D > t| X_{fi}^D, \omega_{fi})
#' P(Y_i(t)|X_{fi}^R(t), \omega_{i}) P(H_{f(-i)}(t)|X_{f(-i)}(t), \omega_{fi})
#' g_{ui} g_{\omega f} } }
#' 
#' \bold{For Joint models for longitudinal data and a terminal event}
#' 
#' The predicted probabilities are calculated in a specific time window given
#' the history of biomarker measurements before the time of prediction t
#' (\eqn{\mathcal{Y}_i(t)}). The probabilities are conditional also on
#' covariates before time t and that the subject was at risk at t.  The
#' marginal predicted probability of the terminal event is
#' 
#' \deqn{P(t,t+w)=P(D_i \le
#' t+w|D_i>t,\mathcal{Y}_i(t))=\frac{\int_0^\infty[S_i^D(t)-S_i^D(t+w)]f(\mathcal{Y}_i(t)|\bold{X}_{Li},\bold{b}_i)f(\bold{b}_i)d\bold{b}_i}{\int_0^\infty
#' S_i^D(t)f(\mathcal{Y}_i(t)|\bold{X}_{Li},\bold{b}_i)f(\bold{b}_i)d\bold{b}_i}}
#' 
#' These probabilities can be calculated in several time points with fixed time
#' of prediction t and varying window w or with fixed window w and varying time
#' of prediction t. See Details for an example of how to construct time
#' windows.
#' 
#' \bold{For Trivariate joint models for longitudinal data, recurrent events
#' and a terminal event}
#' 
#' The predicted probabilities are calculated in a specific time window given
#' the history of biomarker measurements \eqn{\mathcal{Y}_i(t)} and recurrences
#' \eqn{H_i^{J,1}} (complete history of recurrences with known \eqn{J} number
#' of observed events) before the time of prediction t. The probabilities are
#' conditional also on covariates before time t and that the subject was at
#' risk at t.  The marginal predicted probability of the terminal event is
#' 
#' \deqn{ \begin{array}{ll} P(t,t+w)&=P(D_i \le
#' t+w|D_i>t,H_i^{J,1},\mathcal{Y}_i(t))\\
#' &=\frac{\int_0^\infty[S_i^D(t)-S_i^D(t+w)])]\exp(J(v_i+g(t)^\top\bold{\eta}_R))S_{i(J+1)}^R(t)f(\mathcal{Y}_i(t)|\bold{X}_{Li},\bold{b}_i)
#' f(\bold{u}_i)d\bold{u}_i}{\int_0^\infty
#' S_i^D(t)\exp(J(v_i+g(t)^\top\bold{\eta}_R))S_{i(J+1)}^R(t)f(\mathcal{Y}_i(t)|\bold{X}_{Li},\bold{b}_i)f(\bold{u}_i)d\bold{u}_i}
#' \end{array}}
#' 
#' The biomarker history can be represented using a linear (\code{trivPenal})
#' or non-linear mixed-effects model (\code{trivPenalNL}).
#' 
#' These probabilities can be calculated in several time points with fixed time
#' of prediction t and varying window w or with fixed window w and varying time
#' of prediction t. See Details for an example of how to construct time
#' windows.
#' }
#' }
#' 
#' @details{
#' To compute predictions with a prediction time t fixed and a variable window:
#' \preformatted{prediction(fit, datapred, t=10, window=seq(1,10,by=1))}
#' Otherwise, you can have a variable prediction time and a fixed window.
#' \preformatted{prediction(fit, datapred, t=seq(10,20,by=1), window=5)} Or fix
#' both prediction time t and window. \preformatted{prediction(fit, datapred,
#' t=10, window=5)}
#' 
#' The data frame building is an important step. It will contain profiles of
#' patient on which you want to do predictions. To make predictions on a Cox
#' proportional hazard or a shared frailty model, only covariates need to be
#' included. You have to distinguish between numerical and categorical
#' variables (factors). If we fit a shared frailty model with two covariates
#' sex (factor) and age (numeric), here is the associated data frame for three
#' profiles of prediction.
#' 
#' \preformatted{ datapred <- data.frame(sex=0,age=0) datapred$sex <-
#' as.factor(datapred$sex) levels(datapred$sex)<- c(1,2) datapred[1,] <-
#' c(1,40) # man, 40 years old datapred[2,] <- c(2,45) # woman, 45 years old
#' datapred[3,] <- c(1,60) # man, 60 years old }
#' 
#' \bold{Time-dependent covariates:} In the context of time-dependent
#' covariate, the last previous value of the covariate is used before the time
#' t of prediction.
#' 
#' It should be noted, that in a data frame for both marginal and conditional
#' prediction on a shared frailty model for clustered data, the group must be
#' specified. In the case of marginal predictions this can be any number as it
#' does not influence predictions. However, for conditional predictions, the
#' group must be also included in the data set used for the model fitting. The
#' conditional predictions apply the empirical Bayes estimate of the frailty
#' from the specified cluster.  Here, three individuals belong to group 5.
#' 
#' \preformatted{ datapred <- data.frame(group=0, sex=0,age=0) datapred$sex <-
#' as.factor(datapred$sex) levels(datapred$sex)<- c(1,2) datapred[1,] <-
#' c(5,1,40) # man, 40 years old (cluster 5) datapred[2,] <- c(5,2,45) # woman,
#' 45 years old (cluster 5) datapred[3,] <- c(5,1,60) # man, 60 years old
#' (cluster 5) }
#' 
#' To use the prediction function on joint frailty models and trivariate joint
#' models, the construction will be a little bit different. In these cases, the
#' prediction for the terminal event takes into account covariates but also
#' history of recurrent event times for a patient. You have to create a data
#' frame with the relapse times, the indicator of event, the cluster variable
#' and the covariates. Relapses occurring after the prediction time may be
#' included but will be ignored for the prediction. A joint model with
#' calendar-timescale need to be fitted with Surv(start,stop,event), relapse
#' times correspond to the "stop" variable and indicators of event correspond
#' to the "event" variable (if event=0, the relapse will not be taken into
#' account). For patients without relapses, all the values of "event" variable
#' should be set to 0. Finally, the same cluster variable name needs to be in
#' the joint model and in the data frame for predictions ("id" in the following
#' example). For instance, we observe relapses of a disease and fit a joint
#' model adjusted for two covariates sex (1:male 2:female) and chemo (treatment
#' by chemotherapy 1:no 2:yes). We describe 3 different profiles of prediction
#' all treated by chemotherapy: 1) a man with four relapses at 100, 200, 300
#' and 400 days, 2) a man with only one relapse at 1000 days, 3) a woman
#' without relapse.
#' 
#' \preformatted{ datapred <- data.frame(time=0,event=0,id=0,sex=0,chemo=0)
#' datapred$sex <- as.factor(datapred$sex) levels(datapred$sex) <- c(1,2)
#' datapred$chemo <- as.factor(datapred$chemo) levels(datapred$chemo) <- c(1,2)
#' datapred[1,] <- c(100,1,1,1,2) # first relapse of the patient 1 datapred[2,]
#' <- c(200,1,1,1,2) # second relapse of the patient 1 datapred[3,] <-
#' c(300,1,1,1,2) # third relapse of the patient 1 datapred[4,] <-
#' c(400,1,1,1,2) # fourth relapse of the patient 1 datapred[5,] <-
#' c(1000,1,2,1,2) # one relapse at 1000 days for patient 2 datapred[6,] <-
#' c(100,0,3,2,2) # patient 3 did not relapse }
#' 
#' The data can also be the dataset used to fit the joint model. In this case,
#' you will obtain as many prediction rows as patients.
#' 
#' Finally, for the predictions using joint models for longitudinal data and a
#' terminal event and trivariate joint models, a data frame with the history of
#' the biomarker measurements must be provided. It must include data on
#' measurements (values and time points), cluster variable and covariates.
#' Measurements taken after the prediction time may be included but will be
#' ignored for the prediction. The same cluster variable name must be in the
#' data frame, in the data frame used for the joint model and in the data frame
#' with the recurrent event and terminal event times. For instance, we observe
#' two patients and each one had 5 tumor size measurements (patient 1 had an
#' increasing tumor size and patient 2, decreasing). The joint model used for
#' the predictions was adjusted on sex (1: male, 2: female), treatment (1:
#' sequential arm, 2: combined arm), WHO baseline performance status (1: 0
#' status, 2: 1 status, 3: 2 status) and previous resection of the primate
#' tumor (0: no, 1: yes). The data frame for the biomarker measurements can be:
#' \preformatted{ datapredj_longi <- data.frame(id = 0, year = 0, tumor.size =
#' 0, treatment = 0, age = 0, who.PS = 0, prev.resection = 0)
#' datapredj_longi$treatment <- as.factor(datapredj_longi$treatment)
#' levels(datapredj_longi$treatment) <- 1:2 datapredj_longi$age <-
#' as.factor(datapredj_longi$age) levels(datapredj_longi$age) <- 1:3
#' datapredj_longi$who.PS <- as.factor(datapredj_longi$who.PS)
#' levels(datapredj_longi$who.PS) <- 1:3 datapredj_longi$prev.resection <-
#' as.factor (datapredj_longi$prev.resection)
#' levels(datapredj_longi$prev.resection) <- 1:2 # patient 1: increasing tumor
#' size datapredj_longi[1,] <- c(1, 0,1.2 ,2,1,1,1) datapredj_longi[2,] <-
#' c(1,0.3,1.4,2,1,1,1) datapredj_longi[3,] <- c(1,0.6,1.9,2,1,1,1)
#' datapredj_longi[4,] <- c(1,0.9,2.5,2,1,1,1) datapredj_longi[5,] <-
#' c(1,1.5,3.9,2,1,1,1)
#' 
#' # patient 2: decreasing tumor size datapredj_longi[6,] <- c(2, 0,1.2
#' ,2,1,1,1) datapredj_longi[7,] <- c(2,0.3,0.7,2,1,1,1) datapredj_longi[8,] <-
#' c(2,0.5,0.3,2,1,1,1) datapredj_longi[9,] <- c(2,0.7,0.1,2,1,1,1)
#' datapredj_longi[10,] <- c(2,0.9,0.1,2,1,1,1) }
#' }
#' 
#' @usage
#' 
#' prediction(fit, data, data.Longi, t, window, event="Both", conditional =
#' FALSE, MC.sample=0, individual)
#' @param fit A frailtyPenal, jointPenal, longiPenal, trivPenal or trivPenalNL
#' object.
#' @param data Data frame for the prediction. See Details.
#' @param data.Longi Data frame for the prediction used for joint models with
#' longitudinal data. See Details.
#' @param event Only for joint and shared models. The type of event you want to
#' predict : "Terminal" for a terminal event, "Recurrent" for a recurrent event
#' or "Both". Default value is "Both". For joint nested model, only 'Terminal'
#' is allowed.  In a shared model, if you want to predict a new recurrent event
#' then the argument "Recurrent" should be use. If you want to predict a new
#' event from clustered data, do not use this option.
#' @param t Time or vector of times for prediction.
#' @param window Window or vector of windows for prediction.
#' @param conditional Only for prediction method applied on shared models.
#' Provides distinction between the conditional and marginal prediction
#' methods. Default is FALSE.
#' @param MC.sample Number of samples used to calculate confidence bands with a
#' Monte-Carlo method (with a maximum of 1000 samples). If MC.sample=0 (default
#' value), no confidence intervals are calculated.
#' @param individual Only for joint nested model. Vector of individuals (of the
#' same family) you want to make prediction.
#' @return
#' 
#' The following components are included in a 'predFrailty' object obtained by
#' using prediction function for Cox proportional hazard and shared frailty
#' model.
#' 
#' \item{npred}{Number of individual predictions} \item{x.time}{A vector of
#' prediction times of interest (used for plotting predictions): vector of
#' prediction times t if fixed window. Otherwise vector of prediction times
#' t+w} \item{window}{Prediction window or vector of prediction windows}
#' \item{pred}{Predictions estimated for each profile} \item{icproba}{Logical
#' value. Were confidence intervals estimated ?} \item{predLow}{Lower limit of
#' Monte-Carlo confidence interval for each prediction} \item{predHigh}{Upper
#' limit of Monte-Carlo confidence interval for each prediction}
#' \item{type}{Type of prediction probability (marginal or conditional)}
#' \item{group}{For conditional probability, the list of group on which you
#' make predictions}
#' 
#' The following components are included in a 'predJoint' object obtained by
#' using prediction function for joint frailty model.
#' 
#' \item{npred}{Number of individual predictions} \item{x.time}{A vector of
#' prediction times of interest (used for plotting predictions): vector of
#' prediction times t if fixed window. Otherwise vector of prediction times
#' t+w} \item{window}{Prediction window or vector of prediction windows}
#' \item{group}{Id of each patient} \item{pred1}{Estimation of probability of
#' type 1: exactly j recurrences} \item{pred2}{Estimation of probability of
#' type 2: at least j recurrences} \item{pred3}{Estimation of probability of
#' type 3} \item{pred1_rec}{Estimation of prediction of relapse}
#' \item{icproba}{Logical value. Were confidence intervals estimated ?}
#' \item{predlow1}{Lower limit of Monte-Carlo confidence interval for
#' probability of type 1} \item{predhigh1}{Upper limit of Monte-Carlo
#' confidence interval for probability of type 1} \item{predlow2}{Lower limit
#' of Monte-Carlo confidence interval for probability of type 2}
#' \item{predhigh2}{Upper limit of Monte-Carlo confidence interval for
#' probability of type 2} \item{predlow3}{Lower limit of Monte-Carlo confidence
#' interval for probability of type 3} \item{predhigh3}{Upper limit of
#' Monte-Carlo confidence interval for probability of type 3}
#' \item{predhigh1_rec}{Upper limit of Monte-Carlo confidence interval for
#' prediction of relapse} \item{predlow1_rec}{Lower limit of Monte-Carlo
#' confidence interval for prediction of relapse}
#' 
#' The following components are included in a 'predLongi' object obtained by
#' using prediction function for joint models with longitudinal data.
#' 
#' \item{npred}{Number of individual predictions} \item{x.time}{A vector of
#' prediction times of interest (used for plotting predictions): vector of
#' prediction times t if fixed window. Otherwise vector of prediction times
#' t+w} \item{window}{Prediction window or vector of prediction windows}
#' \item{group}{Id of each patient} \item{pred}{Estimation of probability}
#' \item{icproba}{Logical value. Were confidence intervals estimated?}
#' \item{predLow}{Lower limit of Monte-Carlo confidence intervals}
#' \item{predHigh}{Upper limit of Monte-Carlo confidence intervals}
#' \item{trivariate}{Logical value. Are the prediction calculated from the
#' trivariate model?}
#' @references A. Krol, L. Ferrer, JP. Pignon, C. Proust-Lima, M. Ducreux, O.
#' Bouche, S. Michiels, V. Rondeau (2016). Joint Model for Left-Censored
#' Longitudinal Data, Recurrent Events and Terminal Event: Predictive Abilities
#' of Tumor Burden for Cancer Evolution with Application to the FFCD 2000-05
#' Trial. \emph{Biometrics} \bold{72}(3) 907-16.
#' 
#' A. Mauguen, B. Rachet, S. Mathoulin-Pelissier, G. MacGrogan, A. Laurent, V.
#' Rondeau (2013). Dynamic prediction of risk of death using history of cancer
#' recurrences in joint frailty models. \emph{Statistics in Medicine},
#' \bold{32(30)}, 5366-80.
#' 
#' V. Rondeau, A. Laurent, A. Mauguen, P. Joly, C. Helmer (2015). Dynamic
#' prediction models for clustered and interval-censored outcomes:
#' investigating the intra-couple correlation in the risk of dementia.
#' \emph{Statistical Methods in Medical Research}
#' @keywords misc
#' @export
#' @examples
#' 
#' 
#' \dontrun{
#' 
#' #####################################################
#' #### prediction on a COX or SHARED frailty model ####
#' #####################################################
#' 
#' data(readmission)
#' #-- here is a generated cluster (31 clusters of 13 subjects)
#' readmission <- transform(readmission,group=id%%31+1)
#' 
#' #-- we compute predictions of death
#' #-- we extract last row of each subject for the time of death
#' readmission <- aggregate(readmission,by=list(readmission$id),
#'                          FUN=function(x){x[length(x)]})[,-1]
#' 
#' ##-- predictions on a Cox proportional hazard model --##
#' cox <- frailtyPenal(Surv(t.stop,death)~sex+dukes,
#' n.knots=10,kappa=10000,data=readmission)
#' 
#' #-- construction of the data frame for predictions
#' datapred <- data.frame(sex=0,dukes=0)
#' datapred$sex <- as.factor(datapred$sex)
#' levels(datapred$sex)<- c(1,2)
#' datapred$dukes <- as.factor(datapred$dukes)
#' levels(datapred$dukes)<- c(1,2,3)
#' datapred[1,] <- c(1,2) # man, dukes 2
#' datapred[2,] <- c(2,3) # woman, dukes 3
#' 
#' #-- prediction of death for two patients between 100 and 100+w,
#' #-- with w in (50,100,...,1900)
#' pred.cox <- prediction(cox,datapred,t=100,window=seq(50,1900,50))
#' plot(pred.cox)
#' 
#' #-- prediction of death for two patients between t and t+400,
#' #-- with t in (100,150,...,1500)
#' pred.cox2 <- prediction(cox,datapred,t=seq(100,1500,50),window=400)
#' plot(pred.cox2)
#' 
#' ##-- predictions on a shared frailty model for clustered data --##
#' sha <- frailtyPenal(Surv(t.stop,death)~cluster(group)+sex+dukes,
#' n.knots=10,kappa=10000,data=readmission)
#' 
#' #-- marginal prediction
#' # a group must be specified but it does not influence the results 
#' # in the marginal predictions setting
#' datapred$group[1:2] <- 1
#' pred.sha.marg <- prediction(sha,datapred,t=100,window=seq(50,1900,50))
#' plot(pred.sha.marg)
#' 
#' #-- conditional prediction, given a specific cluster (group=5)
#' datapred$group[1:2] <- 5
#' pred.sha.cond <- prediction(sha,datapred,t=100,window=seq(50,1900,50),
#'                             conditional = TRUE)
#' plot(pred.sha.cond)
#' 
#' ##-- marginal prediction of a recurrent event, on a shared frailty model
#' data(readmission)
#' 
#' datapred <- data.frame(t.stop=0,event=0,id=0,sex=0,dukes=0)
#' datapred$sex <- as.factor(datapred$sex)
#' levels(datapred$sex)<- c(1,2)
#' datapred$dukes <- as.factor(datapred$dukes)
#' levels(datapred$dukes)<- c(1,2,3)
#' 
#' datapred[1,] <- c(100,1,1,1,2) #man, dukes 2, 3 recurrent events
#' datapred[2,] <- c(200,1,1,1,2) 
#' datapred[3,] <- c(300,1,1,1,2) 
#' datapred[4,] <- c(350,0,2,1,2) #man, dukes 2  0 recurrent event
#' 
#' #-- Shared frailty model with gamma distribution
#' sha <- frailtyPenal(Surv(t.stop,event)~cluster(id)+sex+dukes,n.knots=10,
#' kappa=10000,data=readmission)
#' pred.sha.rec.marg <- prediction(sha,datapred,t=200,window=seq(50,1900,50),
#' event='Recurrent',MC.sample=100)
#' 
#' plot(pred.sha.rec.marg,conf.bands=TRUE)
#' 
#' ##-- conditional prediction of a recurrent event, on a shared frailty model
#' pred.sha.rec.cond <- prediction(sha,datapred,t=200,window=seq(50,1900,50),
#' event='Recurrent',conditional = TRUE,MC.sample=100)
#' 
#' plot(pred.sha.rec.cond,conf.bands=TRUE)
#' #####################################################
#' ######## prediction on a JOINT frailty model ########
#' #####################################################
#' 
#' data(readmission)
#' 
#' ##-- predictions of death on a joint model --##
#' joi <- frailtyPenal(Surv(t.start,t.stop,event)~cluster(id)
#' +sex+dukes+terminal(death),formula.terminalEvent=~sex
#' +dukes,data=readmission,n.knots=10,kappa=c(100,100),recurrentAG=TRUE)
#' 
#' #-- construction of the data frame for predictions
#' datapredj <- data.frame(t.stop=0,event=0,id=0,sex=0,dukes=0)
#' datapredj$sex <- as.factor(datapredj$sex)
#' levels(datapredj$sex) <- c(1,2)
#' datapredj$dukes <- as.factor(datapredj$dukes)
#' levels(datapredj$dukes) <- c(1,2,3)
#' datapredj[1,] <- c(100,1,1,1,2)
#' datapredj[2,] <- c(200,1,1,1,2)
#' datapredj[3,] <- c(300,1,1,1,2)
#' datapredj[4,] <- c(400,1,1,1,2)
#' datapredj[5,] <- c(380,1,2,1,2)
#' 
#' #-- prediction of death between 100 and 100+500 given relapses
#' pred.joint0 <- prediction(joi,datapredj,t=100,window=500,event = "Terminal")
#' print(pred.joint0)
#' 
#' #-- prediction of death between 100 and 100+w given relapses 
#' # (with confidence intervals)
#' pred.joint <- prediction(joi,datapredj,t=100,window=seq(50,1500,50),
#' event = "Terminal",MC.sample=100)
#' plot(pred.joint,conf.bands=TRUE)
#' # each y-value of the plot corresponds to the prediction between [100,x]
#' 
#' #-- prediction of death between t and t+500 given relapses
#' pred.joint2 <- prediction(joi,datapredj,t=seq(100,1000,50),
#' window=500,event = "Terminal")
#' plot(pred.joint2)
#' # each y-value of the plot corresponds to the prediction between [x,x+500], 
#' #or in the next 500
#' 
#' #-- prediction of relapse between 100 and 100+w given relapses 
#' # (with confidence intervals)
#' pred.joint <- prediction(joi,datapredj,t=100,window=seq(50,1500,50),
#' event = "Recurrent",MC.sample=100)
#' plot(pred.joint,conf.bands=TRUE)
#' # each y-value of the plot corresponds to the prediction between [100,x]
#' 
#' #-- prediction of relapse and death between 100 and 100+w given relapses 
#' # (with confidence intervals)
#' pred.joint <- prediction(joi,datapredj,t=100,window=seq(50,1500,50),
#' event = "Both",MC.sample=100)
#' plot(pred.joint,conf.bands=TRUE)
#' # each y-value of the plot corresponds to the prediction between [100,x]
#' 
#' #############################################################################
#' ### prediction on a JOINT model for longitudinal data and a terminal event ####
#' #############################################################################
#' 
#' 
#' data(colorectal)
#' data(colorectalLongi)
#' 
#' # Survival data preparation - only terminal events 
#' colorectalSurv <- subset(colorectal, new.lesions == 0)
#' 
#' #-- construction of the data-frame for predictions
#' #-- biomarker observations
#' datapredj_longi <- data.frame(id = 0, year = 0, tumor.size = 0, treatment = 0,
#'  age = 0, who.PS = 0, prev.resection = 0)
#' datapredj_longi$treatment <- as.factor(datapredj_longi$treatment)
#' levels(datapredj_longi$treatment) <- 1:2
#' datapredj_longi$age <- as.factor(datapredj_longi$age)
#' levels(datapredj_longi$age) <- 1:3
#' datapredj_longi$who.PS <- as.factor(datapredj_longi$who.PS)
#' levels(datapredj_longi$who.PS) <- 1:3
#' datapredj_longi$prev.resection <- as.factor(datapredj_longi$prev.resection)
#' levels(datapredj_longi$prev.resection) <- 1:2
#' 
#' # patient 1: increasing tumor size
#' datapredj_longi[1,] <- c(1, 0,1.2 ,2,1,1,1)
#' datapredj_longi[2,] <- c(1,0.3,1.4,2,1,1,1)
#' datapredj_longi[3,] <- c(1,0.6,1.9,2,1,1,1)
#' datapredj_longi[4,] <- c(1,0.9,2.5,2,1,1,1)
#' datapredj_longi[5,] <- c(1,1.5,3.9,2,1,1,1)
#' 
#' # patient 2: decreasing tumor size
#' datapredj_longi[6,] <- c(2, 0,1.2 ,2,1,1,1)
#' datapredj_longi[7,] <- c(2,0.3,0.7,2,1,1,1)
#' datapredj_longi[8,] <- c(2,0.5,0.3,2,1,1,1)
#' datapredj_longi[9,] <- c(2,0.7,0.1,2,1,1,1)
#' datapredj_longi[10,] <- c(2,0.9,0.1,2,1,1,1)
#' 
#' #-- terminal event
#' datapredj <- data.frame(id = 0, treatment = 0, age = 0, who.PS = 0,
#' prev.resection = 0)
#' datapredj$treatment <- as.factor(datapredj$treatment)
#' levels(datapredj$treatment) <- 1:2
#' datapredj$age <- as.factor(datapredj$age)
#' levels(datapredj$age) <- 1:3
#' datapredj$who.PS <- as.factor(datapredj$who.PS)
#' datapredj$prev.resection <- as.factor(datapredj$prev.resection)
#' levels(datapredj$prev.resection) <- 1:2
#' levels(datapredj$who.PS) <- 1:3
#' datapredj[1,] <- c(1,2,1,1,1)
#' datapredj[2,] <- c(2,2,1,1,1)
#' 
#' model.spli.CL <- longiPenal(Surv(time1, state) ~ age + treatment + who.PS
#' + prev.resection, tumor.size ~  year * treatment + age + who.PS , 
#' colorectalSurv, data.Longi = colorectalLongi, random = c("1", "year"),
#' id = "id", link = "Current-level", left.censoring = -3.33, n.knots = 6, 
#' kappa = 1)
#' 
#' #-- prediction of death between 1 year and 1+2 given history of the biomarker
#' pred.jointLongi0 <- prediction(model.spli.CL, datapredj, datapredj_longi,
#' t = 1, window = 2)
#' print(pred.jointLongi0)
#' 
#' #-- prediction of death between 1 year and 1+w given history of the biomarker
#' pred.jointLongi <- prediction(model.spli.CL, datapredj, datapredj_longi,
#' t = 1, window = seq(0.5, 2.5, 0.2), MC.sample = 100)
#' plot(pred.jointLongi, conf.bands = TRUE)
#' # each y-value of the plot corresponds to the prediction between [1,x]
#' 
#' #-- prediction of death between t and t+0.5 given history of the biomarker
#' pred.jointLongi2 <- prediction(model.spli.CL, datapredj, datapredj_longi,
#' t = seq(1, 2.5, 0.5), window = 0.5, MC.sample = 100)
#' plot(pred.jointLongi2, conf.bands = TRUE)
#' # each y-value of the plot corresponds to the prediction between [x,x+0.5], 
#' #or in the next 0.5
#' 
#' #############################################################################
#' ##### marginal prediction on a JOINT NESTED model for a terminal event ######
#' #############################################################################
#' #*--Warning! You can compute this prediction method with ONLY ONE family 
#' #*--by dataset of prediction. 
#' #*--Please make sure your data frame contains a column for individuals AND a 
#' #*--column for the reference number of the family chosen.
#' 
#' data(readmission)
#' readmissionNested <- transform(readmission,group=id%%30+1)
#' 
#' #-- construction of the data frame for predictions : 
#' #-- family 5 was selected for the prediction
#' 
#' DataPred <- readmissionNested[which(readmissionNested$group==5),]
#' 
#' #-- Fitting the model
#' modJointNested_Splines <- 
#' frailtyPenal(formula = Surv(t.start, t.stop, event)~subcluster(id)+ 
#' cluster(group) + dukes + terminal(death),formula.terminalEvent
#' =~dukes, data = readmissionNested, recurrentAG = TRUE,n.knots = 8, 
#' kappa = c(9.55e+9, 1.41e+12), initialize = TRUE)
#' 
#' #-- Compute prediction over the individuals 274 and 4 of the family 5
#' predRead <- prediction(modJointNested_Splines, data=DataPred,t=500,
#' window=seq(100,1500,200), conditional=FALSE, individual = c(274, 4))
#' 
#' 
#' #########################################################################
#' ##### prediction on TRIVARIATE JOINT model (linear and non-linear) ######
#' #########################################################################
#' 
#' data(colorectal)
#' data(colorectalLongi)
#' 
#' #-- construction of the data frame for predictions
#' #-- history of recurrences and terminal event
#' datapredj <- data.frame(time0 = 0, time1 = 0, new.lesions = 0, id = 0, 
#' treatment = 0, age = 0, who.PS = 0, prev.resection =0)
#' datapredj$treatment <- as.factor(datapredj$treatment)
#' levels(datapredj$treatment) <- 1:2
#' datapredj$age <- as.factor(datapredj$age)
#' levels(datapredj$age) <- 1:3
#' datapredj$who.PS <- as.factor(datapredj$who.PS)
#' levels(datapredj$who.PS) <- 1:3
#' datapredj$prev.resection <- as.factor(datapredj$prev.resection)
#' levels(datapredj$prev.resection) <- 1:2
#' 
#' datapredj[1,] <- c(0,0.4,1,1,2,1,1,1)
#' datapredj[2,] <- c(0.4,1.2,1,1,2,1,1,1)
#' datapredj[3,] <- c(0,0.5,1,2,2,1,1,1)
#' 
#' # Linear trivariate joint model
#' # (computation takes around 40 minutes)
#' model.trivPenal <-trivPenal(Surv(time0, time1, new.lesions) ~ cluster(id)
#' + age + treatment + who.PS +  terminal(state),
#' formula.terminalEvent =~ age + treatment + who.PS + prev.resection, 
#' tumor.size ~ year * treatment + age + who.PS, data = colorectal, 
#' data.Longi = colorectalLongi, random = c("1", "year"), id = "id", 
#' link = "Random-effects", left.censoring = -3.33, recurrentAG = TRUE,
#' n.knots = 6, kappa=c(0.01, 2), method.GH="Pseudo-adaptive",
#' n.nodes=7, init.B = c(-0.07, -0.13, -0.16, -0.17, 0.42, #recurrent events covarates
#' -0.23, -0.1, -0.09, -0.12, 0.8, -0.23, #terminal event covariates
#' 3.02, -0.30, 0.05, -0.63, -0.02, -0.29, 0.11, 0.74)) #biomarker covariates
#' 
#' #-- prediction of death between 1 year and 1+2
#' pred.jointTri0 <- prediction(model.trivPenal, datapredj, 
#' datapredj_longi, t = 1, window = 2)
#' print(pred.jointTri0)
#' 
#' #-- prediction of death between 1 year and 1+w
#' pred.jointTri <- prediction(model.trivPenal, datapredj, 
#' datapredj_longi, t = 1, window = seq(0.5, 2.5, 0.2), MC.sample = 100)
#' plot(pred.jointTri, conf.bands = TRUE)
#' 
#' #-- prediction of death between t and t+0.5
#' pred.jointTri2 <- prediction(model.trivPenal, datapredj, 
#' datapredj_longi, t = seq(1, 2.5, 0.5), window = 0.5, MC.sample = 100)
#' plot(pred.jointTri2, conf.bands = TRUE)
#' 
#' 
#' ###############################
#' 
#' # No information on dose - creation of a dummy variable 
#' colorectalLongi$dose <- 1
#' 
#' # (computation can take around 40 minutes)
#' model.trivPenalNL <- trivPenalNL(Surv(time0, time1, new.lesions) ~ cluster(id) + age + treatment
#'  + terminal(state), formula.terminalEvent =~ age + treatment, biomarker = "tumor.size",
#'  formula.KG ~ 1, formula.KD ~ treatment, dose = "dose", time.biomarker = "year", 
#'  data = colorectal, data.Longi =colorectalLongi, random = c("y0", "KG"), id = "id", 
#'  init.B = c(-0.22, -0.16, -0.35, -0.19, 0.04, -0.41, 0.23), init.Alpha = 1.86,
#'  init.Eta = c(0.5, 0.57, 0.5, 2.34), init.Biomarker = c(1.24, 0.81, 1.07, -1.53),
#'  recurrentAG = TRUE, n.knots = 5, kappa = c(0.01, 2), method.GH = "Pseudo-adaptive")
#' 
#' #-- prediction of death between 1 year and 1+2
#' pred.jointTriNL0 <- prediction(model.trivPenalNL, datapredj, 
#' datapredj_longi, t = 1, window = 2)
#' print(pred.jointTriNL0)
#' 
#' #-- prediction of death between 1 year and 1+w 
#' pred.jointTriNL <- prediction(model.trivPenalNL, datapredj, 
#' datapredj_longi, t = 1, window = seq(0.5, 2.5, 0.2), MC.sample = 100)
#' plot(pred.jointTriNL, conf.bands = TRUE)
#' 
#' #-- prediction of death between t and t+0.5
#' pred.jointTriNL2 <- prediction(model.trivPenalNL, datapredj, 
#' datapredj_longi, t = seq(2, 3, 0.2), window = 0.5, MC.sample = 100)
#' plot(pred.jointTriNL2, conf.bands = TRUE)
#' 
#' }
#' 
#' 
prediction <- function(fit, data, data.Longi, t, window, event = "Both", conditional=FALSE, MC.sample=0, individual){
  # set.global <- function (x, value) { # Combinee a la fonction aggregate personnalisee, permet de tenir compte de variable dependantes du temps
  # x <- deparse(substitute(x))
  # assign(x, value, pos=.GlobalEnv)
  # }
  
  if (missing(fit)) stop("Need a fit")
  if ((class(fit)!="frailtyPenal")&(class(fit)!="jointPenal")&(class(fit)!='longiPenal')&(class(fit)!='trivPenal')&(class(fit)!='jointNestedPenal')&class(fit)!='trivPenalNL') stop("The argument fit must be one of these objects : frailtyPenal; jointPenal; longiPenal; trivPenal or jointNestedPenal")
  if (fit$istop != 1) stop("Attempting to do predictions with a wrong model")
  
  if (class(fit) == "jointPenal"){
    if (fit$joint.clust == 0) stop("Prediction method is not available for joint model for clustered data")
    else if (fit$joint.clust == 2) stop("Prediction method is not available for joint general model")
  }
  
  if ((class(fit) != "jointNestedPenal") && (!missing(individual))) warning("No need for 'individual' option to predict anything other than a joint nested model.")
  
  if (missing(data)) stop("Need data to do some predictions")
  
  if (missing(data.Longi) & (class(fit)=="longiPenal") & (class(fit)=="trivPenal") & (class(fit) == "trivPenalNL")) stop("Need data.Longi to do predictions")
  
  if (missing(t) | missing(window)) stop("Need times and a window to do predictions")
  if (length(t)!=1 & length(window)!=1) stop("t and window can not be vector both at the same time")
  if (is.unsorted(t)) stop("Last time of predictions must be greater than first one")
  if (any(t < 0)) stop("Be careful, negative time input")
  
  if ((class(fit) == "frailtyPenal") && ((fit$Frailty==TRUE) & (!missing (event)) && (event != "Recurrent"))){ # Prediction modele shared pour evenement repete
    stop("Only 'Recurrent' event is allowed for a shared frailty modeling of parameters")
  }
  
  if (any(window <= 0)) stop("Window must be positive")
  
  event.type <- charmatch(event, c("Both", "Terminal", "Recurrent"), nomatch = 0)
  if (event.type == 0) {
    stop("event must be 'Both', 'Terminal' or 'Recurrent'")
  }
  if (class(fit) == "jointNestedPenal"){
    if(!missing(event) && (event.type != 2)) stop ("Only 'Terminal' event is allowed for a joint nested frailty modeling of parameters")
  }	
  
  if ((MC.sample < 0) | (MC.sample > 1000))  stop("MC.sample needs to be positive integer up to 1000")
  
  if ((class(fit)=="jointPenal" | class(fit)=='longiPenal' | class(fit)=='trivPenal' | class(fit)=='trivPenalNL') & (conditional)) stop("No conditional prediction available on a joint model")
  
  if(class(fit)=='jointPenal' | class(fit)=='trivPenal' | class(fit)=='trivPenalNL'){
    if (max(t+window) > max(fit$xR)) stop("Prediction times cannot exceed maximum time of observation")
    if (max(t+window) > max(fit$xD)) stop("Prediction times cannot exceed maximum time of observation")
  }
  
  if(class(fit)=='frailtyPenal'){
    if (max(t+window) > max(fit$x)) stop("Prediction times cannot exceed maximum time of observation")
  }
  
  if(class(fit)=='longiPenal'){
    if (max(t+window) > max(fit$xD)) stop("Prediction times cannot exceed maximum time of observation")
  }
  # seulement dans le cas du shared
  # if (missing(group)) type <- "marginal"
  # else type <- "conditional"
  
  # if (!(predTime >= min(fit$x1))) stop("predtime must be in the right range")
  # mettre un warning quand une variable est un factor dans le fit et pas dans le datapred => source d'erreur
  
  if (MC.sample==0) ICproba <- FALSE
  else ICproba <- TRUE
  
  # Ajout Julien pour eviter que des variables numeriques ne soient considerees comme "character"
  for (i in 1:length(data)) {
    if (!is.factor(data[,i])) data[,i] <- as.numeric(data[,i])
  }
  if (!missing(data.Longi)) {
    for (i in 1:length(data.Longi)) {
      if (!is.factor(data.Longi[,i])) data.Longi[,i] <- as.numeric(data.Longi[,i])
    }
  }
  
  np <- fit$npar
  b <- fit$b
  typeof <- fit$typeof
  nva1 <- fit$nvarRec
  nva2 <- fit$nvarEnd
  nva3 <- fit$nvarY
  if(class(fit) == "trivPenalNL"){
    nva3 <- fit$nvarKG
    nva4 <- fit$nvarKD
  }
  ng <- fit$groups
  nst <- 2
  HIH <- fit$varHIHtotal
  
  # a definir meme si non utilise
  nz <- 1
  zi <- 0
  nbintervR <- 1
  nbintervDC <- 1
  time <- 0
  timedc <- 0
  
  if(typeof == 0){
    nz <- fit$n.knots.temp
    zi <- fit$zi
  }
  if(typeof == 1){
    nbintervR <- fit$nbintervR
    nbintervDC <- fit$nbintervDC
    time <- fit$time
    timedc <- fit$timedc
  }
  
  # nombre de predictions a faire pour chaque individu
  moving.window <- FALSE
  if (length(t)==1) moving.window <- TRUE
  
  if (moving.window){
    predTime <- t
    timeAll <- t+window #seq(predTime+window,predMax,by=window)
    if (class(fit) == "jointPenal" | class(fit)== "trivPenal" | class(fit)== "trivPenalNL") window <- 0
  }else{
    predTime <- t[1]
    timeAll <- t+window
  }
  ntimeAll <- length(timeAll)
  formula_fit <- fit$formula
  
  # recuperation des profils d'individus pour la prediction
  m <- fit$call
  m2 <- match.call()
  
  m$formula.terminalEvent <- m$formula.LongitudinalData <- m$data.Longi <- m$random <- m$id  <- m$link <- m$left.censoring <- m$n.knots <- m$recurrentAG <- m$cross.validation <- m$kappa <- m$maxit <- m$hazard <- m$nb.int <- m$RandDist <- m$betaorder <- m$betaknots <- m$init.B <- m$LIMparam <- m$LIMlogl <- m$LIMderiv <- m$print.times <- m$init.Theta <- m$init.Alpha <- m$init.Random <- m$init.Eta <- m$Alpha <- m$method.GH <- m$intercept <- m$n.nodes <- m$jointGeneral <- m$initialize <- m$biomarker <- m$formula.KG <- m$formula.KD <- m$dose <- m$time.biomarker <- m$init.Biomarker <- m$BoxCox  <- m$Ksi <- m$init.Ksi <- m$... <- NULL
  
  m[[1]] <- as.name("model.frame")
  m3 <- m # pour recuperer les donnees du dataset initial en plus
  m3$formula <- fit$formula
  m[[3]] <- as.name(m2$data)
  
  if (class(fit) == "jointPenal" | class(fit)=="trivPenal" | class(fit)=="trivPenalNL" | class(fit) == "jointNestedPenal"){
    temp <- as.character(fit$formula[[2]])
    if (temp[1]=="Surv"){
      if (length(temp) == 4) fit$formula[[2]] <- paste(c("cbind(",temp[3],",",temp[4],")"),collapse=" ")
      else if (length(temp) == 3) fit$formula[[2]] <- paste(c("cbind(",temp[2],",",temp[3],")"),collapse=" ")
      else stop("Wrong Surv() function")
    }else{ # SurvIC
      if (length(temp) == 4) fit$formula[[2]] <- paste(c("cbind(",temp[2],",",temp[3],",",temp[4],")"),collapse=" ")
      else if (length(temp) == 5) fit$formula[[2]] <- paste(c("cbind(",temp[2],",",temp[3],",",temp[4],",",temp[5],")"),collapse=" ")
      else stop("Wrong SurvIC() function")
    }
    fit$formula <- unlist(strsplit(deparse(fit$formula)," "))
    fit$formula <- gsub("\"","",fit$formula)
    
    ter <- grep("terminal",fit$formula)
    if (ter==length(fit$formula)) m$formula <- as.formula(paste(fit$formula[-c(ter,max(which(fit$formula=="+")))],collapse=""))
    else m$formula <- as.formula(paste(fit$formula[-ter],collapse=""))
    
    if (fit$joint.clust==0){
      fit$formula <- unlist(strsplit(deparse(fit$formula)," "))
      clus <- grep("cluster",fit$formula)
      if (clus==length(fit$formula)) m$formula <- as.formula(paste(fit$formula[-c(clus,max(which(fit$formula=="+")))],collapse=""))
      else m$formula <- as.formula(paste(fit$formula[-clus],collapse=""))
    }
  }else{
    if (fit$Frailty==TRUE ){
      if(event == 'Recurrent'){
        if ((ICproba) && (conditional)){
          fitIC <- fit
          mIC <- m					
          fitIC$formula <- unlist(strsplit(deparse(fitIC$formula)," "))				
          clus <- grep("cluster",fitIC$formula)			
          mIC$formula <- as.formula(paste(fitIC$formula,collapse=""))	
        }
        temp <- as.character(fit$formula[[2]])
        if (temp[1]=="Surv"){
          # if (length(temp) == 4) 
          fit$formula[[2]] <- paste(c("cbind(",temp[length(temp)-1],",",temp[length(temp)],")"),collapse=" ")
          # else if (length(temp) == 3) fit$formula[[2]] <- paste(c("cbind(",temp[length(temp)-1],",",temp[length(temp)],")"),collapse=" ")
          # if (all(length(temp) != c(3,4))) stop("Wrong Surv() function")
        }
        fit$formula <- unlist(strsplit(deparse(fit$formula)," "))
        fit$formula <- gsub("\"","",fit$formula)
        clus <- grep("cluster",fit$formula)
        m$formula <- as.formula(paste(fit$formula,collapse=""))					
      }else{
        
        fit$formula <- unlist(strsplit(deparse(fit$formula)," "))				
        clus <- grep("cluster",fit$formula)
        # if (clus==length(fit$formula)) m$formula <- as.formula(paste(fit$formula[-c(clus,max(which(fit$formula=="+")))],collapse=""))
        # else m$formula <- as.formula(paste(fit$formula[-clus],collapse=""))				
        m$formula <- as.formula(paste(fit$formula,collapse=""))
        
      }
    }else if(class(fit)=='longiPenal'){
      fit$formula <- unlist(strsplit(deparse(fit$formula)," "))
      m$formula <- as.formula(paste(fit$formula,collapse=""))
    }else{#Cox
      m$formula <- fit$formula
    }
    if (!(fit$Frailty==TRUE && event == 'Recurrent')) m$formula[[2]] <- NULL # pas besoin du Surv dans formula, sauf si prediction pour evenement recurrent 
    # m$formula[[2]] <- NULL # pas besoin du Surv dans formula		
  }
  
  if(class(fit) == "jointPenal" | class(fit) == "jointNestedPenal"){
    if(is.null(fit$alpha))indic.alpha <- 0
    else indic.alpha <- 1 
    
  }
  indic.ksi <- 1
  if(class(fit) == "jointNestedPenal"&& is.null(fit$ksi))indic.ksi <- 0
  dataset <- eval(m, sys.parent())
  dataset3 <- eval(m3, sys.parent())
  
  typeofY <- attr(model.extract(dataset3, "response"),"type")
  Y <- model.extract(dataset3, "response")
  
  if (typeofY=="right") tt1 <- Y[,1]
  else tt1 <- Y[,2]
  class(m$formula) <- "formula"
  if ((ICproba) && (conditional) && (event == "Recurrent")) class(mIC$formula) <- "formula"
  special <- c("strata", "cluster", "subcluster", "terminal", "num.id", "timedep")		
  Terms <- terms(m$formula, special, data = data)
  
  if((event == "Recurrent") && (ICproba) && (conditional)) Terms2 <- terms(mIC$formula, special, data = data)	
  fit$formula <- Terms
  dropx <- NULL	
  
  if (class(fit) == 'jointPenal' | class(fit) == 'trivPenal' | class(fit) == 'trivPenalNL' | class(fit) == "jointNestedPenal"){	
    if (fit$joint.clust==1){ # joint classique	
      tempc <- untangle.specials(Terms, "cluster", 1:10)		
      nameGrup <- substr(tempc$vars,9,nchar(tempc$vars)-1)
      
      dropx <- c(dropx,tempc$terms)
      cluster <- strata(dataset[, tempc$vars], shortlabel = TRUE)
      uni.cluster <- unique(cluster)
      
      ic <- model.extract(dataset, "response")[,2]	
      
      npred <- length(uni.cluster)
      nrec <- max(table(cluster[ic==1]))
      
      if (temp[1]=="Surv"){
        Y <- NULL
        for (i in uni.cluster) {
          temp <- model.extract(dataset, "response")[,1]
          temp <- temp[cluster==i & ic==1]
          Y <- c(Y,c(temp,rep(NA,nrec-length(temp))))
          #Y <- c(Y,c(temp,rep(0,nrec-length(temp))))
          # NA car permet de differencier les t=0 donnes par l utilisateur de ceux qui sont remplis par defaut
        }
        predtimerec <- matrix(Y,nrow=npred,byrow=TRUE)
        trunctime <- rep(0,npred)
        lowertime <- rep(0,npred)
        uppertime <- rep(0,npred)			
        
      }else{
        stop("Predictions not allowed for interval-censored yet...") # a enlever plus tard
        predtimerec <- matrix(0,nrow=npred)
        if (length(temp) == 4){ # pas troncature
          temp <- model.extract(dataset, "response")
          trunctime <- rep(0,npred)
          lowertime <- temp[,1]
          uppertime <- temp[,2]
        }
        if (length(temp) == 5){ # troncature
          temp <- model.extract(dataset, "response")
          trunctime <- temp[,1]
          lowertime <- temp[,2]
          uppertime <- temp[,3]
        }
      }
    }else if (fit$joint.clust==3){ #Joint nested
      tempc <- untangle.specials(Terms, "cluster", 1:10)
      tempsbc <- untangle.specials(Terms, "subcluster", 1:10)	
      
      nameGrup <- substr(tempc$vars,9,nchar(tempc$vars)-1)
      nameSbGrup <- substr(tempsbc$vars,12,nchar(tempsbc$vars)-1)
      
      dropx1 <- c(dropx,tempc$terms)
      dropx2 <- c(dropx,tempsbc$terms)
      cluster <- strata(dataset[, tempc$vars], shortlabel = TRUE)
      uni.cluster <- unique(cluster)			
      subcluster <- strata(dataset[, tempsbc$vars], shortlabel = TRUE)
      uni.subcluster <- unique(subcluster)
      
      ic <- model.extract(dataset, "response")[,2]	
      Time_fam <- model.extract(dataset, "response")[,1]
      
      npred <- length(uni.subcluster)
      nrec <- max(table(subcluster[ic==1]))
      nrecList <- table(subcluster)
      
      if (temp[1]!="Surv") stop("Predictions not allowed for interval-censored yet...")
      Y <- NULL
      temp0 <- model.extract(dataset, "response")[,1]
      for (i in uni.subcluster) {	
        temp <- temp0[subcluster==i & ic==1]
        Y <- c(Y,c(temp,rep(NA,nrec-length(temp))))
      }
      predtimerec <- matrix(Y,nrow=npred,byrow=TRUE)
      
      trunctime <- rep(0,npred)
      lowertime <- rep(0,npred)
      uppertime <- rep(0,npred)
      
    }else{ # Conjoint donnees groupees
      tempnum <- untangle.specials(Terms, "num.id", 1:10)
      dropx <- c(dropx,tempnum$terms)
      num.id <- strata(dataset[, tempnum$vars], shortlabel = TRUE)
      uni.num.id <- unique(num.id)
      
      ic <- model.extract(dataset, "response")[,2]
      npred <- length(uni.num.id)
      nrec <- max(table(num.id[ic==1]))
      
      if (temp[1]=="Surv"){
        Y <- NULL
        for (i in uni.num.id){
          temp <- model.extract(dataset, "response")[,1]
          temp <- temp[num.id==i & ic==1]
          Y <- c(Y,c(temp,rep(0,nrec-length(temp))))
        }
        predtimerec <- matrix(Y,nrow=npred,byrow=TRUE)
        trunctime <- rep(0,npred)
        lowertime <- rep(0,npred)
        uppertime <- rep(0,npred)
      }else{
        stop("Predictions not allowed for interval-censored yet...") # a enlever plus tard
        predtimerec <- matrix(0,nrow=npred)
        if (length(temp) == 4){ # pas troncature
          temp <- model.extract(dataset, "response")
          trunctime <- rep(0,npred)
          lowertime <- temp[,1]
          uppertime <- temp[,2]
        }
        if (length(temp) == 5){ # troncature
          temp <- model.extract(dataset, "response")
          trunctime <- temp[,1]
          lowertime <- temp[,2]
          uppertime <- temp[,3]
        }
      }
    }
  }else{
    if (fit$Frailty){		
      
      #--Traitement donnees prediction
      tempc <- untangle.specials(Terms, "cluster", 1:10)
      dropx <- c(dropx,tempc$terms)
      cluster <- strata(dataset[, tempc$vars], shortlabel = TRUE)		
      uni.cluster <- unique(cluster)
      
      #--Traitement donnees fit$
      class(m3$formula) <- "formula"
      Terms3 <- terms(m3$formula, special, data = data)
      m3$formula <- Terms3
      tempc3 <- untangle.specials(Terms3, "cluster", 1:10)
      clusterfit <- strata(dataset3[, tempc3$vars], shortlabel = TRUE)
      uni.clusterfit <- unique(clusterfit)
      
      #--Pour labelliser les lignes de l output
      nameGrup <- substr(tempc3$vars,9,nchar(tempc3$vars)-1)
      
      if(event == 'Recurrent'){
        # class(m3$formula) <- "formula"
        # Terms3 <- terms(m3$formula, special, data = data)
        # tempc3 <- untangle.specials(Terms3, "cluster", 1:10)
        # clusterfit <- strata(dataset3[, tempc3$vars], shortlabel = TRUE)
        # uni.clusterfit <- unique(clusterfit)
        # nameGrup <- substr(tempc$vars,9,nchar(tempc$vars)-1)
        
        ic <- model.extract(dataset, "response")[,2]				
        npred <- length(uni.cluster)
        nrec <- max(table(cluster[ic==1]))
        
        if (temp[1]=="Surv"){
          Y <- NULL
          for (i in uni.cluster) {
            temp <- model.extract(dataset, "response")[,1]
            temp <- temp[cluster==i & ic==1]
            Y <- c(Y,c(temp,rep(NA,nrec-length(temp))))
            # NA car permet de differencier les t=0 donnes par l utilisateur de ceux qui sont remplis par defaut
          }
          predtimerec <- matrix(Y,nrow=npred,byrow=TRUE)
          predtimerectmp <- predtimerec 
          predtimerectmp[which(is.na(predtimerectmp))] <- 0
        }
      } #else{
      # # class(m3$formula) <- "formula"
      # # Terms3 <- terms(m3$formula, special, data = data)
      # # m3$formula <- Terms3
      # # tempc3 <- untangle.specials(Terms3, "cluster", 1:10)
      # # cluster <- strata(dataset3[, tempc3$vars], shortlabel = TRUE)
      # # uni.cluster <- unique(cluster)
      # nameGrup <- substr(tempc3$vars,9,nchar(tempc3$vars)-1)
      # }
    }
  }
  if (!is.null(dropx)){ 
    newTerms <- Terms[-dropx]
    if ((ICproba) && (event == 'Recurrent') && (conditional)) newTerms2 <- Terms2[-dropx]
  }
  else {
    newTerms <- Terms
    if ((ICproba) && (event == 'Recurrent') && (conditional)) newTerms2 <- Terms2		
  }
  
  X <- model.matrix(newTerms, dataset)
  
  if (ncol(X) > 1) X <- X[, -1, drop = FALSE]
  
  #####-----------------------------------------------------------------------------#####
  #####-&-&-&-&-&-&-&-&-&-&-&-Prediction for joint frailty model-&-&-&-&-&-&-&-&-&-&#####
  #####-----------------------------------------------------------------------------#####	
  
  if (class(fit) == "jointPenal"){	
    #--------ML 30-11-16...
    
    predtimerec <- predtimerec[order(unique(cluster)),]
    
    # listPrec <- NULL	
    # for (k in 1:nrow(predtimerec)){
    # tPrec <- which(predtimerec[k,] < predTime)
    # tPrec <- tPrec[length(tPrec)]
    # if (length(tPrec) == 0) tPrec <- 1 
    # listPrec <- c(listPrec,tPrec)
    # }
    
    taille = 0
    listPrec <- NULL  
    if(!is.matrix(predtimerec)) predtimerec <- matrix(predtimerec, nrow = npred)
    
    for (k in 1:nrow(predtimerec)){
      tPrec <- which(predtimerec[k,] < predTime)   
      if (length(tPrec) == 0) tPrec <- taille + 1 
      tPrec <- taille + length(tPrec)  
      
      rowTimes <- predtimerec[k,][which(!is.na(predtimerec[k,]))]
      if (length(rowTimes)==0) rowTimes <- 1
      taille = length(rowTimes)+taille
      listPrec <- c(listPrec,tPrec)                				
    }
    
    predtimerec <- replace(predtimerec, is.na(predtimerec),0) #30-11-16	
    # listPrec <- rep(listPrec,ncol(X))		
    
    if (fit$joint.clust==1){#vaxpred <- aggregate(X,by=list(cluster),FUN=function(x) {x[1]})[,-1]		
      X <- X[order(cluster),]
      if(!is.matrix(X)) X <- matrix(X, nrow = 1)
      
      vaxpred <- X[listPrec,]			
    }else{ #vaxpred <- aggregate(X,by=list(num.id),FUN=function(x) {x[1]})[,-1]	
      X <- X[order(num.id),]
      if(!is.matrix(X)) X <- matrix(X, nrow = 1)
      vaxpred <- X[listPrec,]
    } 		
    #... ----ML 30-11-16		
    
    # recuperation des variables partie deces
    m3 <- fit$call
    m2 <- match.call()
    
    m3$formula.LongitudinalData <- m3$data.Longi <- m3$random <- m3$id <- m3$link <- m3$left.censoring <- m3$n.knots <- m3$recurrentAG <- m3$cross.validation <- m3$kappa <- m3$maxit <- m3$hazard <- m3$nb.int <- m3$RandDist <- m3$betaorder <- m3$betaknots <- m3$init.B <- m3$LIMparam <- m3$LIMlogl <- m3$LIMderiv <- m3$print.times <- m3$init.Theta <- m3$init.Alpha <- m3$init.Ksi <- m3$Alpha <- m3$init.Random <- m3$init.Eta <- m3$method.GH <- m3$intercept <- m3$n.nodes <- m3$jointGeneral <- m3$initialize <-  m3$init.Biomarker <- m3$... <- NULL
    
    m3$formula <- formula_fit
    m3$formula[[3]] <- fit$formula.terminalEvent[[2]]		
    m3$formula.terminalEvent <- NULL
    m3[[1]] <- as.name("model.frame")
    m3[[3]] <- as.name(m2$data)
    
    temp <- as.character(m3$formula[[2]])
    if (temp[1]=="Surv"){
      if (length(temp) == 4) m3$formula[[2]] <- as.name(temp[3])
      else if (length(temp) == 3) m3$formula[[2]] <- as.name(temp[2])
      else stop("Wrong Surv function")
    }else{ # SurvIC
      if (length(temp) == 4) m3$formula[[2]] <- as.name(temp[2])
      else if (length(temp) == 5) m3$formula[[2]] <- as.name(temp[3])
      else stop("Wrong SurvIC function")
    }
    datasetdc <- eval(m3, sys.parent())
    class(m3$formula) <- "formula"
    special2 <- c("strata", "timedep")
    Terms2 <- terms(m3$formula, special2, data = data)
    
    X2 <- model.matrix(Terms2, datasetdc)
    if (ncol(X2) > 1) X2 <- X2[, -1, drop = FALSE]
    
    #------------ML: 01-12-16
    # listPrec <- rep(listPrec,ncol(X))
    if(!is.matrix(X2)) X2 <- matrix(X2, nrow = 1)
    if (fit$joint.clust==1){ #vaxdcpred <- aggregate(X2,by=list(cluster),FUN=function(x) {x[1]})[,-1]			
      X2 <- X2[order(cluster),]
    }else { # vaxdcpred <- aggregate(X2,by=list(num.id),FUN=function(x) {x[1]})[,-1]
      X2 <- X2[order(num.id),]
    }
    vaxdcpred <- as.matrix(X2)[listPrec,]
    cat("\n")
    cat("Calculating the probabilities ... \n")
    #if(fit$logNormal==0){	#Myriam modifie le 18-08-16	
    ans <- .Fortran(C_predict,
                    as.integer(np),
                    as.double(b),
                    as.integer(nz),
                    as.integer(nbintervR),
                    as.integer(nbintervDC),
                    as.integer(nva1),
                    as.integer(nva2),
                    as.integer(nst),
                    as.integer(typeof),
                    as.integer(event.type),
                    as.double(zi),
                    as.double(HIH),
                    as.double(time),
                    as.double(timedc),
                    as.integer(ntimeAll),
                    as.integer(npred),
                    as.double(predTime),
                    as.double(window),
                    as.double(predtimerec),
                    as.integer(nrec),
                    as.double(as.matrix(vaxpred)),
                    as.double(as.matrix(vaxdcpred)),
                    pred1=as.double(matrix(0,nrow=npred,ncol=ntimeAll)),
                    pred2=as.double(matrix(0,nrow=npred,ncol=ntimeAll)),
                    pred3=as.double(matrix(0,nrow=npred,ncol=ntimeAll)),
                    pred1_rec=as.double(matrix(0,nrow=npred,ncol=ntimeAll)),
                    predlow1=as.double(matrix(0,nrow=npred,ncol=ntimeAll)),
                    predhigh1=as.double(matrix(0,nrow=npred,ncol=ntimeAll)),
                    predlow2=as.double(matrix(0,nrow=npred,ncol=ntimeAll)),
                    predhigh2=as.double(matrix(0,nrow=npred,ncol=ntimeAll)),
                    predlow3=as.double(matrix(0,nrow=npred,ncol=ntimeAll)),
                    predhigh3=as.double(matrix(0,nrow=npred,ncol=ntimeAll)),
                    predlow1_rec=as.double(matrix(0,nrow=npred,ncol=ntimeAll)),
                    predhigh1_rec=as.double(matrix(0,nrow=npred,ncol=ntimeAll)),
                    icproba=as.integer(ICproba),
                    as.integer(MC.sample),
                    as.integer(fit$intcens),
                    as.double(trunctime),
                    as.double(lowertime),
                    as.double(uppertime),
                    as.integer(moving.window),
                    as.double(timeAll),
                    as.integer(fit$logNormal),
                    as.integer(indic.alpha)
    )
    
    # Myriam 18-08-2016 Fusion des fichiers predict et predict_logN
    
    # }else{ #AK: joint log-normal
    # ans <- .Fortran("predict",
    # as.integer(np),
    # as.double(b),
    # as.integer(nz),
    # as.integer(nbintervR),
    # as.integer(nbintervDC),
    # as.integer(nva1),
    # as.integer(nva2),
    # as.integer(nst),
    # as.integer(typeof),
    # as.double(zi),
    # as.double(HIH),
    # as.double(time),
    # as.double(timedc),
    # as.integer(ntimeAll),
    # as.integer(npred),
    # as.double(predTime),
    # as.double(window),
    # #as.integer(event.type),
    # as.double(predtimerec),
    # as.integer(nrec),
    # as.double(as.matrix(vaxpred)),
    # as.double(as.matrix(vaxdcpred)),
    # pred1=as.double(matrix(0,nrow=npred,ncol=ntimeAll)),
    # pred2=as.double(matrix(0,nrow=npred,ncol=ntimeAll)),
    # pred3=as.double(matrix(0,nrow=npred,ncol=ntimeAll)),
    # pred1_rec=as.double(matrix(0,nrow=npred,ncol=ntimeAll)),
    # predlow1=as.double(matrix(0,nrow=npred,ncol=ntimeAll)),
    # predhigh1=as.double(matrix(0,nrow=npred,ncol=ntimeAll)),
    # predlow2=as.double(matrix(0,nrow=npred,ncol=ntimeAll)),
    # predhigh2=as.double(matrix(0,nrow=npred,ncol=ntimeAll)),
    # predlow3=as.double(matrix(0,nrow=npred,ncol=ntimeAll)),
    # predhigh3=as.double(matrix(0,nrow=npred,ncol=ntimeAll)),
    # predlow1_rec=as.double(matrix(0,nrow=npred,ncol=ntimeAll)),
    # predhigh1_rec=as.double(matrix(0,nrow=npred,ncol=ntimeAll)),
    # icproba=as.integer(ICproba),
    # as.integer(MC.sample),
    # as.integer(fit$intcens),
    # as.double(trunctime),
    # as.double(lowertime),
    # as.double(uppertime),
    # as.integer(moving.window),
    # as.double(timeAll), # 38
    # PACKAGE = "frailtypack")
    # }
    out <- NULL
    out$call <- match.call()
    out$name.fit <- match.call()[[2]]
    out$npred <- npred
    out$window <- window
    out$predtimerec <- predtimerec
    out$moving.window <- moving.window
    out$event <- event.type
    if (moving.window){
      out$x.time <- timeAll
      out$t <- predTime
    }else{
      out$x.time <- timeAll - window
    }
    if (fit$joint.clust==1) out$group <- uni.cluster[order(uni.cluster)]
    else out$group <- uni.num.id
    
    if (!fit$intcens){
      if ((event.type == 1) || (event.type == 2)){
        out$pred1 <- matrix(ans$pred1,nrow=npred,ncol=ntimeAll)
        rownames(out$pred1) <- paste(nameGrup,out$group)
        colnames(out$pred1) <- paste("time=", out$x.time)
        
        out$pred3 <- matrix(ans$pred3,nrow=npred,ncol=ntimeAll)
        rownames(out$pred3) <- paste(nameGrup,out$group)
        colnames(out$pred3) <- paste("time=", out$x.time)
      }
      
      if ((event.type == 1) || (event.type == 3)){
        out$pred1_rec <- matrix(ans$pred1_rec,nrow=npred,ncol=ntimeAll)
        rownames(out$pred1_rec) <- paste(nameGrup,out$group)
        colnames(out$pred1_rec) <- paste("time=", out$x.time)	
      }
    }
    if ((event.type == 1) || (event.type == 2)){
      out$pred2 <- matrix(ans$pred2,nrow=npred,ncol=ntimeAll)
      rownames(out$pred2) <- paste(nameGrup,out$group)
      colnames(out$pred2) <- paste("time=", out$x.time)
    }
    # Myriam : Modification de l'affichage des resultats (ajout des temps de prediction)
    out$icproba <- ICproba
    if (ICproba){
      if (!fit$intcens){ 
        if ((event.type == 1) || (event.type == 2)){
          out$predlow1 <- matrix(ans$predlow1,nrow=npred,ncol=ntimeAll)
          out$predhigh1 <- matrix(ans$predhigh1,nrow=npred,ncol=ntimeAll)
          rownames(out$predlow1) <- paste(nameGrup,out$group)
          colnames(out$predlow1) <- paste("time=", out$x.time)
          rownames(out$predhigh1) <- paste(nameGrup,out$group)
          colnames(out$predhigh1) <- paste("time=", out$x.time)
          
          out$predlow3 <- matrix(ans$predlow3,nrow=npred,ncol=ntimeAll)
          out$predhigh3 <- matrix(ans$predhigh3,nrow=npred,ncol=ntimeAll)
          rownames(out$predlow3) <- paste(nameGrup,out$group)
          colnames(out$predlow3) <- paste("time=", out$x.time)
          rownames(out$predhigh3) <- paste(nameGrup,out$group)
          colnames(out$predhigh3) <- paste("time=", out$x.time)
        }
        if ((event.type == 1) || (event.type == 3)){
          out$predlow1_rec <- matrix(ans$predlow1_rec,nrow=npred,ncol=ntimeAll)
          out$predhigh1_rec <- matrix(ans$predhigh1_rec,nrow=npred,ncol=ntimeAll)
          rownames(out$predlow1_rec) <- paste(nameGrup,out$group)
          colnames(out$predlow1_rec) <- paste("time=", out$x.time)
          rownames(out$predhigh1_rec) <- paste(nameGrup,out$group)
          colnames(out$predhigh1_rec) <- paste("time=", out$x.time)
        }
      }
      if ((event.type == 1) || (event.type == 2)){
        out$predlow2 <- matrix(ans$predlow2,nrow=npred,ncol=ntimeAll)
        out$predhigh2 <- matrix(ans$predhigh2,nrow=npred,ncol=ntimeAll)
        rownames(out$predlow2) <- paste(nameGrup,out$group)
        colnames(out$predlow2) <- paste("time=", out$x.time)
        rownames(out$predhigh2) <- paste(nameGrup,out$group)
        colnames(out$predhigh2) <- paste("time=", out$x.time)
      }
    }
    out$joint.clust <- fit$joint.clust
    out$intcens <- fit$intcens
    
    cat("Predictions done for",npred,"subjects and",ntimeAll,"times \n")
    class(out) <- c("predJoint")
    
    #####----------------------------------------------------------------------------#####
    #####-&-&-&-&-&-&-Prediction joint for longitudinal data and terminal event-&-&-&#####
    #####-&-&-&-&-&-&-		or longitudinal datas, recurrent events 		   -&-&-&#####
    #####-&-&-&-&-&-&-				and a terminal event					   -&-&-&#####
    #####----------------------------------------------------------------------------#####
  }else if(class(fit)=="longiPenal" | class(fit)=="trivPenal" | class(fit)=="trivPenalNL"){
    cat("\n")
    cat("Calculating the probabilities ... \n")
    
    if(class(fit)=="longiPenal"){	
      expBX <- exp(X %*% fit$coef[1:fit$nvarEnd])
    }else{
      #-----------ML:07-12-16
      taille = 0
      listPrec <- NULL  
      if(!is.matrix(predtimerec)) predtimerec <- matrix(predtimerec, nrow = npred)
      for (k in 1:nrow(predtimerec)){
        tPrec <- which(predtimerec[k,] < predTime)   
        if (length(tPrec) == 0) tPrec <- taille + 1 
        tPrec <- taille + length(tPrec) 
        
        rowTimes <- predtimerec[k,][which(!is.na(predtimerec[k,]))]
        if (length(rowTimes)==0) rowTimes <- 1
        taille = length(rowTimes)+taille
        listPrec <- c(listPrec,tPrec)
      }			
      predtimerec <- replace(predtimerec, is.na(predtimerec),0) #30-11-16
      listPrec <- listPrec[order(unique(cluster))]
      vaxpred <- X[listPrec,]
      X <- X[order(cluster),]
      #-----------ML:07-12-16
      
      # recuperation des variables partie deces
      m3 <- fit$call
      m2 <- match.call()
      
      
      m3$formula.LongitudinalData <- m3$data.Longi <- m3$random <- m3$id <- m3$link <- m3$left.censoring <- m3$n.knots <- m3$recurrentAG <- m3$cross.validation <- m3$kappa <- m3$maxit <- m3$hazard <- m3$nb.int <- m3$RandDist <- m3$betaorder <- m3$betaknots <- m3$init.B <- m3$LIMparam <- m3$LIMlogl <- m3$LIMderiv <- m3$print.times <- m3$init.Theta <- m3$init.Alpha <- m3$Alpha <- m3$init.Random <- m3$init.Eta <- m3$method.GH <- m3$intercept <- m3$n.nodes <- m3$biomarker <- m3$formula.KG <- m3$formula.KD <- m3$dose <- m3$time.biomarker <- m3$BoxCox <- m3$init.Biomarker <- m3$init.Ksi<- m3$... <- NULL
      
      m3$formula <- formula_fit
      m3$formula[[3]] <- fit$formula.terminalEvent[[2]]
      m3$formula.terminalEvent <- NULL
      m3[[1]] <- as.name("model.frame")
      m3[[3]] <- as.name(m2$data)
      
      temp <- as.character(m3$formula[[2]])
      
      if (length(temp) == 4) m3$formula[[2]] <- as.name(temp[3])
      else if (length(temp) == 3) m3$formula[[2]] <- as.name(temp[2])
      else stop("Wrong Surv function")
      
      datasetdc <- eval(m3, sys.parent())
      class(m3$formula) <- "formula"
      special2 <- c("strata", "timedep")
      Terms2 <- terms(m3$formula, special2, data = data)
      X2 <- model.matrix(Terms2, datasetdc)
      if (ncol(X2) > 1) X2 <- X2[, -1, drop = FALSE]			
      
      vaxdcpred <- aggregate(X2,by=list(cluster),FUN=function(x) {x[1]})[,-1]			
    }
    # nombre de predictions a faire pour chaque individu
    if (moving.window){
      sequence2 <- t+window
      sequence <- rep(predTime,times=length(sequence2))
    }else{
      sequence <- t 
      sequence2 <- t+window
    }
    predMat <- NULL
    
    if(class(fit)!= "trivPenalNL"){
      m2 <- fit$call
      
      m2$formula <-  m2$data <- m2$random <- m2$id <- m2$link <- m2$n.knots <- m2$kappa <- m2$maxit <- m2$hazard <- m2$nb.int <- m2$betaorder <- m2$betaknots <- m2$init.B <- m2$LIMparam <- m2$LIMlogl <- m2$LIMderiv <- m2$print.times <- m2$left.censoring <- m2$init.Random <- m2$init.Eta <- m2$method.GH <- m2$init.Alpha <- m2$intercept <- m2$n.nodes  <- m2$... <- NULL
      
      special <- c("strata", "cluster", "subcluster", "terminal","num.id","timedep")
      
      #========= Longitudinal Data preparation =========================
      class(fit$formula.LongitudinalData) <- "formula"
      TermsY <- terms(fit$formula.LongitudinalData, special, data = data.Longi)
      llY <- attr(TermsY, "term.labels")#liste des variables explicatives
      ord <- attr(TermsY, "order")
      
      #=========================================================>
      name.Y <- as.character(attr(TermsY, "variables")[[2]])	       
      
      if(class(fit)=="longiPenal") X <- X[order(unique(data.Longi$id)),]
      
      data.Longi <- data.Longi[order(data.Longi$id),]
      yy <- data.Longi[,which(names(data.Longi)==name.Y)]
      
      # on identifie les variables explicatives facteurs avec nombre de niveau plus que 2
      ind.placeY <- which(llY%in%names(which(lapply(data.Longi[,which(names(data.Longi)%in%llY)],function(x) length(levels(x)))>2)))
      vec.factorY <- NULL
      vec.factorY <- c(vec.factorY,llY[ind.placeY])
      
      mat.factorY <- matrix(vec.factorY,ncol=1,nrow=length(vec.factorY))
      
      # Fonction servant a prendre les termes entre "as.factor"
      vec.factorY <-apply(mat.factorY,MARGIN=1,FUN=function(x){
        if (length(grep("as.factor",x))>0){
          pos1 <- grep("r",unlist(strsplit(x,split="")))[1]+2
          pos2 <- length(unlist(strsplit(x,split="")))-1
          return(substr(x,start=pos1,stop=pos2))
        }else{
          return(x)
        }
      })
      
      ind.placeY <- grep(paste(vec.factorY,collapse="|"),llY)
      if(is.factor(data.Longi[,names(data.Longi)==llY[1]])) X_L<- as.numeric(data.Longi[,names(data.Longi)==llY[1]])-1
      else X_L<- data.Longi[,names(data.Longi)==llY[1]]
      if(length(llY)>1){
        for(i in 2:length(llY)){
          if(is.factor(data.Longi[,names(data.Longi)==llY[i]]))X_L<- cbind(X_L,as.numeric(data.Longi[,names(data.Longi)==llY[i]])-1)
          else X_L<- cbind(X_L,data.Longi[,names(data.Longi)==llY[i]])
        }
      }
      #X_L<- data.Longi[,names(data.Longi)%in%(llY)]
      
      if(sum(ord)>length(ord)){
        for(i in 1:length(ord)){
          if(ord[i]>1){
            v1 <- strsplit(as.character(llY[i]),":")[[1]][1]
            v2 <- strsplit(as.character(llY[i]),":")[[1]][2]
            if(is.factor(data.Longi[,names(data.Longi)==v1]) && length(levels(data.Longi[,names(data.Longi)==v1]))>2) stop("Interactions not allowed for factors with 3 or more levels (yet)")
            if(is.factor(data.Longi[,names(data.Longi)==v2]) && length(levels(data.Longi[,names(data.Longi)==v2]))>2) stop("Interactions not allowed for factors with 3 or more levels (yet)")
            if(is.factor(data.Longi[,names(data.Longi)==v1]) || !is.factor(data.Longi[,names(data.Longi)==v2])) {
              X_L <- cbind(X_L,(as.numeric(data.Longi[,names(data.Longi)==v1])-1)*data.Longi[,names(data.Longi)==v2])
              llY[i]<-paste(llY[i],levels(data.Longi[,names(data.Longi)==v1])[2],sep="")
            }else if (!is.factor(data.Longi[,names(data.Longi)==v1]) || is.factor(data.Longi[,names(data.Longi)==v2])) {
              X_L <- cbind(X_L,data.Longi[,names(data.Longi)==v1]*(as.numeric(data.Longi[,names(data.Longi)==v2])-1))
              llY[i]<-paste(llY[i],levels(data.Longi[,names(data.Longi)==v2])[2],sep="")
            }else{
              X_L <- cbind(X_L,data.Longi[,names(data.Longi)==v1]*data.Longi[,names(data.Longi)==v2])
            }
          }
        }
      }
      
      if(dim(X_L)[2]!=length(llY))stop("The variables in the longitudinal part must be in the data.Longi")
      X_L <- as.data.frame(X_L)
      names(X_L) <- llY
      Intercept <- rep(1,dim(X_L)[1])
      if(fit$intercept)X_L <- cbind(Intercept,X_L)
      X_Lall<- X_L
      "%+%"<- function(x,y) paste(x,y,sep="")
      if(length(vec.factorY) > 0){
        for(i in 1:length(vec.factorY)){
          X_L <- cbind(X_L[,-(which(names(X_L)==vec.factorY[i]))],model.matrix(as.formula("~"%+%0%+%"+"%+%paste(vec.factorY[i], collapse= "+")), data.Longi)[,-1])
        }
        vect.factY<-names(X_L)[which(!(names(X_L)%in%llY))]
        occurY <- rep(0,length(vec.factorY))
        for(i in 1:length(vec.factorY)){
          #occur[i] <- sum(vec.factor[i] == vect.fact)
          occurY[i] <- length(grep(vec.factorY[i],vect.factY))
        }
      }
      
      if (ncol(X_L) == 0){
        noVarY <- 1
      }else{
        noVarY <- 0
      }
      #=========================================================>
      clusterY <- data.Longi$id
      maxy_rep <- max(table(clusterY))
      uni.cluster<-as.factor(unique(clusterY))
      npred <- length(uni.cluster)
      nvarY<-ncol(X_L) #nvar==1 correspond a 2 situations:
      varY <- as.matrix(sapply(X_L, as.numeric))
      
      #=======================================>
      #======= Construction du vecteur des indicatrice
      if(length(vec.factorY) > 0){
        #ind.place <- ind.place -1	
        k <- 0
        for(i in 1:length(vec.factorY)){
          ind.placeY[i] <- ind.placeY[i]+k
          k <- k + occurY[i]-1		
        }
      }
      
      
    }else{#trviPenalNL
      m3 <- fit$call # longitudinal (KG)
      m3$formula <- m3$formula.terminalEvent <- m3$biomarker <- m3$formula.KD <- m3$dose <- m3$data <- m3$recurrentAG <- m3$random <- m3$id <- m3$link <- m3$n.knots <- m3$kappa <- m3$maxit <- m3$hazard <- m3$init.B <- m3$LIMparam <- m3$LIMlogl <- m3$LIMderiv <- m3$print.times <- m3$left.censoring <- m3$init.Random <- m3$init.Eta <- m3$init.Alpha <- m3$method.GH <- m3$n.nodes  <- m3$init.GH <- m3$time.biomarker <- m3$BoxCox <- m3$init.Biomarker <- m3$... <- NULL
      Names.data.Longi <- m3$data.Longi
      formula.KG <- fit$formula.KG
      
      m4 <- fit$call # longitudinal (KD)
      m4$formula <- m4$formula.terminalEvent <- m4$biomarker <- m4$formula.KG <- m4$dose <- m4$data <- m4$recurrentAG <- m4$random <- m4$id <- m4$link <- m4$n.knots <- m4$kappa <- m4$maxit <- m4$hazard <- m4$init.B <- m4$LIMparam <- m4$LIMlogl <- m4$LIMderiv <- m4$print.times <- m4$left.censoring <- m4$init.Random <- m4$init.Eta <- m4$init.Alpha <- m4$method.GH <- m4$n.nodes <- m4$init.GH <- m4$time.biomarker <- m4$BoxCox <- m4$init.Biomarker <- m4$... <- NULL
      
      Y <- data.Longi[,which(names(data.Longi)==fit$biomarker)]
      
      if(!is.null(formula.KG[3]) && formula.KG[3] != "1()"){
        
        TermsKG <- if (missing(data.Longi)){
          terms(formula.KG, special)
        }else{
          terms(formula.KG, special, data = data.Longi)
        }
        
        ord <- attr(TermsKG, "order") # longueur de ord=nbre de var.expli
        
        #si pas vide tous si il ya au moins un qui vaut 1 on arrete
        
        m2$formula.KG <- TermsKG
        
        
        m2[[1]] <- as.name("model.frame") # m[[1]]=frailtypenal, il le remplace par model.frame en fait
        
        
        if (NROW(m3) == 0)stop("No (non-missing) observations") #nombre ligne different de 0
        
        llKG <- attr(TermsKG, "term.labels")#liste des variables explicatives
        
        
        #=========================================================>
        
        name.KG <- as.character(attr(TermsKG, "variables")[[2]])
        KG <- data.Longi[,which(names(data.Longi)==name.KG)]
        
        
        # on identifie les variables explicatives facteurs avec nombre de niveau plus que 2
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
        
        
        
        if(is.factor(data.Longi[,names(data.Longi)==llKG.real.names[1]]))X_L<- as.numeric(data.Longi[,names(data.Longi)==llKG.real.names[1]])-1
        else X_L<- data.Longi[,names(data.Longi)==llKG.real.names[1]]
        
        if(length(llKG) == 1)X_L <- as.data.frame(X_L)
        
        
        if(length(llKG)>1){
          
          
          for(i in 2:length(llKG.real.names)){
            
            if(is.factor(data.Longi[,names(data.Longi)==llKG.real.names[i]]))X_L<- cbind(X_L,as.numeric(data.Longi[,names(data.Longi)==llKG.real.names[i]])-1)
            else X_L<- cbind(X_L,data.Longi[,names(data.Longi)==llKG.real.names[i]])
          }
        }
        #X_L<- data.Longi[,names(data.Longi)%in%(llY)]
        
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
                  X_L <- cbind(X_L,dummy[,j]*v2)
                  if(i>1 && i<length(llKG.fin))llKG.fin <- c(llKG.fin[1:(i-1+j-2)],paste(name_v1,".",levels(v1)[j],":",name_v2,sep=""),llKG.fin[(i+1+j-2):length(llKG.fin)])
                  else if(i==length(llKG.fin))llKG.fin <- c(llKG.fin[1:(i-1+j-2)],paste(name_v1,".",levels(v1)[j],":",name_v2,sep=""))
                  else llKG.fin <- c(paste(name_v1,".",levels(v1)[j],":",name_v2,sep=""),llKG.fin[(2+j-2):length(llKG.fin)])
                }
                
              }else if(!is.factor(v1) && is.factor(v2)){
                
                dummy <- model.matrix( ~ v2 - 1)
                #  if(length(levels(v2)>2))vec.factorY <- c(vec.factorY,paste(name_v1,":",name_v2,sep=""))
                for(j in 2:length(levels(v2))){
                  
                  X_L <- cbind(X_L,dummy[,j]*v1)
                  
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
                    
                    X_L <- cbind(X_L,dummy1[,j]*dummy2[,k])
                    if(i>1 && i<length(llKG.fin))llKG.fin <- c(llKG.fin[1:(i-1+j-2+k-2)],paste(name_v1,levels(v1)[j],":",name_v2,levels(v2)[k],sep=""),llKG.fin[(i+1+j-2+k-2):length(llKG.fin)])
                    else if(i==length(llKG.fin))llKG.fin <- c(llKG.fin[1:(i-1+j-2+k-2)],paste(name_v1,levels(v1)[j],":",name_v2,levels(v2)[k],sep=""))
                    else llKG.fin <- c(paste(name_v1,levels(v1)[j],":",name_v2,levels(v2)[k],sep=""),llKG.fin[(2+j-2+k-2):length(llKG.fin)])
                  }
                } 
              }else{
                
                X_L <- cbind(X_L,v1*v2)
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
        X_L <- as.data.frame(X_L)
        if(dim(X_L)[2]!=length(llKG.fin))stop("The variables in the longitudinal part must be in the data.Longi")
        
        names(X_L) <- llKG.fin
        
        
        X_Lall<- X_L
        "%+%"<- function(x,y) paste(x,y,sep="")
        if(length(vec.factorKG) > 0){
          for(i in 1:length(vec.factorKG)){
            if(length(grep(":",vec.factorKG[i]))==0){
              
              factor.spot <- which(names(X_L)==vec.factorKG[i])
              
              if(factor.spot<ncol(X_L))  X_L <- cbind(X_L[1:(factor.spot-1)],model.matrix(as.formula("~"%+%0%+%"+"%+%paste(vec.factorKG[i], collapse= "+")), model.frame(~.,data.Longi,na.action=na.pass))[,-1],X_L[(factor.spot+1):ncol(X_L)])
              else X_L <- cbind(X_L[1:(factor.spot-1)],model.matrix(as.formula("~"%+%0%+%"+"%+%paste(vec.factorKG[i], collapse= "+")), model.frame(~.,data.Longi,na.action=na.pass))[,-1])
              
            } }
          
          
        }
        
        varKG <- as.matrix(sapply(X_L, as.numeric))
        
        nsujety<-nrow(X_L)
        
      }
      
      
      #=========================================================>
      
      clusterY <- data.Longi$id
      
      max_rep <- max(table(clusterY))
      uni.clusterY<-as.factor(unique(clusterY))
      
      
      if(is.null(formula.KG[3]) | formula.KG[3] == "1()"){
        varKG <- c()#rep(0, dim(data.Longi)[1])
        nsujety <- length(Y)
      }
      
      formula.KD <- fit$formula.KD
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
        
        
        if(is.factor(data.Longi[,names(data.Longi)==llKD.real.names[1]]))X_L<- as.numeric(data.Longi[,names(data.Longi)==llKD.real.names[1]])-1
        else X_L<- data.Longi[,names(data.Longi)==llKD.real.names[1]]
        
        
        
        
        if(length(llKD)>1){
          for(i in 2:length(llKD.real.names)){
            
            if(is.factor(data.Longi[,names(data.Longi)==llKD.real.names[i]]))X_L<- cbind(X_L,as.numeric(data.Longi[,names(data.Longi)==llKD.real.names[i]])-1)
            else X_L<- cbind(X_L,data.Longi[,names(data.Longi)==llKD.real.names[i]])
          }}
        
        #X_L<- data.Longi[,names(data.Longi)%in%(llY)]
        
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
                  X_L <- cbind(X_L,dummy[,j]*v2)
                  if(i>1 && i<length(llKD.fin))llKD.fin <- c(llKD.fin[1:(i-1+j-2)],paste(name_v1,".",levels(v1)[j],":",name_v2,sep=""),llKD.fin[(i+1+j-2):length(llKD.fin)])
                  else if(i==length(llKD.fin))llKD.fin <- c(llKD.fin[1:(i-1+j-2)],paste(name_v1,".",levels(v1)[j],":",name_v2,sep=""))
                  else llKD.fin <- c(paste(name_v1,".",levels(v1)[j],":",name_v2,sep=""),llKD.fin[(2+j-2):length(llKD.fin)])
                }
                
              }else if(!is.factor(v1) && is.factor(v2)){
                
                dummy <- model.matrix( ~ v2 - 1)
                #  if(length(levels(v2)>2))vec.factorY <- c(vec.factorY,paste(name_v1,":",name_v2,sep=""))
                for(j in 2:length(levels(v2))){
                  
                  X_L <- cbind(X_L,dummy[,j]*v1)
                  
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
                    
                    X_L <- cbind(X_L,dummy1[,j]*dummy2[,k])
                    if(i>1 && i<length(llKD.fin))llKD.fin <- c(llKD.fin[1:(i-1+j-2+k-2)],paste(name_v1,levels(v1)[j],":",name_v2,levels(v2)[k],sep=""),llKG.fin[(i+1+j-2+k-2):length(llKG.fin)])
                    else if(i==length(llKD.fin))llKD.fin <- c(llKD.fin[1:(i-1+j-2+k-2)],paste(name_v1,levels(v1)[j],":",name_v2,levels(v2)[k],sep=""))
                    else llKD.fin <- c(paste(name_v1,levels(v1)[j],":",name_v2,levels(v2)[k],sep=""),llKD.fin[(2+j-2+k-2):length(llKD.fin)])
                  }
                } 
              }else{
                
                X_L <- cbind(X_L,v1*v2)
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
        X_L <- as.data.frame(X_L)
        if(dim(X_L)[2]!=length(llKD.fin))stop("The variables in the longitudinal part must be in the data.Longi")
        
        names(X_L) <- llKD.fin
        
        
        X_Lall<- X_L
        "%+%"<- function(x,y) paste(x,y,sep="")
        if(length(vec.factorKD) > 0){
          for(i in 1:length(vec.factorKD)){
            if(length(grep(":",vec.factorKD[i]))==0){
              
              factor.spot <- which(names(X_L)==vec.factorKD[i])
              
              if(factor.spot<ncol(X_L))  X_L <- cbind(X_L[1:(factor.spot-1)],model.matrix(as.formula("~"%+%0%+%"+"%+%paste(vec.factorKD[i], collapse= "+")), model.frame(~.,data.Longi,na.action=na.pass))[,-1],X_L[(factor.spot+1):ncol(X_L)])
              else X_L <- cbind(X_L[1:(factor.spot-1)],model.matrix(as.formula("~"%+%0%+%"+"%+%paste(vec.factorKD[i], collapse= "+")), model.frame(~.,data.Longi,na.action=na.pass))[,-1])
              
            } 
          }
          
        }
        
        
        varKD <- as.matrix(sapply(X_L, as.numeric))
        
        
      }else{
        varKD <- c()#rep(0, dim(data.Longi)[1])
      }
      
      matzy <- NULL       # here matzy is for time and dose
      matzy <- cbind(data.Longi[,which(colnames(data.Longi)==fit$time.biomarker)],data.Longi[,which(colnames(data.Longi)==fit$dose)])
      
      
      
      matzy <- as.matrix(matzy)
      
    }
    
    
    if(fit$link=="Random-effects")link <- 1
    if(fit$link=="Current-level") link <- 2
    if(fit$leftCensoring==FALSE){
      s_cag_id = 0
      s_cag = 0
    }else{
      s_cag_id = 1
      s_cag = fit$leftCensoring.threshold
    }				
    
    if(class(fit)=="longiPenal"){
      ans <- .Fortran(C_predict_biv,
                      as.integer(np),
                      as.double(b),
                      as.integer(nz),
                      as.integer(nva2),
                      as.integer(nva3),
                      as.integer(fit$ne_re),
                      as.integer(fit$netadc),
                      as.integer(link),
                      as.integer(nst),
                      as.integer(typeof),
                      as.double(zi),
                      as.double(HIH),
                      as.integer(ntimeAll),
                      as.integer(npred),
                      as.double(predTime),
                      as.double(window),
                      as.integer(fit$max_rep),
                      as.double(yy),
                      as.double(as.matrix(X)),
                      as.double(as.matrix(varY)),
                      as.integer(clusterY),
                      as.integer(unique(clusterY)),
                      as.integer(length(clusterY)),
                      pred=as.double(matrix(0,nrow=npred,ncol=ntimeAll)),
                      predlow=as.double(matrix(0,nrow=npred,ncol=ntimeAll)),
                      predhigh=as.double(matrix(0,nrow=npred,ncol=ntimeAll)),
                      icproba=as.integer(ICproba),
                      as.integer(MC.sample),
                      as.integer(moving.window),
                      as.double(timeAll),
                      as.integer(s_cag_id),
                      as.double(s_cag)
      )
      
      predMat <- matrix(ans$pred,nrow=nrow(data),ncol=ntimeAll)
      predMatLow <- matrix(ans$predlow,nrow=nrow(data),ncol=ntimeAll)
      predMatHigh <- matrix(ans$predhigh,nrow=nrow(data),ncol=ntimeAll)
      
      out <- NULL
      out$call <- match.call()
      out$name.fit <- match.call()[[2]]
      out$npred <- npred
      out$moving.window <- moving.window
      if (moving.window){
        out$x.time <- sequence2
        out$t <- predTime
      }else{
        out$x.time <- sequence
      }
      out$group <- uni.cluster
      out$pred <- predMat
      colnames(out$pred) <- c("times",rep(" ",dim(out$pred)[2]-1))
      #rownames(out$pred) <- paste("ind",1:out$npred)
      rownames(out$pred) <- paste("ind",unique(clusterY))
      
      out$icproba <- ICproba
      if (ICproba){
        out$predLow <- predMatLow
        out$predHigh <- predMatHigh
        colnames(out$predLow) <- c("times",rep(" ",dim(out$predLow)[2]-1))
        # rownames(out$predLow) <- paste("ind",1:out$npred)
        rownames(out$predLow) <- paste("ind",unique(clusterY))
        colnames(out$predHigh) <- c("times",rep(" ",dim(out$predHigh)[2]-1))
        # rownames(out$predHigh) <- paste("ind",1:out$npred)
        rownames(out$predHigh) <- paste("ind",unique(clusterY))
      }
      out$window <- window
      out$trivariate <- FALSE	
      
    }else if(class(fit)=="trivPenal"){             
      predtimerec <- predtimerec[order(unique(cluster)),]	
      if(!is.matrix(predtimerec)) predtimerec <- matrix(predtimerec, nrow = npred)
      ans <- .Fortran(C_predict_tri,
                      as.integer(np),
                      as.double(b),
                      as.integer(nz),
                      as.integer(nva1),
                      as.integer(nva2),
                      as.integer(nva3),
                      as.integer(fit$ne_re),
                      as.integer(fit$netar),
                      as.integer(fit$netadc),
                      as.integer(link),
                      as.integer(nst),
                      as.integer(typeof),
                      as.double(zi),
                      as.double(HIH),
                      as.integer(ntimeAll),
                      as.integer(npred),
                      as.double(predTime),
                      as.double(window),
                      as.double(predtimerec),
                      as.integer(nrec),
                      as.integer(fit$max_rep),
                      as.double(yy),
                      as.double(as.matrix(vaxpred)),
                      as.double(as.matrix(X)),
                      as.double(as.matrix(varY)),
                      as.integer(clusterY),
                      as.integer(unique(clusterY)),
                      as.integer(length(clusterY)),
                      as.integer(npred),
                      pred=as.double(matrix(0,nrow=npred,ncol=ntimeAll)),
                      predlow=as.double(matrix(0,nrow=npred,ncol=ntimeAll)),
                      predhigh=as.double(matrix(0,nrow=npred,ncol=ntimeAll)),
                      icproba=as.integer(ICproba),
                      as.integer(MC.sample),
                      as.integer(moving.window),
                      as.double(timeAll),
                      as.integer(s_cag_id),
                      as.double(s_cag)
      )
      
      out <- NULL
      out$call <- match.call()
      out$name.fit <- match.call()[[2]]
      out$npred <- npred
      out$window <- window
      out$predtimerec <- predtimerec
      out$moving.window <- moving.window
      if (moving.window){
        out$x.time <- timeAll
        out$t <- predTime
      }else{
        out$x.time <- timeAll - window
      } 
      out$group <- uni.cluster
      out$pred <- matrix(ans$pred,nrow=npred,ncol=ntimeAll)
      rownames(out$pred) <- paste("ind",out$group)
      colnames(out$pred) <- c("times",rep(" ",ntimeAll-1))
      
      out$icproba <- ICproba
      if (ICproba){
        out$predLow <- matrix(ans$predlow,nrow=npred,ncol=ntimeAll)
        out$predHigh <- matrix(ans$predhigh,nrow=npred,ncol=ntimeAll)
        rownames(out$predLow) <- paste("ind",out$group)
        colnames(out$predLow) <- c("times",rep(" ",ntimeAll-1))
        rownames(out$predHigh) <- paste("ind",out$group)
        colnames(out$predHigh) <- c("times",rep(" ",ntimeAll-1))
      }
      out$trivariate <- TRUE
    }else if(class(fit)=="trivPenalNL"){             
      predtimerec <- predtimerec[order(unique(cluster)),]	
      if(!is.matrix(predtimerec)) predtimerec <- matrix(predtimerec, nrow = npred)
      box_cox <- c(0,1)
      box_cox[1] <- ifelse(fit$BoxCox == TRUE, 1, 0)
      if(!is.null(fit$BoxCox_parameter))box_cox[2] <- fit$BoxCox_parameter
      
      
      ans <- .Fortran(C_predicttrinl,
                      as.integer(np),
                      as.double(b),
                      as.integer(nz),
                      as.integer(nva1),
                      as.integer(nva2),
                      as.integer(nva3),
                      as.integer(nva4),
                      as.integer(fit$ne_re),
                      as.integer(fit$random.which),
                      as.double(box_cox),
                      as.integer(fit$netar),
                      as.integer(fit$netadc),
                      as.integer(link),
                      as.integer(nst),
                      as.integer(typeof),
                      as.double(zi),
                      as.double(HIH),
                      as.integer(ntimeAll),
                      as.integer(npred),
                      as.double(predTime),
                      as.double(window),
                      as.double(predtimerec),
                      as.integer(nrec),
                      as.integer(fit$max_rep),
                      as.double(Y),
                      as.double(matzy),
                      as.double(as.matrix(vaxpred)),
                      as.double(as.matrix(X)),
                      as.double(cbind(varKG,varKD)),
                      as.integer(clusterY),
                      as.integer(unique(clusterY)),
                      as.integer(length(clusterY)),
                      as.integer(npred),
                      pred=as.double(matrix(0,nrow=npred,ncol=ntimeAll)),
                      predlow=as.double(matrix(0,nrow=npred,ncol=ntimeAll)),
                      predhigh=as.double(matrix(0,nrow=npred,ncol=ntimeAll)),
                      icproba=as.integer(ICproba),
                      as.integer(MC.sample),
                      as.integer(moving.window),
                      as.double(timeAll),
                      as.integer(s_cag_id),
                      as.double(s_cag)
                      
      )
      
      out <- NULL
      out$call <- match.call()
      out$name.fit <- match.call()[[2]]
      out$npred <- npred
      out$window <- window
      out$predtimerec <- predtimerec
      out$moving.window <- moving.window
      if (moving.window){
        out$x.time <- timeAll
        out$t <- predTime
      }else{
        out$x.time <- timeAll - window
      } 
      out$group <- uni.cluster
      #	  print(npred)
      out$pred <- matrix(ans$pred,nrow=npred,ncol=ntimeAll)
      rownames(out$pred) <- paste("ind",out$group)
      colnames(out$pred) <- c("times",rep(" ",ntimeAll-1))
      
      out$icproba <- ICproba
      if (ICproba){
        out$predLow <- matrix(ans$predlow,nrow=npred,ncol=ntimeAll)
        out$predHigh <- matrix(ans$predhigh,nrow=npred,ncol=ntimeAll)
        rownames(out$predLow) <- paste("ind",out$group)
        colnames(out$predLow) <- c("times",rep(" ",ntimeAll-1))
        rownames(out$predHigh) <- paste("ind",out$group)
        colnames(out$predHigh) <- c("times",rep(" ",ntimeAll-1))
      }
      out$trivariate <- TRUE
    }
    cat("Predictions done for",npred,"subjects and",ntimeAll,"times \n")
    class(out) <- "predLongi"
    
    ####*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*####
    ####-*-*-*-*-Prediction pour un modele Shared -*-*-*-*-####
    ####*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*####
  }else if(class(fit)=="frailtyPenal"){	
    cat("\n")
    cat("Calculating the probabilities ... \n")		
    
    if (event == 'Recurrent'){	
      taille = 0
      listPrec <- NULL  
      if(!is.matrix(predtimerec)) predtimerec <- matrix(predtimerec, nrow = npred)
      for (k in 1:nrow(predtimerec)){
        tPrec <- which(predtimerec[k,] < predTime)   
        if (length(tPrec) == 0) tPrec <- taille + 1 
        tPrec <- taille + length(tPrec)    
        
        rowTimes <- predtimerec[k,][which(!is.na(predtimerec[k,]))]
        if (length(rowTimes)==0) rowTimes <- 1
        taille = length(rowTimes)+taille
        listPrec <- c(listPrec,tPrec)                				
      }			
      # X <- as.matrix(aggregate(X,by=list(cluster),FUN=function(x){x[1]})[,-1])
      X <- X[listPrec,]
      if (length(unique(cluster)) > 1) X <- X[order(unique(cluster)),]
    }
    expBX <- exp(X %*% fit$coef)	
    
    # nombre de predictions a faire pour chaque individu
    if (moving.window){ # 2 facons differentes de faire des predictions, soit h evolue, soit t evolue
      sequence2 <- t+window #seq(predTime+window,predMax,by=window)
      sequence <- rep(predTime,times=length(sequence2))
    }else{
      sequence <- t #seq(predTime,predMax,length=50)
      sequence2 <- t+window #sequence+window
    }
    predMat <- NULL
    
    if (fit$Frailty){
      ###############################
      ###   Prediction marginale  ###
      ###############################
      if (!conditional){	
        #======================================================#
        #************ML: Prediction pour evenement recurrent****
        #======================================================#
        if(event == 'Recurrent'){				
          mat.survival.X <- NULL
          mat.survival.X.horizon <- NULL
          mat.survival.LastRec <- NULL
          if(!is.matrix(predtimerec)) predtimerec <- matrix(predtimerec, nrow = 1)
          npred0 <- nrow(predtimerec) #nb subjects
          nbrec <- rep(0,npred0)
          
          for (k in 1:npred0){					
            # vect.survival.X <- sapply(sequence,FUN=survival,ObjFrailty=fit)**expBX[k]
            vect.survival.X <- sapply(sequence,FUN=survival,ObjFrailty=fit)
            # vect.survival.X.horizon <- sapply(sequence2,FUN=survival,ObjFrailty=fit)**expBX[k]
            vect.survival.X.horizon <- sapply(sequence2,FUN=survival,ObjFrailty=fit)
            mat.survival.X <- rbind(mat.survival.X,vect.survival.X)
            mat.survival.X.horizon <- rbind(mat.survival.X.horizon,vect.survival.X.horizon)					
            
            recurr <- which(!is.na(predtimerec[k,which(predtimerec[k,] <= predTime)]))
            nbrec[k] <- length(recurr)						
            
            if(length(recurr) == 0) LastRec <- 0
            else if (length(recurr) > 1) LastRec <- predtimerec[k, recurr[-c(1:length(recurr)-1)]]
            else LastRec <- predtimerec[k,recurr]
            
            if((length(recurr) == 1)&&(LastRec == 0)) nbrec[k] <- 0
            
            # vect.survival.LastRec <- survival(LastRec,fit)**expBX[k]
            if (LastRec == 0) vect.survival.LastRec <- 1
            else vect.survival.LastRec <- survival(LastRec,fit)
            mat.survival.LastRec <- rbind(mat.survival.LastRec,vect.survival.LastRec)
          }					
          ######### Distribution gamma #########			
          if(fit$logNormal==0) variance <- fit$theta
          ######### Distribution LogNormale #########
          else variance <- fit$sigma2	
          
          nbrec <- nbrec[order(uni.cluster)]
          mat.survival.LastRec[,1] <- mat.survival.LastRec[order(uni.cluster)]	
          
          ans <- .Fortran(C_predict_recurr_sha,
                          as.integer(fit$logNormal),
                          as.integer(npred0),
                          as.double(mat.survival.X),
                          as.double(mat.survival.X.horizon),
                          as.double(mat.survival.LastRec),
                          as.double(expBX),
                          as.double(variance),
                          pred=as.double(matrix(0,nrow=npred0,ncol=ntimeAll)),
                          as.integer(nbrec),
                          as.integer(ntimeAll),
                          as.integer(0),						
                          as.integer(MC.sample),
                          as.double(rep(0,MC.sample)),
                          as.double(matrix(0,nrow=npred0*MC.sample,ncol=ntimeAll)),
                          as.double(matrix(0,nrow=npred0*MC.sample,ncol=ntimeAll)),
                          as.double(rep(0,nrow=npred0*MC.sample)),
                          as.double(rep(0,nrow=npred0*MC.sample)),
                          predlow1=as.double(matrix(0,nrow=npred0,ncol=ntimeAll)),
                          predhigh1=as.double(matrix(0,nrow=npred0,ncol=ntimeAll))
          )
          
          predMat <- matrix(ans$pred,nrow=npred0,ncol=ntimeAll)
          
          #=================================================#
          #*         Prediction pour donnees groupees       *
          #=================================================#
        }else{				    
          ######### Distribution gamma #########			
          if(fit$logNormal==0){	
            for (k in 1:nrow(data)){
              vect.survival.X <- sapply(sequence,FUN=survival,ObjFrailty=fit)**expBX[k]
              vect.survival.X.horizon <- sapply(sequence2,FUN=survival,ObjFrailty=fit)**expBX[k]
              pred <- 1-((1+fit$theta*(-log(vect.survival.X)))/(1+fit$theta*(-log(vect.survival.X.horizon))))**(1/fit$theta)
              predMat <- rbind(predMat,pred)
            }			
          }else{
            ######### AK: Distribution LogNormale #########				
            mat.survival.X <- NULL
            mat.survival.X.horizon <- NULL
            for (k in 1:nrow(data)){						
              #vect.survival.X <- sapply(sequence,FUN=survival,ObjFrailty=fit)**expBX[k] #ML : Survie doit tenir compte du terme de fragilite 
              vect.survival.X <- sapply(sequence,FUN=survival,ObjFrailty=fit)
              #vect.survival.X.horizon <- sapply(sequence2,FUN=survival,ObjFrailty=fit)**expBX[k]	#ML
              vect.survival.X.horizon <- sapply(sequence2,FUN=survival,ObjFrailty=fit)
              mat.survival.X <- rbind(mat.survival.X,vect.survival.X)
              mat.survival.X.horizon <- rbind(mat.survival.X.horizon,vect.survival.X.horizon)							
            }								
            ans <- .Fortran(C_predict_logn_sha,
                            as.integer(nrow(data)),
                            as.double(mat.survival.X),
                            as.double(mat.survival.X.horizon),
                            as.double(expBX),
                            as.double(fit$sigma2),
                            pred=as.double(matrix(0,nrow=nrow(data),ncol=ntimeAll)),
                            as.integer(0),
                            as.integer(ntimeAll),
                            as.integer(MC.sample),
                            as.double(rep(0,MC.sample)),
                            as.double(matrix(0,nrow=nrow(data)*MC.sample,ncol=ntimeAll)),
                            as.double(matrix(0,nrow=nrow(data)*MC.sample,ncol=ntimeAll)),
                            as.double(matrix(0,nrow=nrow(data),ncol=ntimeAll)),
                            predlow1=as.double(matrix(0,nrow=nrow(data),ncol=ntimeAll)),
                            predhigh1=as.double(matrix(0,nrow=nrow(data),ncol=ntimeAll))
            )					
            predMat <- matrix(ans$pred,nrow=nrow(data),ncol=ntimeAll)
          }
        }				
        ###############################
        ###Prediction conditionnelle###
        ###############################
      }else{
        if (event =="Recurrent"){
          if (is.null(nrow(X))){	
            if (!(unique(cluster) %in% uni.clusterfit)) stop("Are you sure that the group ",unique(cluster)," is present in your cluster variable ?")			
            vect.survival.X <- sapply(sequence,FUN=survival,ObjFrailty=fit)**expBX
            vect.survival.X.horizon <- sapply(sequence2,FUN=survival,ObjFrailty=fit)**expBX
            pred <- 1-(vect.survival.X.horizon/vect.survival.X)**fit$frailty.pred[as.integer(uni.clusterfit)==unique(as.integer(cluster))]						
            predMat <- rbind(predMat,pred)
          }else{
            uni.cluster <- uni.cluster[order(uni.cluster)]
            for (k in 1:nrow(X)){					
              if (!(uni.cluster[k] %in% uni.clusterfit)) stop("Are you sure that the group ",uni.cluster[k]," is present in your cluster variable ?")
              vect.survival.X <- sapply(sequence,FUN=survival,ObjFrailty=fit)**expBX[k]
              vect.survival.X.horizon <- sapply(sequence2,FUN=survival,ObjFrailty=fit)**expBX[k]
              pred <- 1-(vect.survival.X.horizon/vect.survival.X)**fit$frailty.pred[as.integer(uni.clusterfit)==unique(as.integer(cluster))[k]]
              predMat <- rbind(predMat,pred)
            }
          }
        }else{
          cluster <- as.integer(as.vector(cluster))														
          if (any(!(uni.cluster %in% uni.clusterfit))) stop("Are you sure that the group ", uni.cluster, "is present in your cluster variable ?")					
          for (i in 1:nrow(data)){										
            vect.survival.X <- sapply(sequence,FUN=survival,ObjFrailty=fit)**expBX[i]
            vect.survival.X.horizon <- sapply(sequence2,FUN=survival,ObjFrailty=fit)**expBX[i]
            if (fit$logNormal==0){ # Gamma distribution
              pred <- 1-(vect.survival.X.horizon/vect.survival.X)**fit$frailty.pred[as.integer(uni.clusterfit)==as.integer(cluster)[i]]
            }else{ #AK: Normal distribution
              pred <- 1-(vect.survival.X.horizon/vect.survival.X)**exp(fit$frailty.pred[as.integer(uni.clusterfit)==as.integer(cluster)[i]])
            }
            predMat <- rbind(predMat,pred)
          }
        }
      }
      #*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*#
      #-*-*-*-*-Pour un modele de Cox-*-*-*-*-#
      #*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*#		
    }else{
      for (k in 1:nrow(data)){
        vect.survival.X <- sapply(sequence,FUN=survival,ObjFrailty=fit)**expBX[k]
        vect.survival.X.horizon <- sapply(sequence2,FUN=survival,ObjFrailty=fit)**expBX[k]				
        pred <- 1-(vect.survival.X.horizon/vect.survival.X)
        predMat <- rbind(predMat,pred)
      }
    }		
    # -------------------------------------------------------- #
    # calcul des bornes de confiances (methode de Monte Carlo) #
    # -------------------------------------------------------- #
    if (ICproba){			
      balea <- mvrnorm(MC.sample,fit$b,fit$varHtotal)
      if (fit$Frailty){ #AK: For Gamma we have variance theta and for Normal we have variance sigma2
        if(fit$logNormal==0)theta.mc <- balea[,fit$np-fit$nvar]^2
        if(fit$logNormal==1)sigma2.mc <- balea[,fit$np-fit$nvar]^2
      }
      aleaCoef <- balea[,(fit$np-fit$nvar+1):(fit$np)]		
      expBX.mc <- exp(X %*% t(aleaCoef))
      
      # recuperation parametres de la fonction de risque/survie (splines,piecewise,weibull)
      if (fit$typeof == 0){ #Splines
        para.mc <- balea[,1:(fit$n.knots+2)]^2
        if(fit$n.strat == 2) para.mc2 <- balea[,(fit$n.knots+3):(2*(fit$n.knots+2))]^2
        else para.mc2 <- matrix(0,nrow=MC.sample,ncol=fit$n.knots+2)				
      }else if (fit$typeof == 1){ #Piecewise
        para.mc <- balea[,1:(fit$nbintervR)] # attention de ne pas elever au carre
        if(fit$n.strat == 2) para.mc2 <- balea[,(fit$nbintervR+1):(2*fit$nbintervR)]
        else para.mc2 <- matrix(0,nrow=MC.sample,ncol=fit$nbintervR)				
      }else{ #Weibull
        para.mc <- balea[,1:2]^2
        if(fit$n.strat == 2) para.mc2 <- balea[,2:4]^2
        else para.mc2 <- matrix(0,nrow=MC.sample,ncol=2)
      }
      
      survival.mc <- function(t,ObjFrailty,para1,para2){ # dans les trois cas para1 et para2 seront traites differemment
        if (ObjFrailty$typeof == 0){ # splines
          nz <- ObjFrailty$n.knots
          zi <- ObjFrailty$zi
          res <- NULL
          nst <- ObjFrailty$n.strat
          
          out <- .Fortran(C_survival_frailty,
                          as.double(t),
                          as.double(para1),
                          as.double(para2),
                          as.integer(nz+2),
                          as.double(zi),
                          survival=as.double(c(0,0)),
                          lam=as.double(c(0,0)),
                          as.integer(nst)#lam ajoute suite aux modif de survival
          )
          
          if(ObjFrailty$n.strat == 2){
            res <- c(res,out$survival)
          }else{
            res <- c(res,out$survival[1])
          }
          return(res)
        }
        
        if (ObjFrailty$typeof == 1){ # piecewise
          res <- NULL
          if (ObjFrailty$n.strat == 2) b <- c(para1,para2)
          else b <- para1
          time <- ObjFrailty$time
          
          out <- .Fortran(C_survival_cpm,
                          as.double(t),
                          as.double(b),
                          as.integer(ObjFrailty$n.strat),
                          as.integer(ObjFrailty$nbintervR),
                          as.double(time),
                          survival=as.double(c(0,0))
          )
          
          if(ObjFrailty$n.strat == 2){
            res <- c(res,out$survival)
          }else{
            res <- c(res,out$survival[1])
          }
          return(res)
        }
        
        if (ObjFrailty$typeof == 2){ # weibull
          res <- NULL
          sh1 <- para1[1]
          sc1 <- para1[2]
          res <- c(res,exp(-(t/sc1)^sh1))
          if(ObjFrailty$n.strat == 2){
            sh1 <- para2[1]
            sc1 <- para2[2]
            res <- c(res,exp(-(t/sc1)^sh1))
          }
          return(res)
        }				
      }# end of survival.mc function
      
      predMatLow <- NULL
      predMatHigh <- NULL
      frailty.mc <- NULL			
      #*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*#
      #-*-*-*-*-Pour un modele Shared-*-*-*-*-#
      #*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*#
      if (fit$Frailty){
        ###############################
        ###   Prediction marginale  ###
        ###############################
        if (!conditional){   
          #======================================================#
          #************ML: Prediction pour evenement recurrent****
          #======================================================#
          if(event == 'Recurrent'){	
            mat.survival.X.mc <- NULL
            mat.survival.X.horizon.mc <- NULL
            mat.survival.LastRec.mc <- NULL	
            
            for(i in 1:MC.sample){
              mat.survival.X.samp <- NULL
              mat.survival.X.horizon.samp <- NULL
              mat.survival.LastRec.samp <- NULL							
              for(k in 1:npred0){
                # vect.survival.X.samp <- sapply(sequence,FUN=survival.mc,ObjFrailty=fit,para1=para.mc[i,],para2=para.mc2[i,])**expBX.mc[k,i]
                vect.survival.X.samp <- sapply(sequence,FUN=survival.mc,ObjFrailty=fit,para1=para.mc[i,],para2=para.mc2[i,])
                # vect.survival.X.horizon.samp <- sapply(sequence2,FUN=survival.mc,ObjFrailty=fit,para1=para.mc[i,],para2=para.mc2[i,])**expBX.mc[k,i]
                vect.survival.X.horizon.samp <- sapply(sequence2,FUN=survival.mc,ObjFrailty=fit,para1=para.mc[i,],para2=para.mc2[i,])
                mat.survival.X.samp <- rbind(mat.survival.X.samp,vect.survival.X.samp)
                mat.survival.X.horizon.samp <- rbind(mat.survival.X.horizon.samp,vect.survival.X.horizon.samp)								
                recurr <- which(!is.na(predtimerec[k,which(predtimerec[k,] <= predTime)]))
                nbrec[k] <- length(recurr)
                
                if(length(recurr) == 0) LastRec <- 0
                else if (length(recurr) > 1) LastRec <- predtimerec[k,recurr[-1]]
                else LastRec <- predtimerec[k,recurr]	
                # vect.survival.LastRec.samp <- sapply(LastRec,FUN=survival.mc,ObjFrailty=fit,para1=para.mc[i,],para2=para.mc2[i,])**expBX.mc[k,i]
                vect.survival.LastRec.samp <- sapply(LastRec,FUN=survival.mc,ObjFrailty=fit,para1=para.mc[i,],para2=para.mc2[i,])
                mat.survival.LastRec.samp <- cbind(mat.survival.LastRec.samp,vect.survival.LastRec.samp)							
              }									
              mat.survival.X.mc <- rbind(mat.survival.X.mc,mat.survival.X.samp)
              mat.survival.X.horizon.mc <- rbind(mat.survival.X.horizon.mc,mat.survival.X.horizon.samp)
              mat.survival.LastRec.mc <- rbind(mat.survival.LastRec.mc,mat.survival.LastRec.samp)
            }
            
            ######### Distribution gamma #########				
            if(fit$logNormal==0){ 
              variance <- fit$theta
              variance.mc <- theta.mc
              ######### Distribution LogNormale #########
            }else{
              variance <- fit$sigma2
              variance.mc <- sigma2.mc
            }	
            
            ans <- .Fortran(C_predict_recurr_sha,
                            as.integer(fit$logNormal),
                            as.integer(npred0),
                            as.double(mat.survival.X),
                            as.double(mat.survival.X.horizon),
                            as.double(mat.survival.LastRec),
                            as.double(expBX),
                            as.double(variance),
                            pred=as.double(matrix(0,nrow=npred0,ncol=ntimeAll)),
                            as.integer(nbrec),
                            as.integer(ntimeAll),
                            as.integer(1),						
                            as.integer(MC.sample),
                            as.double(variance.mc),
                            as.double(mat.survival.X.mc),
                            as.double(mat.survival.X.horizon.mc),
                            as.double(mat.survival.LastRec.mc),
                            as.double(expBX.mc),
                            predlow1=as.double(matrix(0,nrow=npred0,ncol=ntimeAll)),
                            predhigh1=as.double(matrix(0,nrow=npred0,ncol=ntimeAll))
            )
            
            predMatLow <- matrix(ans$predlow1,nrow=npred0,ncol=ntimeAll)
            predMatHigh <- matrix(ans$predhigh1,nrow=npred0,ncol=ntimeAll)
            
            #=================================================#
            #************Prediction pour evenement de deces****
            #=================================================#
          }else{
            ######### Distribution gamma #########
            if(fit$logNormal==0){
              for (k in 1:nrow(data)){
                realisations <- NULL
                for (i in 1:MC.sample){
                  vect.survival.X <- sapply(sequence,FUN=survival.mc,ObjFrailty=fit,para1=para.mc[i,],para2=para.mc2[i,])**expBX.mc[k,i]
                  vect.survival.X.horizon <- sapply(sequence2,FUN=survival.mc,ObjFrailty=fit,para1=para.mc[i,],para2=para.mc2[i,])**expBX.mc[k,i]
                  pred <- 1-((1+theta.mc[i]*(-log(vect.survival.X)))/(1+theta.mc[i]*(-log(vect.survival.X.horizon))))**(1/theta.mc[i])
                  realisations <- cbind(realisations,pred)
                }
                predMatLow <- rbind(predMatLow,apply(realisations,1,quantile,probs=0.025))
                predMatHigh <- rbind(predMatHigh,apply(realisations,1,quantile,probs=0.975))
              }
              
              ######### AK: Distribution LogNormale #########
            }else{
              mat.survival.X.mc <- NULL
              mat.survival.X.horizon.mc <- NULL
              
              for(i in 1:MC.sample){
                mat.survival.X.samp <- NULL
                mat.survival.X.horizon.samp <- NULL
                
                for(k in 1:nrow(data)){
                  # vect.survival.X.samp <- sapply(sequence,FUN=survival.mc,ObjFrailty=fit,para1=para.mc[i,],para2=para.mc2[i,])**expBX.mc[k,i]
                  vect.survival.X.samp <- sapply(sequence,FUN=survival.mc,ObjFrailty=fit,para1=para.mc[i,],para2=para.mc2[i,])
                  # vect.survival.X.horizon.samp <- sapply(sequence2,FUN=survival.mc,ObjFrailty=fit,para1=para.mc[i,],para2=para.mc2[i,])**expBX.mc[k,i]
                  vect.survival.X.horizon.samp <- sapply(sequence2,FUN=survival.mc,ObjFrailty=fit,para1=para.mc[i,],para2=para.mc2[i,])
                  mat.survival.X.samp <- rbind(mat.survival.X.samp,vect.survival.X.samp)
                  mat.survival.X.horizon.samp <- rbind(mat.survival.X.horizon.samp,vect.survival.X.horizon.samp)
                }
                
                mat.survival.X.mc <- rbind(mat.survival.X.mc,mat.survival.X.samp)
                mat.survival.X.horizon.mc <- rbind(mat.survival.X.horizon.mc,mat.survival.X.horizon.samp)
              }		
              
              ans <- .Fortran(C_predict_logn_sha,
                              as.integer(nrow(data)),
                              as.double(mat.survival.X),
                              as.double(mat.survival.X.horizon),
                              as.double(expBX),
                              as.double(fit$sigma2),
                              pred=as.double(matrix(0,nrow=nrow(data),ncol=ntimeAll)),
                              as.integer(1),
                              as.integer(ntimeAll),
                              as.integer(MC.sample),
                              as.double(sigma2.mc),
                              as.double(mat.survival.X.mc),
                              as.double(mat.survival.X.horizon.mc),
                              as.double(expBX.mc),
                              predlow1=as.double(matrix(0,nrow=nrow(data),ncol=ntimeAll)),
                              predhigh1=as.double(matrix(0,nrow=nrow(data),ncol=ntimeAll))
              )
              
              predMatLow <- matrix(ans$predlow1,nrow=nrow(data),ncol=ntimeAll)
              predMatHigh <- matrix(ans$predhigh1,nrow=nrow(data),ncol=ntimeAll)
            }
          }
          
          ###############################
          ###Prediction conditionnelle###
          ###############################
        }else{
          if (event == 'Recurrent') X3 <- model.matrix(newTerms2, dataset3)	
          else  X3 <- model.matrix(newTerms, dataset3)
          
          if (ncol(X3) > 1) X3 <- X3[, -1, drop = FALSE]
          expBX3 <- exp(X3 %*% fit$coef)
          
          cluster <- as.integer(cluster)
          # cluster <- as.integer(as.vector(cluster))
          
          if (event == "Recurrent") nbInd <- length(unique(cluster))
          else nbInd <- nrow(data)
          
          for (k in 1:nbInd){								
            realisations <- NULL
            frailty.mc <- NULL						
            mi <- fit$n.eventsbygrp[cluster[k]]
            
            # calcul de la somme des risques cumules juste pour le groupe defini			
            res1 <- sum((-log(sapply(tt1[which(clusterfit==cluster[k])],survival,ObjFrailty=fit))) %*% expBX3[which(clusterfit==cluster[k])])
            
            for (i in 1:MC.sample){
              vect.survival.X <- sapply(sequence,FUN=survival.mc,ObjFrailty=fit,para1=para.mc[i,],para2=para.mc2[i,])**expBX.mc[k,i]
              vect.survival.X.horizon <- sapply(sequence2,FUN=survival.mc,ObjFrailty=fit,para1=para.mc[i,],para2=para.mc2[i,])**expBX.mc[k,i]
              
              ######### Distribution gamma #########
              if(fit$logNormal==0){
                # if (k == 1) 
                frailty.mc <- c(frailty.mc,rgamma(1,shape=mi+1/theta.mc[i],scale=1/(res1+1/theta.mc[i])))
                pred <- 1-(vect.survival.X.horizon/vect.survival.X)**frailty.mc[i]
                ######### AK: Distribution LogNormale #########
              }else{
                # if (k==1){
                res<-.Fortran(C_frailpred_sha_nor_mc,
                              as.integer(fit$npar),
                              frail.out=as.double(0),
                              as.double(sigma2.mc[i]),
                              as.double(res1),
                              as.integer(mi)
                )
                frailty.mc[i] <- res$frail.out
                # }
                pred <- 1-(vect.survival.X.horizon/vect.survival.X)**exp(frailty.mc[i])
              }
              realisations <- cbind(realisations,pred)
            }
            predMatLow <- rbind(predMatLow,apply(realisations,1,quantile,probs=0.025))
            predMatHigh <- rbind(predMatHigh,apply(realisations,1,quantile,probs=0.975))
          }
        }
        
        #*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*#
        #-*-*-*-*-Pour un modele de Cox-*-*-*-*-#
        #*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*#
      }else{
        for (k in 1:nrow(data)){
          realisations <- NULL
          for (i in 1:MC.sample){
            vect.survival.X <- sapply(sequence,FUN=survival.mc,ObjFrailty=fit,para1=para.mc[i,],para2=para.mc2[i,])**expBX.mc[k,i]
            vect.survival.X.horizon <- sapply(sequence2,FUN=survival.mc,ObjFrailty=fit,para1=para.mc[i,],para2=para.mc2[i,])**expBX.mc[k,i]
            pred <- 1-(vect.survival.X.horizon/vect.survival.X)
            realisations <- rbind(realisations,pred)
          }
          predMatLow <- rbind(predMatLow,apply(realisations,1,quantile,probs=0.025))
          predMatHigh <- rbind(predMatHigh,apply(realisations,1,quantile,probs=0.975))
        }
      }
    } # Fin du calcul des bornes de confiance	
    
    out <- NULL
    out$call <- match.call()
    out$name.fit <- match.call()[[2]]
    
    if (event == 'Recurrent'){
      if(conditional) out$npred <- length(uni.cluster)
      else out$npred <- npred0
      out$event <- "Recurrent"
      out$predtimerec <- predtimerectmp
    }else{
      out$npred <- nrow(data)	
      out$event <- "Terminal"
    }
    
    out$moving.window <- moving.window
    if (moving.window){
      out$x.time <- sequence2
      out$t <- predTime
    }else{
      out$x.time <- sequence
    }
    out$pred <- predMat
    # colnames(out$pred) <- c("times",rep(" ",dim(out$pred)[2]-1))		
    colnames(out$pred) <- paste("time=", out$x.time)
    
    if (out$event == "Terminal") rownames(out$pred) <- paste("ind",1:out$npred)
    else{
      if (conditional) rownames(out$pred) <- paste("ind",unique(cluster))
      else rownames(out$pred) <- paste(nameGrup,unique(cluster)[order(unique(cluster))])
    }		
    out$icproba <- ICproba
    if (ICproba){
      out$predLow <- predMatLow
      out$predHigh <- predMatHigh
      # colnames(out$predLow) <- c("times",rep(" ",dim(out$predLow)[2]-1))
      colnames(out$predLow) <- paste("time=", out$x.time)
      
      rownames(out$predLow) <- paste("ind",1:out$npred)
      # colnames(out$predHigh) <- c("times",rep(" ",dim(out$predHigh)[2]-1))
      colnames(out$predHigh) <- paste("time=", out$x.time)
      rownames(out$predHigh) <- paste("ind",1:out$npred)
    }
    if (fit$Frailty) {
      if (conditional) out$type <- 'conditional'
      else out$type <- 'marginal'
    }
    out$window <- window
    if (conditional) out$group <- unique(cluster)
    cat("Predictions done for",out$npred,"subjects \n")
    class(out) <- "predFrailty"
    
    
    ####*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*####
    ####-*-*-*-*-Prediction pour un Joint Nested Model-*-*-*-*-*-####
    ####*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*####
    
  }else if (class(fit) == "jointNestedPenal"){
    if(!all(individual%in%subcluster))stop("The individual for whom you want to calculate predictions must be in the data for predictions")
    if(moving.window == FALSE){
      window <- rep(window, length(t))
      Time <- t
    }else{
      Time <- rep(t, length(window))
    }
    nst <- 2
    indID <- 2
    
    # predtimerec <- predtimerec[order(unique(cluster)),]
    taille = 0
    listPrec <- NULL  
    
    if(!is.matrix(predtimerec)) predtimerec <- matrix(predtimerec, nrow = 1)
    
    for (k in 1:npred){
      tPrec <- which(predtimerec[k,] < predTime)   
      if (length(tPrec) == 0) tPrec <- taille + 1 
      tPrec <- taille + length(tPrec)  
      
      rowTimes <- predtimerec[k,][which(!is.na(predtimerec[k,]))]
      if (length(rowTimes)==0) rowTimes <- 1
      taille = length(rowTimes)+taille
      listPrec <- c(listPrec,tPrec)                				
    }		
    predtimerec <- replace(predtimerec, is.na(predtimerec),0)
    
    vaxpred <- X[,-c(1,2), drop=FALSE]		
    
    m3 <- fit$call
    m2 <- match.call() # formule appelee pour prediction()
    
    m3$formula.LongitudinalData <- m3$data.Longi <- m3$random <- m3$id <- m3$link <- m3$left.censoring <- m3$n.knots <- m3$recurrentAG <- m3$cross.validation <- m3$kappa <- m3$maxit <- m3$hazard <- m3$nb.int <- m3$RandDist <- m3$betaorder <- m3$betaknots <- m3$init.B <- m3$LIMparam <- m3$LIMlogl <- m3$LIMderiv <- m3$print.times <- m3$init.Theta <- m3$init.Alpha <- m3$Alpha <- m3$init.Random <- m3$init.Eta <- m3$method.GH <- m3$intercept <- m3$n.nodes <- m3$jointGeneral <- m3$initialize <- m3$Ksi <- m3$init.Ksi <- m3$... <- NULL
    
    mPred <- m3
    m3$formula <- formula_fit
    m3$formula[[3]][[2]] <- fit$formula.terminalEvent[[2]]				
    m3$formula.terminalEvent <- NULL
    m3[[1]] <- as.name("model.frame")
    m3[[3]] <- as.name(m2$data)
    
    temp <- as.character(m3$formula[[2]])
    
    if (temp[1]=="Surv"){
      if (length(temp) == 4) m3$formula[[2]] <- as.name(temp[3])
      else if (length(temp) == 3) m3$formula[[2]] <- as.name(temp[2])
      else stop("Wrong Surv function")
    }else{ # SurvIC
      if (length(temp) == 4) m3$formula[[2]] <- as.name(temp[2])
      else if (length(temp) == 5) m3$formula[[2]] <- as.name(temp[3])
      else stop("Wrong SurvIC function")
    }
    datasetdc <- eval(m3, sys.parent())		
    
    class(m3$formula) <- "formula"
    special2 <- c("strata", "timedep")				
    Terms2 <- terms(m3$formula, special2, data = data)		
    X2 <- model.matrix(Terms2, datasetdc)			
    
    icdc <- X2[,ncol(X2)]	
    if(fit$AG == TRUE){
      tt1dc <- matrix(rep(aggregate(Time_fam, by = list(subcluster), function(x) x[length(x)])[-1][,1], length(Time)), ncol = length(Time))
    }else{
      tt1dc <- matrix(rep(aggregate(Time_fam, by = list(subcluster), sum)[-1][,1], length(Time)), ncol = length(Time))
    }
    icdcT <- matrix(rep(aggregate(icdc, by = list(subcluster), FUN = function(x){x[length(x)]})[-1][[1]],length(Time)),ncol = length(Time))
    
    if(length(Time_fam) == 1)tt1dc <- matrix(rep(Time_fam, length(Time)),  ncol = length(Time))
    if(length(Time)>1){
      for(k in 1:length(Time)){
        for (i in 1:dim(icdcT)[1]){
          if (tt1dc[i,k] > Time[k]) icdcT[i,k] <- 0
        }
        tt1dc[,k] <- sapply(tt1dc[,k], FUN = function(x,tpred=Time[k]){
          if(x > tpred){tpred} 
          else{x}
        }, simplify=TRUE)
      }
      tt1dc[icdcT[,k]==0,k] <- Time[k]
    }else{
      for (i in 1:dim(icdcT)[1]){
        if (tt1dc[i,1] > t) icdcT[i,1] <- 0
      }
      tt1dc[,1] <- sapply(tt1dc[,1], FUN = function(x,tpred=Time){
        if(x > tpred){tpred} 
        else{x}
      }, simplify=TRUE)
      tt1dc[icdcT[,1]==0,1] <- Time
    }
    
    #	tt1dc <- sapply(tt1dc, FUN = function(x,tpred=t){
    #                                 if(x > tpred){tpred} 
    #                                else{x}
    #                               }, simplify=TRUE)
    
    if (ncol(X2) > 1) X2 <- X2[, -1, drop = FALSE]	
    X2 <- X2[, -(ncol(X2)), drop = FALSE]		
    
    vaxdcpred <- aggregate(X2, by = list(subcluster), FUN = function(x){x[1]})[-1] 
    # X2 <- X2[order(cluster),]
    # vaxdcpred <- X2[listPrec,]	
    
    tt1 <- Time_fam		
    indiv <- subcluster[1]
    if(fit$AG == FALSE && length(subcluster)>1){
      for (i in 2:length(tt1)){
        if (subcluster[i] == indiv){
          tt1[i] <- tt1[i-1]+tt1[i]					
        }else{
          indiv = subcluster[i]
        }
      }	
    }
    
    
    icT <- matrix(rep(ic, length(Time)), ncol = length(Time))	
    
    if(length(Time) > 1){
      for(k in 1:length(Time)){
        for (i in 1:dim(icT)[1]){
          if (tt1[i] > Time[k]) icT[i,k] <- 0
        }
      }
    }else{
      for (i in 1:dim(icT)[1]){
        if (tt1[i] > Time) icT[i,1] <- 0
      }
    }
    
    tt1T <- matrix(rep(Time_fam, length(Time)), ncol = length(Time)) 
    indiv <- subcluster[1]
    cpt <- rep(1, length(Time))
    nrecT <- matrix(rep(0, length(unique(subcluster))*length(Time)), ncol = length(Time))
    if(length(subcluster)>1){
      if(length(Time) >1){
        for(k in 1:length(Time)){
          for (i in 1: dim(tt1T)[1]){
            if (indiv == subcluster[i]){
              if(i == 1){
                if(tt1[i] > Time[k]) tt1T[i,k] <- Time[k]
              }else{
                if(tt1[i] > Time[k]){
                  if ((icT[i-1,k] == 1)&&(icT[i,k] == 0)) tt1T[i,k] <- Time[k]-tt1[i-1] #AK: on n'a pas besoin car dans la formule les survies se simplifient (voir l'article Audrey, appendix)
                  else{
                    tt1T[i,k] <- 0
                    vaxpred[i,] <- 0
                  }
                } 
              }
              cpt[k] <- cpt[k]+1
            }else{
              if(tt1[i] > Time[k]) tt1T[i,k] <- Time[k]
              indiv <- subcluster[i]
            }
          }
          nrecT[,k] <- unlist(aggregate(icT[,k], by = list(subcluster), FUN = sum)[,-1])	
        }
      }else{
        for (i in 1: dim(tt1T)[1]){
          if (indiv == subcluster[i]){
            if(i == 1){
              if(tt1[i] > Time) tt1T[i,1] <- t 
            }else{
              if(tt1[i] > Time){
                if ((icT[i-1,1] == 1)&&(icT[i,1] == 0)) tt1T[i,1] <- Time-tt1[i-1]  
                else{
                  tt1T[i,1] <- 0
                  vaxpred[i,] <- 0
                }
              } 
            }
            cpt <- cpt+1
          }else{
            if(tt1[i] >Time) tt1T[i,1] <- t
            indiv <- subcluster[i]
          }
        }
        nrecT[,1] <- unlist(aggregate(icT[,1], by = list(subcluster), FUN = sum)[,-1])	
      }
    }
    out <- NULL	
    out$pred <- NULL
    out$predLow <- NULL
    out$predHigh <- NULL
    VecIND <- NULL
    
    for (i in individual){	
      
      indiceID <- which(as.integer(names(nrecList))==i)
      ans <- .Fortran(C_predictfam,
                      as.integer(fit$npar),
                      as.double(fit$b),
                      as.integer(fit$n.knots.temp),
                      as.integer(fit$nvarRec),
                      as.integer(fit$nvarEnd),
                      as.integer(nst),
                      as.integer(fit$typeof),
                      as.double(fit$zi),
                      as.double(fit$varHIHtotal),
                      as.integer(indiceID), 
                      as.double(tt1T), 
                      as.double(unlist(tt1dc)), 
                      as.integer(unlist(icdcT)),
                      as.integer(ntimeAll),
                      as.integer(nrow(data)),			
                      as.integer(npred), 
                      as.double(window), 
                      as.integer(max(nrecList)), 
                      as.integer(nrecList), 
                      as.integer(nrecT), 
                      as.double(as.matrix(vaxpred)),
                      as.double(as.matrix(vaxdcpred)), 
                      as.integer(ICproba), 
                      as.integer(MC.sample),	
                      predAll=as.double(matrix(0,nrow=1,ncol=ntimeAll)), 
                      predAlllow=as.double(matrix(0,nrow=1,ncol=ntimeAll)), 
                      predAllhigh=as.double(matrix(0,nrow=1,ncol=ntimeAll)), 
                      as.double(0), 
                      as.double(0),
                      pred = as.double(matrix(0,nrow=1,ncol=ntimeAll)),
                      as.integer(c(indic.alpha, indic.ksi))
      )
      # PACKAGE = "frailtypack" )#31 arguments
      
      out$pred <- rbind(out$pred,ans$predAll)
      out$predLow <- rbind(out$predLow,ans$predAlllow)
      out$predHigh <- rbind(out$predHigh, ans$predAllhigh)
      
      VecIND <- c(VecIND, indiceID)
      
    }
    
    out$call <- match.call()
    out$name.fit <- match.call()[[2]]
    out$npred <- length(individual)
    out$window <- window
    out$moving.window <- moving.window
    
    out$nrecList <- nrecList
    out$nrecT <- nrecT
    
    if (moving.window){
      out$x.time <- timeAll
      out$t <- predTime
    }else{
      out$x.time <- timeAll - window
    } 
    
    colnames(out$pred) <- paste("time=", out$x.time)
    rownames(out$pred) <- paste("ind", individual)
    out$icproba <- ICproba
    
    if (ICproba){
      colnames(out$predLow) <- paste("time=", out$x.time)
      rownames(out$predLow) <- paste("ind",individual)
      colnames(out$predHigh) <- paste("time=", out$x.time)
      rownames(out$predHigh) <- paste("ind",individual)
    }
    
    
    out$predtimerec <- rbind(NULL, predtimerec[VecIND,])		
    out$group <- individual		
    cat("Predictions done for",out$npred,"subjects and",ntimeAll,"times \n")
    class(out) <- "predJointNested"		
  }	
  out
} 
