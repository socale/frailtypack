#' Breast Cosmesis Data
#' 
#' The often used data set for interval-censored data, described and given in
#' full in Finkelstein and Wolfe (1985). It involves 94 breast cancer patients
#' who were randomized to either radiation therapy with chemotherapy or
#' radiation therapy alone. The outcome is time until the onset of breast
#' retraction which is interval-censored between the last clinic visit before
#' the event was observed and the first visit when the event was observed.
#' Patients without breast retraction were right-censored.
#' 
#' 
#' @name bcos
#' @docType data
#' @usage data(bcos)
#' @format A data frame with 94 observations and 3 variables: \describe{
#' \item{left}{left end point of the breast retraction interval}
#' \item{right}{right end point of the breast retraction interval}
#' \item{treatment}{type of treatment received} }
#' @source Finkelstein, D.M. and Wolfe, R.A. (1985). A semiparametric model for
#' regression analysis of interval-censored failure time data.
#' \emph{Biometrics} \bold{41}, 731-740.
#' @keywords datasets
NULL





#' Follow-up of metastatic colorectal cancer patients: times of new lesions
#' appearance and death
#' 
#' Randomly chosen 150 patients from the follow-up of the FFCD 2000-05
#' multicenter phase III clinical trial originally including 410 patients with
#' metastatic colorectal cancer randomized into two therapeutic strategies:
#' combination and sequential. The dataset contains times of observed
#' appearances of new lesions censored by a terminal event (death or
#' right-censoring) with baseline characteristics (treatment arm, age, WHO
#' performance status and previous resection).
#' 
#' 
#' @name colorectal
#' @docType data
#' @usage data(colorectal)
#' @format This data frame contains the following columns: \describe{
#' \item{id}{identification of each subject. Repeated for each recurrence}
#' \item{time0}{start of interval (0 or previous recurrence time)}
#' \item{time1}{recurrence or censoring time} \item{new.lesions}{Appearance of
#' new lesions status. 0: censsored or no event, 1: new lesions}
#' \item{treatment}{To which treatment arm a patient was allocated? 1:
#' sequential (S); 2: combination (C)} \item{age}{Age at baseline: 1: <50
#' years, 2: 50-69 years, 3: >69 years} \item{who.PS}{WHO performance status at
#' baseline: 1: status 0, 2: status 1, 3: status 2}
#' \item{prev.resection}{Previous resection of the primate tumor?  0: No, 1:
#' Yes} \item{state}{death indicator. 0: alive, 1: dead}
#' \item{gap.time}{interocurrence time or censoring time} }
#' @note We thank the Federation Francophone de Cancerologie Digestive and
#' Gustave Roussy for sharing the data of the FFCD 2000-05 trial supported by
#' an unrestricted Grant from Sanofi.
#' @references M. Ducreux, D. Malka, J. Mendiboure, P.-L. Etienne, P. Texereau,
#' D. Auby, P. Rougier, M. Gasmi, M. Castaing, M. Abbas, P. Michel, D. Gargot,
#' A. Azzedine, C. Lombard- Bohas, P. Geoffroy, B. Denis, J.-P., Pignon,
#' L.,Bedenne, and O.  Bouche (2011). Sequential versus combination
#' chemotherapy for the treatment of advanced colorectal cancer (FFCD 2000-05):
#' an open-label, randomised, phase 3 trial.  \emph{The Lancet Oncology}
#' \bold{12}, 1032-44.
#' @keywords datasets
NULL





#' Follow-up of metastatic colorectal cancer patients : longitudinal
#' measurements of tumor size
#' 
#' Randomly chosen 150 patients from the follow-up of the FFCD 2000-05
#' multicenter phase III clinical trial originally including 410 patients with
#' metastatic colorectal cancer randomized into two therapeutic strategies:
#' combination and sequential. The dataset contains measurements of tumor size
#' (left-censored sums of the longest diameters of target lesions; transformed
#' using Box-Cox) with baseline characteristics(treatment arm, age, WHO
#' performance status and previous resection).
#' 
#' 
#' @name colorectalLongi
#' @docType data
#' @usage data(colorectalLongi)
#' @format This data frame contains the following columns: \describe{
#' \item{id}{identification of each subject. Repeated for each recurrence}
#' \item{year}{time of visit counted in years from baseline}
#' \item{tumor.size}{Individual longitudinal measurement of transformed
#' (Box-Cox with parameter 0.3) sums of the longest diameters, left-censored
#' due to a detection limit (threshold \eqn{s=-3.33}). } \item{treatment}{To
#' which treatment arm a patient was allocated? 1: sequential (S); 2:
#' combination (C)} \item{age}{Age at baseline: 1: <50 years, 2: 50-69 years,
#' 3: >69 years} \item{who.PS}{WHO performance status at baseline: 1: status 0,
#' 2: status 1, 3: status 2} \item{prev.resection}{Previous resection of the
#' primate tumor?  0: No, 1: Yes} }
#' @note We thank the Federation Francophone de Cancerologie Digestive and
#' Gustave Roussy for sharing the data of the FFCD 2000-05 trial supported by
#' an unrestricted Grant from Sanofi.
#' @references Ducreux, M., Malka, D., Mendiboure, J., Etienne, P.-L.,
#' Texereau, P., Auby, D., Rougier, P., Gasmi, M., Castaing, M., Abbas, M.,
#' Michel, P., Gargot, D., Azzedine, A., Lombard- Bohas, C., Geoffroy, P.,
#' Denis, B., Pignon, J.-P., Bedenne, L., and Bouche, O. (2011). Sequential
#' versus combination chemotherapy for the treatment of advanced colorectal
#' cancer (FFCD 2000-05): an open-label, randomised, phase 3 trial.  \emph{The
#' Lancet Oncology} \bold{12}, 1032-44.
#' @keywords datasets
NULL





#' Simulated data as a gathering of clinical trials databases
#' 
#' This contains simulated samples of 100 clusters with 100 subjects in each
#' cluster, like a gathering of clinical trials databases. Two correlated
#' centred gaussian random effects are generated with the same variance fixed
#' at 0.3 and the covariance at -0.2. The regression coefficient \eqn{\beta} is
#' fixed at -0.11. The percentage of right-censored data is around 30 percent
#' which are generated from a uniform distribution on [1,150]. The baseline
#' hazard function is considered as a simple Weibull.
#' 
#' 
#' @name dataAdditive
#' @docType data
#' @usage data(dataAdditive)
#' @format This data frame contains the following columns: \describe{
#' \item{group}{identification variable} \item{t1}{start of interval (=0,
#' because left-truncated data are not allowed)} \item{t2}{end of interval
#' (death or censoring time)} \item{event}{censoring status (0:alive, 1:death,
#' as acensoring indicator} \item{var1}{dichotomous covariate (=0 or 1,as a
#' treatment variable)} \item{var2}{dichotomous covariate (=0 or 1,as a
#' treatment variable)} }
#' @source
#' 
#' V. Rondeau, S. Michiels, B.Liquet, and J.P. Pignon (2008). Investigating
#' trial and treatment heterogeneity in an individual patient data
#' meta-analysis of survival data by mean of the maximum penalized likelihood
#' approach. \emph{Statistics in Medecine}, \bold{27}, 1894-1910.
#' @keywords datasets
NULL





#' Simulated data for two types of recurrent events and a terminal event
#' 
#' This contains a simulated sample of of 800 subjects and 1652 observations.
#' This dataset can be used to illustrate how to fit a joint multivariate
#' frailty model. Two gaussian correlated random effects were generated with
#' mean 0, variances 0.5 and a correlation coefficient equals to 0.5. The
#' coefficients \eqn{\alpha_1} and \eqn{\alpha_2} were fixed to 1. The three
#' baseline hazard functions followed a Weibull distribution and right
#' censoring was fixed at 5.
#' 
#' 
#' @name dataMultiv
#' @docType data
#' @usage data(dataMultiv)
#' @format This data frame contains the following columns: \describe{
#' \item{PATIENT}{identification of patient} \item{obs}{number of observation
#' for a patient} \item{TIME0}{start of interval} \item{TIME1}{end of interval
#' (death or censoring time)} \item{INDICREC}{recurrent of type 1 status (0:no,
#' 1:yes)} \item{INDICMETA}{recurrent of type 2 status (0:no, 1:yes)}
#' \item{INDICDEATH}{censoring status (0:alive, 1:death)} \item{v1}{dichotomous
#' covariate (0,1)} \item{v2}{dichotomous covariate (0,1)}
#' \item{v3}{dichotomous covariate (0,1)} \item{TIMEGAP}{time to event} }
#' @keywords datasets
NULL





#' Simulated data for recurrent events and a terminal event with weigths using
#' nested case-control design
#' 
#' This contains a simulated sample of of 819 subjects and 1510 observations.
#' This dataset can be used to illustrate how to fit a joint frailty model for
#' data from nested case-control studies.
#' 
#' 
#' @name dataNCC
#' @docType data
#' @usage data(dataNCC)
#' @format This data frame contains the following columns: \describe{
#' \item{id}{identification of patient} \item{cov1}{dichotomous covariate
#' (0,1)} \item{cov2}{dichotomous covariate (0,1)} \item{t.start}{start of
#' interval} \item{t.stop}{end of interval (death or censoring time)}
#' \item{gaptime}{time to event} \item{event}{recurrent event status (0:no,
#' 1:yes)} \item{deathdays}{time of terminal event (death or right-censoring)}
#' \item{death}{censoring status (0:alive, 1:death)} \item{ncc.wts}{weights for
#' NCC design} }
#' @keywords datasets
NULL


#' Simulated data with two levels of grouping
#' 
#' This contains a simulated sample of 400 observations which allow
#' establishing 20 clusters with 4 subgroups and 5 subjects in each subgroup,
#' in order to obtain two levels of grouping. This data set is useful to
#' illustrate how to fit a nested model. Two independent gamma frailty
#' parameters with a variance fixed at 0.1 for the cluster effect and at 0.5
#' for the subcluster effect were generated. Independent survival times were
#' generated from a simple Weibull baseline risk function. The percentage of
#' censoring data was around 30 per cent. The right-censoring variables were
#' generated from a uniform distribution on [1,36] and a left-truncating
#' variable was generated with a uniform distribution on [0,10]. Observations
#' were included only if the survival time is greater than the truncated time.
#' 
#' 
#' @name dataNested
#' @docType data
#' @usage data(dataNested)
#' @format This data frame contains the following columns: \describe{
#' \item{group}{group identification variable} \item{subgroup}{subgroup
#' identification variable} \item{t1}{start of interval (0 or truncated time)}
#' \item{t2}{end of interval (death or censoring time)} \item{event}{censoring
#' status (0: alive, 1: death)} \item{cov1}{dichotomous covariate (0,1)}
#' \item{cov2}{dichotomous covariate (0,1)} }
#' @source
#' 
#' V. Rondeau, L. Filleul, P. Joly (2006). Nested frailty models using maximum
#' penalized likelihood estimation. \emph{Statistics in Medecine}, \bold{25},
#' 4036-4052.
#' @keywords datasets
NULL



#' Rehospitalization colorectal cancer
#' 
#' This contains rehospitalization times after surgery in patients diagnosed
#' with colorectal cancer
#' 
#' 
#' @name readmission
#' @docType data
#' @usage data(readmission)
#' @format This data frame contains the following columns: \describe{
#' \item{id}{identification of each subject. Repeated for each recurrence}
#' \item{enum}{which readmission} \item{t.start}{start of interval (0 or
#' previous recurrence time)} \item{t.stop}{recurrence or censoring time}
#' \item{time}{interocurrence or censoring time} \item{event}{rehospitalization
#' status. All event are 1 for each subject excepting last one that it is 0}
#' \item{chemo}{Did patient receive chemotherapy? 1: No; 2:Yes}
#' \item{sex}{gender: 1:Males 2:Females} \item{dukes}{Dukes' tumoral stage:
#' 1:A-B; 2:C 3:D} \item{charlson}{Comorbidity Charlson's index. Time-dependent
#' covariate.  0: Index 0; 1: Index 1-2; 3: Index >=3 } \item{death}{death
#' indicator. 1:dead and 0:alive } }
#' @source
#' 
#' Gonzalez, JR., Fernandez, E., Moreno, V., Ribes, J., Peris, M., Navarro, M.,
#' Cambray, M. and Borras, JM (2005). Sex differences in hospital readmission
#' among colorectal cancer patients. \emph{Journal of Epidemiology and
#' Community Health}, \bold{59}, 6, 506-511.
#' @keywords datasets
NULL

##' Advanced Ovarian Cancer dataset
##' 
##' This dataset combines the data  that were collected in four double-blind randomized
##' clinical trials in advanced ovarian cancer. In these trials, the objective was to 
##' examine the efficacy of cyclophosphamide plus cisplatin (CP) versus cyclophosphamide 
##' plus adriamycin plus cisplatin (CAP) to treat advanced ovarian cancer. The candidate 
##' surrogate endpoint \bold{S} is progression-free survival time, defined as the time (in years)
##' from randomization to clinical progression of the disease or death. The true endpoint
##' \bold{T} is survival time, defined as the time (in years) from randomization to death of any 
##' cause
##' 
##' 
##' @name dataOvarian
##' @docType data
##' @usage data(dataOvarian)
##' @format This data frame contains the following columns: 
##' \describe{
##' \item{patientID}{The identification number of a patient} 
##' \item{trialID}{The center in which a patient was treated}
##' \item{trt}{The treatment indicator, coded as 0 = cyclophosphamide plus cisplatin (CP)
##' and 1 = cyclophosphamide plus adriamycin plus cisplatin(CAP)}
##' \item{timeS}{The candidate surrogate (progression-free survival)}
##' \item{statusS}{Censoring indicator for for Progression-free survival}
##' \item{timeT}{The true endpoint (survival time)}
##' \item{statusT}{Censoring indicator for survival time}
##' }
##' @source
##' 
##' Ovarian cancer Meta-Analysis Project (1991). Cyclophosphamide plus cisplatin plus adriamycin
##' versus Cyclophosphamide, doxorubicin, and cisplatin chemotherapy of ovarian carcinoma: 
##' A meta-analysis. \emph{Classic Papers and Current Comments}, \bold{3}, 237-234.
##' @keywords datasets
NULL

##' Advanced Gastric Cancer dataset
##' 
##' This meta-analysis was carried out by the GASTRIC (Global Advanced/Adjuvant Stomach 
##' Tumor Research international Collaboration) group, using individual data on patients 
##' with curatively resected gastric cancer. Data from all published randomized trials, 
##' with a patient recruitment end date before 2004, and comparing adjuvant chemotherapy 
##' with surgery alone for resectable gastric cancers, were searched electronically. 
##' The candidate surrogate endpoint \bold{S} was Desease-free survival time, defined 
##' as the time (in days) to relapse, second cancer or dead from any cause. The true 
##' endpoint \bold{T} was the overall survival time, defined as the time (in days) from 
##' randomization to death of any cause or to the last follow-up.
##' 
##' 
##' @name gastadj
##' @docType data
##' @usage data(gastadj)
##' @format This data frame contains the following columns: 
##' \describe{
##' \item{trialID}{The trial in which the patient was treated}
##' \item{patientID}{The identification number of a patient} 
##' \item{trt}{The treatment indicator, coded as 0 = Control and 1 = Experimental}
##' \item{timeS}{The candidate surrogate (progression-free survival in days) }
##' \item{statusS}{Censoring indicator for for Progression-free survival 
##' (0 = alive and progression-free, 1 = with progression or dead)}
##' \item{timeT}{The true endpoint (overall survival time in days)}
##' \item{statusT}{Censoring indicator for survival time (0 = alive, 1 = dead)}
##' }
##' @source
##' 
##' Oba K, Paoletti X, Alberts S, Bang YJ, Benedetti J, Bleiberg H, Catalona P, 
##' Lordick F, Michiels S, Morita A, Okashi Y, Pignon JP, Rougier P, Sasako M, 
##' Sakamoto J, Sargent D, Shitara K, Van Cutsem E, Buyse M, Burzykowski T on 
##' behalf of the GASTRIC group (2013). Disease-Free Survival as a Surrogate 
##' for Overall Survival in Adjuvant Trials of Gastric Cancer: A Meta-Analysis. 
##' \emph{JNCI: Journal of the National Cancer Institute};\bold{105(21)}:1600-1607
##' @keywords datasets
NULL

##' Longitudinal semicontinuous biomarker dataset (TPJM)
##' 
##' This is a simulated dataset used to illustrate the two-part 
##' joint model included in the longiPenal function.
##' 
##' @name longDat
##' @docType data
##' @usage data(longDat)
##' @format This data frame contains the following columns: 
##' \describe{
##' \item{id}{The identification number of a patient}
##' \item{timej}{The measurement times of the biomarker}
##' \item{trtY}{Treatment covariate}
##' \item{Y}{Biomarker value}
##' }
NULL

##' Survival dataset (TPJM)
##' 
##' This is a simulated dataset used to illustrate the two-part 
##' joint model included in the longiPenal function.
##' 
##' @name survDat
##' @docType data
##' @usage data(survDat)
##' @format This data frame contains the following columns: 
##' \describe{
##' \item{id}{The identification number of a patient}
##' \item{deathTimes}{The event times (death or censoring)}
##' \item{d}{Censoring indicator}
##' \item{trt}{Treatment covariate}
##' }
NULL