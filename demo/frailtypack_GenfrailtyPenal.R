# Jocelyn CHAUVET 2021-04-01


################################################################################
###                           loading PACKAGES                               ###
################################################################################
options(digits=5)
if(!require("frailtypack"))stop("this test requires package frailtypack.")
if(!require("survival"))stop("this test requires survival.")
if(!require("boot"))stop("this test requires boot.")
if(!require("MASS"))stop("this test requires MASS.")
if(!require("survC1"))stop("this test requires survC1.")
if(!require("timereg"))stop("this test requires timereg.")
########################################################################


################################################################################
###                             loading DATA                                 ###
################################################################################
adult.retino <- retinopathy[retinopathy$type == "adult", ]
adult.retino[adult.retino$futime >= 50, "status"] <- 0
adult.retino[adult.retino$futime >= 50, "futime"] <- 50
data("readmission")
readmission[, 3:5] <- readmission[, 3:5]/365.25
################################################################################





################################################################################
###                     PARAMETRIC SHARED FRAILTY GSMs                       ###
###                        (PHM, AHM, POM and PROM)                          ###
################################################################################
sf.param.ph <- GenfrailtyPenal(
  formula = Surv(futime,status)~trt+cluster(id), 
  data = adult.retino, print.times = FALSE,
  hazard = "parametric", family = "PH")
summary(sf.param.ph)
plot(sf.param.ph, ylim = c(0,0.06), main = "Proportional Hazards Model")


sf.param.po <- GenfrailtyPenal(
  formula = Surv(futime,status)~trt+cluster(id), 
  data = adult.retino, print.times = FALSE,
  hazard ="parametric", family = "PO")
summary(sf.param.po)
plot(sf.param.po, ylim = c(0,0.06), main = "Proportional Odds Model")


sf.param.pr <- GenfrailtyPenal(
  formula = Surv(futime,status)~trt+cluster(id), 
  data = adult.retino, print.times = FALSE, 
  hazard = "parametric", family = "probit")
summary(sf.param.pr) 
plot(sf.param.pr, ylim = c(0,0.06), main = "Probit Model")


sf.param.ah <- GenfrailtyPenal(
  formula = Surv(futime,status)~trt+cluster(id), 
  data = adult.retino, print.times = FALSE,
  hazard = "parametric", family = "AH")
summary(sf.param.ah)
plot(sf.param.ah, ylim = c(0,0.1), main = "Additive Hazards Model")
################################################################################


################################################################################
###    SEMI-PARAMETRIC SHARED FRAILTY GSMs with time-varying coefficients    ###
###                               (PHM and AHM)                              ###
################################################################################
sf.semiparam.ph <- GenfrailtyPenal(
  formula = Surv(futime,status)~cluster(id)+timedep(trt), 
  data = adult.retino, print.times = FALSE,
  n.knots = 4, kappa = 10^6, betaknots = 1, betaorder = 3, 
  hazard = "Splines", family = "PH")
print(sf.semiparam.ph)
plot(sf.semiparam.ph, ylim = c(0,0.05), main = "Proportional Hazards Model")


sf.semiparam.ah <- GenfrailtyPenal(
  formula = Surv(futime,status)~cluster(id)+timedep(trt), 
  data = adult.retino, print.times = FALSE,
  n.knots = 4, kappa = 10^6, betaknots = 1, betaorder = 3, 
  hazard = "Splines", family = "AH")
print(sf.semiparam.ah)
plot(sf.semiparam.ah, ylim = c(0,0.05), main = "Additive Hazards Model")
################################################################################







##################################################################################
###                         PARAMETRIC JOINT FRAILTY GSMs                      ###
###                         (dual- PHM, AHM, POM and PROM)                     ###
##################################################################################
param.dualph.cst <- GenfrailtyPenal(
  formula = Surv(t.start,t.stop,event)~cluster(id)+terminal(death)+sex+dukes+chemo,
  formula.terminalEvent = ~sex+dukes+chemo, 
  data = readmission, recurrentAG = TRUE, print.times = FALSE, 
  hazard = "parametric", family = c("PH","PH"))
summary(param.dualph.cst)
par(mfrow=c(1,1))
plot(param.dualph.cst, conf.bands = TRUE, ylim = c(0,1), main = "Proportional Hazards Models")


param.dualpo.cst <- GenfrailtyPenal(
  formula = Surv(t.start,t.stop,event)~cluster(id)+terminal(death)+sex+dukes+chemo,
  formula.terminalEvent = ~sex+dukes+chemo, 
  data = readmission, recurrentAG = TRUE, print.times = FALSE, 
  init.B = c(-0.77,0.61,2.10,-0.33,-0.39,2.02,5.36,1.27), 
  hazard = "parametric", family = c("PO","PO"))
summary(param.dualpo.cst)
plot(param.dualpo.cst, conf.bands = TRUE, ylim = c(0,1), main = "Proportional Odds Models")


param.dualpr.cst <- GenfrailtyPenal(
  formula = Surv(t.start,t.stop,event)~cluster(id)+terminal(death)+sex+dukes+chemo,
  formula.terminalEvent = ~sex+dukes+chemo, 
  data = readmission, recurrentAG = TRUE, print.times = FALSE, 
  hazard = "parametric", family = c("probit","probit"))
summary(param.dualpr.cst)
plot(param.dualpr.cst, conf.bands = TRUE, ylim = c(0,1), main = "Probit Models")


param.dualah.cst <- GenfrailtyPenal(
  formula = Surv(t.start,t.stop,event)~cluster(id)+terminal(death)+sex+dukes+chemo,
  formula.terminalEvent = ~sex+dukes+chemo, 
  data = readmission, recurrentAG = TRUE, print.times = FALSE, 
  hazard = "parametric", family = c("AH","AH"))
summary(param.dualah.cst)
plot(param.dualah.cst, conf.bands = TRUE, ylim = c(0,1), main = "Additive Hazards Models")
##################################################################################




##################################################################################
###      Semi-PARAMETRIC JOINT FRAILTY GSMs with time-varying coefficients     ###
###                         (dual- PHM, AHM, POM and PROM)                     ###
##################################################################################
semiparam.dualph.timedep = GenfrailtyPenal(
  formula = Surv(t.start,t.stop,event)~cluster(id)+terminal(death)+sex+timedep(chemo),
  formula.terminalEvent = ~sex+chemo,
  data = readmission, recurrentAG = TRUE,
  n.knots = 4, kappa = c(5,10), betaknots = 1, betaorder = 3, 
  hazard = "Splines", family=c("PH","PH"))
print(semiparam.dualph.timedep)
par(mfrow=c(1,1))
plot(semiparam.dualph.timedep, conf.band = TRUE, ylim = c(0,1.1), main = "Proportional Hazards Models")


semiparam.dualah.timedep = GenfrailtyPenal(
  formula = Surv(t.start,t.stop,event)~cluster(id)+terminal(death)+sex+timedep(chemo),
  formula.terminalEvent = ~sex+chemo,
  data = readmission, recurrentAG = TRUE,
  n.knots = 4, kappa = c(5,10), betaknots = 1, betaorder = 3, 
  hazard = "Splines", family=c("AH","AH"))
print(semiparam.dualah.timedep)
plot(semiparam.dualah.timedep, conf.band = TRUE, ylim = c(0,1.1), main = "Additive Hazards Models")
##################################################################################
