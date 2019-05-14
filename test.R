library(frailtypack)

# test de la fonction de generation des donnees avec les copules
data.sim <- jointSurrCopSimul()
summary(data.sim)
head(data.sim)
tail(data.sim)

#evaluation de la generation a partir des statistiques empirique 
result <- frailtypack:::param.empirique(nsim = 400, cens.adm=549, ver =2, 
                lambda.S = 1.3, nu.S = 0.0025,lambda.T = 1.1, nu.T = 0.0025, 
                seed = 0, dec = 3, betas = c(-1.25, 0.5),betat = c(-1.25, 0.5), 
                filter.surr = c(1,1), filter.true = c(1,1), variatio.seed = 0,
                random.generator = 1)
result
# censure a 349 jours: remps en jours
result <- frailtypack:::param.empirique(nsim = 100, n.obs = 1000, cens.adm=349, ver =1, 
                                        lambda.S = 1.3, nu.S = 0.0025,lambda.T = 1.1, nu.T = 0.0025, 
                                        seed = 0, dec = 3, betas = c(-1.25),betat = c(-1.25), 
                                        filter.surr = c(1), filter.true = c(1), variatio.seed = 0,
                                        cor = 0.9, random.generator = 1)
result
# censure a 349 jours: remps en jours
result2 <- frailtypack:::param.empirique(nsim = 100, cens.adm=349, ver =1, 
                                        lambda.S = 1.3, nu.S = 0.0025,lambda.T = 1.1, nu.T = 0.0025, 
                                        seed = 0, dec = 3, betas = c(-1.25),betat = c(-1.25), 
                                        filter.surr = c(1), filter.true = c(1), variatio.seed = 0,
                                        cor = 0.9, random.generator = 2)
result2

# censure a 349 jours: remps en jours
result <- frailtypack:::param.empirique(nsim = 500, cens.adm=349, ver =1, 
                                        lambda.S = 1.3, nu.S = 0.0025,lambda.T = 1.1, nu.T = 0.0025, 
                                        seed = 0, dec = 3, betas = c(-1.25),betat = c(-1.25), 
                                        filter.surr = c(1), filter.true = c(1), variatio.seed = 0,
                                        random.generator = 1)
result

# censure a 8 ans: lambda =1 ie une exponentielle
result <- frailtypack:::param.empirique(nsim = 100, cens.adm = 10, ver = 1, 
                                        lambda.S = 1,lambda.T = 1, nu.S = 1, nu.T = 0.5, 
                                        seed = 0, dec = 3, betas = c(-1.25),betat = c(-1.25), 
                                        filter.surr = c(1), filter.true = c(1))
result

data.sim.cop <- jointSurrCopSimul(n.obs=100, n.trial = 5,cens.adm = 8, 
                                  alpha = 1, gamma = 2.5, lambda.S = 1, nu.S = 1, lambda.T = 1, nu.T = 0.5, 
                                  cor = 0.9, betas = c(-1.25), betat = c(-1.25), 
                                  full.data = 0, random.generator = 1,ver = 1, covar.names = "trt", 
                                  nb.reject.data = 0, random.nb.sim = 1, thetacopule = 6, filter.surr = c(1), 
                                  filter.true = c(1), seed = 0)

data.sim.cop$trt <- as.factor(data.sim.cop$trt) 
data.sim.cop$statusS <- as.factor(data.sim.cop$statusS) 
data.sim.cop$statusT <- as.factor(data.sim.cop$statusT)
tapply(data.sim.cop[, "trt"], data.sim.cop$trialID, summary)
tapply(data.sim.cop[, "statusS"], data.sim.cop$trialID, summary)
tapply(data.sim.cop[, "statusT"], data.sim.cop$trialID, summary)

describe(data.sim.cop)

# simulation study
joint.simul <- jointSurroPenalSimul(nb.dataset = 1, nb.mc = 50, nb.gh = 5, true.init.val = 1, 
                nbSubSimul=300, ntrialSimul=30, type.joint.simul = 3, maxit = 2,
                  LIMparam = 0.001, LIMlogl = 0.001, LIMderiv = 0.001,nb.gh2 = 9, 
                print.iter=T, kappa.use = 1
                 )

joint.simul1 <- jointSurroPenalSimul(nb.dataset = 1, nbSubSimul=100, 
                                    ntrialSimul=10, nb.mc = 50, nb.gh = 5, 
                                    nb.gh2 = 32, print.iter=T, kappa.use = 1)

# results
summary(joint.simul1, d = 3, R2boot = 1) # bootstrap

joint.simul2 <- jointSurroPenalSimul(nb.dataset = 1, nbSubSimul=100, 
                                     ntrialSimul = 10, nb.mc = 50, nb.gh = 5, 
                                     nb.gh2 = 32, print.iter = T, kappa.use = 1, type.joint.estim = 1,
                                     type.joint.simul = 3, time.cens = 549, 
                                     lambdas = 1.3, nus = 0.0025, lambdat = 1.1, nut = 0.0025, 
                                     seed = 0, betas = c(-1.25), betat = c(-1.25), 
                                     init.kappa = c(10,20), maxit = 4)

# results
summary(joint.simul2, d = 3, R2boot = 1) # bootstrap

data.sim.cop <- jointSurrCopSimul(n.obs=100, n.trial = 10,cens.adm=549, 
                                  alpha = 1, gamma = 2.5, lambda.S = 1.3, nu.S = 0.0025, lambda.T = 1.1, nu.T = 0.0025, 
                                  cor = 0.9, betas = c(-1.25, 0.5), betat = c(-1.25, 0.5), 
                                  full.data = 0, random.generator = 1,ver = 2, covar.names = "trt", 
                                  nb.reject.data = 0, random.nb.sim = 1, thetacopule = 6, filter.surr = c(1,1), 
                                  filter.true = c(1,1), seed = 0)

head(data.sim.cop)

# integration MC
joint.simul2 <- frailtypack::jointSurroPenalSimul(nb.dataset = 1, nbSubSimul=100, int.method = 0, 
                                     ntrialSimul = 10, nb.mc = 100,
                                     print.iter = T, kappa.use = 1, type.joint.estim = 3,
                                     type.joint.simul = 3, time.cens = 549, 
                                     lambdas = 1.3, nus = 0.0025, lambdat = 1.1, nut = 0.0025, 
                                     seed = 0, betas = c(-1.25, 0.5), betat = c(-1.25, 0.5),
                                     betas.init = c(-0.25, 0.25), betat.init = c(-0.25, 0.25),
                                     init.kappa = c(1000,2000), maxit = 10, true.init.val = 2, 
                                     typecopula = 1)

# results
summary(joint.simul2, d = 3, R2boot = 1) # bootstrap

# integration GH
joint.simul2 <- frailtypack::jointSurroPenalSimul(nb.dataset = 1, nbSubSimul=100, 
   int.method = 0, nb.mc = 50, ntrialSimul = 10, nb.gh = 5, nb.gh2 = 9, adaptatif = 1,
   nb.iterPGH = 0,
   print.iter = T, kappa.use = 0, type.joint.estim = 3,
   type.joint.simul = 3, time.cens = 549, 
   lambdas = 1.3, nus = 0.0025, lambdat = 1.1, nut = 0.0025, 
   seed = 0, betas = c(-1.25, 0.5), betat = c(-1.25, 0.5),
   betas.init = c(-0.25, 0.25), betat.init = c(-0.25, 0.25),
   init.kappa = c(1000,2000), maxit = 10, true.init.val = 1, 
   typecopula = 1)   

# one covariate

joint.simul2 <- frailtypack::jointSurroPenalSimul(nb.dataset = 1, nbSubSimul=100, 
      int.method = 0, nb.mc = 50, ntrialSimul = 10, nb.gh = 5, nb.gh2 = 9, adaptatif = 1,
      nb.iterPGH = 0,
      print.iter = T, kappa.use = 0, type.joint.estim = 3,
      type.joint.simul = 3, time.cens = 549, 
      lambdas = 1.3, nus = 0.0025, lambdat = 1.1, nut = 0.0025, 
      seed = 0, betas = c(-1.25), betat = c(-1.25), filter.surr = 1,
      betas.init = c(-0.25), betat.init = c(-0.25), filter.true = 1,
      init.kappa = c(1000,2000), maxit = 10, true.init.val = 1, 
      typecopula = 1)  

# one covariate and simulation using joint-frailty model

joint.simul2 <- frailtypack::jointSurroPenalSimul(nb.dataset = 4, nbSubSimul=100, 
      int.method = 0, nb.mc = 5, ntrialSimul = 10, nb.gh = 5, nb.gh2 = 9, adaptatif = 1,
      nb.iterPGH = 0, maxit = 20,
      print.iter = T, kappa.use = 0, type.joint.estim = 3,
      type.joint.simul = 1, time.cens = 549, 
      lambdas = 1.3, nus = 0.0025, lambdat = 1.1, nut = 0.0025, 
      seed = 0, betas = c(-1.25), betat = c(-1.25), filter.surr = 1,
      betas.init = c(-0.25), betat.init = c(-0.25), filter.true = 1,
      init.kappa = NULL, true.init.val = 0, 
      typecopula = 2, nb.MC.kendall = 10, nboot.kendall = 10, LIMparam = 2.0, 
      LIMlogl = 1.01, LIMderiv = 1.000)
                                                  
joint.simul2 <- frailtypack::jointSurroPenalSimul(nb.dataset = 1, nbSubSimul=100, 
      int.method = 0, nb.mc = 5, ntrialSimul = 10, nb.gh = 5, nb.gh2 = 9, adaptatif = 1,
      nb.iterPGH = 0, maxit = 20,
      print.iter = T, kappa.use = 0, type.joint.estim = 3,
      type.joint.simul = 1, time.cens = 549, 
      lambdas = 1.3, nus = 0.0025, lambdat = 1.1, nut = 0.0025, 
      seed = 0, betas = c(-1.25), betat = c(-1.25), filter.surr = 1,
      betas.init = c(-0.25), betat.init = c(-0.25), filter.true = 1,
      init.kappa = NULL, ckappa = c(100,100), true.init.val = 0, 
      typecopula = 2, nb.MC.kendall = 10, nboot.kendall = 10, LIMparam = 2.0, 
      LIMlogl = 1.01, LIMderiv = 1.000)

joint.simul2 <- frailtypack::jointSurroPenalSimul(nb.dataset = 1, nbSubSimul=600, ntrialSimul=30, 
              int.method = 0, nb.mc = 300, maxit = 1, 
              #nb.gh = 5, nb.gh2 = 9, adaptatif = 1, nb.iterPGH = 0,
              print.iter = T, kappa.use = 0, type.joint.estim = 3,
              type.joint.simul = 3, time.cens = 349, n.knots =  6,
              lambdas = 1.3, nus = 0.0025, lambdat = 1.1, nut = 0.0025, 
              seed = 0, betas = c(-1.25), betat = c(-1.25), filter.surr = 1,
              betas.init = c(-0.25), betat.init = c(-0.25), filter.true = 1,
              init.kappa = NULL, ckappa = c(-900000,-900000), true.init.val = 0, 
              typecopula = 1, theta.copula = 1)#, LIMparam = 2.0, LIMlogl = 1.01, LIMderiv = 1.000)                                                  



# results
summary(joint.simul2, d = 5) 
summary(joint.simul2, d = 5, R2boot = 1) # bootstrap

#integer



# Displaying of the baseline hazards functions
par(mfrow = c(2,2))
n.s <- 0.0025
gamma.s <- 1.3
x <- seq(from = 1, to = 100, length.out = 1000)
y <- n.s * gamma.s * x**(gamma.s-1)

plot(x,y, main = "Surrogate endpoint")

n.s <- 0.0025
gamma.s <- 1.1
x <- seq(from = 1, to = 100, length.out = 1000)
y <- n.s * gamma.s * x**(gamma.s-1)

plot(x,y, main = "True endpoint", col = "red")

#Censure 5 ans
n.s <- 3.25
gamma.s <- 3.3
x <- seq(from = 1, to = 100, length.out = 1000)
y <- n.s * gamma.s * x**(gamma.s-1)

plot(x,y, main = "Surrogate endpoint")

n.s <- 0.45
gamma.s <- 0.8
x <- seq(from = 1, to = 100, length.out = 1000)
y <- n.s * gamma.s * x**(gamma.s-1)

plot(x,y, main = "True endpoint", col = "red")





taux_couverture_tauk
vrai_tau_copula

# transformer en vecteur


# revoir taille vecteur
varcov <- matrix(c(.7, .63, 0, .63, .7, 0, 0, 0, 2.5), 3, 3)
a = rmnorm(n = 300, mean = rep(0, 3), varcov, sqrt=NULL) 
var(a)


# calsul des taux de kendall
setwd("G:/socale/PHD-Thesis/programmes/Creation_Package/package_CRAN/Curta/exploitation_result")
load("joint.simul2_15.RData")
rangthetacop <- 22
h_theta_chap <- joint.simul2$dataHessian[rangthetacop,rangthetacop]
theta_chap <- joint.simul2$datab[1,rangthetacop]
# Clayton
theta <- exp(theta_chap)
KTau = theta/(theta+2)
se_ktau <- 2*exp(theta_chap) *sqrt(h_theta_chap)/(exp(theta_chap) + 2)^2
ic = c(KTau - 1.96 * se_ktau , KTau + 1.96 * se_ktau)

# Hougard
theta <- theta_chap^2
KTau = theta/(theta+1)
se_ktau <- 2*theta_chap *sqrt(h_theta_chap)/(theta_chap^2 + 1)^2
ic = c(KTau - 1.96 * se_ktau , KTau + 1.96 * se_ktau)
