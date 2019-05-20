library(frailtypack)

# test de la fonction de generation des donnees avec les copules
n.sim = 1
cens = 349 # DOIT ETRE FIXE EXTREMEMENT GRAND POUR TESTER LE TAU DE KENDALL
n.obs = 300
n.trial = 10
lambdas = 1.3
lambdat = 1.1
nus = 0.0025
nut = 0.0025
nb.mc = 200
R2 = 0.81
n.knots =  8
cor = sqrt(R2)
true.init.val = 0
theta.copula = 3
typecopula = 1
type.joint.estim = 3
type.joint.simul = 3
kappa.use = 0
maxit = 100
print.iter = T

result <- frailtypack:::param.empirique(nsim = n.sim, cens.adm = cens, ver = 1, n.obs = n.obs, n.trial = n.trial, 
             lambda.S = lambdas,lambda.T = lambdat, nu.S = nus, nu.T = nut, 
             seed = 0, dec = 3, betas = c(-1.25),betat = c(-1.25), cor = cor, 
             typecopula = typecopula, filter.surr = c(1), filter.true = c(1))
result                                        

# estimation

joint.simul4 <- frailtypack::jointSurroPenalSimul(nb.dataset = n.sim, nbSubSimul=n.obs, ntrialSimul=n.trial, 
              int.method = 0, nb.mc = nb.mc, maxit = maxit, R2 = R2,
              #nb.gh = 5, nb.gh2 = 9, adaptatif = 1, nb.iterPGH = 0,
              print.iter = print.iter, kappa.use = kappa.use, type.joint.estim = type.joint.estim,
              type.joint.simul = type.joint.simul, time.cens = cens, n.knots =  n.knots,
              lambdas = lambdas, nus = nus, lambdat = lambdat, nut = nut, 
              seed = 0, betas = c(-1.25), betat = c(-1.25), filter.surr = 1,
              betas.init = c(0.25), betat.init = c(0.25), filter.true = 1,
              init.kappa = NULL, ckappa = c(0,0), true.init.val = true.init.val, 
              typecopula = typecopula, theta.copula = theta.copula)#, LIMparam = 2.0, LIMlogl = 1.01, LIMderiv = 1.000)                                                  



# resultsw
summary(joint.simul4, d = 5)

summary(joint.simul2, d = 5) 
summary(joint.simul2, d = 5, R2boot = 1) # bootstrap

#integer



# Displaying of the baseline hazards functions
par(mfrow = c(1,2))
n.s <- nus
gamma.s <- lambdas
x <- seq(from = 1, to = 100, length.out = 1000)
y <- n.s * gamma.s * x**(gamma.s-1)

plot(x,y, main = "Surrogate endpoint", xlab = "Time", ylab = "Baseline hazard")

n.s <- nut
gamma.s <- lambdat
x <- seq(from = 1, to = 100, length.out = 1000)
y <- n.s * gamma.s * x**(gamma.s-1)

plot(x,y, main = "True endpoint", col = "red", xlab = "Time", ylab = "Baseline hazard")

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
