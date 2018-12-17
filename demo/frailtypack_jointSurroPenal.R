# Casimir Ledoux SOFEU 2018-11-02

options(digits=12)
if(!require("frailtypack"))stop("this test requires package frailtypack.")
if(!require("survival"))stop("this test requires survival.")
if(!require("boot"))stop("this test requires boot.")
if(!require("MASS"))stop("this test requires MASS.")
if(!require("survC1"))stop("this test requires survC1.")
if(!require("doBy"))stop("this test requires doBy")

cat("frailtypack test for the one-step joint surrogate model ...")


###########################################################################
### ONE-STEP JOINT SURROGATE MODEL FOR THE VALIDATION OF SURROGATE ENDPOINT
###########################################################################

# Generation of data to use

data.sim <- jointSurrSimul(n.obs=600, n.trial = 30,cens.adm=549.24, 
            alpha = 1.5, theta = 3.5, gamma = 2.5, zeta = 1, sigma.s = 0.7, 
            sigma.t = 0.7, rsqrt = 0.8, betas = -1.25, betat = -1.25, 
            full.data = 0, random.generator = 1, seed = 0, nb.reject.data = 0)


# model estimation
joint.surro.sim <- jointSurroPenal(data = data.sim, n.knots = 6, 
          int.method = 2, nb.mc = 300, nb.gh = 20, true.init.val = 0,
          print.itter=FALSE)


# Estimation using real data
data(dataOvarian)

joint.surro.ovar <- jointSurroPenal(data = dataOvarian, indicator.alpha = 0, 
        n.knots = 8, int.method = 2, nb.mc = 200, nb.gh = 20, true.init.val = 0, 
           init.kappa = c(2000,1000), scale = 1/365)


# ========================Results summaries=============================

summary(joint.surro.sim)
summary(joint.surro.ovar)
summary(joint.surro.ovar, int.method.kt = 1)

# ==============================Figures=================================

pdf(file="fig_jointSurroPenal.pdf")
plot(joint.surro.sim, type.plot = "Su", main = "Generated data", conf.bands=TRUE, 
     pos.legend="topright", cex.legend=0.7, main, color=2, endpoint = 2)
plot(joint.surro.sim, type.plot = "Haz", main = "Generated data", conf.bands=TRUE, 
     pos.legend="topleft", cex.legend=0.7, main, color=2, endpoint = 2)
plot(joint.surro.ovar, type.plot = "Su", main = "Ovarian Cancer data", conf.bands=TRUE, 
     pos.legend="topright", cex.legend=0.7, main, color=2, endpoint = 2)
plot(joint.surro.ovar, type.plot = "Haz", main = "Ovarian Cancer data", conf.bands=TRUE, 
     pos.legend="topleft", cex.legend=0.7, main, color=2, endpoint = 2)

#dev.off()
#dev.off()