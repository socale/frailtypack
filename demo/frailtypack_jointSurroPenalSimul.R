# Casimir Ledoux SOFEU 2018-11-02

options(digits=12)
if(!require("frailtypack"))stop("this test requires package frailtypack.")
if(!require("survival"))stop("this test requires survival.")
if(!require("boot"))stop("this test requires boot.")
if(!require("MASS"))stop("this test requires MASS.")
if(!require("survC1"))stop("this test requires survC1.")
if(!require("doBy"))stop("this test requires doBy")

cat("frailtypack test for the simulation studies of the one-step joint surrogate model ...")


###########################################################################
### SIMULATION STUDIES FOR THE ONE-STEP JOINT SURROGATE MODEL FOR THE 
### VALIDATION OF SURROGATE ENDPOINT
###########################################################################

 joint.simul <- jointSurroPenalSimul(nb.dataset = 10, nbSubSimul=600, 
                    ntrialSimul=30, LIMparam = 0.001, LIMlogl = 0.001, 
                    LIMderiv = 0.001, nb.mc = 200, nb.gh = 20, 
                    nb.gh2 = 32, true.init.val = 1, print.itter=F)

# ========================Results summaries=============================

summary(joint.simul, d = 3, R2boot = 1) # bootstrap
summary(joint.simul, d = 3, R2boot = 0) # Delta-method
