General Frailty Models: Shared, Joint and Nested Frailty Models with Prediction, Evaluation of Failure-Time Surrogate Endpoints
================
Virginie Rondeau, Juan R. Gonzalez, Yassin Mazroui, Audrey Mauguen, Amadou Diakite, Alexandre Laurent, Myriam Lopez, Agnieszka Krol and Casimir L. Sofeu
07 mars 2019

Description:
============

The following several classes of frailty models using a penalized likelihood estimation on the hazard function but also a parametric estimation can be fit using this R package:

    1) A shared frailty model (with gamma or log-normal frailty distribution) and Cox proportional hazard model. Clustered and recurrent survival times can be studied.
    2) Additive frailty models for proportional hazard models with two correlated random effects (intercept random effect with random slope).
    3) Nested frailty models for hierarchically clustered data (with 2 levels of clustering) by including two iid gamma random effects.
    4) Joint frailty models in the context of the joint modelling for recurrent events with terminal event for clustered data or not. A joint frailty model for two semi-competing risks and clustered data is also proposed.
    5) Joint general frailty models in the context of the joint modelling for recurrent events with terminal event data with two independent frailty terms.
    6) Joint Nested frailty models in the context of the joint modelling for recurrent events with terminal event, for hierarchically clustered data (with two levels of clustering) by including two iid gamma random effects.
    7) Multivariate joint frailty models for two types of recurrent events and a terminal event.
    8) Joint models for longitudinal data and a terminal event.
    9) Trivariate joint models for longitudinal data, recurrent events and a terminal event. 
    10) Joint frailty models for the validation of surrogate endpoints in multiple randomized clinical trials with failure-time endpoints 
        

Prediction values are available (for a terminal event or for a new recurrent event). Left-truncated (not for Joint model), right-censored data, interval-censored data (only for Cox proportional hazard and shared frailty model) and strata are allowed. In each model, the random effects have the gamma or normal distribution. Now, you can also consider time-varying covariates effects in Cox, shared and joint frailty models (1-5). The package includes concordance measures for Cox proportional hazards models and for shared frailty models.

Installation
------------

The R package frailtypack is available on cran at <https://cran.r-project.org/web/packages/frailtypack/index.html>

You can install frailtypack from github with:

``` r
# install.packages("frailtypack")
devtools::install_github("socale/frailtypack", ref = "master")
```

Example
-------

This is a basic example which shows you how to use frailtypack in the context of the validation of a candidate surrogate endpoint. We use genarated data.

``` r
## Data Generation
data.sim <- jointSurrSimul(n.obs=600, n.trial = 30,cens.adm=549.24, 
    alpha = 1.5, theta = 3.5, gamma = 2.5, zeta = 1, sigma.s = 0.7, 
    sigma.t = 0.7, rsqrt = 0.8, betas = -1.25, betat = -1.25, full.data = 0, 
    random.generator = 1, seed = 0, nb.reject.data = 0)
```

``` r
# Estimation
joint.surro.sim.MCGH <- jointSurroPenal(data = data.sim, int.method = 2,
    nb.mc = 300, nb.gh = 20)
```

``` r
# Result
summary(joint.surro.sim.MCGH)
```

Check resullt on branch Master
------------------------------

Travis [![Travis-CI Build Status](https://travis-ci.org/socale/frailtypack.svg?branch=master)](https://travis-ci.org/socale/frailtypack)

AppVeyor [![AppVeyor Build Status](https://ci.appveyor.com/api/projects/status/github/socale/frailtypack?branch=master&svg=true)](https://ci.appveyor.com/project/socale/frailtypack)
