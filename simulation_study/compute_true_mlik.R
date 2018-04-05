#
#Compute marg. lik. with true groups and ONE group
#

library(INLA)
library(MCMCpack)

#Gaussian data
load("gaussian_data.RData")

#One group
inla.hyperpar(inla(y ~ 1, data = d,
  control.fixed = list(prec.intercept = 0.001),
  control.inla = list(strategy = "laplace")
))$mlik
attr(MCMCregress(y ~ 1, data = d, B0 = 0.001, c0 = 1, d0 = 5e-05,
  marginal.likelihood = "Chib95",
  mcmc = 50000, thin = 10),
   "logmarglike")

#Three groups
inla.hyperpar(inla(y ~ -1 +  as.factor(grp), data = d, 
   control.fixed = list(prec = 0.001), control.inla = list(strategy = "laplace")))$mlik
attr(MCMCregress(y ~ as.factor(grp), data = d, B0 = rep(0.001, 3),  
  c0 = 1, d0 = 5e-05,
     marginal.likelihood = "Chib95"),
   "logmarglike")

#Poisson data
load("poisson_data.RData")

#One group
inla.hyperpar(inla(y ~ 1, data = d, family = "poisson", 
  control.fixed = list(prec.intercept = 0.001)))$mlik
inla.hyperpar(inla(y ~ 1, data = d, family = "poisson", 
  control.fixed = list(prec.intercept = 0.001),
  control.inla = list(strategy = "laplace")))$mlik
attr(MCMCpoisson(y ~ 1, data = d, B0 = 0.001, marginal.likelihood = "Laplace"),
   "logmarglike")

#Three groups
inla.hyperpar(inla(y ~ -1 +  as.factor(grp), data = d, family = "poisson"))$mlik
attr(MCMCpoisson(y ~ as.factor(grp), data = d, B0 = 0.001,
  marginal.likelihood = "Laplace"), "logmarglike")


