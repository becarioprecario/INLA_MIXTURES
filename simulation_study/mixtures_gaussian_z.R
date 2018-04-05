## Example: Fitting mixtures with R-INLA
#
# Samples from z using a EM-like approach

library(MASS)#For 'galaxies' dataset
library(MCMCpack)#For dirichlet distribution
library(INLA)
library(parallel)
options(mc.cores = 2)

#Load code to fit mixtures
source("../Rcode/INLAmixture.R")
source("../Rcode/compute_marglik.R")
source("../Rcode/z_sample.R")
source("../Rcode/fitmix_inla.R")


#Get data
load("gaussian_data.RData")

#Number of groups
n.grp <- 5 

#Model with 1 component
d.mix0 <- data.mixture(d$y, n.grp = 1)
m0 <- fitmix.inla.gaussian(d.mix0, rep(1, length(d$y)), 1)
summary(m0$model)

#Create data for model fitting
set.seed(1)
d.mix <- data.mixture(d$y, n.grp = n.grp)


#Check
m1 <- fitmix.inla.gaussian(d.mix, d.mix$z, n.grp)
summary(m1$model)


res <- INLAmixture(d.mix, fitmix.inla.gaussian, d.mix$z, z.sample.gaussian, 
  n.sim = 1000, n.burnin = 200, n.thin = 10, verbose = TRUE)

#Compute marginal likelihood
mixture.mlik(d.mix, fitmix.inla.gaussian, res$z.sim, n.grp)


#Subset
zz <- res$z.sim
mm <- lapply(res$model.sim, function(X){ X$model})
ww <- NULL
mliks <- lapply(res$model.sim, function(X){X$mlik})

iter.tot <- length(zz)

zz.tab <- do.call(cbind, zz)

probs <- apply(zz.tab, 1, function(X){
  sum(X[] == 1) / iter.tot
})
probs

plot(d$y, probs)

#Display marginals
marg.mean <- INLABMA:::fitmargBMA2(mm, rep(1/iter.tot, iter.tot),
  "marginals.fixed")
marg.prec <- INLABMA:::fitmargBMA2(mm, rep(1/iter.tot, iter.tot),
  "marginals.hyperpar")

par(mfcol = c(2,n.grp))
for(i in 1:n.grp) {
  plot(marg.mean[[i]], type ="l", main = paste0("mu", i))
  plot(marg.prec[[i]], type ="l", main = paste0("prec", i))
}

#Probabilities of each group weight
margz <- do.call(rbind, lapply(zz, table)) / length(d$y)
#margw <- do.call(rbind, ww)


#Summary statistics
lapply(marg.mean, function(X) { unlist(inla.zmarginal(X, TRUE))} )
lapply(marg.prec, function(X) { unlist(inla.zmarginal(X, TRUE))} )
#apply(margw, 2, mean)
#apply(margw, 2, sd)

#Compute marginal likelihood
mixture.mlik(d.mix, fitmix.inla.gaussian, zz, n.grp)

save(file = paste0("results-gaussian-ngrp_z-", n.grp, ".RData"), list = ls())

