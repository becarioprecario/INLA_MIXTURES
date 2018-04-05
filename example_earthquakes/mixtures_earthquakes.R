## Example: Fitting mixtures with R-INLA

library(MCMCpack)#For dirichlet distribution
library(MixtureInf)#For dataset
library(INLA)
library(parallel)
options(mc.cores = 4)


#Load code to fit mixtures
source("../Rcode/INLAmixture.R")
source("../Rcode/compute_marglik.R")
source("../Rcode/z_sample.R")
source("../Rcode/fitmix_inla.R")



#Get data
data(earthquake)

#Display data
years <- 1900:2006
yy <- earthquake$number

pdf(file = "earthquakes.pdf", width = 10, height = 5)
par(mfrow = c(1, 2))
plot(years, yy, xlab = "Year", ylab = "Count", type = "l")
points(years, yy, pch = 19)

hist(yy, xlab = "Count", main = "")
dev.off()




#Number of groups
n.grp <- 4

#Model with 1 component
d.mix0 <- data.mixture(earthquake$number, n.grp = 1)
m0 <- fitmix.inla.poisson(d.mix0, rep(1, length(d.mix0$y)), 1)
summary(m0$model)

#Create data for model fitting
set.seed(1)
d.mix <- data.mixture(earthquake$number, n.grp = n.grp)


#Check
m1 <- fitmix.inla.poisson(d.mix, d.mix$z, n.grp)
summary(m1$model)


res <- INLAmixture(d.mix, fitmix.inla.poisson, d.mix$z, z.sample.poisson,
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

plot(earthquake$number, probs)


#Display marginals
marg.mean <- INLABMA:::fitmargBMA2(mm, rep(1/iter.tot, iter.tot),
  "marginals.fixed")
marg.mean <- lapply(marg.mean, function(X) {
  inla.tmarginal(function(x) exp(x), X)
})

par(mfcol = c(1,n.grp))
for(i in 1:n.grp) {
  plot(marg.mean[[i]], type ="l", main = paste0("log-mu", i))
}

#Probabilities of each group weight
margz <- do.call(rbind, lapply(zz, table)) / length(earthquake$number)
#margw <- do.call(rbind, ww)


#Summary statistics
lapply(marg.mean, function(X) { unlist(inla.zmarginal(X, TRUE))} )
#apply(margw, 2, mean)
#apply(margw, 2, sd)


save(file = paste0("results-poisson-ngrp_z-", n.grp, ".RData"), list = ls())
