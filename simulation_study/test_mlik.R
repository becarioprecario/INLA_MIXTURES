#Code to test how the marg. lik. is computed using the
#alternative methods of sum pi(y|z) pi(z)

#Load results
dataset <- "poisson"
n.grp <- 3
load(paste0("results-", dataset, "-ngrp_z-", n.grp, ".RData"))

#Source code
source("../Rcode/z_unique.R")
source("../Rcode/prior_z.R")
source("../Rcode/compute_mliks3.R")

#Get mliks as vector
mliks <- unlist(mliks)
names(mliks) <- NULL

#Get unique values of the configuration
zz.unique <- z_unique(zz, mliks, rep(2, n.grp))


compute_mliks3(zz.unique)


#Check agreement
plot(zz.unique$pprob.mliks, zz.unique$pprob.samp); abline(0, 1)

#Compare posterior probabilities between MCMC and MGS
load(paste0("jags-", dataset, "-ngrp-", n.grp, ".RData"))

#Z from JAGS
zz.jags <- apply(res.jags$Z[,,1], 2, list)
zz.jags <- zz.jags <- lapply(zz.jags, function(X){X[[1]]})

zz.unique.jags <- z_unique(zz.jags, NA, rep(2, n.grp))

#Compare number of configurations
length(zz.unique$unique)
length(zz.unique.jags$unique)

#Find index in MCMC
idx.mcmc <- lapply(zz.unique$unique, function(X) {
  which(unlist(lapply(zz.unique.jags$unique, function(Y) {
    isTRUE(all.equal(X, Y))
  })))
})

idx.mgs <- which(unlist(lapply(idx.mcmc, function(X) {length(X)>0})))
idx.mcmc <- unlist(idx.mcmc)

#Compare probabilities
par(mfrow = c(1, 2))
plot(log(zz.unique$pprob.samp[idx.mgs]), log(zz.unique.jags$pprob.samp[idx.mcmc]), main = "Post. Prob with MGS")
abline(0, 1)
plot(log(zz.unique$pprob.mliks[idx.mgs]), log(zz.unique.jags$pprob.samp[idx.mcmc]), main = "Post. Prob with MLIKS")
abline(0, 1)

