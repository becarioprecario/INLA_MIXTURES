## Example: Fitting mixtures with R-INLA
#
# Samples from z using a EM-like approach
#
#USing queal variances in all groups

library(MASS)#For 'galaxies' dataset
library(MCMCpack)#For dirichlet distribution
library(INLA)
library(parallel)
options(mc.cores = 2)


#Get data
yy <- galaxies/1000

#Number of data
n <- length(yy)

#Display data
pdf(file = "galaxies.pdf")
plot(density(yy, bw = 1.5), main = "", xlab = "")
dev.off()


#Number of groups
n.grp <- 5 

# Using 3 as a cut-off point
#grp <- 1 + (yy >=3)

# 50-50 random assignment
set.seed(1)
grp <- sample(rep(1:n.grp, length.out = n))

#Data in two-column format
y <- matrix(NA, ncol = n.grp, nrow = n)
for(i in 1:n.grp) {
  idx <- which(grp == i)
  y[idx, i] <- yy[idx]
}

#X stores the intercept in the model
x <- y
x[!is.na(x)] <- 1

d <- list(y = yy, x = x)

#Precision priors
prec.prior1 <- list(list(hyper = list(prec = list(prior = "loggamma", 
#  param = c(0.38, 0.5)))))
  param = c(0.5, 0.5)))))
#  param = c(6/2, 40/2))))) #Chi (1995)
prec.prior <- replicate(n.grp, prec.prior1)

#Prior on the mean effects 
#c.fixed <- list(mean = 20, prec = 1/100)#(Chib, 1995)
c.fixed <- list(mean.intercept = 0, prec.intercept = 0.001, 
  mean = 0, prec = 0.001)

#One group
m0 <- inla(y ~ 1, data = d,
  control.family = prec.prior1[[1]], control.fixed = c.fixed
)

#Initial group
m1 <- inla(y ~ -1 + x, data = d,
  control.family = prec.prior1[[1]], control.fixed = c.fixed
)

iter <- 1

n.burnin <- 200
n.iter <- 1000
n.thin <- 10
total.iter <- n.burnin + n.iter  * n.thin

print(paste0("Total number of iterations: ", total.iter))

zz <- as.list (rep(NA, total.iter))
mm <- as.list (rep(NA, total.iter))
ww <- as.list(rep(NA, total.iter))
mliks <- rep(NA, total.iter)

d.aux <- d
m.aux <- m1


while(iter <= total.iter) {

#Define group index
#z.grp <- rep(1:2, each = n)
z.grp <- apply(d.aux$x, 1, function(X) {which(!is.na(X))})

#Number of observations per group
tab.grp <- table(z.grp) 

#Fix groups with zero observations
if(length(tab.grp) < n.grp) {
  tab.grp <- tab.grp[as.character(1:n.grp)]
  names(tab.grp) <- as.character(1:n.grp)
  tab.grp[is.na(tab.grp)] <- 0

}

d.param <- rep(1, n.grp)
d.param[as.numeric(names(tab.grp))] <- 1 + tab.grp

w.grp <- (tab.grp + 1)/(n + n.grp)  #rdirichlet(1, d.param)

#Get means of mixtures
means <- m.aux$summary.fixed[, "mean"]
precs <- rep(m.aux$summary.hyperpar[, "mean"], n.grp)

#Reorder means
idx <- order(means)
means <- means[idx]
precs <- precs[idx]

#Probability of being in each group
prob.z <- sapply(1:n.grp, function(X) {
  w.grp[X] * dnorm(yy, mean = means[X], sd = 1/sqrt(precs[X]))
})

prob.z <- t(apply(prob.z, 1, function(X) { X / sum(X)}))

set.seed(iter + 1000)
grp <- apply(prob.z, 1, function(X) {sample(1:n.grp, 1, prob = X)})


d.aux$x[,] <- NA
for(i in 1:n.grp) {
  idx <- which(grp == i)
  d.aux$x[idx, i] <- 1
}

m.aux <- inla(y ~ -1 + x, data = d.aux, family = "gaussian",
  control.family = prec.prior1, control.fixed = c.fixed
)
mliks[iter] <- m.aux$mlik[1,1]
print(paste0("Iter: ", iter, " mlik: ", mliks[iter]))

zz[[iter]] <- grp
mm[[iter]] <- m.aux
ww[[iter]] <- w.grp

iter <- iter + 1

}


#Burnin and thinning
iter.idx <- n.burnin + seq(n.thin, n.iter * n.thin, by = n.thin)
iter.tot <- length(iter.idx)

#Subset
zz <- zz[iter.idx]
mm <- mm[iter.idx]
ww <- ww[iter.idx]
mliks <- mliks[iter.idx]

zz.tab <- do.call(cbind, zz)

probs <- apply(zz.tab, 1, function(X){
  sum(X[] == 1) / iter.tot
})
probs

plot(yy, probs)

#Display marginals
marg.mean <- INLABMA:::fitmargBMA2(mm, rep(1/iter.tot, iter.tot),
  "marginals.fixed")
marg.prec <- INLABMA:::fitmargBMA2(mm, rep(1/iter.tot, iter.tot),
  "marginals.hyperpar")

par(mfcol = c(2,n.grp))
for(i in 1:n.grp) {
  plot(marg.mean[[i]], type ="l", main = paste0("mu", i))
}
  plot(marg.prec[[1]], type ="l", main = "prec")

#Probabilities of each group weight
margz <- do.call(rbind, lapply(zz, table))/n
margw <- do.call(rbind, ww)


#Summary statistics
lapply(marg.mean, function(X) { unlist(inla.zmarginal(X, TRUE))} )
lapply(marg.prec, function(X) { unlist(inla.zmarginal(X, TRUE))} )
apply(margw, 2, mean)
apply(margw, 2, sd)


save(file = paste0("results-ngrp_z_equalvars-", n.grp, ".RData"), list = ls())


