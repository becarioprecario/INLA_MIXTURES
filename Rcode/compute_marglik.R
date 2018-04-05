## Compute marginal likelihood using Bayes' rule

#d: Data set.
#fit.inla: Function to fit INLA, as in INLAmixture().
#z.sample: A list of vectors with the sampled indices from MCMC.
#n.grp: Number of groups in the mixture

mixture.mlik <- function(d, fit.inla, z.mcmc, n.grp) {

n.mcmc <- length(z.mcmc)

#Create vector of factors with samples
zz.factor <- factor(unlist(lapply(z.mcmc, function(X) {paste0(X, collapse = ".")})))

zz.prob <- table(zz.factor)/n.mcmc


#Compute mode of z
zz.mode <- which.max(zz.prob)
zz.mode.post <- zz.prob[which.max(zz.prob)]

#Recover z
z <- as.integer(unlist(strsplit(names(zz.mode), "\\.")))

# Compute model at mode
m.mode <- fit.inla(d, z, n.grp)

#Compute posterior probabilites at mode
log.zpost <- log(zz.mode.post)


## Compute prior at mode (using a multinomial-Dirichlet distribution)
# Assuming Dirichlet prior with alpha_i = 1
# https://en.wikipedia.org/wiki/Dirichlet-multinomial_distribution#Conditional_distribution

#Compute number of observations of each group at the mode
res.mode <- rep(0, n.grp)
names(res.mode) <- as.character(1:n.grp)

tab.mode <- table(z)
res.mode[names(tab.mode)] <- tab.mode

#Parameters in the m-D distribution
A <- n.grp * 1
N <- length(z) #Number of observations
#Original equation
#log.zprior <- lgamma(A) - lgamma(A +  N) + sum(lgamma(res.mode + 1)) - n.grp * lgamma(1)
#Corrected equation
log.zprior <- lgamma(A) + lgamma(N) - lgamma(A +  N)


#Compute marginal likelihood
mlik <- m.mode$model$mlik[1, 1] + log.zprior - log.zpost
return(mlik)

}


