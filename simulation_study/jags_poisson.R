## Example: Fitting mixtures with R-INLA
#
# Samples from z using a EM-like approach

library(rjags)

#Get data
load("poisson_data.RData")

#Number of groups
#n.grp <- 3 

for(n.grp in 2:5) {


#Create jags data
jags.data <- list(y = d$y, N = length(d$y),
  K = n.grp, alpha = rep(1, n.grp),
  mean.grp = 0, prec.grp = 0.001)

jags.inits <- list(Z = sample(1:n.grp, jags.data$N, TRUE)
  #, sigma = rep(1, n.grp)
)

# This model uses Gamma priors on the precisions and
# the results match thoser obtained with INLA
m.jags <- jags.model("../Rcode/mixture_poisson.bug", jags.data, jags.inits)
update(m.jags, 200)
res.jags <- jags.samples(m.jags, c("mu", "Z", "w"),
  n.iter = 10000, thin = 10 )

res.jags

save(file = paste0("jags-poisson-ngrp-", n.grp, ".RData"), list = ls())
}

