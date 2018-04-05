#Code to test how the marg. lik. is computed using the
#alternative methods of sum pi(y|z) pi(z)

#Load results
n.grp <- 4
load(paste0("results-ngrp_z_equalvars-", n.grp, ".RData"))

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

