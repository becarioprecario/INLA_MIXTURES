#Show results for a certain model

n.grp <- 3

#Load MCMC results
load(paste0("jags-poisson-ngrp-", n.grp, ".RData"))
#Load INLA results
load(paste0("results-poisson-ngrp_z-", n.grp, ".RData"))

pdf(file = paste0("poisson-ngrp-", i, ".pdf" ), height = 3, width = 9)
par(mfcol = c(1, n.grp))
for(i in 1:n.grp) {
  plot(marg.mean[[i]], type = "l", main = bquote(lambda[.(i)]), xlab = "",
    ylab = "density" )
  lines(density(res.jags$mu[i,,]), lty = 2)
  legend("topleft", legend = c("INLA", "MCMC"), lty = 1:2, bty = "n")
}
dev.off()


