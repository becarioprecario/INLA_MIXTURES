#Show results for a certain model

n.grp <- 3

#Load MCMC results
load(paste0("jags-gaussian-ngrp-", n.grp, ".RData"))
#Load INLA results
load(paste0("results-gaussian-ngrp_z-", n.grp, ".RData"))

pdf(file = paste0("gaussian-ngrp-", i, ".pdf" ))
par(mfcol = c(2, n.grp))
for(i in 1:n.grp) {
  plot(marg.mean[[i]], type = "l", main = bquote(mu[.(i)]), xlab = "",
    ylab = "density" )
  lines(density(res.jags$mu[i,,]), lty = 2)
  legend("topleft", legend = c("INLA", "MCMC"), lty = 1:2, bty = "n")
  plot(marg.prec[[i]], type = "l", main = bquote(tau[.(i)]), xlab = "",
    ylab = "density" )
  lines(density(res.jags$prec[i,,]), lty = 2)
  legend("topright", legend = c("INLA", "MCMC"), lty = 1:2, bty = "n")
}
dev.off()


