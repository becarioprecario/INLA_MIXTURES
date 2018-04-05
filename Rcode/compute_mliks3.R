#Estimate the log-marginal likelihood using 3 different criteria:
# pi_I, pi_G, pi_M

#zz.unique: Unique structure of auxiliary variables; output from z_unique().

compute_mliks3 <- function(zz.unique) {

  mliks <- zz.unique$mliks
  prior <- zz.unique$prior

  max.mliks <- max(mliks)
  max.prior <- max(prior)

  #Compute log(pi(y)) using sum_ pi(y|z) * pi(z)
  pi_I <- log(sum(exp(mliks - max.mliks) * exp(prior - max.prior))) + 
    max.mliks + max.prior

  #Reported log(pi(y)) using Bayes rule and MGS probabilities
  idx <- which.max(zz.unique$pprob.samp)
  pi_G <- mliks[idx] + prior[idx] - zz.unique$pprob.samp[idx]


  #Compute log(pi(y)) with Bayes Rule with corrected posterior probabilities
  idx2 <- which.max(zz.unique$pprob.mliks)
  pi_M <- mliks[idx2] + (prior[idx2]) - zz.unique$pprob.mliks[idx2]

  #Compute criterion by Raftery et al. from MGS
  aux <- sum(exp(-(mliks - max.mliks)) * zz.unique$pprob.samp)
  pi_RG <- -log(aux) + max.mliks

  #Compute criterion by Raftery et al. from marg. liks.
  aux <- sum(exp(-(mliks - max.mliks)) * zz.unique$pprob.mliks)
  pi_RM <- -log(aux) + max.mliks

  return (list(log_pi_I = pi_I, log_pi_G = pi_G, log_pi_M = pi_M,
    log_pi_RG = pi_RG, log_pi_RM = pi_RM))
}
