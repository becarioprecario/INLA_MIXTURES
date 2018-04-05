#Select unique configurations, index in original list and proportion
# zz: list of vector with configurations
# mliks: vector of marginal likelihoods extracted from models
# alpha: Parameters for Dirichlet prior
z_unique <- function (zz, mliks, alpha) {
  if(!is.list(zz)) {
    stop("zz must be a list.")
  }

  n <- length(zz)
  zz.u <- unique(zz)

  #Index in original vector
  zz.idx <- lapply(zz.u, function(X) {
    which(unlist(lapply(zz, function(Y) { (as.vector(identical(X, Y)))})))
  })

  #Posterior probabilities using sampling from MGS
  zz.prob <- unlist(lapply(zz.idx, length)) / n

  #Unique mliks
  zz.mliks <- mliks[unlist(lapply(zz.idx, function(X){X[1]}))]

  #Compute pi(z)
  zz.prior <- unlist(lapply(zz.u, function(X) {
    prior.z(X, alpha)
  }))

  #Posterior probabilities using marg. lik.
  mliks.pprob <- exp(zz.mliks - max(zz.mliks) + zz.prior - max(zz.prior))
  mliks.pprob <- mliks.pprob / sum(mliks.pprob)

  return(list(unique = zz.u, index = zz.idx, pprob.samp = zz.prob, 
    prior = zz.prior, pprob.mliks = mliks.pprob, mliks = zz.mliks
    ))
}

