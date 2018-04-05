#Functions to sample the indices in a mixture for different
#  types of mixtures

# ------ GAUSSIAN MIXTURE
z.sample.gaussian <- function(...) {
  res <- z.sample(..., likelihood = "gaussian")
  return(res)
}


# ---------- POISSON mixture
z.sample.poisson <- function(...) {
  res <- z.sample(..., likelihood = "poisson")
  return(res)
}

# -------- GENERIC z.sample
#z: Vector of group indices
#model: Model fited using 'z'
#n.grp: Number of components in the mixture
#d.param: Parameters on the Dirichlet prior. Must be of length n.grp
z.sample <- function(z, model, n.grp, d.param = NULL, likelihood = "gaussian") {

  #Get dataset
  d <- model$.args$data

  #Get response
  y.orig <- apply(d$y, 1, function(X) {X[!is.na(X)]})
  n <- length(y.orig)

  #Get groups
  z.grp <- apply(d$y, 1, function(X) {which(!is.na(X))})
  tab.grp <- table(z.grp)

  #Fix groups with zero observations
  if(length(tab.grp) < n.grp) {
    tab.grp <- tab.grp[as.character(1:n.grp)]
    names(tab.grp) <- as.character(1:n.grp)
    tab.grp[is.na(tab.grp)] <- 0
  }

  #Parameters of the Dirchlet prior
  if(is.null(d.param)) {
    d.param <- rep(1, n.grp)
  }
  d.param[as.numeric(names(tab.grp))] <- 
    d.param[as.numeric(names(tab.grp))] +  tab.grp

  w.grp <- d.param/sum(d.param)  #tab.grp/n  #rdirichlet(1, d.param)

  #Get means of mixtures
  means <- model$summary.fixed[, "mean"]
  #Reorder means
  idx <- order(means)
  means <- means[idx]

  #If Gaussian keep precissions
  if(likelihood %in% c("gaussian")) {
    precs <- model$summary.hyperpar[, "mean"]
    precs <- precs[idx]
  }


  #Probability of being in each group
  prob.z <- switch(likelihood,
    gaussian = sapply(1:n.grp, function(X) {
    w.grp[X] * dnorm(y.orig, mean = means[X], sd = 1/sqrt(precs[X]))
  }),
    poisson = sapply(1:n.grp, function(X) {
      w.grp[X] * dpois(y.orig, lambda = exp(means[X])) })
  )

  prob.z <- t(apply(prob.z, 1, function(X) { X / sum(X)}))

  grp <- apply(prob.z, 1, function(X) {sample(1:n.grp, 1, prob = X)})

  return(grp)
}

