#Fits a mixture model given a vector of indices

# ------ GAUSSIAN mixture
fitmix.inla.gaussian <- function(...) {
  res <- fitmix.inla(..., likelihood = "gaussian")
  return(res)
}

# ------ POISSON mixture
fitmix.inla.poisson <- function(...) {
  res <- fitmix.inla(..., likelihood = "poisson")
  return(res)
}

#d: Data in 'mixture format'
#z: Vector of indices
#n.grp: Number of groups in the mixture
fitmix.inla <- function(d, z, n.grp, likelihood = "gaussian") {


  #FIXME: Pass these as extra arguments
  #Precision priors
  if(likelihood %in% c("gaussian")) {
    prec.prior1 <- list(list(hyper = list(prec = list(prior = "loggamma",
    #  param = c(0.38, 0.5)))))
    param = c(0.5, 0.5)))))
    #  param = c(6/2, 40/2))))) #Chi (1995)
    prec.prior <- replicate(n.grp, prec.prior1)
  }

  #Prior on the mean effects 
  c.fixed <- list(mean.intercept = 0, prec.intercept = 0.001,
    mean = 0, prec = 0.001)

  yy <- apply(d$y, 1, function(X){X[!is.na(X)]})

  d$y[,] <- NA
  for(i in 1:n.grp) {
    idx <- which(z == i)
    d$y[idx, i] <- yy[idx]
  }

  d$x[,] <- d$y
  d$x[!is.na(d$x)] <- 1

  m <- switch(likelihood, 
    gaussian =  inla(y ~ -1 + x, data = d, family = rep(likelihood, n.grp),
    control.family = prec.prior, control.fixed = c.fixed),
    poisson =  inla(y ~ -1 + x, data = d, family = rep(likelihood, n.grp),
    control.fixed = c.fixed)
  )
  return(list(mlik = m$mlik[1,1], model = m))
}

