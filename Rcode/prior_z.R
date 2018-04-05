#Compute the prior of a configuration of the assignement assuming a
# Dirichlet prior
#
#
#z: Configuration (vector of values from 1 to K)
#alpha: Paramters of Dirichlet prior (vector of length K)
#log: Whether compute this in the log scale

prior.z <- function(z, alpha, log = TRUE) {
  #Number of components
  K <- length(alpha)
  n <- length(z)

  z.tab <- sapply(1:K, function(X) {sum(z == X)})

  #Compute table with number of observatios per group

  val <- lgamma(sum(alpha)) - sum(lgamma(alpha))
  val <- val + sum(lgamma(z.tab + alpha)) - lgamma(n + sum(alpha))

  if(!log) {
   val <- exp(val)
  } 

  return(val)
}

