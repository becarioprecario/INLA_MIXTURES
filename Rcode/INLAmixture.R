#d: data to be passed to inla() in a suitable format
#fit.inla: function that fits a model given d and a. Takes two arguments: 'data'
#  for the data and 'z' for the group indices. Returns a list with 'mlik' and
#  'model'.
#z.init: Initial values of the group indices.
#z.sample. function to sample new values of the group indices
#dirichlet.prior: Paramters for the Dirichlet prior of group weights. It should
#be a vector of length the number of components in the mixture.
#n.sim: Total number of simulations
#n.burin: Simulations for burn in
#n.thin: Thinning AFTER burnin
#n.errors: Number of total INLA errors allowed
#verbose: Whether to plot some information messages (default to FALSE).

INLAmixture <- function(d, fit.inla, z.init, z.sample, dirichlet.prior, 
  n.sim = 200, 
  n.burnin = 100, n.thin = 1, n.errors = 20, verbose = FALSE) {


z.sim <- as.list(rep(NA, n.sim))
model.sim <- as.list(rep(NA, n.sim))

#Initial values
#b.sim[[1]] <- b.init
#model.sim[[1]] <- fit.inla(d, b.sim[[1]])

#Number of groups
n.grp <- length(unique(z.init))


z.cur <- z.init
model.cur <- fit.inla(d, z.cur, n.grp)

#Index for saved iterations
save.idx <- 0

#Total number of simulations
n.sim.tot <- n.burnin + n.thin * n.sim

#Set number of errors to zero
n.err.idx <- 0

i <- 1
while (i <= n.sim.tot) {

   #Sample new proposal
   z.new <- z.sample(z.cur, model.cur$model, n.grp)
   #Fit model using try to handle possible errors in INLA
   model.new <- try(fit.inla(d, z.new, n.grp))

   if(class(model.new) == "try-error") {
     #Update number of errors
     n.err.idx <- n.err.idx + 1

     if(verbose) {
        cat("INLA error number ", n.err.idx, " at iteration ", i, ".\n")
     }

     if(n.err.idx > n.errors) {
       stop("INLA failed too many times to fit a model.")
     }


     #Set current to random previous state (it is from posterior)
     if(i > n.burnin & save.idx > 100) {
       idx.aux <- sample(save.idx - 1, 1) #i - 1 to avoid previous state
       b.cur <- b.sim[[idx.aux]]
       model.cur <- model.sim[[idx.aux]]
     } else {
       stop("INLA failed early in the sampling process.")
     }

     #DO NOT incrase 'i'

  } else {

      z.cur <- z.new
      model.cur <- model.new
   }

   #Save results (if needed)
   if(i > n.burnin & (i - n.burnin + 1 ) %% n.thin == 0) {
     save.idx <- save.idx + 1
     z.sim[[save.idx]] <- z.cur
     model.sim[[save.idx]] <- model.cur
   }

   if(verbose) {
    if(i %% 100 ==0) {cat("Iteration ", i, "completed. Time:", 
      as.character(Sys.time()), "\n")}
   }

   #INCREASE 'i'
   i <- i + 1 
  }#Try-error

  res <- list(model.sim = model.sim,
    z.sim = z.sim)

  return(res)
}


#Utility function to create a data in the INLA format to use several
#likelihoods
#y: Vector of data
#z: Vector of indices. If NULL observations are allocated at random.
#n.grp: Number of components
data.mixture <- function(y, z = NULL, n.grp) {

  #Number of data
  n <- length(y)

  #Group indices
  if(is.null(z)) {
    z <- sample(rep(1:n.grp, length.out = n))
  }

  #Data in n.grp-column format
  yy <- matrix(NA, ncol = n.grp, nrow = n)
  for(i in 1:n.grp) {
    idx <- which(z == i)
    yy[idx, i] <- y[idx]
  }

  #x stores the intercept in the model
  x <- yy
  x[!is.na(x)] <- 1

  d <- list(y = yy, x = x, z = z)
}
