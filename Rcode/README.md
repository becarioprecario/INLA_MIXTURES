# `R` code to fit mixture models

This directory contains some utility functions to fit mixture models with
`INLA`. Files included in this directory are:

* `mixture.bug`

  BUGS code to fit a mixture model of Gaussian distributions.

* `mixture_poisson.bug`

  BUGS code to fit a mixture model of Poisson distributions.

* `fitmix_inla.R`

  Fits a mixture model given a vector of indices for Gaussian and
  Poisson mixtures.

* `z_sample.R`

  `R` code to sample new values of the latent variable `z` that assigns
  observations to components in the mixture.

* `z_unique.R`

  Given a list of observations of latent variable `z`, this function
  creates a list of unique values of `z` and associated information
  (frecuency, marg. lik. of associated model, etc.).

* `compute_mliks3.R`

  Estimate the log-marginal likelihood using 3 different criteria discussed
  in the paper: pi_I, pi_G, pi_M

* `prior_z.R`

  Prior distribution for latent variable `z` which assigns observations to
  components.

* `compute_marglik.R`

  Compute marginal likelihood using Bayes' rule.
