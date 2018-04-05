# Simulation study

This directory contains the data to reproduce the simulation study described
in the paper. Files included here are:

* `gaussian_data.R`

  Simulation of Gaussian data.

* `poisson_data.R`

  Simulation of Poisson data.

* `mixtures_gaussian_z.R`

  Code to fit a mixture of Gaussian distributions. The number of components
  in the mixture must be set in variable `n.grp`.

* `mixtures_poisson_z.R`

  Code to fit a mixture of Poisson distributions. The number of components
  in the mixture must be set in variable `n.grp`.

* `jags_gaussian.R`

  Code to fit a mixture of Gaussian distributions with JAGS.

* `jags_poisson.R`

  Code to fit a mixture of Poisson distributions with JAGS.

* `show_results.R`

  Display marginals of parameters in the mixture for the Gaussian data example.
  The number of components in the mixture must be set in variable `n.grp`.

* `show_results-pois.R`

  Display marginals of parameters in the mixture for the Poisson data example.
  The number of components in the mixture must be set in variable `n.grp`.

* `test_mlik.R`

  Compute marginal likelihood for mixture models using the different methods
  described in the paper.

* `compute_true_mlik.R `

  Compute the marginal likelihood for the model with a single component and the
  actual model (i.e., three components and correct assignment of observations to
  components) for the Gaussian and Poisson example. This uses `R` package
  `MCMCpack`.

**Warning:** This code may take a while to run. In particular, mixture models
with `INLA` can take some time. Values `n.sim` and `n.burnin` can be reduced in
the call to function `INLAmixture` in order to test the `R` code.

Also, this code can generate **huge** `RData` files. In particular, the size
of the output files can be of several Gigabytes.


