# Example on the earthquakes data

This directory includes the `R` code to fit a mixture model to the earthquakes
data. Files included are:

* `mixtures_earthquakes.R`

  Fit mixture model (with Poisson components) to the data. Variable `n.grp`
  sets the number of Poisson components in the mixture.

* `test_mlik.R`

  Compute marg. lik. of the fitted models using the methods described in the
  paper.  Variable `n.grp` sets the number of components in the mixture.

