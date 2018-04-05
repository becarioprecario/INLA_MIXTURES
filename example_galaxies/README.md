# Example on the galaxies data

This directory includes the `R` code to fit a mixture model to the galaxy
data. FIles included are:

* `mixtures_galaxies_z_equalvars.R`

  Fit mixture model (with Gaussian components) to the data. Variable `n.grp`
  sets the number of Gaussian components in the mixture.

* `test_mlik.R`

  Compute marg. lik. of the fitted models using the methods described in the
  paper.  Variable `n.grp` sets the number of components in the mixture.
