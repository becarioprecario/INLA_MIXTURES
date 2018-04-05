# Install required packages to run the R code

# Set CRAN mirror (0-Cloud [https])
chooseCRANmirror(ind = 1)

# Main packages
pkgs <- c("MCMCpack", "MixtureInf", "INLABMA")
install.packages(pkgs, dep = TRUE)


# Install JAGS
warning("This may require the installation of the main JAGS software")
install.packages("rjags", dep = TRUE)

# Install INLA (stable version)
install.packages("INLA", repos = c(getOption("repos"),
  INLA = "https://inla.r-inla-download.org/R/stable"), dep = TRUE)
