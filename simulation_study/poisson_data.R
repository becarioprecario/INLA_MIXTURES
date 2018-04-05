#Simulate Gaussian data

n.grp <- 3
n.means <- c(1, 15, 45)
s.size <- 50

#Sampel data
set.seed(1)
yy <- lapply(n.means, function(X) {rpois(s.size, lambda = X)})
d <- data.frame(y = unlist(yy), grp = rep(1:n.grp, each = s.size))

pdf(file = "poisson_data.pdf")
hist(d$y, main = "Poisson data", xlab = "")
dev.off()

save(file = "poisson_data.RData", list = ls())

