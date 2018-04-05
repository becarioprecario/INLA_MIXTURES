#Simulate Gaussian data

n.grp <- 3
n.means <- c(0, 5, 10)
s.size <- 50

#Sampel data
set.seed(1)
yy <- lapply(n.means, function(X) {rnorm(s.size, mean = X, sd = 1)})
d <- data.frame(y = unlist(yy), grp = rep(1:n.grp, each = s.size))

pdf(file = "gaussian_data.pdf")
plot(density(d$y), main = "Gaussian data", xlab = "")
dev.off()

save(file = "gaussian_data.RData", list = ls())

