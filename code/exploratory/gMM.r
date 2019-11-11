# Fit generalized Michaelis-Menten function to data

gMM <- function(x, a, b, k) (a * x ^ b) / (k + x ^b)

log_gMM <- function(x, a, b, k) a + b * log(x) - log(k + x ^ b)

