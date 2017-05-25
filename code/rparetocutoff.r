# Create the cutoff distributions from the pareto_cutoff density function 
# Sample the same number of individuals from gap, intermediate and shade distributions.
# See whether we get as many big shade trees and few big gap trees when we sample the same large number of individuals from each of the distributions.

pareto_cutoff <- function(x, alpha, xmin, L) {
  C <- (1/L) / (expint::gammainc(1-alpha, xmin/L))
  fx <- ( (x/L)^ -alpha ) * exp(-x/L)
  px <- C * fx
  return(px)
}

rparetocutoff <- function(n, alpha, xmin, L) {
  require(distr)
  mydist <- AbscontDistribution(d = function(x) pareto_cutoff(x, alpha=alpha, xmin=xmin, L=L), low1=xmin)
  r(mydist)(n)
}

