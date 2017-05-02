# Maximum likelihood estimation of power law coefficients
# First try to recreate pareto, then try to do it with a cutoff

library(actuar)
library(stats4)

# Must constrain parameter space
nll <- function(shape_par, min_par) {
  #if (min_par < min(x)) 
    -sum(dpareto1(x, shape_par, min_par, log = TRUE)) 
  #else 
  #  NA
  
}

true_shape <- 2
true_min <- 1

set.seed(30303)
x <- rpareto1(n=100000, shape=true_shape, min=true_min)

mle(minuslogl=nll, start=list(shape_par=3), fixed = list(min_par = 1), method='BFGS')


######################################
# Now try to do an exponential piecewise function!


# Load some of the density data from bci to test.
source('~/GitHub/forestlight/code/load50ha.r')
alltreedat <- subset(bcicensusdat, !is.na(dbh) & production67 > 0 & comp_idx > 0)
x <- alltreedat$dbh

nll <- function(shape_par, min_par) {
  -sum(dpareto1(x, shape_par-1, min_par, log = TRUE)) 
}  

# Try to create the pareto manually to make sure it works.
nll_pareto_manual <- function(shape_par, min_par)
  -sum(log(shape_par * min_par ^ shape_par / (x ^ (shape_par + 1))))

fit1 <- mle(nll, start=list(shape_par=3), fixed=list(min_par=min(x)), method='BFGS') # this gives the same result as poweRlaw.
fit1 <- mle(nll_pareto_manual, start=list(shape_par=3), fixed=list(min_par=min(x)), method='BFGS') # this gives the same result as poweRlaw.


bigtreedist <- function(x, breakpoint, alpha, beta, min_par) {
  (exp(-x * beta)/exp(-breakpoint * beta)) * (alpha * min_par ^ alpha * breakpoint ^ -(alpha + 1))
  #(alpha * min_par ^ alpha / (breakpoint ^ (alpha + 1)))/exp(-breakpoint*beta)) * exp(-x * beta)
}

bigtreedist <- function(x, breakpoint, alpha, beta, min_par)
  (x ^ (-alpha - 1)) * (exp(-x * beta)/exp(-breakpoint * beta))

nll_cutoff <- function(alpha, beta, min_par, breakpoint) {
  if (breakpoint <= min(x) | breakpoint >= max(x) | beta > 0) return(1e9)
  smalltrees <- x[x < breakpoint]
  bigtrees <- x[x >= breakpoint]
  
  nll_smalltrees <- -sum(log(alpha * min_par ^ alpha / (x ^ (alpha + 1))))
  nll_bigtrees <- -sum(log(bigtreedist(bigtrees, breakpoint, alpha, beta, min_par)))
  nll_smalltrees + nll_bigtrees
}

nll_powerlaw <- function(alpha, xmin) {
  C <- (alpha - 1) * ( xmin ^ (alpha - 1) )
  fx <- x ^ ( -alpha )
  px <- C * fx
  -sum(log(px))
} 

nll_powerlaw_cutoff <- function(alpha, xmin, lambda) {
  C <- lambda ^ (1 - alpha) / ( gamma(1 - alpha, lambda * xmin) )
  fx <- x ^ ( -alpha ) * exp(-lambda * x)
  px <- C * fx
  -sum(log(px))
} 

# See: http://www.stat.cmu.edu/~cshalizi/2010-10-18-Meetup.pdf
nll_powerlaw_cutoff2 <- function(alpha, xmin, L) {
  C <- (1/L) / (expint::gammainc(1-alpha, xmin/L))
  fx <- ( (x/L)^ -alpha ) * exp(-x/L)
  px <- C * fx
  -sum(log(px))
}

nll_powerlaw_breakpoint <- function(alpha, xmin, beta, L) {
  if (L > max(x) | L < min(x)) return(1e9)
  smalltrees <- x[x < L]
  bigtrees <- x[x >= L]
  
  #C <- (alpha - 1) * ( xmin ^ (alpha - 1) )
  C <- (1/L) / (expint::gammainc(1-alpha, xmin/L))
  fx_small <- x ^ ( -alpha )
  fx_big <- ((L^-alpha)/exp(-L*beta)) * exp(-x*beta)
  px <- C * c(fx_small, fx_big)
  -sum(log(px))
}

x <- alltreedat$dbh
fit1 <- mle(nll_powerlaw, start = list(alpha=3), fixed = list(xmin=min(x)), method='BFGS')
#fit2 <- mle(nll_powerlaw_cutoff, start = list(alpha=3, lambda=1), fixed = list(xmin=min(x)), method='BFGS')
fit2 <- mle(nll_powerlaw_cutoff2, start = list(alpha=3, L=1), fixed = list(xmin=min(x)), method='BFGS')
fit3 <- mle(nll_powerlaw_breakpoint, start = list(L=2, beta=1), fixed = list(alpha=as.numeric(fit1@coef[1]), xmin=min(x)), method='BFGS')

x <- subset(alltreedat, tol_wright=='S')$dbh
fit1shade <- mle(nll_powerlaw, start = list(alpha=3), fixed = list(xmin=min(x)), method='BFGS')
fit2shade <- mle(nll_powerlaw_cutoff2, start = list(alpha=3, L=1), fixed = list(xmin=min(x)), method='BFGS')
fit3shade <- mle(nll_powerlaw_breakpoint, start = list(L=2, beta=1), fixed = list(alpha=as.numeric(fit1shade@coef[1]), xmin=min(x)), method='BFGS')


x <- subset(alltreedat, tol_wright=='I')$dbh
fit1int <- mle(nll_powerlaw, start = list(alpha=3), fixed = list(xmin=min(x)), method='BFGS')
fit2int <- mle(nll_powerlaw_cutoff2, start = list(alpha=3, L=1), fixed = list(xmin=min(x)), method='BFGS')

x <- subset(alltreedat, tol_wright=='G')$dbh
fit1gap <- mle(nll_powerlaw, start = list(alpha=3), fixed = list(xmin=min(x)), method='BFGS')
fit2gap <- mle(nll_powerlaw_cutoff2, start = list(alpha=3, L=1), fixed = list(xmin=min(x)), method='BFGS')

# Output coefficients
coefdat <- data.frame(guild = c('all','shade-tolerant','intermediate','gap'),
                      slope1 = c(fit1@coef, fit1shade@coef, fit1int@coef, fit1gap@coef),
                      slope2 = c(fit2@coef[1], fit2shade@coef[1], fit2int@coef[1], fit2gap@coef[1]),
                      cutoff = c(fit2@coef[2], fit2shade@coef[2], fit2int@coef[2], fit2gap@coef[2]))
coefdat$logcutoff <- log10(coefdat$cutoff)

##########################################
# crap? below
fit1 <- mle(nll_pareto_manual, start=list(shape_par=3), fixed=list(min_par=min(x)), method='BFGS') # this gives the same result as poweRlaw.
fit2 <- mle(nll_cutoff, start=list(beta=-1, breakpoint = 1.5), fixed=list(alpha=fit1@coef, min_par=min(x)), method='BFGS')
fit2b <- mle(nll_cutoff2, start=list(beta=-1, alpha=3), fixed = list(min_par=min(x)), method='BFGS')


x <- subset(alltreedat, tol_wright=='S')$dbh
fit1shade <- mle(nll, start=list(shape_par=3), fixed=list(min_par=min(x)), method='BFGS') # this gives the same result as poweRlaw.
fit2shade <- mle(nll_cutoff, start=list(beta=2, breakpoint = 5), fixed=list(alpha=fit1shade@coef, min_par=1), method='BFGS')


x <- subset(alltreedat, tol_wright=='G')$dbh
fit1gap <- mle(nll, start=list(shape_par=3), fixed=list(min_par=min(x)), method='BFGS') # this gives the same result as poweRlaw.
fit2gap <- mle(nll_cutoff, start=list(beta=1, breakpoint = 5), fixed=list(alpha=fit1gap@coef, min_par=1), method='BFGS')
