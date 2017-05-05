# Maximum likelihood estimation of power law coefficients
# First try to recreate pareto, then try to do it with a cutoff

library(actuar)
library(stats4)

# Load some of the density data from bci to test.
source('~/GitHub/forestlight/code/load50ha.r')
alltreedat <- subset(bcicensusdat, !is.na(dbh) & production67 > 0 & comp_idx > 0)

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

# Try Farrior's breakpoint analysis. Ignore normalization constant since it is not needed to minimize log likelihood.
nll_powerlaw_breakpoint <- function(alpha, beta, breakpoint) {
  if (breakpoint < min(x) | breakpoint > max(x)) return(1e9) # Penalize if breakpoint is invalid.
  smalltrees <- x[x < breakpoint]
  bigtrees <- x[x >= breakpoint]
  
  #C <- (alpha - 1) * ( xmin ^ (alpha - 1) )
  #C <- (1/L) / (expint::gammainc(1-alpha, xmin/L))
  C <- 1
  fx_small <- smalltrees ^ ( -alpha )
  fx_big <- (breakpoint ^ -alpha) * (exp(-bigtrees * beta) / exp(-breakpoint * beta))
  px <- C * c(fx_small, fx_big)
  -sum(log(px))
}


x <- alltreedat$dbh
fit1 <- mle(nll_powerlaw, start = list(alpha=3), fixed = list(xmin=min(x)), method='BFGS')
fit2 <- mle(nll_powerlaw_cutoff2, start = list(alpha=3, L=1), fixed = list(xmin=min(x)), method='BFGS')
# Modify fit3 to let alpha vary.
fit3 <- mle(nll_powerlaw_breakpoint, start = list(beta=1, breakpoint=10), fixed = list(alpha=as.numeric(fit1@coef[1])), method='BFGS')
fit3 <- mle(nll_powerlaw_breakpoint, start = list(beta=2, breakpoint=10, alpha=2), method='BFGS')


x <- subset(alltreedat, tol_wright=='S')$dbh
fit1shade <- mle(nll_powerlaw, start = list(alpha=3), fixed = list(xmin=min(x)), method='BFGS')
fit2shade <- mle(nll_powerlaw_cutoff2, start = list(alpha=3, L=1), fixed = list(xmin=min(x)), method='BFGS')
fit3shade <- mle(nll_powerlaw_breakpoint, start = list(beta=1, breakpoint=10), fixed = list(alpha=as.numeric(fit1shade@coef[1])), method='BFGS')

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

# Bootstrap CIs

boot_mle <- function(xboot, nboot) {
  boot_stats <- list()
  pb <- txtProgressBar(0, nboot, style = 3)
  for (i in 1:nboot) {
    setTxtProgressBar(pb, i)
    x <<- sample(xboot, size = length(xboot), replace = TRUE) # Sets x in global environment.
    fit1_i <- mle(nll_powerlaw, start = list(alpha=3), fixed = list(xmin=min(x)), method='BFGS')
    fit2_i <- mle(nll_powerlaw_cutoff2, start = list(alpha=3, L=1), fixed = list(xmin=min(x)), method='BFGS')
    boot_stats[[i]] <- c(slope1=fit1_i@coef[1], slope2=fit2_i@coef[1], cutoff=fit2_i@coef[2])
  }
  close(pb)
  do.call('rbind', boot_stats)
}

boot_shade <- boot_mle(xboot = subset(alltreedat, tol_wright=='S')$dbh, nboot = 99)
boot_int <- boot_mle(xboot = subset(alltreedat, tol_wright=='I')$dbh, nboot = 99)
boot_gap <- boot_mle(xboot = subset(alltreedat, tol_wright=='G')$dbh, nboot = 99)
boot_all <- boot_mle(xboot = alltreedat$dbh, nboot = 99)

shade_ci <- apply(boot_shade, 2, quantile, probs = c(0.025, 0.975))
int_ci <- apply(boot_int, 2, quantile, probs = c(0.025, 0.975))
gap_ci <- apply(boot_gap, 2, quantile, probs = c(0.025, 0.975))
all_ci <- apply(boot_all, 2, quantile, probs = c(0.025, 0.975))

# Make sure confidence interval has converged
nboot <- 99
bootgapci <- matrix(0, ncol=2, nrow=length(10:nboot))
for (i in 10:nboot) {
  bootgapci[i-9, ] <- quantile(boot_gap[1:i, 3], probs=c(0.025,0.975))
}

plot(1:90, bootgapci[,2], type='l', ylim=c(0,70))
lines(1:90, bootgapci[,1]) # it converges.

save(boot_shade, boot_int, boot_gap, boot_all, file = 'C:/Users/Q/Dropbox/projects/forestlight/bootout.r')

# Plot the new slopes.
slopedat <- expand.grid(parameter = c('slope1','slope2','cutoff'), guild = c('all','shade','intermediate','gap'))
slopedat$cimin <-  c(all_ci[1,], shade_ci[1,], int_ci[1,], gap_ci[1,])
slopedat$cimax <-  c(all_ci[2,], shade_ci[2,], int_ci[2,], gap_ci[2,])
slopedat$est <- c(fit1@coef, fit2@coef, fit1shade@coef, fit2shade@coef, fit1int@coef, fit2int@coef, fit1gap@coef, fit2gap@coef)

library(cowplot)

# Plot parameter estimates and confidence intervals
dummy <- data.frame(expand.grid(parameter = c('slope1','slope2','cutoff'), guild = c('all','shade','intermediate','gap')),
                    est=c(0,2),cimin=0,cimax=0)

ggplot(slopedat, aes(x = guild, y = est, ymin = cimin, ymax = cimax)) +
  geom_pointrange() +
  geom_blank(data=dummy) +
  facet_wrap(~ parameter, scales='free_y', nrow=2) +
  panel_border(colour = 'black') +
  theme(strip.background = element_blank()) +
  labs(y = 'parameter estimate')

ggsave('C:/Users/Q/Google Drive/ForestLight/figs/bootstrapci_cutoffs.png', height=7, width=7, dpi=300)
