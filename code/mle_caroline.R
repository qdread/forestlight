# MLE with Caroline's code.

source('~/GitHub/FunctionFitting/PowerLawFit_etc_160509.r')

# Try with all trees.

# Load some of the density data from bci to test.
source('~/GitHub/forestlight/code/load50ha.r')
alltreedat <- subset(bcicensusdat, !is.na(dbh) & production67 > 0 & comp_idx > 0)

x <- alltreedat$dbh
# This does not seem to work. 
powerlaw_fit(x = round(x), xmin = 1, interval = c(1,3), plot = TRUE, lines = TRUE) 

x <- subset(alltreedat, tol_wright=='S')$dbh
fit1shade <- mle(nll_powerlaw, start = list(alpha=3), fixed = list(xmin=min(x)), method='BFGS')

powerexp_fit(x, xmin1=1, alpha1=fit1shade@coef['alpha'], plotting=T)
