# BOOTSTRAP CI IN PARALLEL (forest light)
# Run locally :)

setwd('~/forestlight')
load('line353.RData')

library(stats4)
library(dplyr)
library(Hmisc)

library(foreach)
library(doParallel)

registerDoParallel(cores = 10)

nb <- 99
qprobs <- c(0.025, 0.975)

allyears_paretoboot <- foreach(i = allyears_names) %dopar% {
	dat <- get(i)
	lapply(dat, function(x) boot_mle(xboot = x$dbh_corr, nboot = nb, L_init = 1))
}

for (i in 1:length(allyears_names)) {
	assign(paste0(allyears_names[i], '_paretoboot'), allyears_paretoboot[[i]])
	assign(paste0(allyears_names[i], '_paretobootci'), lapply(get(paste0(allyears_names[i], '_paretoboot')), function(x) apply(x, 2, quantile, probs = qprobs)))
} 

save.image('tempsave.RData')

ninety_paretoboot <- foreach(i = names1990) %dopar% {
	dat <- get(i)
	boot_mle(xboot = dat$dbh_corr, nboot = nb, L_init = 1)
}

for (i in 1:length(names1990)) {
	assign(paste0(names1990[i], '_paretoboot'), ninety_paretoboot[[i]])
	assign(paste0(names1990[i], '_paretobootci'), apply(ninety_paretoboot[[i]], 2, quantile, probs = qprobs))
} 

save.image('tempsave.RData')

ninetyfive_paretoboot <- foreach(i = names1995) %dopar% {
	dat <- get(i)
	boot_mle(xboot = dat$dbh_corr, nboot = nb, L_init = 1)
}

for (i in 1:length(names1995)) {
	assign(paste0(names1995[i], '_paretoboot'), ninetyfive_paretoboot[[i]])
	assign(paste0(names1995[i], '_paretobootci'), apply(ninetyfive_paretoboot[[i]], 2, quantile, probs = qprobs))
} 

save.image('tempsave.RData')



# density scaling by light

ninety_light <- foreach(i = names1990) %dopar% {
	dat <- get(i)
	pareto_cutoff_fits(dat, 'light_received', L_init = 1000)
}

coef_tables_9095light <- list()

for (i in 1:length(names1990)) {
	assign(paste0(names1990[i], '_lightparetofits'), ninety_light[[i]])
	coef_tables_9095light[[length(coef_tables_9095light) + 1]] <- with(get(paste0(names1990[i],'_lightparetofits')), extractcoeffs(fit_pareto, fit_cutoff, aic_pareto, aic_cutoff))
}

ninetyfive_light <- foreach(i = names1995) %dopar% {
	dat <- get(i)
	pareto_cutoff_fits(dat, 'light_received', L_init = 3000)
}

for (i in 1:length(names1995)) {
	assign(paste0(names1995[i], '_lightparetofits'), ninetyfive_light[[i]])
	coef_tables_9095light[[length(coef_tables_9095light) + 1]] <- with(get(paste0(names1995[i],'_lightparetofits')), extractcoeffs(fit_pareto, fit_cutoff, aic_pareto, aic_cutoff))
}

# bootstrap light ci

ninety_lightparetoboot <- foreach (i = names1990) %dopar% {
	dat <- get(i)
	boot_mle(xboot = dat$light_received, nboot = nb, L_init = 1000)
}

for (i in 1:length(names1990)) {
	assign(paste0(names1990[i], '_lightparetoboot'), ninety_lightparetoboot[[i]])
	assign(paste0(names1990[i], '_lightparetobootci'), apply(ninety_lightparetoboot[[i]], 2, quantile, probs = qprobs))
} 

save.image('tempsave.RData')

ninetyfive_lightparetoboot <- foreach (i = names1995) %dopar% {
	dat <- get(i)
	boot_mle(xboot = dat$light_received, nboot = nb, L_init = 3000)
}

for (i in 1:length(names1995)) {
	assign(paste0(names1995[i], '_lightparetoboot'), ninetyfive_lightparetoboot[[i]])
	assign(paste0(names1995[i], '_lightparetobootci'), apply(ninetyfive_lightparetoboot[[i]], 2, quantile, probs = qprobs))
} 

save.image('allscalings_30aug.RData')

