stanmodel_paretoxpowerll <- stan_model(file = 'stan/model_ppow_withlik.stan', model_name = 'paretoxpow')
stanmodel_weibullxpowerll <- stan_model(file = 'stan/model_wpow_withlik.stan', model_name = 'weibullxpow')

NC <- 3
NI <- 6000
NW <- 5000

fit_ppow_all <- sampling(stanmodel_paretoxpowerll, data = data1995_alltree, chains = NC, iter = NI, warmup = NW)

fit_wpow_all <- sampling(stanmodel_weibullxpowerll, data = data1995_alltree, chains = NC, iter = NI, warmup = NW)

library(loo)

ppow_lldens <- extract_log_lik(fit_ppow_all, 'log_lik_dens')
waic(ppow_lldens)
loo(ppow_lldens)

ppow_llprod <- extract_log_lik(fit_ppow_all, 'log_lik_prod')
waic(ppow_llprod)
loo(ppow_llprod)

wpow_lldens <- extract_log_lik(fit_wpow_all, 'log_lik_dens')
waic(wpow_lldens)
loo(wpow_lldens)

wpow_llprod <- extract_log_lik(fit_wpow_all, 'log_lik_prod')
waic(wpow_llprod)
loo(wpow_llprod)
