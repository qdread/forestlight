stanmodel_paretoxpowerll <- stan_model(file = 'stan/model_ppow_withlik.stan', model_name = 'paretoxpow')
stanmodel_weibullxpowerll <- stan_model(file = 'stan/model_wpow_withlik.stan', model_name = 'weibullxpow')
stanmodel_paretoxexpll <- stan_model(file = 'stan/model_pexp_withlik.stan', model_name = 'paretoxexp')
stanmodel_weibullxexpll <- stan_model(file = 'stan/model_wexp_withlik.stan', model_name = 'weibullxexp')

mod_par <- stan_model(file = 'stan/model_pareto_withlik.stan')
mod_wei <- stan_model(file = 'stan/model_weibull_withlik.stan')

NC <- 2
NI <- 3000
NW <- 2000
options(mc.cores = NC)

xpred <- seq(1.1, 316, 0.1)
Npred <- length(xpred)

fit_par <- sampling(mod_par, data = c(data1995_alltree, list(N_pred=Npred, x_pred=xpred)), chains = NC, iter = NI, warmup = NW)
fit_wei <- sampling(mod_wei, data = c(data1995_alltree, list(N_pred=Npred, x_pred=xpred)), chains = NC, iter = NI, warmup = NW)

fit_ppow_all <- sampling(stanmodel_paretoxpowerll, data = data1995_alltree, chains = NC, iter = NI, warmup = NW)
fit_wpow_all <- sampling(stanmodel_weibullxpowerll, data = data1995_alltree, chains = NC, iter = NI, warmup = NW)
fit_pexp_all <- sampling(stanmodel_paretoxexpll, data = data1995_alltree, chains = NC, iter = NI, warmup = NW)
fit_wexp_all <- sampling(stanmodel_weibullxexpll, data = data1995_alltree, chains = NC, iter = NI, warmup = NW)

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

pexp_lldens <- extract_log_lik(fit_pexp_all, 'log_lik_dens')
waic(pexp_lldens)
loo(pexp_lldens)

pexp_llprod <- extract_log_lik(fit_pexp_all, 'log_lik_prod')
waic(pexp_llprod)
loo(pexp_llprod)

wexp_lldens <- extract_log_lik(fit_wexp_all, 'log_lik_dens')
waic(wexp_lldens)
loo(wexp_lldens)

wexp_llprod <- extract_log_lik(fit_wexp_all, 'log_lik_prod')
waic(wexp_llprod)
loo(wexp_llprod)

par_ll <- extract_log_lik(fit_par, 'log_lik_dens')
wei_ll <- extract_log_lik(fit_wei, 'log_lik_dens')

waic(par_ll)
waic(wei_ll)
loo(par_ll)
loo(wei_ll)
