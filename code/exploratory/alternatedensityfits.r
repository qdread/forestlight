# Fit only Weibull distribution to density in Stan. Do not mess with the production

# Load data (changing file path if necessary)
fpdata <- 'C:/Users/Q/google_drive/ForestLight/data/data_22jan2018'
load(file.path(fpdata, 'rawdataobj_22jan.r'))

library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = 3) 

prod_dump <- function(dat, to_file = FALSE, fn = NULL, subsample = NULL) {
  require(rstan)
  if (!is.null(subsample) && nrow(dat) > subsample) {
    dat <- dat[sample(nrow(dat), subsample, replace = FALSE), ]
  }
  x <- dat$dbh_corr
  y <- dat$production
  xdat <- list(N = length(x), x = x, y = y, x_min = min(x))
  if (to_file) {
    with(xdat, stan_rdump(names(xdat), file = fn))
  } else {
    return(xdat)
  }
}

data1995_alltree <- prod_dump(alltreedat[[3]], subsample = NULL)

mod_weib <- '
data {
	int<lower=0> N;
  vector<lower=0>[N] x;
}

parameters {
  // Weibull density
  real<lower=0> shape;
  real<lower=0> scale;
}

model {
  // Priors: Weibull density
  shape ~ lognormal(1, 1);
  scale ~ lognormal(1, 1);
  
  // Likelihood: Weibull density
  x ~ weibull(shape, scale);
}
'

stanmod_weib <- stan_model(model_name = 'weibull', model_code = mod_weib)

NC <- 2
NI <- 1000
NW <- 500

fit_weib <- sampling(stanmod_weib, data = data1995_alltree, chains = NC, iter = NI, warmup = NW)

summ_weib <- summary(fit_weib)

library(bayesplot)
mcmc_trace(as.array(fit_weib))

weibull_ci <- function(fit, dbh_pred, n_indiv = 1) {
  pdf_weibull <- function(x, shape, scale) (shape/scale) * (x/scale)^(shape-1) * exp(-(x/scale)^shape)

  pars <- do.call('cbind', extract(fit))
  dens_pred <- sapply(dbh_pred, pdf_weibull, shape = pars[,'shape'], scale = pars[,'scale'])
  dens_pred <- dens_pred * n_indiv

  # Get some quantiles from each of these.
  qprobs <- c(0.025, 0.05, 0.25, 0.5, 0.75, 0.95, 0.975)
  dens_pred_quant <- apply(dens_pred, 2, quantile, probs = qprobs)

  out <- data.frame(dbh = dbh_pred,
                    t(dens_pred_quant))
  
  setNames(out, c('dbh', 'q025', 'q05', 'q25', 'q50', 'q75', 'q95', 'q975'))
  
}

densitybin_5census <- read.csv('C:/Users/Q/google_drive/ForestLight/data/data_22jan2018/densitybin_5census.csv', stringsAsFactors = FALSE)
dbh_pred_bins <- densitybin_5census[densitybin_5census$fg == 'all', 'bin_midpoint']

# Minimum x values for Pareto.
min_n <- read.csv('C:/Users/Q/Dropbox/projects/forestlight/stanoutput/min_n.csv', stringsAsFactors = FALSE)
n1995 <- min_n$n[min_n$fg == 'alltree' & min_n$year == 1995]

weibull_ci(fit = fit_weib, dbh_pred = dbh_pred_bins, n_indiv = n1995) # Gives same very low number as before

pdf_weibull <- function(x, shape, scale) (shape/scale) * (x/scale)^(shape-1) * exp(-(x/scale)^shape)
pdf_weibull(277, .93, 5.71) * n1995


###############################################
# Generalized extreme value

mod_gev <- '
functions {

  real gev_log (real y, real mu, real sigma, real xi){
    real z;
    z = 1 + (y - mu) * xi / sigma;
    return -log(sigma) - (1 + 1/xi)*log(z) - pow(z,-1/xi);  
  }

}

data {
  int<lower=0> N;
  real y[N];
}

transformed data {
  real min_y;
  real max_y;
  min_y = min(y);
  max_y = max(y);
}
parameters {
  real<lower=-0.5, upper=0.5> xi;
  real<lower=0> sigma;
  // location has upper/lower bounds depending on the value of xi
  real<lower=if_else( xi < 0, min_y + sigma / xi, negative_infinity() ),
    upper=if_else( xi > 0, positive_infinity(), max_y + sigma / xi )> mu;
}
model {
  // priors on component parameters
  sigma ~ gamma(.0001,.0001);
  // mu ~ uniform
  xi ~ uniform(-.5,.5);
  
  for(i in 1:N) {
    y[i] ~ gev(mu, sigma, xi);  
  }
}
'
stanmod_gev <- stan_model(model_name = 'gev', model_code = mod_gev)

NC <- 2
NI <- 1000
NW <- 500

dump_gev <- function(dat) {
  list(N = nrow(dat), y = sort(dat$dbh_corr))
}

dat_gev_alltree95 <- dump_gev(alltreedat[[3]])

fit_gev <- sampling(stanmod_gev, data = dat_gev_alltree95, chains = NC, iter = NI, warmup = NW)

summ_gev <- summary(fit_gev)

gev_ci <- function(fit, dbh_pred, n_indiv = 1) {
  pdf_gev <- function(x, mu, sigma, xi) {
    tx <- function(x, mu, sigma, xi) {
      ifelse(xi == 0, exp(-(x - mu)/sigma), (1 + xi * ((x - mu)/sigma))^(-1/xi))
    }
    txx <- tx(x, mu, sigma, xi)
    (1/sigma) * txx^(xi + 1) * exp(-txx)  
  }
  
  pars <- do.call('cbind', extract(fit))
  dens_pred <- sapply(dbh_pred, pdf_gev, mu = pars[,'mu'], sigma = pars[,'sigma'], xi = pars[,'xi'])
  dens_pred <- dens_pred * n_indiv
  
  # Get some quantiles from each of these.
  qprobs <- c(0.025, 0.05, 0.25, 0.5, 0.75, 0.95, 0.975)
  dens_pred_quant <- apply(dens_pred, 2, quantile, probs = qprobs)
  
  out <- data.frame(dbh = dbh_pred,
                    t(dens_pred_quant))
  
  setNames(out, c('dbh', 'q025', 'q05', 'q25', 'q50', 'q75', 'q95', 'q975'))
  
}

gev_ci_df <- gev_ci(fit_gev, dbh_pred_bins, n1995)

###########################################

# generate data for plotting fit of gev

library(dplyr)

load('C:/Users/Q/google_drive/ForestLight/data/data_22jan2018/rawdataobj_22jan.r')
group_names <- c('all','all_classified','fg1','fg2','fg3','fg4','fg5','unclassified')
source('code/allfunctions27july.r')

fakebin_across_years <- function(dat_values, dat_classes, edges, mean_type = 'geometric', n_census = 5) {
  qprobs <- c(0.025, 0.25, 0.5, 0.75, 0.975)
  # add some padding just in case
  mins <- edges$bin_min
  mins[1] <- 0
  maxes <- edges$bin_max
  maxes[length(maxes)] <- Inf
  
  binstats <- t(sapply(1:length(mins), function(x) {
    indivs <- dat_values[dat_classes >= mins[x] & dat_classes < maxes[x]]
    if (mean_type == 'geometric') {
      mean_n <- exp(mean(log(indivs)))
      sd_n <- exp(sd(log(indivs)))
      ci_width <- qnorm(0.975) * sd(log(indivs)) / sqrt(length(indivs))
      ci_min <- exp(mean(log(indivs)) - ci_width)
      ci_max <- exp(mean(log(indivs)) + ci_width)
      
    } else {
      mean_n <- mean(indivs)
      sd_n <- sd(indivs)
      ci_width <- qnorm(0.975) * sd(indivs) / sqrt(length(indivs))
      ci_min <- mean_n - ci_width
      ci_max <- mean_n + ci_width
    }
    c(mean = mean_n, 
      sd = sd_n,
      quantile(indivs, probs = qprobs),
      ci_min = ci_min,
      ci_max = ci_max)
  }))
  dimnames(binstats)[[2]] <- c('mean', 'sd', 'q025', 'q25', 'median', 'q75', 'q975', 'ci_min', 'ci_max')
  data.frame(bin_midpoint = edges$bin_midpoint,
             bin_min = edges$bin_min,
             bin_max = edges$bin_max,
             mean_n_individuals = edges$bin_count / n_census,
             binstats)
}



# Set number of bins       
numbins <- 20

# Bin all trees including unclassified
allyeardbh <- unlist(lapply(alltreedat[2:6], '[', , 'dbh_corr'))
dbhbin_all <- logbin(x = allyeardbh, y = NULL, n = numbins)
dbhbin_all_byyear <- lapply(alltreedat[2:6], function(z) logbin_setedges(x = z$dbh_corr, y = NULL, edges = dbhbin_all))

all_dens_1995 <- dbhbin_all_byyear[[2]]


fpfig <- 'C:/Users/Q/google_drive/ForestLight/figs/credible_interval_plots'

area_core <- 42.84 # Area of the plot without the edge strip and without the "young" forest area.
global_diam_xlimits <- c(1, 316) # Maximum and minimum for x axis if it's diameter.

library(cowplot)


p_all_dens <- all_dens_1995 %>%
  filter(!is.na(bin_value), bin_value > 0) %>%
  mutate(bin_min = ifelse(bin_min == 0, bin_value, bin_min)) %>%
  mutate(bin_yvalue = bin_value/area_core, bin_ymin = bin_min/area_core, bin_ymax = bin_max/area_core) %>%
  ggplot() +
  geom_point(aes(x = bin_midpoint, y = bin_yvalue), color = 'black') +
  scale_x_log10(name = 'Diameter (cm)', limits = global_diam_xlimits) + 
  scale_y_log10(name = expression(paste('Density (individuals ha'^-1,')'))) +
  panel_border(colour = 'black')


p_all_dens +
  geom_ribbon(data = gev_ci_df, aes(ymin = q025/area_core, ymax = q975/area_core, x = dbh), fill = 'red', alpha = 0.5) +
  geom_line(data = gev_ci_df, aes(y = q50/area_core, x = dbh), color = 'red') +
  theme(legend.position = 'bottom', plot.title = element_text(hjust = 0)) + ggtitle('Density', 'GEV distribution fit shown')
ggsave(file.path(fpfig, 'alltree_1995_density_GEV.png'), height=5, width=5, dpi=400)

############################################
# Muller Landau's so-called Quasi Weibull

# Regular (not quasi) Weibull
mod_weib <- '
functions {
  // Wrote this by taking the log of right hand side of equation 4b in Muller-Landau 2006
  // Must include the normalizing constant (same as for other Weibull)
  real myweib_log (real y, real mu, real nu){
    return log(mu / nu) + (mu - 1) * log(y / nu) - (mu * y / nu);
  }

}

data {
  int<lower=0> N;
  real y[N];
}

parameters {
  real<lower=0> mu;
  real<lower=0> nu;
}

model {
  // priors on component parameters
  mu ~ gamma(.0001, .0001);
  nu ~ gamma(.0001, .0001);

  for(i in 1:N) {
    y[i] ~ myweib(mu, nu);  
  }
}
'

stanmod_weib <- stan_model(model_name = 'MLweibull', model_code = mod_weib)

NC <- 3
NI <- 3000
NW <- 2000

fit_weib <- sampling(stanmod_weib, data = dat_gev_alltree95, chains = NC, iter = NI, warmup = NW)

summ_weib <- summary(fit_weib)
summ_weib$summary
mcmc_trace(as.array(fit_weib))

myweib_ci <- function(fit, dbh_pred, n_indiv = 1) {
  pdf_myweib <- function(D, mu, nu) {
    # based on Equation 4b
    (mu/nu) * (D/nu)^(mu-1) * exp(-(D/nu)^mu)
  }
  
  pars <- do.call('cbind', extract(fit))
  dens_pred <- sapply(dbh_pred, pdf_myweib, mu = pars[,'mu'], nu = pars[,'nu'])
  dens_pred <- dens_pred * n_indiv
  
  # Get some quantiles from each of these.
  qprobs <- c(0.025, 0.05, 0.25, 0.5, 0.75, 0.95, 0.975)
  dens_pred_quant <- apply(dens_pred, 2, quantile, probs = qprobs)
  
  out <- data.frame(dbh = dbh_pred,
                    t(dens_pred_quant))
  
  setNames(out, c('dbh', 'q025', 'q05', 'q25', 'q50', 'q75', 'q95', 'q975'))
  
}

myweib_ci_df <- myweib_ci(fit_weib, dbh_pred_bins, n1995)

p_all_dens <- all_dens_1995 %>%
  filter(!is.na(bin_value), bin_value > 0) %>%
  mutate(bin_min = ifelse(bin_min == 0, bin_value, bin_min)) %>%
  mutate(bin_yvalue = bin_value/area_core, bin_ymin = bin_min/area_core, bin_ymax = bin_max/area_core) %>%
  ggplot() +
  geom_point(aes(x = bin_midpoint, y = bin_yvalue), color = 'black') +
  scale_x_log10(name = 'Diameter (cm)', limits = global_diam_xlimits) + 
  scale_y_log10(name = expression(paste('Density (individuals ha'^-1,')'))) +
  panel_border(colour = 'black')


p_all_dens +
  geom_ribbon(data = myweib_ci_df, aes(ymin = q025/area_core, ymax = q975/area_core, x = dbh), fill = 'red', alpha = 0.5) +
  geom_line(data = myweib_ci_df, aes(y = q50/area_core, x = dbh), color = 'red') +
  theme(legend.position = 'bottom', plot.title = element_text(hjust = 0)) + ggtitle('Density', 'Weibull distribution fit shown')

#############################
# Quasi Weibull
mod_qweib <- '
functions {
  // Wrote this by taking the log of right hand side of equation 7b in Muller-Landau 2006
  // Included shape over scale as the normalizing constant
  real qweib_log (real y, real alpha, real beta, real gamma) {
    return log(alpha / beta) + (alpha - 1) * log(y / beta) - (alpha + gamma) * (y / beta);
  }
}

data {
  int<lower=0> N;
  real y[N];
}

parameters {
  real<lower=0> alpha;
  real<lower=0> beta;
  real<lower=0, upper=2> gamma;
}

model {
  // priors on component parameters
  alpha ~ gamma(.0001, .0001);
  beta ~ gamma(.0001, .0001);
  gamma ~ uniform(0, 2);
  
  for(i in 1:N) {
    y[i] ~ qweib(alpha, beta, gamma);  
  }
}
'

stanmod_qweib <- stan_model(model_name = 'Quasi-Weibull', model_code = mod_qweib)

fit_qweib <- sampling(stanmod_qweib, data = dat_gev_alltree95, chains = NC, iter = NI, warmup = NW)

summ_qweib <- summary(fit_qweib)
summ_qweib$summary
mcmc_trace(as.array(fit_qweib))

# The quasi Weibull just becomes Weibull.