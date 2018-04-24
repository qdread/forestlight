# Production versus light: Fit all functional groups.
# Do with CMDstan 

# Edited 24 Apr 2018: Extract the "corrected" slope in log space also.

# Data dumps --------------------------------------------------------------

# Load data (changing file path if necessary)
fpdata <- '~/forestlight/stanrdump'
load(file.path(fpdata, 'rawdataobj_22jan.r'))

library(dplyr)
library(rstan)


dat90 <- alltree_light_90 %>%
  select(dbh_corr, production, light_received, crownarea, fg) %>%
  mutate(production_area = production/crownarea, light_area = light_received/crownarea) %>%
  mutate(fg = if_else(is.na(fg), 'unclassified', paste0('fg', fg)))

dat95 <- alltree_light_95 %>%
  select(dbh_corr, production, light_received, crownarea, fg) %>%
  mutate(production_area = production/crownarea, light_area = light_received/crownarea) %>%
  mutate(fg = if_else(is.na(fg), 'unclassified', paste0('fg', fg)))

make_logistic_data <- function(x, n_sub) {
  if (n_sub < nrow(x)) {
    sample_rows <- sample(nrow(x), n_sub)
  } else {
    sample_rows <- 1:nrow(x)
    n_sub <- nrow(x)
  }
  with(x[sample_rows, ],
       list(N = n_sub, x = light_received/crownarea, y = production/crownarea))
}

# Do with no subsampling.
dat90_fg <- dat90 %>%
  group_by(fg) %>%
  do(data = make_logistic_data(., n_sub = nrow(.)))
dat90_all <- make_logistic_data(dat90, n_sub = nrow(dat90))

dat95_fg <- dat95 %>%
  group_by(fg) %>%
  do(data = make_logistic_data(., n_sub = nrow(.)))
dat95_all <- make_logistic_data(dat90, n_sub = nrow(dat95))

# Create stan rdumps
fgnames <- c('fg1', 'fg2', 'fg3', 'fg4', 'fg5', 'unclassified')

library(purrr)
pmap(dat90_fg, function(fg, data) with(data, stan_rdump(names(data), file=file.path(fpdata, paste('dump_light', fg, '1990.r', sep = '_')))))
pmap(dat95_fg, function(fg, data) with(data, stan_rdump(names(data), file=file.path(fpdata, paste('dump_light', fg, '1995.r', sep = '_')))))

with(dat90_all, stan_rdump(names(dat90_all), file=file.path(fpdata, 'dump_light_alltree_1990.r')))
with(dat95_all, stan_rdump(names(dat90_all), file=file.path(fpdata, 'dump_light_alltree_1995.r')))


# Load fits and get CIs ---------------------------------------------------

fp <- '~/forestlight/stanoutput'

library(purrr)
library(tidyr)
library(dplyr)
library(rstan)

fgnames <- c('fg1', 'fg2', 'fg3', 'fg4', 'fg5', 'unclassified', 'alltree')

z <- expand.grid(fg=fgnames, year = c(1990, 1995), stringsAsFactors = FALSE)

# Define function to get credible intervals around the parameters and to get the interval around the fitted values

extract_light_ci <- function(fg, year) {
  getpar <- function(fit) as.data.frame(do.call(cbind, extract(fit)))
  fn_logistic3 <- function(x, G, b1, k, ...) G * (1 - b1 * exp(-k * x)) ^ 3
  light_pred <- exp(seq(log(1.1), log(412), length.out = 101))
  
  qprob <- c(0.025, 0.25, 0.5, 0.75, 0.975)
  qname <- c('q025', 'q25', 'q50', 'q75', 'q975')
  
  filenames <- paste0('fit_light_', fg, '_', year, '_', 1:3, '.csv')
  fit <- read_stan_csv(file.path(fp, filenames))
  pars <- getpar(fit)
  
  summ_fit <- summary(fit)
  
  get_pars <- c('G','b1','k','max_slope', 'x_max', 'y_max', 'log_slope')
  
  param_cis <- cbind(data.frame(fg = fg, year = year), summ_fit[[1]][get_pars, ])
  
  param_cis <- cbind(parameter = dimnames(param_cis)[[1]], param_cis)
  names(param_cis)[7:11] <- c('q025', 'q25', 'q50', 'q75', 'q975')
  
  pred <- do.call(rbind, pmap(pars, fn_logistic3, x = light_pred))
  pred_quant <- pred %>%
    as.data.frame %>%
    map(~ quantile(., probs = qprob)) %>%
    do.call(rbind, .) %>%
    as.data.frame %>%
    setNames(nm = qname)
  pred_quant <- data.frame(fg = fg, year = year, light_area = light_pred, pred_quant)
  
  list(pars = param_cis, preds = pred_quant)
}

all_output <- pmap(z, extract_light_ci)

# Combine the parameter confidence intervals and the predicted values into two single data frames.
all_pars <- do.call(rbind, map(all_output, 'pars'))
all_preds <- do.call(rbind, map(all_output, 'preds'))

write.csv(all_pars, file = '~/forestlight/lightbyarea_paramci_by_fg.csv', row.names = FALSE)
write.csv(all_preds, file = '~/forestlight/lightbyarea_predci_by_fg.csv', row.names = FALSE)