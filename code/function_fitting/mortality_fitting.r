# fit logistic regression trendlines to mortality data in Stan
# Do this locally with RStan unless it is too slow, in which case run CmdStan on cluster
# QDR / Forestlight / 07 Oct 2019

# Load data ---------------------------------------------------------------

library(tidyverse)
library(rstan)

# Check out this slick trick so that we no longer need to comment out any paths - just run this and it sees if it's Quentin or not.
user <- Sys.info()['user']
gdrive_path <- ifelse(user == 'qread', '~/google_drive/ForestLight/', file.path('/Users',user,'Google Drive/ForestLight'))
github_path <- ifelse(user == 'qread', '~/Documents/GitHub/forestlight', file.path('/Users',user,'Documents/GitHub/forestlight'))

mort <- read_csv(file.path(gdrive_path, 'data/data_forplotting_aug2018/obs_mortalityindividuals.csv'))

# Compile stan model
logreg <- stan_model(file.path(github_path, 'stan/mortreg.stan'), model_name = 'logreg') # Model for each FG
logreg_mixed <- stan_model(file.path(github_path, 'stan/mortreg_fg_v3.stan'), model_name = 'logreg_fg') # Mixed model with random slopes and intercepts for FGs

# Set stan options
options(mc.cores = 2)
rstan_options(auto_write = TRUE)



# Test fit model on subset of data ----------------------------------------


set.seed(711)
mort_subset <- mort %>%
  filter(!fg %in% 'unclassified') %>%
  sample_n(10000) %>%
  mutate(died = alive == 0) %>%
  select(fg, died, light_received_byarea)

with(mort_subset, table(fg,died)) # Good numbers for all.

mort_subset_data <- with(mort_subset, list(N = nrow(mort_subset), M = 5, x = light_received_byarea, y = as.numeric(died), fg = as.numeric(factor(fg))))

mort_subset_mixed_fit <- sampling(logreg_mixed, data = mort_subset_data, chains = 2, iter = 1000, warmup = 500, thin = 1, seed = 333)
summary(mort_subset_mixed_fit)
# This is working but it will be too slow to run locally. Move to cluster to run.

# Generate a plot of the fitted lines.
# Just show median right now
sfit <- summary(mort_subset_mixed_fit)
coef_data <- data.frame(fg = paste0('fg', 1:5),
                        intercept = sfit$summary[grep('intercept', row.names(sfit$summary)), '50%'],
                        slope = sfit$summary[grep('slope', row.names(sfit$summary)), '50%'])

x_pred <- seq(log10(1), log10(412), length.out = 50)

fitted_values <- coef_data %>%
  group_by(fg) %>%
  group_modify(~ data.frame(x = x_pred, y = plogis(.x$intercept + .x$slope * x_pred)))

ggplot(fitted_values, aes(x = 10^x, y = y, group = fg, color = fg)) +
  geom_line() +
  scale_x_log10() +
  theme_bw()
# Looks good!


# Write rdump for all data ------------------------------------------------

mort_data <- mort %>%
  filter(!fg %in% 'unclassified') %>%
  mutate(died = alive == 0) %>%
  select(fg, died, light_received_byarea)

with(mort_data, table(fg,died))

mort_data_dump <- with(mort_data, list(N = nrow(mort_data), M = 5, x = light_received_byarea, y = as.numeric(died), fg = as.numeric(factor(fg))))

with(mort_data_dump, stan_rdump(list = names(mort_data_dump), file = '~/Dropbox/projects/forestlight/stanrdump_final/mortalitydump.r'))
