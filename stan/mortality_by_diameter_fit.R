# New mortality fit and plot
# QDR 25 Jan 2022

mort <- read_csv('data/data_forplotting/obs_mortalityindividuals.csv')

library(cmdstanr)
library(tidyverse)

mort_mod <- cmdstan_model('~/GitHub/old_projects/forestscalingworkflow/model_scripts/mortreg_fg_v3.stan')

mort_data <- mort %>%
  filter(!fg %in% 'unclassified') %>%
  mutate(died = alive == 0) %>%
  select(fg, died, dbh)

mort_data_dump <- with(mort_data, list(N = nrow(mort_data), M = 5, x = dbh, y = as.numeric(died), fg = as.numeric(factor(fg))))

mort_fit <- mort_mod$sample(
  data = mort_data_dump,
  seed = 919,
  chains = 3,
  parallel_chains = 3,
  iter_warmup = 5000,
  iter_sampling = 1000,
  adapt_delta = 0.9,
  max_depth = 20
)
