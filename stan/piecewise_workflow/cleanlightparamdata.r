# Clean data frame of parameter CIs to remove duplicates
# QDR Forestlight 18 Feb 2019

params <- read.csv('~/google_drive/ForestLight/data/data_piecewisefits/lightpiecewise/lightpiecewise_paramci_by_fg.csv', stringsAsFactors = FALSE)
density_params <- c('alpha','alpha_high','alpha_low','alpha_mid','m','n','mu_logn','sigma_logn','tau','tau_high','tau_low')
production_params <- c('beta0','beta1','beta1_high','beta1_low','x0','delta')

library(dplyr)

params <- params %>%
  filter((prod_model == 1 & parameter %in% density_params) | (dens_model == '1' & parameter %in% production_params)) %>%
  mutate(variable = if_else(parameter %in% density_params, 'density', 'production'),
         model = if_else(parameter %in% density_params, as.character(dens_model), as.character(prod_model))) %>%
  select(year, variable, model, everything()) %>% 
  select(-dens_model, -prod_model) %>%
  arrange(year, fg, variable, model)

write.csv(params, '~/google_drive/ForestLight/data/data_piecewisefits/lightpiecewise/lightpiecewise_paramci_by_fg_noduplicates.csv', row.names = FALSE)
