library(tidyverse)
density <- read_csv('/Users/jgradym/Google Drive/ForestLight/data/data_piecewisefits/clean_parameters_density.csv')
growth <- read_csv('/Users/jgradym/Google Drive/ForestLight/data/data_piecewisefits/clean_parameters_production.csv')
r2_growth <- read_csv('/Users/jgradym/Google Drive/ForestLight/data/data_piecewisefits/clean_rsquared_production.csv')
N_fg <- read_csv('/Users/jgradym/Google Drive/ForestLight/data/data_piecewisefits/clean_N_by_functionalgroup.csv')

density_short <- density %>%
  select(year:mean, q025, q975) %>%
  filter(model == "three segment") %>%
  filter(parameter_description == "slope middle trees")
density_short$fg <- recode(density_short$fg, fg1 = "fast", fg2 = "large pioneer", fg3 = "slow", fg4 = "small breeder", fg5 = "medium")
density_short

growth_short <- growth %>%
  select(year:mean, q025, q975) %>%
  filter(model == "one segment", parameter_description == "slope") 
growth_short$fg <- recode(growth_short$fg, fg1 = "fast", fg2 = "large pioneer", fg3 = "slow", fg4 = "small breeder", fg5 = "medium")
growth_short 

r2_growth_short <- r2_growth %>%
  filter(model == "one segment") %>%
  select(fg, r2)
r2_growth_short$fg <- recode(r2_growth_short$fg, "all trees" = "alltree")

growth_short <- left_join(growth_short, r2_growth_short, by = "fg")
growth_short

N_fg_short <- N_fg %>%
  select(fg, n_indiv, n_spp)

growth_short <- left_join(growth_short, N_fg_short, by = "fg")
growth_short$n_indiv[growth_short$fg == "alltree"] <- sum(growth_short$n_indiv, na.rm = T)
growth_short$n_spp[growth_short$fg == "alltree"] <- sum(growth_short$n_spp, na.rm = T)
growth_short

growth_combo <- growth
