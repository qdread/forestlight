library(tidyverse)
gdrive_path <- file.path('/Users/jgradym/Google Drive/ForestLight/')

pred_density <- read_csv(file.path(gdrive_path,'/data/data_piecewisefits/pred_dens.csv'))
pred_density_all <- pred_density %>%
  filter(dens_model == "3", fg == "all")
  
pred_indiv_light <- read_csv(file.path(gdrive_path,'data/data_piecewisefits/pred_indivlight.csv'))
pred_indiv_light_all <- pred_indiv_light %>%
  filter(prod_model == "2", fg == "all")

pred_tot_light <- read_csv(file.path(gdrive_path, 'data/data_piecewisefits/pred_totallight.csv'))
pred_tot_light_all <- pred_tot_light %>%
  filter(dens_model == "3", fg == "all", prod_model == "2")

calc_tot_light_all <- pred_indiv_light_all$q50 * pred_density_all$q50

#compare calculated with predicted totals
calc_tot_light_all/pred_tot_light_all$q50

# calculated is 1.57 time times bigger than predicted
