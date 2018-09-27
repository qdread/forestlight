# Create piecewise data for plotting.

library(dplyr)

ci_df <- read.csv('~/google_drive/ForestLight/data/data_piecewisefits/piecewise_ci_by_fg.csv', stringsAsFactors = FALSE)
area_core <- 42.84

ci_df$fg[ci_df$fg == 'alltree'] <- 'all'

pred_dens <- ci_df %>%
  filter(variable == 'density') %>%
  select(-variable) %>%
  mutate_at(vars(starts_with('q')), funs(./area_core)) 

fitted_indivprod <- ci_df %>%
  filter(variable == 'production_fitted') %>%
  select(-variable)

fitted_totalprod <- ci_df %>%
  filter(variable == 'total_production_fitted') %>%
  select(-variable) %>%
  mutate_at(vars(starts_with('q')), funs(./area_core)) 

pred_indivprod <- ci_df %>%
  filter(variable == 'production') %>%
  select(-variable)

pred_totalprod <- ci_df %>%
  filter(variable == 'total_production') %>%
  select(-variable) %>%
  mutate_at(vars(starts_with('q')), funs(./area_core)) 

fp <- '~/google_drive/ForestLight/data/data_piecewisefits'

write.csv(pred_dens, file.path(fp, 'pred_dens.csv'), row.names = FALSE)
write.csv(pred_indivprod, file.path(fp, 'pred_indivprod.csv'), row.names = FALSE)
write.csv(pred_totalprod, file.path(fp, 'pred_totalprod.csv'), row.names = FALSE)
write.csv(fitted_indivprod, file.path(fp, 'fitted_indivprod.csv'), row.names = FALSE)
write.csv(fitted_totalprod, file.path(fp, 'fitted_totalprod.csv'), row.names = FALSE)
