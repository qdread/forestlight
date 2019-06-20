# Create piecewise data for plotting.

library(dplyr)

ci_df <- read.csv('~/google_drive/ForestLight/data/data_piecewisefits/newpiecewise_ci_by_fg.csv', stringsAsFactors = FALSE)
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



#################
# For individual + total light scaling and total volume scaling
# Added 20 June 2019

library(dplyr)

ci_df <- read.csv('~/google_drive/ForestLight/data/data_piecewisefits/totallightscaling/light_piecewise_ci_by_fg.csv', stringsAsFactors = FALSE)
area_core <- 42.84

ci_df$fg[ci_df$fg == 'alltree'] <- 'all'

fitted_indivlight <- ci_df %>%
  filter(variable == 'incoming_light_fitted') %>%
  select(-variable)

fitted_totallight <- ci_df %>%
  filter(variable == 'total_incoming_light_fitted') %>%
  select(-variable) %>%
  mutate_at(vars(starts_with('q')), funs(./area_core)) 

pred_indivlight <- ci_df %>%
  filter(variable == 'incoming_light') %>%
  select(-variable)

pred_totallight <- ci_df %>%
  filter(variable == 'total_incoming_light') %>%
  select(-variable) %>%
  mutate_at(vars(starts_with('q')), funs(./area_core)) 


fp <- '~/google_drive/ForestLight/data/data_piecewisefits/totallightscaling'

write.csv(pred_indivlight, file.path(fp, 'pred_indivlight.csv'), row.names = FALSE)
write.csv(pred_totallight, file.path(fp, 'pred_totallight.csv'), row.names = FALSE)
write.csv(fitted_indivlight, file.path(fp, 'fitted_indivlight.csv'), row.names = FALSE)
write.csv(fitted_totallight, file.path(fp, 'fitted_totallight.csv'), row.names = FALSE)

## Volume
ci_df <- read.csv('~/google_drive/ForestLight/data/data_piecewisefits/volume_piecewise_ci_by_fg.csv', stringsAsFactors = FALSE)
area_core <- 42.84

ci_df$fg[ci_df$fg == 'alltree'] <- 'all'

fitted_totalvol <- ci_df %>%
  select(-variable) %>%
  mutate_at(vars(starts_with('q')), funs(./area_core)) 

fp <- '~/google_drive/ForestLight/data/data_piecewisefits/totallightscaling'

write.csv(fitted_totalvol, file.path(fp, 'fitted_totalvol.csv'), row.names = FALSE)
