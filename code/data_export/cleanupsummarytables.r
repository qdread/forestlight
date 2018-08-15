# Clean up parameter table so that there are no duplicated rows
# Also clean up the information criteria table so that it's separated by variable.
# QDR / Forestlight / 15 Aug 2018

library(dplyr)

# Parameter table

params <- read.csv('~/google_drive/ForestLight/data/summarytables_june2018/paramci_by_fg.csv', stringsAsFactors = FALSE)

params <- params %>%
  filter((dens_model == 'pareto' & prod_model == 'powerlaw') | (dens_model == 'weibull' & prod_model == 'powerlawexp')) %>%
  mutate(response_variable = if_else(dens_model == 'pareto', 'density', 'growth'),
         model = if_else(response_variable == 'density', dens_model, prod_model)) %>%
  select(-dens_model, -prod_model) %>%
  select(year, response_variable, model, everything()) %>%
  arrange(year, response_variable, model, fg)

write.csv(params, '~/google_drive/ForestLight/data/summarytables_june2018/paramci_by_fg_noduplicates.csv', row.names = FALSE)

# Information criteria table

ics <- read.csv('~/google_drive/ForestLight/data/data_forplotting_aug2018/ics_by_fg.csv', stringsAsFactors = FALSE)

# Summarize the ics
ic_production <- ics %>%
  filter(criterion == 'LOOIC', variable == 'production') %>%
  group_by(fg, year, prod_model) %>%
  summarize(LOOIC = mean(ic))
ic_density <- ics %>%
  filter(criterion == 'LOOIC', variable == 'density') %>%
  group_by(fg, year, dens_model) %>%
  summarize(LOOIC = mean(ic))    

library(reshape2)
ic_production_cast <- dcast(ic_production, fg + year ~ prod_model) %>%
  mutate(deltaLOOIC = power - exp)
ic_density_cast <- dcast(ic_density, fg + year ~ dens_model) %>%
  mutate(deltaLOOIC = pareto - weibull)

write.csv(ic_production_cast, file = '~/google_drive/ForestLight/data/summarytables_june2018/LOOIC_production.csv', row.names = FALSE)
write.csv(ic_density_cast, file = '~/google_drive/ForestLight/data/summarytables_june2018/LOOIC_density.csv', row.names = FALSE)
