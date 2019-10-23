# Visualize the information criteria
# Edit 15 Aug: use new data and include new file paths

ics <- read.csv('~/google_drive/ForestLight/data/data_forplotting_aug2018/ics_by_fg.csv', stringsAsFactors = FALSE)

# Midsize
#ics <- read.csv('C:/Users/Q/Dropbox/projects/forestlight/ics_by_fg_midsizetrees.csv', stringsAsFactors = FALSE)

library(dplyr)
library(ggplot2)


ics %>%
  filter(criterion == 'LOOIC', year == 1995, variable == 'production') %>%
  ggplot(aes(x = interaction(dens_model, prod_model), y = ic)) +
    geom_point() +
    facet_wrap(~ fg, scales='free')

ics %>%
  filter(criterion == 'LOOIC', year == 1995, variable == 'density') %>%
  ggplot(aes(x = interaction(dens_model, prod_model), y = ic)) +
  geom_point() +
  facet_wrap(~ fg, scales='free')

ics %>%
  filter(criterion == 'WAIC', year == 1995, variable == 'production') %>%
  ggplot(aes(x = interaction(dens_model, prod_model), y = ic)) +
  geom_point() +
  facet_wrap(~ fg, scales='free')

ics %>%
  filter(criterion == 'WAIC', year == 1995, variable == 'density') %>%
  ggplot(aes(x = interaction(dens_model, prod_model), y = ic)) +
  geom_point() +
  facet_wrap(~ fg, scales='free')


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

write.csv(ic_production_cast, file = 'C:/Users/Q/google_drive/ForestLight/data/summarytables_12apr2018/LOOIC_production.csv', row.names = FALSE)
write.csv(ic_density_cast, file = 'C:/Users/Q/google_drive/ForestLight/data/summarytables_12apr2018/LOOIC_density.csv', row.names = FALSE)

# Midsize
write.csv(ic_production_cast, file = 'C:/Users/Q/google_drive/ForestLight/data/summarytables_12apr2018/LOOIC_production_midsize.csv', row.names = FALSE)
write.csv(ic_density_cast, file = 'C:/Users/Q/google_drive/ForestLight/data/summarytables_12apr2018/LOOIC_density_midsize.csv', row.names = FALSE)
