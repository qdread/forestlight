# Tally species by functional group

gdrive_path <- ifelse(Sys.info()['user'] == 'qread', '~/google_drive/ForestLight/', file.path('/Users',Sys.info()['user'],'Google Drive/ForestLight'))

# Load data (changing file path if necessary)
load(file.path(gdrive_path, 'data/rawdataobj_alternativecluster.r'))

# Tally the functional groups' species and individuals.

library(dplyr)
library(purrr)

tallies <- map(alltreedat, function(x) {
  x %>%
    mutate(fg = if_else(is.na(fg), 'unclassified', paste0('fg',fg))) %>%
    group_by(fg) %>%
    summarize(n_indiv = n())
})

tallies <- cbind(year = rep(seq(1985,2010,5),each=6), do.call(rbind, tallies))

# all year sp
alldat <- data.frame(fg = unlist(map(alltreedat, 'fg')),
                     sp = unlist(map(alltreedat, 'sp')))

n_spp <- alldat %>%
  mutate(fg = if_else(is.na(fg), 'unclassified', paste0('fg',fg))) %>%
  group_by(fg) %>%
  summarize(n_spp = length(unique(sp))) %>%
  mutate(fg = fg_full_names[match(fg, fgs)])

n_ind <- tallies %>% 
  mutate(fg = fg_full_names[match(fg, fgs)]) 


write.csv(n_ind, file.path(gdrive_path, 'data/clean_summary_tables/clean_N_individuals_by_functionalgroup.csv'), row.names = FALSE)
write.csv(n_spp, file.path(gdrive_path, 'data/clean_summary_tables/clean_N_species_by_functionalgroup.csv'), row.names = FALSE)
