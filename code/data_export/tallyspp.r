# Load data (changing file path if necessary)
fpdata <- 'C:/Users/Q/google_drive/ForestLight/data/data_22jan2018'
load(file.path(fpdata, 'rawdataobj_22jan.r'))

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
  summarize(n_spp = length(unique(sp)))

write.csv(tallies, 'C:/Users/Q/google_drive/ForestLight/data/summarytables_12apr2018/tally_indiv_by_fg.csv', row.names = FALSE)
write.csv(n_spp, 'C:/Users/Q/google_drive/ForestLight/data/summarytables_12apr2018/tally_spp_by_fg.csv', row.names = FALSE)
