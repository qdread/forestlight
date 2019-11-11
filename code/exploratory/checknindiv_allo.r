# How many of the individuals in BCI are covered by Martinez Cano et al coefficients

load('~/google_drive/ForestLight/data/rawdataobj_alternativecluster.r')
martinezS2 <- read_csv('~/google_drive/ForestLight/data/BCI_raw/allometry_martinez-cano_table_s2.csv')

# by individual
table(alltreedat[[3]]$sp %in% martinezS2$species) # >90% of individuals

# by functional group
table(alltreedat[[3]]$fg, alltreedat[[3]]$sp %in% martinezS2$species, useNA = 'alw') # All are well represented except that fg 1 and fg 4 are relatively less so

# by number of species

uniquespp <- unique(alltreedat[[3]][,c('fg','sp')])


table(uniquespp$fg, uniquespp$sp %in% martinezS2$species, useNA = 'alw') # Again fg2 and fg3 have best representation.
