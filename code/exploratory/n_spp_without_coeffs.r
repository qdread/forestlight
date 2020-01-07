# Tallies of species that don't have species-specific coefficients from Martinez Cano et al. and from Bohlman et al.
# Included in methods section of MS
# QDR / Forestlight / 07 Jan 2020

library(tidyverse)
library(forestscaling)

gdrive_path <- ifelse(Sys.info()['user'] == 'qread', '~/google_drive/ForestLight/', file.path('/Users',user,'Google Drive/ForestLight'))

bohlman_raw <- read_csv(file.path(gdrive_path, 'data/BCI_raw/allometry_bci_trees_bohlman.csv'))
martinezS2 <- read_csv('~/google_drive/ForestLight/data/BCI_raw/allometry_martinez-cano_table_s2.csv')

load(file.path(gdrive_path, 'data/rawdataobj_alternativecluster.R'))

table(bohlman_raw$species)
length(unique(bohlman_raw$species)) # 80 species

# Cross reference allometry with FG lookup.
fgbci <- read.table(file.path(gdrive_path, 'data/Ruger/fgroups_dynamics_new.txt'), stringsAsFactors = FALSE)

# Correct functional groups so that: 1 fast, 2 pioneer, 3 slow, 4 breeder, 5 intermediate
# Old 1,2,3,4,5 --> New 2,3,1,4,5
fgbci$fg5 <- match(fgbci$fg5, c(2,3,1,4,5))

fg_lookup <- data.frame(species = fgbci$splower, fg = fgbci$fg5)

bohlman_raw <- bohlman_raw %>%
  select(-X9) %>%
  left_join(fg_lookup)


# Tallies -----------------------------------------------------------------

# Species and individuals without coefficients for tree height (Martinez)

table(martinezS2$species %in% alltreedat[[3]]$sp)

length(unique(alltreedat[[3]]$sp[!alltreedat[[3]]$sp %in% martinezS2$species]))

spwithmeas <- sum(alltreedat[[3]]$sp %in% martinezS2$species)

nrow(alltreedat[[3]])-spwithmeas
1 - spwithmeas/nrow(alltreedat[[3]])

# Species and individuals without coefficients for crown depth (Bohlman raw data)

has_coef <- unique(bohlman_raw$species)
has_fg <- tolower(unique(fgbci$sp))
all_spp <- unique(alltreedat[[3]]$sp)

no_coef <- setdiff(all_spp, has_coef)
no_coef_yes_fg <- intersect(no_coef, has_fg)
no_coef_no_fg <- setdiff(no_coef, has_fg)

# Number of individuals and percentage
n_indiv_all <- nrow(alltreedat[[3]])
n_indiv_no_coef_yes_fg <- sum(alltreedat[[3]]$sp %in% no_coef_yes_fg)
n_indiv_no_coef_no_fg <- sum(alltreedat[[3]]$sp %in% no_coef_no_fg)

c(n_indiv_no_coef_yes_fg, n_indiv_no_coef_no_fg) / n_indiv_all
