# New workflow to do all analysis on BCI 50-hectare plot, including light data
# QDR 14 July 2017

# 1. Load BCI data and do the quality controls suggested by Meakem et al (except we are using the pre-made bci.full dataset rather than coming up with our own stem tracking algorithm)
# Run in bciqc.r script
fp <- 'C:/Users/Q/Dropbox/projects/forestlight/'
load(file.path(fp, 'bcidata/bciqcrun.R')) # Removes strangler figs, applies taper correction on dbh, and flags the secondary forest area and the 20-m edge strip.

# Wright groups by 2 different methods
source('code/wright_groups.r')

# 2. Convert all census data frames to the correct units (dbh in cm, not mm; agb in kg, not Mg)
# Simultaneously, calculate biomass increments for each stem from the previous census to the current one. Combine everything into a single data frame if possible.

for (i in 2:7) {
  dat_i <- get(paste0('bci.full',i))
  dat_iminus1 <- get(paste0('bci.full',i-1))
  
  
}

bcicensusdat <- bci.full4 %>%
  filter(DFstatus == 'alive') %>%
  mutate(dbh = dbh/10,
         agb = agb * 1000

# 2. Calculate biomass increments for each stem from the previous census to the current one. Combine everything into a single data frame if possible.

# 3. Join with shade tolerance (Wright) groups.