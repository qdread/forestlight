# Load canopy photo analysis output

fp <- 'C:/Users/Q/Google Drive/ForestLight/data/Canopy Pictures/Canopy_analysis'
bci_canopy <- read.csv(file.path(fp, 'BCI_no_mask_openess_par.csv'), stringsAsFactors = FALSE)
harv_canopy <- read.csv(file.path(fp, 'HFR_openess_par.csv'), stringsAsFactors = FALSE)

# Calculate medians by class
bci_byclass <- bci_canopy %>% group_by(Canopy_class) %>% 
  filter(Canopy_class != '4_3') %>%
  summarize(Openness = median(Openness),
            PPFD = median(PPFDTotalUnderPerDay.MJorMol.m2day.)) %>%
  rename(CII = Canopy_class) %>%
  mutate(Site = 'Barro_Colorado_Island', CII = as.numeric(CII))

harv_byclass <- harv_canopy %>% group_by(Canopy_class) %>% 
  summarize(Openness = median(Openness),
            PPFD = median(PPFDTotalUnderPerDay.MJorMol.m2day.)) %>%
  rename(CII = Canopy_class) %>%
  mutate(Site = 'Harvard_Forest_LTER')

canopyclass <- rbind(bci_byclass, harv_byclass)

# throw out the 

# Load all other data

source('~/GitHub/forestlight/code/mergeshade.r')

pdat <- dat %>% 
  mutate(diameter1 = pmax(Year1_DGH_1, Year1_DBH, 0, na.rm=TRUE),
         diameter2 = pmax(Year2_DGH, Year2_DBH, 0, na.rm=TRUE),
         diameter3 = pmax(Year3_DGH, Year3_DBH, na.rm=TRUE),
         biomass1 = pmax(BiomassDGH_year1_kg, BiomassDBH_year1_kg, 0, na.rm=TRUE),
         biomass2 = pmax(BiomassDGH_year2_kg, BiomassDBH_year2_kg, 0, na.rm=TRUE),
         biomass3 = pmax(BiomassDGH_year3_kg, BiomassDBH_year3_kg, na.rm=TRUE),
         massprod12 = pmax(biomass2 - biomass1, 0, na.rm=TRUE),
         massprod23 = pmax(biomass3 - biomass2, 0, na.rm=TRUE),
         massprod13 = pmax((biomass3 - biomass1)/2, 0, na.rm=TRUE)) %>% 
  select(Site, Genus, Species, diameter1, diameter2, diameter3, biomass1, biomass2, biomass3, massprod12, massprod23, massprod13, CII, tolerance, pca_scores)

pdat <- left_join(pdat, canopyclass)

# Create a new column with tolerance pooled to add intermediate and gap together
pdat$pooledtolerance <- pdat$tolerance
pdat$pooledtolerance[pdat$pooledtolerance %in% c('I','G')] <- 'G'

write.csv(pdat, file = 'C:/Users/Q/Google Drive/ForestLight/data/allforestdata19April.csv', row.names = FALSE)



bcidat <- pdat %>% filter(Site == 'Barro_Colorado_Island', !is.na(CII), !is.na(biomass3))
harvdat <- pdat %>% filter(Site == 'Harvard_Forest_LTER', !is.na(CII), !is.na(biomass3))

# Save the data as R object
save(bci_canopy, harv_canopy, bcidat, harvdat, file = 'C:/Users/Q/Google Drive/ForestLight/data/allforestdata19April.RData')
