# Load Wright's shade tolerance and mortality data, and split into 3 groups using different methods
# Join with the BCI taxon lookup table

library(dplyr)

wright <- read.csv('C:/Users/Q/Dropbox/projects/forestlight/wright2010.csv', stringsAsFactors = FALSE)
wright[wright == -99] <- NA # Unknown values were given a value of -99
wright$SPECIES.[109:110] <- c('simplex_var1', 'simplex_var2') # Correct duplicate named subspecies.
wright$Taxon <- with(wright, paste(GENUS., SPECIES.))

wright_df <- with(wright, data.frame(Taxon, mrt = MRT25SAP/100, rgr = RGR95SAP, stringsAsFactors = FALSE))
wright_df <- subset(wright_df, !is.na(mrt))

wright_pca <- with(wright_df, prcomp(data.frame(qlogis(mrt), log10(rgr)), scale=TRUE, center=TRUE)) # 90% of variation on the growth-mortality single axis. Nice.
pca_scores <- wright_pca$x[,1]

# Grouping by even spacing of the pca axis (best)
pca_groups <- cut(pca_scores, breaks = 3)
pca_groupcodes <- factor(pca_groups, labels = c('G','I','S'))

# Grouping by splitting into even numbers of species in each group
pca_groups_quantile <- Hmisc::cut2(pca_scores, g = 3)
pca_groupcodes_quantile <- factor(pca_groups_quantile, labels = c('G','I','S'))

wright_df <- data.frame(wright_df, pca_scores, tol_wright = as.character(pca_groupcodes), tol_wright_quantile = as.character(pca_groupcodes_quantile), stringsAsFactors = FALSE)

# Correct wright_df entries that are not correct.
wright_df$Taxon[grep('Beilsc',wright_df$Taxon)] <- 'Beilschmiedia pendula'
wright_df$Taxon[grep('Cestrum',wright_df$Taxon)] <- 'Cestrum megalophyllum'
wright_df$Taxon[grep('phyllu arg',wright_df$Taxon)] <- 'Chrysophyllum argenteum'
wright_df$Taxon[grep('Coccol',wright_df$Taxon)] <- 'Coccoloba manzinellensis'
wright_df$Taxon[grep('Tabern',wright_df$Taxon)] <- 'Tabernaemontana arborea'
wright_df$Taxon[grep('var1',wright_df$Taxon)] <- 'Swartzia simplex_var.grandiflora'
wright_df$Taxon[grep('var2',wright_df$Taxon)] <- 'Swartzia simplex_var.ochnacea'
wright_df$Taxon[grep('colorado',wright_df$Taxon)] <- 'Eugenia coloradoensis'
wright_df$Taxon[grep('Croton',wright_df$Taxon)] <- 'Croton billbergianus'

bci_lookup <- read.delim('C:/Users/Q/Dropbox/projects/forestlight/bcidata/ViewTax.txt', stringsAsFactors = FALSE)

bci_lookup <- bci_lookup %>%
  mutate(Taxon = paste(Genus, SpeciesName)) 

taxmatch <- bci_lookup$Taxon %in% wright_df$Taxon

bci_lookup_wright <- left_join(bci_lookup, wright_df) %>%
  rename(sp = Mnemonic) %>%
  select(sp, mrt, rgr, pca_scores, tol_wright, tol_wright_quantile)
