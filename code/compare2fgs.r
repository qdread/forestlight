# Comparison of Nadja's functional groups with the two groups clustered based on Wright's measurements

library(dplyr)

wright <- read.csv('C:/Users/Q/Dropbox/projects/forestlight/wright2010.csv', stringsAsFactors = FALSE)
wright[wright == -99] <- NA # Unknown values were given a value of -99
wright$SPECIES.[109:110] <- c('simplex_var1', 'simplex_var2') # Correct duplicate named subspecies.
wright$Taxon <- with(wright, paste(GENUS., SPECIES.))

wright_df <- with(wright, data.frame(Taxon, mrt = MRT25SAP/100, rgr = RGR95SAP, logit_mrt = qlogis(MRT25SAP/100), log_rgr = log10(RGR95SAP), stringsAsFactors = FALSE))
wright_df <- subset(wright_df, !is.na(mrt))

wright_pca <- with(wright_df, prcomp(data.frame(qlogis(mrt), log10(rgr)), scale=TRUE, center=TRUE)) # 90% of variation on the growth-mortality single axis. Nice.
pca_scores <- wright_pca$x[,1]

set.seed(4930116)
### Amendment 11 Oct.
# Change the clustering method to clara()

library(cluster)
clara_trans <- clara(x = wright_df[,c('logit_mrt','log_rgr')], k = 2)
#k_bothtrans <- kmeans(x = wright_df[,c('logit_mrt','log_rgr')], centers=2, nstart=25)
wright_df$cluster <- clara_trans$clustering
wright_df$pca <- pca_scores

wright_df$tol_wright <- c('G','S')[wright_df$cluster]

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
  select(sp, mrt, rgr, pca, tol_wright)

# Load Nadja's data
fg_both <- full_join(fgbci_less, bci_lookup_wright)

# Make a figure to show the FGs
ggplot(fg_both, aes(x = X1, y = X2, color = factor(fg))) +
  geom_point(data = subset(fg_both, tol_wright %in% 'G'), color = 'black', size = 3) +
  geom_point(data = subset(fg_both, tol_wright %in% 'S'), color = 'gray', size = 3) +
  geom_point() +
  labs(x = 'X1 slow to fast', y = 'X2 pioneers to breeders') +
  scale_color_manual(values = guild_colors)
ggsave('C:/Users/Q/google_drive/ForestLight/figs/figures_22jan2018/funcgroups_withGandS.pdf', height=6, width=6)
