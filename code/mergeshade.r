# Combine the previously merged data frame with additional information on shade tolerance class.
# Make some figures as in John's schematics
# Modified 21 April: retain raw score from wright pca in the bci data.

fp <- 'C:/Users/Q/Google Drive/ForestLight/'
dat <- read.csv(file.path(fp, 'data/Data To Merge/GrowthLightMerged.csv'), stringsAsFactors = FALSE)

# Load shade tolerance.
library(XLConnect)
comita <- readWorksheetFromFile(file.path(fp, 'data/ComitaAppendix1.xlsx'), sheet=1)
laselva <- read.csv(file.path(fp, 'data/shade_tolerance_la_selva.csv'), stringsAsFactors = FALSE)
harv <- read.csv(file.path(fp, 'data/shade_tolerance_harvard.csv'), stringsAsFactors = FALSE)

comita$Shade.tolerance.guild[comita$Shade.tolerance.guild == '-'] <- NA

laselva$shadetol[laselva$shadetol == ''] <- NA
laselva$shadetol[laselva$shadetol == 'intolerant'] <- 'G'
laselva$shadetol[laselva$shadetol == 'intermediate'] <- 'I'
laselva$shadetol[laselva$shadetol == 'tolerant'] <- 'S'

harv$Tolerance <- substr(harv$Tolerance, 1, 1)

# Load Wright's growth data and back calculate the quantiles for shade tolerance
wright <- readWorksheetFromFile(file = file.path(fp, 'data/Shade Tolerance/Demographic/Wright et al 2010, growth mortality tradeoffs.xlsx'), sheet = 1, startRow = 26) # Get rid of the lines above header.
wright[wright == -99] <- NA # Unknown values were given a value of -99
wright$SPECIES.[109:110] <- c('simplex_var1', 'simplex_var2') # Correct duplicate named subspecies.
wright$Taxon <- with(wright, paste(GENUS., SPECIES.))

wright_df <- with(wright, data.frame(Taxon, mrt = MRT25SAP/100, rgr = RGR95SAP, stringsAsFactors = FALSE))
wright_df <- subset(wright_df, !is.na(mrt))

wright_pca <- with(wright_df, prcomp(data.frame(qlogis(mrt), log10(rgr)), scale=TRUE, center=TRUE)) # 90% of variation on the growth-mortality single axis. Nice.
pca_scores <- wright_pca$x[,1]
pca_groups <- cut(pca_scores, breaks = 3)
pca_groupcodes <- factor(pca_groups, labels = c('G','I','S'))

wright_df <- data.frame(wright_df, pca_scores, tol_wright = as.character(pca_groupcodes), stringsAsFactors = FALSE)

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






# Merge shade tolerance with data.

library(dplyr)
dat <- left_join(dat, comita %>% select(-Family) %>% rename(Taxon = Species, tol_comita = Shade.tolerance.guild))
dat <- left_join(dat, laselva %>% mutate(sppcode = gsub('_', ' ', sppcode)) %>% rename(Taxon = sppcode, tol_laselva = shadetol))
dat <- left_join(dat, harv %>% mutate(Species = gsub('_', ' ', Species)) %>% rename(Taxon = Species, tol_harvard = Tolerance) %>% select(-Source))
dat <- left_join(dat, wright_df %>% select(Taxon, pca_scores, tol_wright))

# Combine into a single shade tolerance score. Wright takes precedence, then Comita, then LaSelva.
tol_all <- pmin(dat$tol_harvard, dat$tol_wright, na.rm=T)
tol_all[!is.na(dat$tol_comita) & is.na(tol_all)] <- dat$tol_comita[!is.na(dat$tol_comita) & is.na(tol_all)]
tol_all[!is.na(dat$tol_laselva) & is.na(tol_all)] <- dat$tol_laselva[!is.na(dat$tol_laselva) & is.na(tol_all)]

# 5005 rows have a shade tolerance score and 2815 do not.
dat$tolerance <- tol_all

toltable <- dat %>% group_by(Site, Taxon) %>% summarize(t = any(!is.na(tolerance)), n = n()) %>% arrange(Site, -n)
#write.csv(toltable, file.path(fp, 'data/hastolerance.csv'), row.names=FALSE)
