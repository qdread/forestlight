# Comparison of Comita's shade tolerance ranking with those produced by different flavors of Wright's ranking.

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

wright_pca <- with(wright_df, prcomp(data.frame(qlogis(mrt), log10(rgr)), scale=TRUE, center=TRUE)) # 90% of variation on the growth-mortality single axis. Nice.
wright_pca_raw <- with(wright_df, prcomp(data.frame(mrt, rgr), scale=TRUE, center=TRUE))

pca_scores <- wright_pca$x[,1]
pca_groups <- cut(pca_scores, breaks = 3)
pca_groupcodes <- factor(pca_groups, labels = c('G','I','S'))

pca_rawscores <- wright_pca_raw$x[,1]
pca_rawgroups <- cut(pca_rawscores, breaks = 3)
pca_rawgroupcodes <- factor(pca_rawgroups, labels = c('G','I','S'))

pca_group_quantile <- Hmisc::cut2(pca_scores, g = 3)
pca_groupcodes_quantile <- factor(pca_group_quantile, labels = c('G','I','S'))

pca_rawgroup_quantile <- Hmisc::cut2(pca_rawscores, g = 3)
pca_rawgroupcodes_quantile <- factor(pca_rawgroup_quantile, labels = c('G','I','S'))

wright_df <- data.frame(wright_df, 
                        pca_scores, 
                        tol_wright_pcatransform = as.character(pca_groupcodes), 
                        tol_wright_pcaraw = as.character(pca_rawgroupcodes),
                        tol_wright_pcatransformquantile = as.character(pca_groupcodes_quantile),
                        tol_wright_pcarawquantile = as.character(pca_rawgroupcodes_quantile),
                        stringsAsFactors = FALSE)

comparison_df <- full_join(wright_df, comita %>% select(-Family) %>% rename(Taxon = Species, tol_comita = Shade.tolerance.guild))

# Agreement between Comita and 4 different ways of grouping Wright

# PCA with transformed variables and grouped by the range of the data
(c1 <- with(comparison_df, table(tol_comita, tol_wright_pcatransform)))
c(agree = sum(diag(c1)), disagree = sum(c1) - sum(diag(c1)))

# PCA with raw variables and grouped by the range of the data
(c2 <- with(comparison_df, table(tol_comita, tol_wright_pcaraw)))
c(agree = sum(diag(c2)), disagree = sum(c2) - sum(diag(c2)))

# PCA with transformed variables and split into three equal size groups
(c3 <- with(comparison_df, table(tol_comita, tol_wright_pcatransformquantile)))
c(agree = sum(diag(c3)), disagree = sum(c3) - sum(diag(c3)))

# PCA with raw variables and split into three equal size groups
(c4 <- with(comparison_df, table(tol_comita, tol_wright_pcarawquantile)))
c(agree = sum(diag(c4)), disagree = sum(c4) - sum(diag(c4)))




