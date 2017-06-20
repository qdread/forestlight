# Load and explore Nadja Rueger's modeled relative irradiance values
# Merge with older BCI census data
# QDR 20 Jun 2017

fp <- 'C:/Users/Q/Dropbox/projects/forestlight/'

growth8590 <- read.delim(file.path(fp, 'BCI_light/growth_final8590.txt'), stringsAsFactors = FALSE)
growth9095 <- read.delim(file.path(fp, 'BCI_light/growth_final9095.txt'), stringsAsFactors = FALSE)

hist(growth8590$light)
hist(growth9095$light)
with(growth8590, plot(dbh, light))
with(growth8590, plot(log10(dbh), light))
with(growth9095, plot(log10(dbh), light))

# Must join these dataframes with the BCI census dataframes to get the biomass values.

# Load older census data.
#load(file.path(fp, 'bcidata/bci.full2.rdata')) # 1985
load(file.path(fp, 'bcidata/bci.full3.rdata')) # 1990
load(file.path(fp, 'bcidata/bci.full4.rdata')) # 1995

# Use 1990-1995 interval.
bci.full4$production34 <- pmax((bci.full4$agb - bci.full3$agb)/5, 0, na.rm = T)

library(dplyr)

bcicensusdat <- bci.full4 %>%
  filter(DFstatus == 'alive') %>%
  mutate(dbh = dbh/10,
         agb = agb * 1000,
         production34 = production34 * 1000)  # mm to cm and tonnes to kg

# Old code to join with the Wright tradeoff data.
library(XLConnect)
wright <- readWorksheetFromFile(file = 'C:/Users/Q/google_drive/ForestLight/data/Shade Tolerance/Demographic/Wright et al 2010, growth mortality tradeoffs.xlsx', sheet = 1, startRow = 26) # Get rid of the lines above header.
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

bci_lookup <- read.delim('C:/Users/Q/Dropbox/projects/forestlight/bcidata/ViewTax.txt', stringsAsFactors = FALSE)

bci_lookup <- bci_lookup %>%
  mutate(Taxon = paste(Genus, SpeciesName)) 

taxmatch <- bci_lookup$Taxon %in% wright_df$Taxon

bcicensusdat <- left_join(bci_lookup, wright_df) %>%
  rename(sp = Mnemonic) %>%
  select(sp, mrt, rgr, pca_scores, tol_wright) %>%
  right_join(bcicensusdat)

####################################################################

# Join bcicensusdat with Nadja's light values.

# Check how many of the tags in the growth9095 data frame are also in the bcicensusdat.
table(growth9095$tag %in% bcicensusdat$tag)
table(growth9095$tag %in% as.numeric(bcicensusdat$tag)) # All are in!

bcicensusdat <- bcicensusdat %>%
  mutate(tag = as.numeric(tag)) %>%
  left_join(growth9095 %>% dplyr::select(tag, light, dinc, interval))

# That was much easier than expected!

#####################################################################

# Correct the two slope plot function.

twoslopeplot <- function(dat, plottitle = 'plot title', xl = 'x label', yl = 'y label', binvar='diameter') {
  if (binvar == 'mass') dat$yv <- dat$agb else dat$yv <- dat$dbh
  lmsize <- lm(log10(production34) ~ log10(yv), data=dat)
  lmsizecomp <- lm(log10(production34) ~ log10(yv) + log10(light), data=dat)
  dbh_slope <- lmsizecomp$coefficients[2] # extract slope from full model
  # refit model, setting slope from full model and estimating intercept
  adjusted_lm <- lm(log10(production34) ~ 1 + offset(dbh_slope * log10(yv)), data=dat) 
  adjusted_intercept <- adjusted_lm$coefficients[1]
  
  ggplot(dat, aes(x = yv, y = production34)) +
    geom_point(alpha = 0.6) +
    geom_abline(slope = lmsize$coefficients[2], intercept = lmsize$coefficients[1], color = 'indianred', size = 2) +
    geom_abline(slope = dbh_slope, intercept = adjusted_intercept, color = 'steelblue1', size = 2) +
    scale_x_log10(name = xl,
                  breaks = scales::trans_breaks("log10", function(x) 10^x),
                  labels = scales::trans_format("log10", scales::math_format(10^.x))) + 
    scale_y_log10(name = yl,
                  breaks = scales::trans_breaks("log10", function(x) 10^x),
                  labels = scales::trans_format("log10", scales::math_format(10^.x))) +
    ggtitle(plottitle) +
    panel_border(colour = 'black')
}

######################################################################

# Plots using light instead of compidx. (log10 light)

xl1 <- 'Diameter (cm)'
yl1 <- expression(paste('Individual production (kg y'^-1,')', sep=''))

library(cowplot)

alltreedat <- subset(bcicensusdat, !is.na(dbh) & production34 > 0 & !is.na(light))
twoslopeplot(dat = alltreedat, 
             plottitle = 'All species', 
             xl = xl1, 
             yl =  yl1)

shadedat <- subset(bcicensusdat, !is.na(dbh) & production34 > 0 & !is.na(light) & tol_wright == 'S')
twoslopeplot(dat = shadedat, 
             plottitle = 'Shade-tolerant species', 
             xl = xl1, 
             yl =  yl1)

intdat <- subset(bcicensusdat, !is.na(dbh) & production34 > 0 & !is.na(light) & tol_wright == 'I')
twoslopeplot(dat = intdat, 
             plottitle = 'Intermediate species', 
             xl = xl1, 
             yl =  yl1)

gapdat <- subset(bcicensusdat, !is.na(dbh) & production34 > 0 & !is.na(light) & tol_wright == 'G')
twoslopeplot(dat = gapdat, 
             plottitle = 'Gap species', 
             xl = xl1, 
             yl =  yl1)

######################################################################

# Hypothesis testing with the new light values.

lmsize_all <- lm(log10(production34) ~ log10(dbh), data=alltreedat)
lmsizecomp_all <- lm(log10(production34) ~ log10(dbh) + log10(light), data=alltreedat)
lmsize_shade <- lm(log10(production34) ~ log10(dbh), data=shadedat)
lmsizecomp_shade <- lm(log10(production34) ~ log10(dbh) + log10(light), data=shadedat)
lmsize_int <- lm(log10(production34) ~ log10(dbh), data=intdat)
lmsizecomp_int <- lm(log10(production34) ~ log10(dbh) + log10(light), data=intdat)
lmsize_gap <- lm(log10(production34) ~ log10(dbh), data=gapdat)
lmsizecomp_gap <- lm(log10(production34) ~ log10(dbh) + log10(light), data=gapdat)

slopedat <- data.frame(withcompidx = rep(c('with light','without light'), each=4),
                       guild = c('all','shade','intermediate','gap'),
                       slope = c(lmsizecomp_all$coef[2], lmsizecomp_shade$coef[2], lmsizecomp_int$coef[2], lmsizecomp_gap$coef[2],
                                 lmsize_all$coef[2], lmsize_shade$coef[2], lmsize_int$coef[2], lmsize_gap$coef[2]),
                       cimin = c(confint(lmsizecomp_all)[2,1], confint(lmsizecomp_shade)[2,1], confint(lmsizecomp_int)[2,1], confint(lmsizecomp_gap)[2,1],
                                 confint(lmsize_all)[2,1], confint(lmsize_shade)[2,1], confint(lmsize_int)[2,1], confint(lmsize_gap)[2,1]),
                       cimax = c(confint(lmsizecomp_all)[2,2], confint(lmsizecomp_shade)[2,2], confint(lmsizecomp_int)[2,2], confint(lmsizecomp_gap)[2,2],
                                 confint(lmsize_all)[2,2], confint(lmsize_shade)[2,2], confint(lmsize_int)[2,2], confint(lmsize_gap)[2,2]))

ggplot(slopedat, aes(x = guild, group = interaction(guild, withcompidx), color = withcompidx, y = slope, ymin = cimin, ymax = cimax)) +
  geom_pointrange() +
  panel_border(colour = 'black') +
  labs(color = 'Regression type') +
  theme(legend.position = 'bottom')