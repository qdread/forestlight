# Additional work on BCI all census data.
# 21 April 2017

# To do list:
# get shade tolerance score or category for the trees in the 50 hectare plot
# if possible, convert competition index into a more direct measure of light availability (use method of Korzhukin, Groot, or Brunner)

# All trees pooled, compare slopes with and without the competition index.

lmsize <- lm(log10(production67) ~ log10(dbh), data=bcicensusdat, subset = !is.na(dbh) & production67 > 0 & comp_idx > 0)
lmsizecomp <- lm(log10(production67) ~ log10(dbh) + log10(comp_idx), data=bcicensusdat, subset = !is.na(dbh) & production67 > 0 & comp_idx > 0)

summary(lmsize)
summary(lmsizecomp)

with(bcicensusdat, plot(dbh, production67, log = 'xy'))
abline(lmsize, col = 'red')
abline(lmsizecomp, col = 'blue')


#################
# load taxon lookup table from 50-hectare plot and join with Wright data frame.

bci_lookup <- read.delim('C:/Users/Q/Dropbox/projects/forestlight/bcidata/ViewTax.txt', stringsAsFactors = FALSE)

bci_lookup <- bci_lookup %>%
  mutate(Taxon = paste(Genus, SpeciesName)) 

taxmatch <- bci_lookup$Taxon %in% wright_df$Taxon
bci_lookup$Taxon[taxmatch]

bcicensusdat <- left_join(bci_lookup, wright_df) %>%
  rename(sp = Mnemonic) %>%
  select(sp, mrt, rgr, pca_scores, tol_wright) %>%
  right_join(bcicensusdat)

table(!is.na(bcicensusdat$pca_scores))

# Gap trees
lmsizegap <- lm(log10(production67) ~ log10(dbh), data=bcicensusdat, subset = !is.na(dbh) & production67 > 0 & comp_idx > 0 & tol_wright == 'G')
lmsizecompgap <- lm(log10(production67) ~ log10(dbh) + log10(comp_idx), data=bcicensusdat, subset = !is.na(dbh) & production67 > 0 & comp_idx > 0 & tol_wright == 'G')
dbh_slope_gap <- lmsizecompgap$coefficients[2] # extract slope from full model

adjusted_lm_gap <- lm(log10(production67) ~ 1 + offset(dbh_slope * log10(dbh)), data=bcicensusdat, subset = !is.na(dbh) & production67 > 0 & comp_idx > 0 & tol_wright == 'G') # refit setting slope from full model and estimating intercept
intercept_gap <- adjusted_lm_gap$coefficients[1]

summary(lmsizegap)
summary(lmsizecompgap)

png('C:/Users/Q/Google Drive/ForestLight/figs/bci50hectare_sensitivity_gap.png', height=5, width=5, res=300, units='in')
with(subset(bcicensusdat, !is.na(dbh) & production67 > 0 & comp_idx > 0 & tol_wright == 'G'), plot(dbh, production67, log = 'xy', main='Gap species', las=1))
abline(lmsizegap, col = 'indianred', lwd=2)
abline(a = intercept_gap, b = dbh_slope_gap, col = 'steelblue1', lwd=2) # plot line with refit intercept and slope from full model
dev.off()

# Intermediate trees
lmsizeint <- lm(log10(production67) ~ log10(dbh), data=bcicensusdat, subset = !is.na(dbh) & production67 > 0 & comp_idx > 0 & tol_wright == 'I')
lmsizecompint <- lm(log10(production67) ~ log10(dbh) + log10(comp_idx), data=bcicensusdat, subset = !is.na(dbh) & production67 > 0 & comp_idx > 0 & tol_wright == 'I')
dbh_slope_int <- lmsizecompint$coefficients[2] # extract slope from full model

adjusted_lm_int <- lm(log10(production67) ~ 1 + offset(dbh_slope * log10(dbh)), data=bcicensusdat, subset = !is.na(dbh) & production67 > 0 & comp_idx > 0 & tol_wright == 'I') # refit setting slope from full model and estimating intercept
intercept_int <- adjusted_lm_int$coefficients[1]

summary(lmsizeint)
summary(lmsizecompint)

png('C:/Users/Q/Google Drive/ForestLight/figs/bci50hectare_sensitivity_int.png', height=5, width=5, res=300, units='in')
with(subset(bcicensusdat, !is.na(dbh) & production67 > 0 & comp_idx > 0 & tol_wright == 'G'), plot(dbh, production67, log = 'xy', main='Intermediate species', las=1))
abline(lmsizeint, col = 'indianred', lwd=2)
abline(a = intercept_int, b = dbh_slope_int, col = 'steelblue1', lwd=2)
dev.off()

# Shade trees
lmsizeshade <- lm(log10(production67) ~ log10(dbh), data=bcicensusdat, subset = !is.na(dbh) & production67 > 0 & comp_idx > 0 & tol_wright == 'S')
lmsizecompshade <- lm(log10(production67) ~ log10(dbh) + log10(comp_idx), data=bcicensusdat, subset = !is.na(dbh) & production67 > 0 & comp_idx > 0 & tol_wright == 'S')
dbh_slope_shade <- lmsizecompshade$coefficients[2] # extract slope from full model

adjusted_lm_shade <- lm(log10(production67) ~ 1 + offset(dbh_slope * log10(dbh)), data=bcicensusdat, subset = !is.na(dbh) & production67 > 0 & comp_idx > 0 & tol_wright == 'S') # refit setting slope from full model and estimating intercept
intercept_shade <- adjusted_lm_shade$coefficients[1]

summary(lmsizeshade)
summary(lmsizecompshade)

png('C:/Users/Q/Google Drive/ForestLight/figs/bci50hectare_sensitivity_shade.png', height=5, width=5, res=300, units='in')
with(subset(bcicensusdat, !is.na(dbh) & production67 > 0 & comp_idx > 0 & tol_wright == 'S'), plot(dbh, production67, log = 'xy', main='Shade-tolerant species', las=1))
abline(lmsizeshade, col = 'indianred', lwd=2)
abline(a = intercept_shade, b = dbh_slope_shade, col = 'steelblue1', lwd=2)
dev.off()

# Tolerance as continuous variable
lmsizecomp_continuous <- lm(log10(production67) ~ log10(dbh) + pca_scores:log10(comp_idx), data=bcicensusdat, subset = !is.na(dbh) & production67 > 0 & comp_idx > 0 & !is.na(pca_scores))

summary(lmsizecomp_continuous)

# Plot sensitivity to competition by PCA score of each species

# Get coefficient from linear model of competition index for each species
# Plot it versus pca score

getslope <- function(x) {
  xlm <- lm(log10(production67) ~ log10(dbh) + log10(comp_idx), data=x, subset = !is.na(dbh) & production67 > 0 & comp_idx > 0)
  data.frame(n_indiv = with(x, sum(!is.na(dbh) & production67 > 0 & comp_idx > 0)),
             comp_sensitivity = as.numeric(xlm$coefficients[3]))
}


compslopes <- bcicensusdat %>% 
  filter(!is.na(pca_scores)) %>%
  group_by(sp, pca_scores) %>%
  do(getslope(.))

sens_lm <- lm(comp_sensitivity ~ pca_scores, data=compslopes, weights = n_indiv)

summary(sens_lm)

png('C:/Users/Q/Google Drive/ForestLight/figs/bci50hectare_sensitivity_continuous.png', height=5, width=5, res=300, units='in')
with(compslopes, plot(pca_scores, comp_sensitivity, cex = log10(n_indiv)))  
abline(sens_lm, col = 'black', lwd=2)
abline(h=0, lty=3, col='blue')
dev.off()