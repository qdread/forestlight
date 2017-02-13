# Attempt to run ordinal category regression on the BCI and Harvard tree data.

source('~/GitHub/forestlight/code/mergeshade.r')

pdat <- dat %>% 
  mutate(biomass3 = pmax(BiomassDGH_year3_kg, BiomassDBH_year3_kg, na.rm=TRUE),
         biomass2 = pmax(BiomassDGH_year2_kg, BiomassDBH_year2_kg, 0, na.rm=TRUE),
         biomass1 = pmax(BiomassDGH_year1_kg, BiomassDBH_year1_kg, 0, na.rm=TRUE),
         massprod23 = pmax(biomass3 - biomass2, 0, na.rm=TRUE),
         massprod12 = pmax(biomass2 - biomass1, 0, na.rm=TRUE),
         massprod13 = pmax((biomass3 - biomass1)/2, 0, na.rm=TRUE),
         diameter3 = pmax(Year3_DGH, Year3_DBH, na.rm=TRUE),
         diameter2 = pmax(Year2_DGH, Year2_DBH, 0, na.rm=TRUE),
         diameter1 = pmax(Year1_DGH_1, Year1_DBH, 0, na.rm=TRUE)) %>% 
  select(Site, biomass3, biomass2, diameter3, diameter2, massprod23, massprod12, massprod13, CII, tolerance)
bcidat <- pdat %>% filter(Site == 'Barro_Colorado_Island', !is.na(CII), !is.na(biomass3))
harvdat <- pdat %>% filter(Site == 'Harvard_Forest_LTER', !is.na(CII), !is.na(biomass3))

library(MASS)

bcidat <- mutate(bcidat, CIIfactor = factor(CII, ordered = TRUE), biomass3scale = scale(biomass3))

# Proportional odds logistic regression, fitting biomass to an ordered-factor response
options(contrasts = c("contr.treatment", "contr.poly"))
polrbci <- polr(CIIfactor ~ I(log10(biomass3)), data = bcidat, Hess = TRUE, method = 'logistic')
resid(polrbci)
studres(polrbci)

#### Note : must log-transform biomass for this to work.
# Works if only some of the rows are used!
polrsmall <- polr(CIIfactor ~ biomass3scale, data = bcidat[1:100,], Hess = TRUE, method = 'logistic')

# Use coefficients from smaller fit to do the fit to all data.
polrbci <- polr(CIIfactor ~ biomass3scale, data = bcidat, Hess = TRUE, method = 'logistic', start = c(polrsmall$coefficients, polrsmall$zeta))

table(predict(polrbci), bcidat$CIIfactor)
