# Modified analysis 02 Oct 2017.

# Do all the function fits using Caroline's MLE code.
# Include 2 clusters rather than 3.
# Don't worry about quantiles.

load('C:/Users/Q/Dropbox/projects/forestlight/bcidata/bciqcrun.R')

# Modification 25 August: get rid of young and edge trees for all datasets. This will reduce the number of hectares but will be most correct.

for (i in 1:7) {
  assign(paste0('bci.full', i), subset(get(paste0('bci.full', i)), !young & !edge))
  assign(paste0('bci.stem', i), subset(get(paste0('bci.stem', i)), !young & !edge))
}

growth8590 <- subset(growth8590, !young & !edge)
growth9095 <- subset(growth9095, !young & !edge)

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
k_bothtrans <- kmeans(x = wright_df[,c('logit_mrt','log_rgr')], centers=2, nstart=25)
wright_df$cluster <- k_bothtrans$cluster
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

# 2. Convert all census data frames to the correct units (dbh in cm, not mm; agb in kg, not Mg)
# Simultaneously, calculate biomass increments for each stem from the previous census to the current one. Combine everything into a single data frame if possible.
# Also get rid of the young (secondary) forest patches.

bci_production <- list()

for (i in 2:7) {
  dat_i <- get(paste0('bci.full',i))
  dat_iminus1 <- get(paste0('bci.full',i-1))
  bci_production[[i-1]] <- dat_i$agb_corr - dat_iminus1$agb_corr
  dat_i$production <- dat_i$agb_corr - dat_iminus1$agb_corr
  assign(paste0('bci.full', i), dat_i)
}

bci_production <- do.call(cbind, bci_production)
bci_production <- as.data.frame(bci_production)
names(bci_production) <- paste0('production', c('8085','8590','9095','9500','0005','0510'))

bcicensusdat <- list()

for (i in 2:7) {
  bcicensusdat[[i-1]] <- get(paste0('bci.full',i)) %>%
    filter(DFstatus == 'alive') %>%
    mutate(dbh_corr = dbh_corr/10,
           agb_corr = agb_corr * 1000,
           production = production * 1000) %>%
    filter(!young)
}


# 3. Join with shade tolerance (Wright) groups.

for (i in 1:6) {
  bcicensusdat[[i]] <- bcicensusdat[[i]] %>%
    mutate(Taxon = paste(Genus, Species)) %>%
    left_join(wright_df)
}

# 4. Join the appropriate years with Nadja's light data.
# bcicensusdat[[2]] is 1985-1990
# bcicensusdat[[3]] is 1990-1995

bcicensusdat[[2]] <- bcicensusdat[[2]] %>%
  mutate(tag = as.numeric(tag)) %>%
  left_join(growth8590 %>% select(tag, light, young, edge))

bcicensusdat[[3]] <- bcicensusdat[[3]] %>%
  mutate(tag = as.numeric(tag)) %>%
  left_join(growth9095 %>% select(tag, light, young, edge))

# 5. Run allometries to get total light received.

# Function to get a rough approximation of insolation by latitude.

insolation <- function(lat) {
  lat <- lat * pi/180 # to radians
  y <- sin(lat)
  0.25 * 1367 * (1 - 0.482 * (3*y^2 - 1)/2)
}

# Insolation at BCI, 9.2 degrees N
(insol_bci <- insolation(9.2))

# Function to get tree height and crown dimensions from dbh
# Use same parameters for all trees, taken from Bohlman and O'Brien

tp <- function(dbh) {
  h <- exp(.438 + .595 * log(dbh))    # Height
  cd <- exp(-.157 + .702 * log(dbh))  # Crown depth
  cr <- exp(-.438 + .658 * log(dbh))  # Crown radius
  cV <- exp(-.681 + 2.02 * log(dbh))  # Crown volume
  data.frame(h=h, cd=cd, cr=cr, cV=cV)
}

for (i in 2:3) {
  
  crowndim <- tp(bcicensusdat[[i]]$dbh_corr) 
  bcicensusdat[[i]]$crownarea <- pi * crowndim$cr^2
  bcicensusdat[[i]] <- transform(bcicensusdat[[i]], light_received = light * crownarea * insol_bci)
  
}

# Classification of light into 3 groups.
# This is percent of full irradiance, not the total light received.

for (i in 2:3) {
  
  light_groups <- cut(bcicensusdat[[i]]$light, breaks = 3)
  table(light_groups)
  light_groupcodes <- factor(light_groups, 
                             labels = c('Low light','Intermediate light','High light'))
  
  bcicensusdat[[i]]$light_group <- light_groupcodes
  
}

alltreedat <- list()
shadedat <- gapdat <- unclassifieddat <- list()

for (i in 1:6) {
  
  alltreedat[[i]] <- subset(bcicensusdat[[i]], 
                            !is.na(dbh_corr) & production > 0)
  shadedat[[i]] <- subset(alltreedat[[i]], tol_wright == 'S')
  gapdat[[i]] <- subset(alltreedat[[i]], tol_wright == 'G')
  unclassifieddat[[i]] <- subset(alltreedat[[i]], is.na(tol_wright))
}

# for 2 and 3, make additional data frames that have only trees with light measurements.
alltree_light_90 <- subset(alltreedat[[2]], !is.na(light))
alltree_light_95 <- subset(alltreedat[[3]], !is.na(light))
shade_light_90 <- subset(shadedat[[2]], !is.na(light))
shade_light_95 <- subset(shadedat[[3]], !is.na(light))
gap_light_90 <- subset(gapdat[[2]], !is.na(light))
gap_light_95 <- subset(gapdat[[3]], !is.na(light))
unclassified_light_90 <- subset(unclassifieddat[[2]], !is.na(light))
unclassified_light_95 <- subset(unclassifieddat[[3]], !is.na(light))




# Scalings ----------------------------------------------------------------

# Required scalings:
# 1. Fig 3a. Individual production scaling.
# 2. Fig 3b. Density scaling.
# 3. Fig 3c. Total production scaling. (binned)
# 4. Fig 4a. Total crown area scaling (binned)
# 5. Fig 4b. Total light received scaling (binned)
# 6. Fig 4c. Production by light received
# 7. Fig 5. Shade to gap ratios, Shade to gap production ratios, Shade tolerance score by light received regression


# Individual production scaling -------------------------------------------

# Straight line fit and curved line fit.

# Problem: Is there even a point to fitting this? It looks like a modeled thing.

# Plot them.
library(cowplot)
ggplot(shadedat[[6]], aes(x=dbh_corr, y=production)) + scale_x_log10() + scale_y_log10() + geom_point()
ggplot(shadedat[[6]], aes(x=dbh_corr, y=production)) + scale_x_log10() + scale_y_log10() + geom_hex()
ggplot(gapdat[[6]], aes(x=dbh_corr, y=production)) + scale_x_log10() + scale_y_log10() + geom_point() +
  stat_smooth(method='lm', color = 'blue') + stat_smooth(method='auto', color = 'red')

nls_fit_gap <- nls(logprod ~ C * (exp(k * logdbh)), 
    data = transform(gapdat[[6]], logprod=log10(production), logdbh=log10(dbh_corr)), 
    start = c(C=100, k=-2), 
    control = list(maxiter=500))

nls_fit_shade <- nls(logprod ~ C * (exp(k * logdbh)), 
    data = transform(shadedat[[6]], logprod=log10(production), logdbh=log10(dbh_corr)), 
    start = c(C=100, k=-2), 
    control = list(maxiter=500))

linear_fit_gap <- lm(logprod ~ logdbh, 
                      data = transform(gapdat[[6]], logprod=log10(production), logdbh=log10(dbh_corr)))

library(mgcv)
gam_fit_gap <- gam(logprod ~ logdbh, 
                   data = transform(gapdat[[6]], logprod=log10(production), logdbh=log10(dbh_corr)), k=3)

plot(gam_fit_gap, pages=1, residuals=T, all.terms=T)
summary(gam_fit_gap)

gam_fit_shade <- gam(logprod ~ s(logdbh, bs='cr', k = 3), 
                   data = transform(shadedat[[6]], logprod=log10(production), logdbh=log10(dbh_corr)))
linear_fit_shade <- lm(logprod ~ logdbh, 
                     data = transform(shadedat[[6]], logprod=log10(production), logdbh=log10(dbh_corr)))

predict(gam_fit_shade, newdata=)

ggplot(shadedat[[6]], aes(x=dbh_corr, y=production)) + 
  geom_point() +
  stat_function(geom='line', fun = function(x,beta) beta[1] + beta[2]*log10(x), args = list(beta=linear_fit_shade$coefficients), color = 'red', size=2) +
  geom_line(data=data.frame(cbind(shadedat[[6]], ypred=predict(gam_fit_shade))), aes(y=ypred), color='skyblue', size=2) +
  scale_x_log10() + scale_y_log10() 

ggplot(shadedat[[6]], aes(x=dbh_corr, y=production)) + 
  geom_point() +
  stat_smooth(method = 'lm', se = FALSE, color = 'red', size = 2) +
  stat_smooth(method = 'lm', formula = y ~ splines::bs(x,3), se = FALSE, color = 'goldenrod', size = 2) +
  stat_smooth(method = 'lm', formula = y ~ splines::bs(x,knots=c(1,2.7)), se = FALSE, color = 'green', size = 2) +
  stat_smooth(method = 'lm', formula = y ~ x + I(x^0.5), se = FALSE, color = 'blue', size = 2) +
  scale_x_log10() + scale_y_log10() 

spline_fit_shade <- lm(logprod ~ splines::bs(logdbh,3), 
                       data = transform(shadedat[[6]], logprod=log10(production), logdbh=log10(dbh_corr)))  
spline_fit2_shade <- lm(logprod ~ splines::bs(logdbh, knots = c(1, 2.7)), 
                       data = transform(shadedat[[6]], logprod=log10(production), logdbh=log10(dbh_corr)))  
sqrt_fit_shade <- lm(logprod ~ logdbh + I(logdbh^0.5),
                        data = transform(shadedat[[6]], logprod=log10(production), logdbh=log10(dbh_corr)))

AIC(sqrt_fit_shade)
AIC(spline_fit_shade)
AIC(spline_fit2_shade)
AIC(linear_fit_shade)

library(lme4)
mixed_sqrt_shade <- lmer(logprod ~ logdbh + I(logdbh^0.5) + (1|sp),
     data = transform(shadedat[[6]], logprod=log10(production), logdbh=log10(dbh_corr)))
mixed_lm_shade <- lmer(logprod ~ logdbh + (1|sp),
                         data = transform(shadedat[[6]], logprod=log10(production), logdbh=log10(dbh_corr)))
mixed_spline_shade <- lmer(logprod ~ splines::bs(logdbh, knots = c(1, 2.7)) + (1|sp),
                       data = transform(shadedat[[6]], logprod=log10(production), logdbh=log10(dbh_corr)))

AIC(mixed_spline_shade)
AIC(mixed_lm_shade)
AIC(mixed_sqrt_shade)
# Plot the fits hopefully with predictive intervals

nls_pars_gap <- nls_fit_gap$m$getPars()

ggplot(gapdat[[6]], aes(x=dbh_corr, y=production)) + 
  stat_function(geom='line', fun = function(x,beta) beta[1] + beta[2]*log10(x), args = list(beta=linear_fit_gap$coefficients), color = 'red', size=2) +
  stat_function(geom='line', fun = function(x,pars) (pars[1] * exp(pars[2] * log10(x))), args = list(pars=nls_pars_gap), color = 'skyblue', size=2) +
  scale_x_log10() + scale_y_log10() + geom_point()
  
  



# Test code below ---------------------------------------------------------



# Check whether there is some artifact in the tree production scalings
ggplot(shadedat[[6]] %>% filter(sp=='hybapr'), aes(x=dbh_corr, y=production)) + scale_x_log10() + scale_y_log10() + geom_point()
ggplot(shadedat[[6]] %>% filter(sp=='caseac'), aes(x=dbh_corr, y=production)) + scale_x_log10() + scale_y_log10() + geom_point()
ggplot(shadedat[[6]] %>% filter(sp=='tachve'), aes(x=dbh_corr, y=production)) + scale_x_log10() + scale_y_log10() + geom_point()
ggplot(shadedat[[6]] %>% filter(sp %in% unique(shadedat[[6]]$sp[1:50])), aes(x=dbh_corr, y=production)) + scale_x_log10() + scale_y_log10() + geom_point() + facet_wrap(~sp)

# These appear to be modeled. Fitting curves is maybe a bad idea?

# get basal area increments.
bci.full7$basalareaincrement <- pi*(bci.full7$dbh_corr/2)^2 - pi*(bci.full6$dbh_corr/2)^2
bci.full7$diameterincrement <- bci.full7$dbh_corr - bci.full6$dbh_corr

bci.full7 %>%
  filter(sp %in% unique(bci.full7$sp[1:100])) %>%
  ggplot(aes(x=dbh_corr, y=basalareaincrement)) +
  scale_x_log10() + scale_y_log10() +
  geom_point() +
  facet_wrap(~ sp)

bci.full7 %>%
  filter(sp %in% unique(bci.full7$sp[1:50])) %>%
  ggplot(aes(x=dbh_corr, y=diameterincrement)) +
  geom_point() +
  facet_wrap(~ sp)
