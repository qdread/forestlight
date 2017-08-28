# New workflow to do all analysis on BCI 50-hectare plot, including light data
# QDR 14 July 2017


# Load and massage data ---------------------------------------------------



# 1. Load BCI data and do the quality controls suggested by Meakem et al (except we are using the pre-made bci.full dataset rather than coming up with our own stem tracking algorithm)
# Run in bciqc.r script
setwd('~/forestlight')
load('bciqcrun.R') # Removes strangler figs, applies taper correction on dbh, and flags the secondary forest area and the 20-m edge strip.

# Modification 25 August: get rid of young and edge trees for all datasets. This will reduce the number of hectares but will be most correct.

for (i in 1:7) {
  assign(paste0('bci.full', i), subset(get(paste0('bci.full', i)), !young & !edge))
  assign(paste0('bci.stem', i), subset(get(paste0('bci.stem', i)), !young & !edge))
}

growth8590 <- subset(growth8590, !young & !edge)
growth9095 <- subset(growth9095, !young & !edge)

# Wright groups by 2 different methods
source('wright_groups.r')

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

# 6. Split each year into shade, intermediate, and gap groups.
# Use both the even axis split and the quantile split.

alltreedat <- list()
shadedat <- intdat <- gapdat <- unclassifieddat <- list()
shadedat_quant <- intdat_quant <- gapdat_quant <- list()

for (i in 1:6) {

alltreedat[[i]] <- subset(bcicensusdat[[i]], 
                     !is.na(dbh_corr) & production > 0)
shadedat[[i]] <- subset(alltreedat[[i]], tol_wright == 'S')
intdat[[i]] <- subset(alltreedat[[i]], tol_wright == 'I')
gapdat[[i]] <- subset(alltreedat[[i]], tol_wright == 'G')
unclassifieddat[[i]] <- subset(alltreedat[[i]], is.na(tol_wright))
shadedat_quant[[i]] <- subset(alltreedat[[i]], tol_wright_quantile == 'S')
intdat_quant[[i]] <- subset(alltreedat[[i]], tol_wright_quantile == 'I')
gapdat_quant[[i]] <- subset(alltreedat[[i]], tol_wright_quantile == 'G')

}

# for 2 and 3, make additional data frames that have only trees with light measurements.
alltree_light_90 <- subset(alltreedat[[2]], !is.na(light))
alltree_light_95 <- subset(alltreedat[[3]], !is.na(light))
shade_light_90 <- subset(shadedat[[2]], !is.na(light))
shade_light_95 <- subset(shadedat[[3]], !is.na(light))
int_light_90 <- subset(intdat[[2]], !is.na(light))
int_light_95 <- subset(intdat[[3]], !is.na(light))
gap_light_90 <- subset(gapdat[[2]], !is.na(light))
gap_light_95 <- subset(gapdat[[3]], !is.na(light))
unclassified_light_90 <- subset(unclassifieddat[[2]], !is.na(light))
unclassified_light_95 <- subset(unclassifieddat[[3]], !is.na(light))
shadequant_light_90 <- subset(shadedat_quant[[2]], !is.na(light))
shadequant_light_95 <- subset(shadedat_quant[[3]], !is.na(light))
intquant_light_90 <- subset(intdat_quant[[2]], !is.na(light))
intquant_light_95 <- subset(intdat_quant[[3]], !is.na(light))
gapquant_light_90 <- subset(gapdat_quant[[2]], !is.na(light))
gapquant_light_95 <- subset(gapdat_quant[[3]], !is.na(light))


# Do scalings -------------------------------------------------------------

# Source functions
source('allfunctions27july.r')

# LIST OF SCALINGS
# ----------------

# 1. Density scaling (dbh as scaling variable)
#   - Do this for all trees, for unclassified trees, for the 3 shade groups made by splitting axis evenly, 
#   and for the 3 shade groups made by splitting species evenly
#   - Do this with cutoff and without cutoff
# 2. Density scaling (light received as scaling variable)
#   - Same as above
# 3. Individual production scaling
# 4. Total production scaling
# 5. Individual light-received scaling
# 6. Total light-received scaling
# 7. 3x3 shade group by light group density scaling
#   - Do this for each type of shade grouping

# Get weighted average of shade, intermediate, and gap proportions at each light level
# Individual light vs individual production plots
# Slope of production~diameter regression, with and without light included as a covariate.

######
# PRODUCTION BY DIAMETER SLOPES
# xl1 <- 'Diameter (cm)'
# yl1 <- expression(paste('Individual production (kg y'^-1,')', sep=''))
# 
# library(cowplot)

# 1990
# Only trees with light measurements. Use light_received!
# twoslopeplot(dat = alltree_light_90, plottitle = 'All species 1990', xl = xl1, yl =  yl1)
# twoslopeplot(dat = unclassified_light_90, plottitle = 'Unclassified species 1990', xl = xl1, yl =  yl1)
# twoslopeplot(dat = shade_light_90, plottitle = 'Shade-tolerant species 1990', xl = xl1, yl =  yl1)
# twoslopeplot(dat = int_light_90, plottitle = 'Intermediate species 1990', xl = xl1, yl =  yl1)
# twoslopeplot(dat = gap_light_90, plottitle = 'Gap species 1990', xl = xl1, yl =  yl1)
# twoslopeplot(dat = shadequant_light_90, plottitle = 'Shade-tolerant species 1990 (by quantile)', xl = xl1, yl =  yl1)
# twoslopeplot(dat = intquant_light_90, plottitle = 'Intermediate species 1990 (by quantile)', xl = xl1, yl =  yl1)
# twoslopeplot(dat = gapquant_light_90, plottitle = 'Gap species 1990 (by quantile)', xl = xl1, yl =  yl1)
# 
# # 1995
# twoslopeplot(dat = alltree_light_95, plottitle = 'All species 1995', xl = xl1, yl =  yl1)
# twoslopeplot(dat = unclassified_light_95, plottitle = 'Unclassified species 1995', xl = xl1, yl =  yl1)
# twoslopeplot(dat = shade_light_95, plottitle = 'Shade-tolerant species 1995', xl = xl1, yl =  yl1)
# twoslopeplot(dat = int_light_95, plottitle = 'Intermediate species 1995', xl = xl1, yl =  yl1)
# twoslopeplot(dat = gap_light_95, plottitle = 'Gap species 1995', xl = xl1, yl =  yl1)
# twoslopeplot(dat = shadequant_light_95, plottitle = 'Shade-tolerant species 1995 (by quantile)', xl = xl1, yl =  yl1)
# twoslopeplot(dat = intquant_light_95, plottitle = 'Intermediate species 1995 (by quantile)', xl = xl1, yl =  yl1)
# twoslopeplot(dat = gapquant_light_95, plottitle = 'Gap species 1995 (by quantile)', xl = xl1, yl =  yl1)

# Get slopes for hypothesis testing and output them to a table.
# Confidence intervals as well.

get_two_slopes <- function(dat, xv = 'dbh_corr', yv = 'production', modv = 'light_received') {
  form1 <- paste0('log10(', yv, ') ~ log10(', xv, ')')
  form2 <- paste0(form1, ' + log10(', modv, ')')
  lm1 <- lm(formula = form1, data = dat)
  lm2 <- lm(formula = form2, data = dat)
  withoutmod <- c(slope = as.numeric(lm1$coef[2]), cimin = confint(lm1)[2,1], cimax = confint(lm1)[2,2])
  withmod <- c(slope = as.numeric(lm2$coef[2]), cimin = confint(lm2)[2,1], cimax = confint(lm2)[2,2])
  return(data.frame(withindex = c('without','with'), rbind(withoutmod, withmod)))
}

allyears_names <- c('alltreedat', 'shadedat', 'intdat', 'gapdat', 'unclassifieddat', 'shadedat_quant', 'intdat_quant', 'gapdat_quant')
names1990 <- c('alltree_light_90', 'shade_light_90', 'int_light_90', 'gap_light_90', 'unclassified_light_90', 'shadequant_light_90', 'intquant_light_90', 'gapquant_light_90')
names1995 <- c('alltree_light_95', 'shade_light_95', 'int_light_95', 'gap_light_95', 'unclassified_light_95', 'shadequant_light_95', 'intquant_light_95', 'gapquant_light_95')


for (i in c(names1990, names1995)) {
  dat <- get(i)
  assign(paste0(i, '_dens_slopes'), get_two_slopes(dat))
}


# Do shade tolerance as continuous variable and plot.
getslope <- function(dat) {
  xlm <- lm(log10(production) ~ log10(dbh_corr) + log10(light_received), data=dat)
  data.frame(n_indiv = nrow(dat),
             comp_sensitivity = as.numeric(xlm$coefficients[3]))
}


compslopes1990 <- alltree_light_90 %>% 
  filter(!is.na(pca_scores)) %>%
  group_by(sp, pca_scores) %>%
  do(getslope(.))

compslopes1995 <- alltree_light_95 %>% 
  filter(!is.na(pca_scores)) %>%
  group_by(sp, pca_scores) %>%
  do(getslope(.))

# sens_lm1990 <- lm(comp_sensitivity ~ pca_scores, data=compslopes1990, weights = n_indiv)
# sens_lm1995 <- lm(comp_sensitivity ~ pca_scores, data=compslopes1995, weights = n_indiv)
# 
# summary(sens_lm1990)
# summary(sens_lm1995)
# 
# ggplot(compslopes1990, aes(x = pca_scores, y = comp_sensitivity)) +
#   geom_hline(yintercept = 0, lty = 3, color = 'gray50', size = 1.5) +
#   geom_point(aes(size = log10(n_indiv)), pch = 1) +
#   geom_abline(slope=sens_lm1990$coef[2], intercept = sens_lm1990$coef[1], color = 'blue', size = 2) +
#   scale_size_continuous(name = 'N', breaks=c(2,3,4), labels=10^c(2,3,4)) +
#   labs(x = 'Shade tolerance', y = 'Production sensitivity') + ggtitle('Competition sensitivity, 1990')
# 
# ggplot(compslopes1995, aes(x = pca_scores, y = comp_sensitivity)) +
#   geom_hline(yintercept = 0, lty = 3, color = 'gray50', size = 1.5) +
#   geom_point(aes(size = log10(n_indiv)), pch = 1) +
#   geom_abline(slope=sens_lm1995$coef[2], intercept = sens_lm1995$coef[1], color = 'blue', size = 2) +
#   scale_size_continuous(name = 'N', breaks=c(2,3,4), labels=10^c(2,3,4)) +
#   labs(x = 'Shade tolerance', y = 'Production sensitivity') + ggtitle('Competition sensitivity, 1990')

#####
# DENSITY SCALINGS BY DIAMETER

# Pareto fits
# All years, all combinations.



for (i in allyears_names) {
  dat <- get(i)
  assign(paste0(i, '_dens_fit'), lapply(dat, function(x) powerlawfit(x$dbh_corr)))
}

for (i in c(names1990, names1995)) {
  dat <- get(i)
  assign(paste0(i, '_dens_fit'), powerlawfit(dat$dbh_corr))
}

# Plot density scalings by diameter

# ***insert code here***

#####
# DENSITY SCALINGS BY LIGHT RECEIVED (only 1990 and 1995)

for (i in c(names1990, names1995)) {
  dat <- get(i)
  assign(paste0(i, '_lightdens_fit'), powerlawfit(dat$light_received))
}

# ***insert plotting code here***

#####
# DENSITY SCALING BY DIAMETER, WITH CUTOFF
# all years

library(stats4)

pareto_cutoff_fits <- function(dataframe, variable, L_init = 1) {
  x <<- dataframe[,variable]
  fit_pareto <- mle(nll_powerlaw, start = list(alpha = 3), fixed = list(xmin = min(x)), method = 'BFGS')
  fit_cutoff <- mle(nll_powerlaw_cutoff2, start = list(alpha = 3, L = L_init), fixed = list(xmin = min(x)), method = 'BFGS')
  aic_pareto <- AICc(n = length(x), k = 1, lhat = fit_pareto@min)
  aic_cutoff <- AICc(n = length(x), k = 2, lhat = fit_cutoff@min)
  return(list(fit_pareto = fit_pareto, fit_cutoff = fit_cutoff, aic_pareto = aic_pareto, aic_cutoff = aic_cutoff))
}

for (i in allyears_names) {
  dat <- get(i)
  assign(paste0(i, '_paretofits'), lapply(dat, pareto_cutoff_fits, variable = 'dbh_corr'))
  
}

for(i in c(names1990, names1995)) {
  dat <- get(i)
  assign(paste0(i, '_paretofits'), pareto_cutoff_fits(dat, 'dbh_corr'))
}

# generate table of coefficients

extractcoeffs <- function(fit1, fit2, aic1, aic2) {
  c(alpha_pareto = fit1@coef, alpha_cutoff = fit2@coef[1], L_cutoff = fit2@coef[2], logL = log10(fit2@coef[2]), deltaaicc = aic1 - aic2)
}

coef_tables <- list()

for (i in allyears_names) {
  coef_tables[[length(coef_tables) + 1]] <- lapply(get(paste0(i,'_paretofits')), function(x) with(x, extractcoeffs(fit_pareto, fit_cutoff, aic_pareto, aic_cutoff)))
}

coef_tables_9095 <- list()

for (i in c(names1990, names1995)) {
  coef_tables_9095[[length(coef_tables_9095) + 1]] <- with(get(paste0(i,'_paretofits')), extractcoeffs(fit_pareto, fit_cutoff, aic_pareto, aic_cutoff))
}

# generate bootstrap confidence interval around coefficients

nb <- 99
qprobs <- c(0.025, 0.975)

for (i in allyears_names) {
  dat <- get(i)
  assign(paste0(i, '_paretoboot'), lapply(dat, function(x) boot_mle(xboot = x$dbh_corr, nboot = nb)))
  assign(paste0(i, '_paretobootci'), lapply(get(paste0(i, '_paretoboot')), function(x) apply(x, 2, quantile, probs = qprobs)))
}

for (i in c(names1990, names1995)) {
  dat <- get(i)
  assign(paste0(i, '_paretoboot'), boot_mle(xboot = dat$dbh_corr, nboot = nb))
  assign(paste0(i, '_paretobootci'), apply(get(paste0(i, '_paretoboot')), 2, quantile, probs = qprobs))
}

# plot results

#####
# DENSITY SCALING BY LIGHT RECEIVED, WITH CUTOFF

# fit parameters with and without cutoff

# L_init must be set to a high value to work.
for(i in c(names1990, names1995)) {
  dat <- get(i)
  assign(paste0(i, '_lightparetofits'), pareto_cutoff_fits(dat, 'light_received', L_init = 10000))
}


# generate table of coefficients

coef_tables_9095light <- list()

for (i in c(names1990, names1995)) {
  coef_tables_9095light[[length(coef_tables_9095light) + 1]] <- with(get(paste0(i,'_lightparetofits')), extractcoeffs(fit_pareto, fit_cutoff, aic_pareto, aic_cutoff))
}

# generate bootstrap confidence interval around coefficients

for (i in c(names1990, names1995)) {
  dat <- get(i)
  assign(paste0(i, '_lightparetoboot'), boot_mle(xboot = dat$light_received, nboot = nb, L_init = 10000))
  assign(paste0(i, '_lightparetobootci'), apply(get(paste0(i, '_paretoboot')), 2, quantile, probs = qprobs))
}

# plot results

#####
# TOTAL (BINNED) PRODUCTION SCALING

numbins <- 20 # Can be edited if desired. ***NOT JUST FOR LOOKS***

# Do logbin and get slopes
for (i in allyears_names) {
  dat <- get(i)
  assign(paste0(i, '_prod_logbin'), lapply(dat, function(z) logbin(x=z$dbh_corr, y=z$production, n=numbins)))
  assign(paste0(i, 'prod_lm'), lapply(get(paste0(i, '_prod_logbin')), function(x) lm(log10(bin_value) ~ log10(bin_midpoint), data=subset(x, bin_value > 0))))
}

for (i in c(names1990, names1995)) {
  dat <- get(i)
  assign(paste0(i, '_prod_logbin'), logbin(x=dat$dbh_corr, y=dat$production, n=numbins))
  assign(paste0(i, 'prod_lm'), lm(log10(bin_value) ~ log10(bin_midpoint), data=subset(get(paste0(i, '_prod_logbin')), bin_value > 0)))
}

# ***insert plotting code here***

#####
# TOTAL (BINNED) LIGHT SCALING

numbins <- 20

# Do logbin and get slopes

for (i in c(names1990, names1995)) {
  dat <- get(i)
  assign(paste0(i, '_light_logbin'), logbin(x=dat$dbh_corr, y=dat$light_received, n=numbins))
  assign(paste0(i, 'light_lm'), lm(log10(bin_value) ~ log10(bin_midpoint), data=subset(get(paste0(i, '_light_logbin')), bin_value > 0)))
}

# ***insert plotting code here



######
# DENSITY SCALING: 2-factor grouping: light grouping vs shade grouping
# 1990 and 1995 only
# Use even-axis shade grouping and quantile shade groupings.

threelevelbins <- function(datanames, nbins=10) {
  dat_shade <- get(datanames[1])
  dat_int <- get(datanames[2])
  dat_gap <- get(datanames[3])
  
  # Get max and min, split into 10 bins
  bin_lims <- c(floor(min(dat_shade$dbh_corr)), ceiling(max(dat_shade$dbh_corr)))
  bin_edges <- seq(log10(bin_lims[1]), log10(bin_lims[2]), length.out=nbins+1)
  
  bin_dims <- data.frame(bin_min = 10^bin_edges[1:nbins], bin_max = 10^bin_edges[2:(nbins+1)])
  bin_dims <- transform(bin_dims, bin_midpoint = (bin_min+bin_max)/2)
  
  shade_lowlight_binset2 <-  with(subset(dat_shade, light_group == 'Low light'), logbin_setedges(x=dbh_corr, y=NULL, edges=bin_dims))
  int_lowlight_binset2 <- with(subset(dat_int, light_group == 'Low light'), logbin_setedges(x=dbh_corr, y=NULL, edges=bin_dims))
  gap_lowlight_binset2 <- with(subset(dat_gap, light_group == 'Low light'), logbin_setedges(x=dbh_corr, y=NULL, edges=bin_dims))
  
  shade_intlight_binset2 <-  with(subset(dat_shade, light_group == 'Intermediate light'), logbin_setedges(x=dbh_corr, y=NULL, edges=bin_dims))
  int_intlight_binset2 <- with(subset(dat_int, light_group == 'Intermediate light'), logbin_setedges(x=dbh_corr, y=NULL, edges=bin_dims))
  gap_intlight_binset2 <- with(subset(dat_gap, light_group == 'Intermediate light'), logbin_setedges(x=dbh_corr, y=NULL, edges=bin_dims))
  
  shade_highlight_binset2 <-  with(subset(dat_shade, light_group == 'High light'), logbin_setedges(x=dbh_corr, y=NULL, edges=bin_dims))
  int_highlight_binset2 <- with(subset(dat_int, light_group == 'High light'), logbin_setedges(x=dbh_corr, y=NULL, edges=bin_dims))
  gap_highlight_binset2 <- with(subset(dat_gap, light_group == 'High light'), logbin_setedges(x=dbh_corr, y=NULL, edges=bin_dims))
  
  # Concatenate each light level into one df.
  lowlight_bins <- rbind(transform(shade_lowlight_binset2, guild = 'shade'),
                         transform(int_lowlight_binset2, guild = 'intermediate'),
                         transform(gap_lowlight_binset2, guild = 'gap'))
  intlight_bins <- rbind(transform(shade_intlight_binset2, guild = 'shade'),
                         transform(int_intlight_binset2, guild = 'intermediate'),
                         transform(gap_intlight_binset2, guild = 'gap'))
  highlight_bins <- rbind(transform(shade_highlight_binset2, guild = 'shade'),
                          transform(int_highlight_binset2, guild = 'intermediate'),
                          transform(gap_highlight_binset2, guild = 'gap'))
  
  lowlight_bins <- subset(lowlight_bins, bin_value > 0)
  intlight_bins <- subset(intlight_bins, bin_value > 0)
  highlight_bins <- subset(highlight_bins, bin_value > 0)
  
  return(list(lowlight_bins, intlight_bins, highlight_bins))
}

threebins_even_1990 <- threelevelbins(datanames = names1990[2:4])
threebins_quantile_1990 <- threelevelbins(datanames = names1990[6:8])
threebins_even_1995 <- threelevelbins(datanames = names1995[2:4])
threebins_quantile_1995 <- threelevelbins(datanames = names1995[6:8])

#####
# weighted averages of shade, intermediate, and gap.

# *** insert code here ***

# Save gigantic workspace for later plotting
save.image(file = 'allscalings.RData')

#####
# individual-level scatterplots of light received by production for both years

# *** insert code here ***
