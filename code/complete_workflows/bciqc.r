# Quality Control run on raw BCI data
# Follows Meakem et al. 2017
# QDR 11 July 2017

# Load Condit's BCI data and Nadja's light data.

fp <- '~/google_drive/ForestLight/data/BCI_raw'

growth8590 <- read.delim(file.path(fp, 'BCI_light/growth_final8590.txt'), stringsAsFactors = FALSE)
growth9095 <- read.delim(file.path(fp, 'BCI_light/growth_final9095.txt'), stringsAsFactors = FALSE)

# Load census data. 1 = 1980, 2 = 1985, 3 = 1990, 4 = 1995, 5 = 2000, 6 = 2005, 7 = 2010.

for (i in 1:7) {
  load(file.path(fp, paste0('bcidata/bci.full', i, '.rdata')))
  load(file.path(fp, paste0('bcidata/bci.stem', i, '.rdata')))
}

# Get rid of tree ferns and strangler figs
# Tree ferns are already removed
species_to_remove <- c('ficubu','ficuc1','ficuc2','ficuci','ficupe','ficupo')
remove_spp <- function(dat) subset(dat, !sp %in% c(species_to_remove, toupper(species_to_remove)))

for (i in 1:7) {
  assign(paste0('bci.full',i), remove_spp(get(paste0('bci.full',i))))
  assign(paste0('bci.stem',i), remove_spp(get(paste0('bci.stem',i))))
}

growth8590 <- remove_spp(growth8590)
growth9095 <- remove_spp(growth9095)

# Apply taper correction for trees not measured at dbh because of their buttress, returning an estimate of what the true dbh would be.
# See Cushman reference for the formula
# Additional data files helpfully provided by KC Cushman, 13 July 2017.
# Details on WSG measurement protocol are linked on this page: http://www.forestgeo.si.edu/group/Plant+Functional+Traits/

wsg.data <- read.csv(file.path(fp, 'bci_taper/wsg.data.csv'), stringsAsFactors = FALSE)
taper.parameters <- read.csv(file.path(fp, 'bci_taper/taper.parameters.csv'), stringsAsFactors = FALSE)

# These are the species that need to be corrected for taper.
taper_spp <- c("alsebl", "anacex", "brosal", "cavapl", "ceibpe", "diptpa", "huracr", "other", "quaras", "tab2ar", "tet2pa")        

#    Make a vector of palm species with individuals measured above 1.3 m height (not 
#    included for taper correction)
palmsp <- c("elaeol","sch1zo","socrex")


# Values of the b1 parameter for each of the species that need to be corrected for taper. See Table 3, Cushman et al. 2014.


#    Define function to calculate equivalent diameter at 1.3 m given a diameter measurement 
#    in cm (d) taken at a nonstandard height in m(h), and a taper parameter (b1). This is 
#    the overall best taper model from Metcalf et el. (2009)
apply.eqn1 <- function(d,h,b1) {d/(exp(-b1*(h-1.3)))}

#    Define functions to estimate total tree height (m) from diameter (cm) using 
#    uncorrected (notaper.H) or taper-corrected (taper.H) allometries parameterized from 
#    data from the BCI plot. The function is the weibull function as described in 
#    Feldpausch et al. 2012
notaper.H <- function(dbh) {43.1427*(1-
                                       exp(-0.04003*dbh^0.82585))
}
taper.H <- function(dbh) {43.4375*(1-
                                     exp(-0.04469*dbh^0.78339))
}

#    Define function (agb.allometry) to estimate aboveground biomass from wood specific 
#    gravity (wsg), tree diameter in centimeters (dbh), and total tree height in meters (H) 
#    using allometries from Chave et al. (2005) for most tropical forests. Also define 
#    function (agb.allometryb) that does not included total tree height. 
agb.allometry <- function(wsg,dbh,H) {0.0509*wsg*((dbh)^2)*H}

library(MASS)
library(nlme)
library(tidyverse)

model1.pars <- taper.parameters[taper.parameters$eqn==1 & taper.parameters$b1 >0,]

#  The best model includes 2010 diameter, 2010 height of measurement, and species group
varmodel <- gls(log(b1)~log(diam.2010)+log(hom.2010)+spb, data=model1.pars, method="ML")
#    Extract coefficients and covariance matrix from the model
varmodel.coefficients <- varmodel$ coefficients
varmodel.cormatrix <- summary(varmodel)$ varBeta
#    Extract coefficients for relationship with 2010 diameter and height of measurement
mdbh <- varmodel$coefficients[2]
mhom <- varmodel$coefficients[3]

# To find b1 for each of the species, we need to calculate it from the equation:
# b1 <- exp(mdbh * log(dbh) + mhom * log(hom) + intercept + species coefficient)

spcoeffs <- rep(0, length(taper_spp))

for (i in 1:length(spcoeffs)) {
  
  # Assign each tree a taper parameter based on its species group, 2010 diameter, and 2010 
  # height of measurement 
     spcoeffs[i] <- ifelse(i == 1,  varmodel.coefficients[1], varmodel.coefficients[1] + varmodel.coefficients[i+2])

}

names(spcoeffs) <- taper_spp

# Vectorized function to quickly calculate the taper-corrected DBH for all individuals in one data frame.
taper_correct_dbh <- function(dat) {
  indiv_coeffs <- spcoeffs[match(dat$sp, names(spcoeffs))] # Match individuals with their species taper coefficients
  indiv_coeffs[is.na(indiv_coeffs)] <- spcoeffs['other']   # All other individuals get "other" coefficient
  b1 <- exp(mdbh * log(dat$dbh/10) + mhom * log(dat$hom) + indiv_coeffs)
  dbh_corr <- case_when(
    is.na(dat$pom) ~ as.numeric(NA),                               # If not measured, return NA
    dat$pom == 1.3 | dat$sp %in% palmsp ~ dat$dbh,                 # Don't correct if already at 1.3, or if it's a palm
    TRUE ~ 10 * apply.eqn1(d = dat$dbh / 10, h = dat$hom, b1 = b1) # Otherwise apply the correction.
  )
  return(data.frame(dbh_corr = dbh_corr))
}

# Apply dbh correction

for (i in 1:7) {
  stem_i <- get(paste0('bci.stem', i))
  full_i <- get(paste0('bci.full', i))
  dbh_corr_stem_i <- taper_correct_dbh(stem_i)
  dbh_corr_full_i <- taper_correct_dbh(full_i)
  assign(paste0('bci.stem', i), cbind(stem_i, dbh_corr_stem_i))
  assign(paste0('bci.full', i), cbind(full_i, dbh_corr_full_i))
}

# Use dbh correction and wsg table to get corrected aboveground biomass values.

for (i in 1:7) {
  assign(paste0('bci.stem', i), 
         get(paste0('bci.stem', i)) %>% 
           left_join(wsg.data) %>%
           mutate(height_corr = taper.H(dbh_corr/10),
                  agb_corr = agb.allometry(wsg, dbh_corr/10, height_corr)/1000))
  assign(paste0('bci.full', i), 
         get(paste0('bci.full', i)) %>% 
           left_join(wsg.data) %>%
           mutate(height_corr = taper.H(dbh_corr/10),
                  agb_corr = agb.allometry(wsg, dbh_corr/10, height_corr)/1000))
}



# Remove the patch of secondary forest in the BCI plot
habitat <- read.table(file.path(fp, 'BCI_light/habitats2.txt'), header = TRUE, stringsAsFactors = FALSE)

# Define area of young forest
young_forest <- subset(habitat, Habitat == 'young')
with(young_forest, plot(GX20, GY20))

# To remove young forest, we need to round all BCI trees to the nearest multiple of 20 lower than their coordinate, and see if it is in one of those.
young_forest_list <- data.frame(t(young_forest[,c('GX20','GY20')]))

# Vectorized function to assign individual stems to young forest area
is_young <- function(x, y) {
  x_quadrat <- plyr::round_any(x, 20, f = floor)
  y_quadrat <- plyr::round_any(y, 20, f = floor)
  quadrat_list <- data.frame(rbind(x_quadrat, y_quadrat))
  !is.na(match(quadrat_list, young_forest_list))
}

for (i in 1:7) {
  assign(paste0('bci.full', i), 
         get(paste0('bci.full', i)) %>% mutate(young = is_young(gx, gy)))
  assign(paste0('bci.stem', i), 
         get(paste0('bci.stem', i)) %>% mutate(young = is_young(gx, gy)))
}

growth8590 <- growth8590 %>% mutate(young = is_young(gx, gy))
growth9095 <- growth9095 %>% mutate(young = is_young(gx, gy))

# Calculate the area that Nadja has light values for, so we can determine the right denominator for the per-hectare values.
# 20-meter strip around the outside: gx < 20 | gx > 980 | gy < 20 | gy > 480

is_edge <- function(gx, gy) gx < 20 | gx > 980 | gy < 20 | gy > 480

for (i in 1:7) {
  assign(paste0('bci.full', i), 
         get(paste0('bci.full', i)) %>% mutate(edge = is_edge(gx, gy)))
  assign(paste0('bci.stem', i), 
         get(paste0('bci.stem', i)) %>% mutate(edge = is_edge(gx, gy)))
}

growth8590 <- mutate(growth8590, edge = is_edge(gx, gy))
growth9095 <- mutate(growth9095, edge = is_edge(gx, gy))

# Calculate number of hectares that is neither edge nor young.
# number of square meters in edge
edge_area <- 2 * 20 * 1000 + 2 * 20 * 460
# number of square meters in young forest
young_area <- nrow(young_forest) * 400
# number of square meters in young forest that isn't edge
young_area_core <- sum(young_forest$GX20 > 0 & young_forest$GX20 < 980 & young_forest$GY20 > 0 & young_forest$GY20 < 480) * 400

notyoung_area <- 5e5 - young_area # 48.08 ha
notyoung_area_core <- 5e5 - edge_area - young_area_core # 42.84 ha
core_area <- 5e5 - edge_area # 44.16 ha

# Save all the edited files
save(list = c('growth8590','growth9095',grep('bci.', ls(), value = TRUE)), file = file.path(fp, 'bcidata/bciqcrun.R'))
