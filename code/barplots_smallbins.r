source('~/GitHub/forestlight/code/mergeshade.r')
n_per_bin <- 10 # minimal number of individuals per bin
bin_width <- 1 # in kilograms


pdat <- dat %>% 
  mutate(biomass = pmax(BiomassDGH_year3_kg, BiomassDBH_year3_kg, na.rm=TRUE),
         biomass2 = pmax(BiomassDGH_year2_kg, BiomassDBH_year2_kg, 0, na.rm=TRUE),
         massprod = pmax(biomass - biomass2, 0, na.rm=TRUE)) %>% 
  select(Site, biomass, massprod, CII)
bcidat <- pdat %>% filter(Site == 'Barro_Colorado_Island', !is.na(CII), !is.na(biomass))
harvdat <- pdat %>% filter(Site == 'Harvard_Forest_LTER', !is.na(CII), !is.na(biomass))

# Biomass binning: using midpoints, and ensuring that at least 10 individuals are in each bin.
# Make a function, and then apply it to both bcidat and harvdat.

biomass_binby10 <- function(x) {
  # First, create 1 kg bins. *bin_width is set to 1*
  biomassrange <- range(x$biomass)
  bin_min <- seq(0, floor(biomassrange[2]), by = bin_width)
  bin_max <- bin_min + bin_width
  bin1kg <- rep(0, nrow(x))
  for (i in 1:nrow(x)) {
    bin1kg[i] <- which(x$biomass[i] > bin_min & x$biomass[i] <= bin_max)
  }
  
  x <- x %>% mutate(bin1kg = bin1kg) %>% arrange(biomass)
  
  # Next, group the bins until at least 10individuals are reached in each bin. *n_per_bin set to 10*
  
  bin10indiv <- rep(1, nrow(x))
  bin_id <- 1
  bin_indiv <- 0
  cumulative_total_indiv <- 1
  
  for (i in 1:nrow(x)) {
    bin_indiv <- bin_indiv + 1 # Increment number of individuals in bin by 1
    # If tree i has crossed the 1 kg boundary, check whether 10 individuals have been reached. If so, make a new bin.
    if (x$bin1kg[i] > bin_max[bin_id] & bin_indiv >= n_per_bin) {
      bin10indiv[(cumulative_total_indiv -  1):(cumulative_total_indiv - 1 + bin_indiv)] <- bin_id
      cumulative_total_indiv <- cumulative_total_indiv + bin_indiv
      bin_id <- bin_id + 1
      bin_indiv <- 0
    }
    # If we've reached the final row without exactly reaching a multiple of 10, put the remaining <10 largest trees into their own bin.
    if (i == nrow(x) & bin_indiv > 0) {
      bin10indiv[(cumulative_total_indiv - 1):nrow(x)] <- bin_id
    }
  }
  
  # For each of these bins containing a minimum of 10 individuals, do the following steps:
  
  x %>% 
    mutate(bin10indiv = bin10indiv) %>% 
    group_by(bin10indiv) %>% 
    summarize(production_sum = sum(massprod),                  # Gets raw biomass sum for the bin.
              median_biomass = median(biomass),                # Finds median biomass of bin.
              mid_biomass = (floor(min(biomass)) + ceiling(max(biomass)))/2, # Find biomass midpoint of the bin.
              n_1kgbins = 1 + max(bin1kg) - min(bin1kg),       # Determines how many 1kg bins are contained within the bin.
              production_perbin = production_sum / n_1kgbins,  # Divides raw sum by number of bins.
              CII_mean = mean(CII),                            # Gets "mean" CII.
              CII_se = sd(CII)/sqrt(length(CII)),              # Gets "standard error" of CII.
              n_individuals = n())                             # Returns number of individual trees in the bin.
  
}

# Evaluate the function for BCI and Harvard.
bci_bybin <- biomass_binby10(bcidat)
harv_bybin <- biomass_binby10(harvdat)

# Plot them, using a log scale for biomass.
library(ggplot2)

sc_y <- scale_y_continuous(limits = c(0,5.5), expand = c(0,0))
sc_x <- scale_x_continuous(name = 'Biomass (kg)')
sc_xlog <- scale_x_log10(name = 'Biomass (kg)', breaks = c(0,1,10,100,1000))
th_bar <- theme_bw() + theme(panel.grid = element_blank())

# This involves fakery with averages. Can't really average the ordinal category but did for figure.

ggplot(bci_bybin, aes(x = mid_biomass, y = CII_mean, ymin = CII_mean - CII_se, ymax = CII_mean + CII_se)) +
  geom_pointrange() +
  sc_xlog + sc_y + th_bar + ggtitle('BCI')

ggplot(harv_bybin, aes(x = mid_biomass, y = CII_mean, ymin = CII_mean - CII_se, ymax = CII_mean + CII_se)) +
  geom_pointrange() +
  sc_xlog + sc_y + th_bar + ggtitle('Harvard Forest')

ggplot(bci_bybin, aes(x = mid_biomass, y = production_perbin)) +
  geom_point() +
  scale_y_log10(name = 'Biomass production (kg/y)', expand = c(0,0), limits = c(0.2, 150)) +
  sc_xlog + th_bar + ggtitle('BCI')

ggplot(harv_bybin, aes(x = mid_biomass, y = production_perbin)) +
  geom_point() +
  scale_y_log10(name = 'Biomass production (kg/y)', limits = c(0.2 ,45), expand = c(0,0)) +
  sc_xlog + th_bar + ggtitle('Harvard Forest')

# Summary info.

print(round(bci_bybin,2), n = nrow(bci_bybin))
print(round(harv_bybin,2), n = nrow(harv_bybin))

