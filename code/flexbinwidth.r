# Improved function for binning in which you can set the width and number per bin

biomass_bin <- function(x, year, bin_width = 1, n_per_bin = 10) {
  # Set the right year columns.
  if (year == 23) {
    x <- transform(x, biomass=biomass3, massprod=massprod23)
  } else if (year == 12) {
    x <- transform(x, biomass=biomass2, massprod=massprod12)
  } else {
    x <- transform(x, biomass=biomass3, massprod=massprod13)
  }
  
  # First, create small bins. Default bin width is 1 kg but it will work for others.
  biomassrange <- range(x$biomass)
  bin_min <- seq(0, plyr::round_any(biomassrange[2], bin_width, floor), by = bin_width)
  bin_max <- bin_min + bin_width
  bin1kg <- rep(0, nrow(x))
  for (i in 1:nrow(x)) {
    # Note: changed this to have the bottom open, not top
    bin1kg[i] <- bin_width * which(x$biomass[i] >= bin_min & x$biomass[i] < bin_max)
  }
  
  x <- x %>% mutate(bin1kg = bin1kg) %>% arrange(biomass)
  
  # Next, group the bins until at least n_per_bin individuals are reached in each bin. 
  # *n_per_bin set to 10 by default*
  
  bin10indiv <- rep(1, nrow(x))
  bin_id <- 1
  bin_indiv <- 0
  cumulative_total_indiv <- 1
  max_of_bin <- bin_width
  
  for (i in 1:nrow(x)) {
    bin_indiv <- bin_indiv + 1 # Increment number of individuals in bin by 1
    # If tree i has crossed into the next bin, 
    # check whether n_per_bin individuals have been reached. If so, make a new bin.
    if (x$bin1kg[i] > max_of_bin & bin_indiv >= n_per_bin) {
      bin10indiv[(cumulative_total_indiv -  1):(cumulative_total_indiv - 1 + bin_indiv)] <- bin_id
      cumulative_total_indiv <- cumulative_total_indiv + bin_indiv
      bin_id <- bin_id + 1
      bin_indiv <- 0
      max_of_bin <- plyr::round_any(x$bin1kg[i], bin_width, ceiling)
    }
    # If we've reached the final row without exactly reaching a multiple of n_per_bin, 
    # put the remaining largest trees into their own bin.
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
              mid_biomass = (plyr::round_any(min(biomass), bin_width, floor) + plyr::round_any(max(biomass), bin_width, ceiling))/2, # Find biomass midpoint of the bin.
              n_smallbins = 1 + max(bin1kg) - min(bin1kg),       # Determines how many virtual bins are contained within the bin.
              production_perbin = production_sum / n_smallbins,  # Divides raw sum by number of bins.
              CII_mean = mean(CII),                            # Gets "mean" CII.
              CII_se = sd(CII)/sqrt(length(CII)),              # Gets "standard error" of CII.
              n_individuals = n())                             # Returns number of individual trees in the bin.
  
 
}


# Run the binning algorithm with bin width decreasing.

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
  select(Site, biomass3, biomass2, diameter3, diameter2, massprod23, massprod12, massprod13, CII)
bcidat <- pdat %>% filter(Site == 'Barro_Colorado_Island', !is.na(CII), !is.na(biomass3))
harvdat <- pdat %>% filter(Site == 'Harvard_Forest_LTER', !is.na(CII), !is.na(biomass3))


bin_widths <- c(10, 5, 1, 0.5, 0.1, 0.05, 0.01)
bin_ns <- c(2, 5, 10, 50, 100)

bci_out <- list()
harv_out <- list()

for (w in bin_widths) {
  for (n in bin_ns) {
    bci_out[[length(bci_out) + 1]] <- cbind(biomass_bin(bcidat, year=23, bin_width=w, n_per_bin = n), bin_width=w, n_per_bin=n)
    harv_out[[length(harv_out) + 1]] <- cbind(biomass_bin(harvdat, year=23, bin_width=w, n_per_bin = n), bin_width=w, n_per_bin=n)
  }
}

bci_outdf <- do.call('rbind', bci_out)
harv_outdf <- do.call('rbind', harv_out)

library(ggplot2)

sc_ylog <- scale_y_log10(name = 'Biomass production (kg/y)')
sc_xlog <- scale_x_log10(name = 'Biomass (kg)', breaks = c(0.1,1,10,100,1000))
th_scatter <- theme_bw() + theme(panel.grid = element_blank())

ggplot(bci_outdf, aes(x=mid_biomass, y=production_perbin)) +
  geom_point() +
  facet_grid(bin_width ~ n_per_bin) +
  sc_ylog + sc_xlog + th_scatter

ggplot(harv_outdf, aes(x=mid_biomass, y=production_perbin)) +
  geom_point() +
  facet_grid(bin_width ~ n_per_bin) +
  sc_ylog + sc_xlog + th_scatter