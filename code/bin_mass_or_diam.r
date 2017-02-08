# Even more flexible binning function that also allows you to select diameter as the binning factor.

# Improved function for binning in which you can set the width and number per bin

biomass_diameter_bin <- function(x, year, bin_width = 1, n_per_bin = 10, grp = 'biomass') {
  # Set the right year columns.
  if (year == 23) {
    x <- transform(x, biomass=biomass3, diameter=diameter3, massprod=massprod23)
  } else if (year == 12) {
    x <- transform(x, biomass=biomass2, diameter=diameter2, massprod=massprod12)
  } else {
    x <- transform(x, biomass=biomass3, diameter=diameter3, massprod=massprod13)
  }
  
  # Set the grouping variable.
  if (grp == 'biomass') x <- transform(x, grpvar = biomass) else x <- transform(x, grpvar = diameter)
  
  # First, create small bins. Default bin width is 1 kg but it will work for others.
  grpvarrange <- range(x$grpvar)
  round_range <- plyr::round_any(grpvarrange, bin_width, floor)
  bin_min <- seq(round_range[1], round_range[2], by = bin_width) # Corrected to deal with diameter cutoff at 1.
  bin_max <- bin_min + bin_width
  bins_small <- rep(0, nrow(x))
  for (i in 1:nrow(x)) {
    # Note: changed this to have the bottom open, not top
    bins_small[i] <- bin_width * which(x$grpvar[i] >= bin_min & x$grpvar[i] < bin_max)
  }
  
  x <- x %>% mutate(bins_small = bins_small) %>% arrange(grpvar)
  
  
  
  # Next, group the bins until at least n_per_bin individuals are reached in each bin. 
  # *n_per_bin set to 10 by default*
  
  bin_n_indiv <- rep(1, nrow(x))
  bin_id <- 1
  bin_indiv <- 0
  cumulative_total_indiv <- 1
  max_of_bin <- bin_width
  
  for (i in 1:nrow(x)) {
    bin_indiv <- bin_indiv + 1 # Increment number of individuals in bin by 1
    # If tree i has crossed into the next bin, 
    # check whether n_per_bin individuals have been reached. If so, make a new bin.
    if (x$bins_small[i] > max_of_bin & bin_indiv >= n_per_bin) {
      bin_n_indiv[(cumulative_total_indiv -  1):(cumulative_total_indiv - 1 + bin_indiv)] <- bin_id
      cumulative_total_indiv <- cumulative_total_indiv + bin_indiv
      bin_id <- bin_id + 1
      bin_indiv <- 0
      max_of_bin <- plyr::round_any(x$bins_small[i], bin_width, ceiling)
    }
    # If we've reached the final row without exactly reaching a multiple of n_per_bin, 
    # put the remaining largest trees into their own bin.
    if (i == nrow(x) & bin_indiv > 0) {
      bin_n_indiv[(cumulative_total_indiv - 1):nrow(x)] <- bin_id
    }
  }
  
  # For each of these bins containing a minimum of n individuals, do the following steps:
  
  x %>% 
    mutate(bin_n_indiv = bin_n_indiv) %>% 
    group_by(bin_n_indiv) %>% 
    summarize(production_sum = sum(massprod),                  # Gets raw biomass production sum for the bin.
              median_grp = median(grpvar),                     # Finds median biomass or diameter of bin.
              mid_grp = (plyr::round_any(min(grpvar), bin_width, floor) + plyr::round_any(max(grpvar), bin_width, ceiling))/2, # Find biomass or diameter midpoint of the bin.
              n_smallbins = 1 + max(bins_small) - min(bins_small),       # Determines how many virtual bins are contained within the bin.
              production_perbin = production_sum / n_smallbins,  # Divides raw sum by number of bins.
              CII_mean = mean(CII),                            # Gets "mean" CII.
              CII_se = sd(CII)/sqrt(length(CII)),              # Gets "standard error" of CII.
              n_individuals = n())                             # Returns number of individual trees in the bin.
  
  
}

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


massbci <- biomass_diameter_bin(bcidat, year=23, bin_width=1, n_per_bin=10, grp='biomass')
diambci <- biomass_diameter_bin(bcidat, year=23, bin_width=0.5, n_per_bin=10, grp='diameter')

library(ggplot2)

sc_ylog <- scale_y_log10(name = 'Biomass production (kg/y)', breaks = c(0.1, 1, 10, 100), limits = c(0.1, 150))
sc_ylinear <- scale_y_continuous(name = 'Biomass production (kg/y)', limits = c(0,150))
sc_xbiomass <- scale_x_log10(name = 'Biomass (kg)', breaks = c(0.1,1,10,100,1000))
sc_xdiam <- scale_x_continuous(name = 'Diameter (cm)')
th_scatter <- theme_bw() + theme(panel.grid = element_blank())

ggplot(massbci, aes(x=mid_grp, y=production_perbin)) +
  geom_point() +
  sc_ylog + sc_xbiomass + th_scatter +
  ggtitle('BCI, binned by biomass, log y axis')

ggplot(diambci, aes(x=mid_grp, y=production_perbin)) +
  geom_point() +
  sc_ylog + sc_xdiam + th_scatter +
  ggtitle('BCI, binned by diameter, log y axis')

ggplot(massbci, aes(x=mid_grp, y=production_perbin)) +
  geom_point() +
  sc_ylinear + sc_xbiomass + th_scatter +
  ggtitle('BCI, binned by biomass, linear y axis')

ggplot(diambci, aes(x=mid_grp, y=production_perbin)) +
  geom_point() +
  sc_ylinear + sc_xdiam + th_scatter +
  ggtitle('BCI, binned by diameter, linear y axis')
