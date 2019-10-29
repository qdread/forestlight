# Bin by size instead of light, corresponds to fig 5 in the MS

# Source code for log binning function
source(file.path(fpdata, 'allfunctions27july.r'))

# Set number of bins       
numbins <- 20

# Log bin density
allyeardbh_alltree <- unlist(lapply(alltreedat[2:6], '[', , 'dbh_corr'))
allyeardbh_shade <- unlist(lapply(shadedat[2:6], '[', , 'dbh_corr'))
allyeardbh_gap <- unlist(lapply(gapdat[2:6], '[', , 'dbh_corr'))
allyeardbh_unclassified <- unlist(lapply(unclassifieddat[2:6], '[', , 'dbh_corr'))

dbhbin_alltree <- logbin(x = allyeardbh_alltree, y = NULL, n = numbins)
dbhbin_shade <- logbin(x = allyeardbh_shade, y = NULL, n = numbins)
dbhbin_gap <- logbin(x = allyeardbh_gap, y = NULL, n = numbins)
dbhbin_unclassified <- logbin(x = allyeardbh_unclassified, y = NULL, n = numbins)



# Have to rebin production and density because we have to have the same cutoffs for shade and gap so the numbers are comparable.
totalprodbin_shade_byyear_samebins <- lapply(shadedat[2:6], function(z) logbin_setedges(x = z$dbh_corr, y = z$production, edges = dbhbin_alltree))
totalprodbin_gap_byyear_samebins <- lapply(gapdat[2:6], function(z) logbin_setedges(x = z$dbh_corr, y = z$production, edges = dbhbin_alltree))

densitybin_shade_byyear_samebins <- lapply(shadedat[2:6], function(z) logbin_setedges(x = z$dbh_corr, y = NULL, edges = dbhbin_alltree))
densitybin_gap_byyear_samebins <- lapply(gapdat[2:6], function(z) logbin_setedges(x = z$dbh_corr, y = NULL, edges = dbhbin_alltree))

shadegap_stats <- list()

for (i in 1:5) {
  shadegap_stats[[i]] <- data.frame(bin = 1:numbins,
                                    year = c(1990,1995,2000,2005,2010)[i],
                                    shade_gap_production_ratio = totalprodbin_shade_byyear_samebins[[i]]$bin_value/totalprodbin_gap_byyear_samebins[[i]]$bin_value,
                                    shade_gap_density_ratio = densitybin_shade_byyear_samebins[[i]]$bin_value/densitybin_gap_byyear_samebins[[i]]$bin_value)  
}

shadegap_stats <- do.call('rbind', shadegap_stats)

shadegap_stats[shadegap_stats == Inf | is.na(shadegap_stats)] <- NA

library(dplyr)
library(cowplot)

shadegap_stats_bin_5census <- shadegap_stats %>% 
  group_by(bin) %>%
  summarize(production_ratio_mean = mean(shade_gap_production_ratio),
            density_ratio_mean = mean(shade_gap_density_ratio),
            production_ratio_min = min(shade_gap_production_ratio),
            production_ratio_max = max(shade_gap_production_ratio),
            density_ratio_min = min(shade_gap_density_ratio),
            density_ratio_max = max(shade_gap_density_ratio)) %>%
  cbind(dbhbin_alltree[,c(1,4,5)])

shadegap_stats_bin_5census %>%
  ggplot(aes(x = bin_midpoint, y = density_ratio_mean, ymin = density_ratio_min, ymax = density_ratio_max)) +
  geom_pointrange() +
  scale_x_log10(name = 'Diameter (cm)') + 
  scale_y_log10(name = 'Shade to gap density ratio') +
  panel_border(colour = 'black')
