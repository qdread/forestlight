# Create production and light bins for use in plotting production vs light.


# Observed data -----------------------------------------------------------


fpdata <- 'C:/Users/Q/google_drive/ForestLight/data/data_22jan2018'
load(file.path(fpdata, 'rawdataobj_22jan.r'))

library(dplyr)

dat90 <- alltree_light_90 %>%
  select(dbh_corr, production, light_received, crownarea, fg) %>%
  mutate(production_area = production/crownarea, light_area = light_received/crownarea) %>%
  mutate(fg = if_else(is.na(fg), 'unclassified', paste0('fg', fg)))

dat95 <- alltree_light_95 %>%
  select(dbh_corr, production, light_received, crownarea, fg) %>%
  mutate(production_area = production/crownarea, light_area = light_received/crownarea) %>%
  mutate(fg = if_else(is.na(fg), 'unclassified', paste0('fg', fg)))

source('code/allfunctions27july.r')

fakebin_across_years <- function(dat_values, dat_classes, edges, mean_type = 'geometric', n_census = 5) {
  qprobs <- c(0.025, 0.25, 0.5, 0.75, 0.975)
  # add some padding just in case
  mins <- edges$bin_min
  mins[1] <- 0
  maxes <- edges$bin_max
  maxes[length(maxes)] <- Inf
  
  binstats <- t(sapply(1:length(mins), function(x) {
    indivs <- dat_values[dat_classes >= mins[x] & dat_classes < maxes[x]]
    if (mean_type == 'geometric') {
      mean_n <- exp(mean(log(indivs)))
      sd_n <- exp(sd(log(indivs)))
      ci_width <- qnorm(0.975) * sd(log(indivs)) / sqrt(length(indivs))
      ci_min <- exp(mean(log(indivs)) - ci_width)
      ci_max <- exp(mean(log(indivs)) + ci_width)
      
    } else {
      mean_n <- mean(indivs)
      sd_n <- sd(indivs)
      ci_width <- qnorm(0.975) * sd(indivs) / sqrt(length(indivs))
      ci_min <- mean_n - ci_width
      ci_max <- mean_n + ci_width
    }
    c(mean_n_individuals = length(indivs) / n_census,
      mean = mean_n, 
      sd = sd_n,
      quantile(indivs, probs = qprobs),
      ci_min = ci_min,
      ci_max = ci_max)
  }))
  dimnames(binstats)[[2]] <- c('mean_n_individuals','mean', 'sd', 'q025', 'q25', 'median', 'q75', 'q975', 'ci_min', 'ci_max')
  data.frame(bin_midpoint = edges$bin_midpoint,
             bin_min = edges$bin_min,
             bin_max = edges$bin_max,
             binstats)
}

numbins <- 20
light_bins_9095 <- logbin(x = c(dat90$light_area, dat95$light_area), n = numbins)

prod_light_bin_90 <- fakebin_across_years(dat_values = dat90$production_area,
                                       dat_classes = dat90$light_area,
                                       edges = light_bins_9095,
                                       n_census = 1)
prod_light_bin_95 <- fakebin_across_years(dat_values = dat95$production_area,
                                          dat_classes = dat95$light_area,
                                          edges = light_bins_9095,
                                          n_census = 1)

prod_light_bin_fg_90 <- dat90 %>%
  group_by(fg) %>%
  do(bin = fakebin_across_years(dat_values = .$production_area,
                                dat_classes = .$light_area,
                                edges = light_bins_9095,
                                n_census = 1))

prod_light_bin_fg_90 <- cbind(fg = rep(c('fg1','fg2','fg3','fg4','fg5','unclassified'), each = numbins),
                           do.call(rbind, prod_light_bin_fg_90$bin)) %>%
  filter(complete.cases(.))

prod_light_bin_fg_95 <- dat95 %>%
  group_by(fg) %>%
  do(bin = fakebin_across_years(dat_values = .$production_area,
                                dat_classes = .$light_area,
                                edges = light_bins_9095,
                                n_census = 1))

prod_light_bin_fg_95 <- cbind(fg = rep(c('fg1','fg2','fg3','fg4','fg5','unclassified'), each = numbins),
                              do.call(rbind, prod_light_bin_fg_95$bin)) %>%
  filter(complete.cases(.))

prod_light_bin_all <- rbind(
  cbind(year = 1990, fg = 'alltree', prod_light_bin_90),
  cbind(year = 1990, prod_light_bin_fg_90),
  cbind(year = 1995, fg = 'alltree', prod_light_bin_95),
  cbind(year = 1995, prod_light_bin_fg_95)
)

write.csv(prod_light_bin_all, file = 'C:/Users/Q/google_drive/ForestLight/data/data_forplotting_light_12apr2018/obs_light_binned.csv', row.names = FALSE)

write.csv(rbind(
  cbind(year = 1990, dat90),
  cbind(year = 1995, dat95)
), 
file = 'C:/Users/Q/google_drive/ForestLight/data/data_forplotting_light_12apr2018/obs_light_raw.csv', 
row.names = FALSE)


# Predicted data ----------------------------------------------------------

# Already created in lightbyarea_predci_fg.csv because there is no need to split things up or divide.
