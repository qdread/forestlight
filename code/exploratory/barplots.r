# Figure 2: does not require grouping by tolerance class

n_bins <- 5

# Explain how the flux was calculated later on.

pdat <- dat %>% 
  mutate(biomass = pmax(BiomassDGH_year3_kg, BiomassDBH_year3_kg, na.rm=TRUE),
         biomass2 = pmax(BiomassDGH_year2_kg, BiomassDBH_year2_kg, 0, na.rm=TRUE),
         massprod = pmax(biomass - biomass2, 0, na.rm=TRUE)) %>% 
  select(Site, biomass, massprod, CII)
bcidat <- pdat %>% filter(Site == 'Barro_Colorado_Island', !is.na(CII), !is.na(biomass))
harvdat <- pdat %>% filter(Site == 'Harvard_Forest_LTER', !is.na(CII), !is.na(biomass))

# Biomass binning (try alternative methods)

hist(bcidat$biomass, breaks = 5)
hist(harvdat$biomass, breaks = 5)

bci_bin <- with(bcidat, cut(biomass, breaks = n_bins))
harv_bin <- with(harvdat, cut(biomass, breaks = n_bins))
bci_logbin <- with(bcidat, cut(log10(biomass), breaks = n_bins))
harv_logbin <- with(harvdat, cut(log10(biomass), breaks = n_bins))

bcidat$biomass_bin <- bci_logbin
harvdat$biomass_bin <- harv_logbin

bcidat$biomass_bin_linear <- bci_bin
harvdat$biomass_bin_linear <- harv_bin


# Pull numbers from bin labels
binlab2n <- function(x, islog) {
  lx <- levels(x)
  lxs <- strsplit(lx, ',')
  mins <- sapply(lxs, '[', 1)
  maxs <- sapply(lxs, '[', 2)
  mins <- as.numeric(substr(mins, 2, nchar(mins)))
  maxs <- as.numeric(substr(maxs, 1, nchar(maxs) - 1))
  if (islog) {
    mins <- round(10^mins,1)
    maxs <- round(10^maxs,1)
  }
  c(paste('<', maxs[1]), paste(mins[2:(length(mins)-1)], 'to', maxs[2:(length(maxs)-1)]), paste('>', mins[length(mins)]))
}


library(ggplot2)

sc_y <- scale_y_continuous(limits = c(0,5.5), expand = c(0,0))
sc_y2 <- scale_y_continuous(name = 'Biomass production (kg/y)', limits = c(0,3500), expand = c(0,0))
sc_x <- scale_x_discrete(name = 'Biomass (kg)')
th_bar <- theme_bw() + theme(panel.grid = element_blank())

# This involves fakery with averages. Can't really average the ordinal category but did for figure.

ggplot(bcidat, aes(x = factor(biomass_bin, labels = binlab2n(biomass_bin, islog=T)), y = CII)) +
  stat_summary(geom = 'bar', fun.y = 'mean') +
  stat_summary( geom = 'errorbar', fun.data = 'mean_se', width = 0) +
  sc_y + sc_x + th_bar + ggtitle('BCI')

ggplot(harvdat, aes(x = factor(biomass_bin, labels = binlab2n(biomass_bin, islog=T)), y = CII)) +
  stat_summary(geom = 'bar', fun.y = 'mean') +
  stat_summary( geom = 'errorbar', fun.data = 'mean_se', width = 0) +
  sc_y + sc_x + th_bar + ggtitle('Harvard Forest')

ggplot(bcidat, aes(x = factor(biomass_bin, labels = binlab2n(biomass_bin, islog=T)), y = massprod)) +
  stat_summary(geom = 'bar', fun.y = 'sum') +
  sc_y2 + sc_x + th_bar + ggtitle('BCI')

ggplot(harvdat, aes(x = factor(biomass_bin, labels = binlab2n(biomass_bin, islog=T)), y = massprod)) +
  stat_summary(geom = 'bar', fun.y = 'sum') +
  sc_y2 + sc_x + th_bar + ggtitle('Harvard Forest')

# Summary info.
bcidat %>% 
  mutate(biomass_bin = factor(biomass_bin, labels = binlab2n(biomass_bin, islog=T))) %>% 
  group_by(biomass_bin) %>% 
  summarize(mean_cii=mean(CII), total_production=sum(massprod), n=n())

harvdat %>% 
  mutate(biomass_bin = factor(biomass_bin, labels = binlab2n(biomass_bin, islog=T))) %>% 
  group_by(biomass_bin) %>% 
  summarize(mean_cii=mean(CII), total_production=sum(massprod), n=n())



# Plots with linear scale
ggplot(bcidat, aes(x = biomass_bin_linear, y = CII)) +
  stat_summary(geom = 'bar', fun.y = 'mean') +
  stat_summary( geom = 'errorbar', fun.data = 'mean_se', width = 0) +
  sc_y + sc_x + th_bar + ggtitle('BCI')

ggplot(harvdat, aes(x = biomass_bin_linear, y = CII)) +
  stat_summary(geom = 'bar', fun.y = 'mean') +
  stat_summary( geom = 'errorbar', fun.data = 'mean_se', width = 0) +
  sc_y + sc_x + th_bar + ggtitle('Harvard Forest')

ggplot(bcidat, aes(x = biomass_bin_linear, y = massprod)) +
  stat_summary(geom = 'bar', fun.y = 'sum') +
  sc_y2 + sc_x + th_bar + ggtitle('BCI')

ggplot(harvdat, aes(x = biomass_bin_linear, y = massprod)) +
  stat_summary(geom = 'bar', fun.y = 'sum') +
  sc_y2 + sc_x + th_bar + ggtitle('Harvard Forest')