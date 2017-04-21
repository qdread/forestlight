# Analyze canopy

pdat <- read.csv('C:/Users/Q/Google Drive/ForestLight/data/allforestdata19April.csv', stringsAsFactors = FALSE)

powerlawfit <- function(dat) {
  library(poweRlaw)
  pl_dat <- conpl$new(dat)
  lognorm_dat <- conlnorm$new(dat)
  xmin_pl <- pl_dat$getXmin()
  xmin_lognorm <- lognorm_dat$getXmin()
  pars_pl <- estimate_pars(pl_dat)
  pars_lognorm <- estimate_pars(lognorm_dat)
  pl_dat$setPars(pars_pl)
  lognorm_dat$setPars(pars_lognorm)
  plotdat <- plot(pl_dat)
  plfit_dat <- lines(pl_dat)
  lognormfit_dat <- lines(lognorm_dat)
  pl_pdf <- dist_pdf(m = pl_dat, q = plfit_dat$x, log = FALSE)
  lognorm_pdf <- dist_pdf(m = lognorm_dat, q = lognormfit_dat$x, log = FALSE)
  
  # bootstrap confidence interval of Pareto fit
  # discard 500 burnin iterations
  # n_boot <- 1499
  # n_burn <- 500
  # pl_boot <- bootstrap(m = pl_dat, xmins = pl_dat$getXmin(), no_of_sims = n_boot)
  # boot_ci <- quantile(pl_boot$bootstraps$pars[-(1:n_burn)], probs = c(0.025, 0.975))
  
  return(list(plotdat = plotdat, 
              plfit = plfit_dat, 
              lognormfit = lognormfit_dat, 
              plpdf = data.frame(x = plfit_dat$x, y = pl_pdf),
              lognormpdf = data.frame(x = lognormfit_dat$x, y = lognorm_pdf),
              xmin = xmin_pl, 
              alpha = pars_pl$pars,
              xmin_lognorm = xmin_lognorm,
              pars_lognorm = pars_lognorm$pars
              #boot_ci = as.numeric(boot_ci)
              ))
}

logbin <- function(x, y = NULL, n) {
  logx <- log10(x)                                           # log transform x value (biomass)
  bin_edges <- seq(min(logx), max(logx), length.out = n + 1) # get edges of bins
  logxbin <- rep(NA, length(logx))                           # create data structure to assign trees to bins
  b <- bin_edges                                             # add a little to the biggest bin temporarily
  b[length(b)] <- b[length(b)] + 1                           # (so that the biggest single tree is put in a bin)
  for (i in 1:length(logx)) {
    logxbin[i] <- sum(logx[i] >= b)                          # assign each tree to a bin
  }
  bin_midpoints <- numeric(n)
  for (i in 1:n) {
    bin_midpoints[i] <- mean(10^(bin_edges[i:(i+1)]))        # backtransform bin edges to linear, and get midpoints
  }
  bin_widths <- diff(10^bin_edges)                           # get linear width of each bin
  bin_factor <- factor(logxbin, levels=1:n)                  # convert bin to factor (required to deal with zeroes if present)
  bin_counts <- table(bin_factor)                            # find number of trees in each bin
  if (!is.null(y)) {
    rawy <- tapply(y, bin_factor, sum)                       # sum y value (production) in each bin
    rawy[is.na(rawy)] <- 0                                   # add zeroes back in if present
    bin_values <- as.numeric(rawy/bin_widths)                # divide production by width for each bin 
  }
  else {
    bin_values <- as.numeric(bin_counts/bin_widths)          # 1-dimensional case.
  }
  
  return(data.frame(bin_midpoint = bin_midpoints,            # return result!
                    bin_value = bin_values,                  # also add bin min and max for bar plot purposes
                    bin_count = as.numeric(bin_counts),
                    bin_min = 10^bin_edges[1:n],
                    bin_max = 10^bin_edges[2:(n+1)]))
  
}

plotlogbin <- function(dat, xl, yl, title, subtitle, reg = FALSE) {
  library(ggplot2)
  th <- theme_bw() + theme(panel.grid = element_blank())
  
  y_min <- 10^floor(log10(min(dat$bin_value, na.rm = TRUE)))
  y_max <- max(dat$bin_value, na.rm = TRUE) * 1.1
  
  p <- ggplot(dat, aes(xmin=bin_min, xmax=bin_max, ymin=0, ymax=bin_value)) + 
    geom_rect(alpha = 0.5) +
    scale_x_log10(name = xl, expand = c(0,0),
                  breaks = scales::trans_breaks("log10", function(x) 10^x),
                  labels = scales::trans_format("log10", scales::math_format(10^.x))) +
    scale_y_log10(name = yl, expand = c(0,0), limits = c(y_min, y_max),
                  breaks = scales::trans_breaks("log10", function(x) 10^x),
                  labels = scales::trans_format("log10", scales::math_format(10^.x))) +
    th + ggtitle(title, subtitle)
  if (reg) {
    p <- p +
      stat_smooth(method = 'lm', se = FALSE, color = 'black', 
                  aes(x = bin_midpoint, y = bin_value)) +
      geom_text(x = -Inf, y = -Inf, 
                label = paste('Slope:', 
                              round(lm(I(log10(bin_value)) ~ I(log10(bin_midpoint)), data=dat)$coef[2], 2)),
                hjust = 0, vjust = -0.25)
  }
  return(p)
}

# Stacked by shade tolerance class, show the total production by biomass or diameter.

# Show the slope of individual production by biomass or diameter. One is predicted from a simple regression, and the other is predicted from a regression that also includes ppfd as a predictor.

#### Regs without ppfd
bci_prod_all <- lm(log10(massprod13) ~ log10(diameter3), 
                   data = bcidat, subset = massprod13 > 0)
bci_prod_shade <- lm(log10(massprod13) ~ log10(diameter3), 
                     data = bcidat, subset = massprod13 > 0 & tolerance == 'S')
bci_prod_inter <- lm(log10(massprod13) ~ log10(diameter3),
                     data = bcidat, subset = massprod13 > 0 & tolerance == 'I')
bci_prod_gap <- lm(log10(massprod13) ~ log10(diameter3), 
                   data = bcidat, subset = massprod13 > 0 & tolerance == 'G')
harv_prod_all <- lm(log10(massprod13) ~ log10(diameter3), 
                    data = harvdat, subset = massprod13 > 0)
harv_prod_shade <- lm(log10(massprod13) ~ log10(diameter3), 
                      data = harvdat, subset = massprod13 > 0 & tolerance == 'S')
harv_prod_inter <- lm(log10(massprod13) ~ log10(diameter3), 
                      data = harvdat, subset = massprod13 > 0 & tolerance == 'I')
harv_prod_gap <- lm(log10(massprod13) ~ log10(diameter3), 
                    data = harvdat, subset = massprod13 > 0 & tolerance == 'G')

#### Regs with ppfd
bci_prod_all_ppfd <- lm(log10(massprod13) ~ log10(diameter3) + log10(PPFD), 
                        data = pdat, subset = massprod13 > 0 & Site == 'Barro_Colorado_Island')
bci_prod_all <- lm(log10(massprod13) ~ log10(diameter3), 
                        data = pdat, subset = massprod13 > 0 & Site == 'Barro_Colorado_Island')
bci_prod_shade_ppfd <- lm(log10(massprod13) ~ log10(diameter3) + log10(PPFD), 
                        data = pdat, subset = massprod13 > 0 & Site == 'Barro_Colorado_Island' & pooledtolerance == 'S')
bci_prod_shade <- lm(log10(massprod13) ~ log10(diameter3), 
                   data = pdat, subset = massprod13 > 0 & Site == 'Barro_Colorado_Island' & pooledtolerance == 'S')
bci_prod_gap_ppfd <- lm(log10(massprod13) ~ log10(diameter3) + log10(PPFD), 
                          data = pdat, subset = massprod13 > 0 & Site == 'Barro_Colorado_Island' & pooledtolerance == 'G')
bci_prod_gap <- lm(log10(massprod13) ~ log10(diameter3), 
                     data = pdat, subset = massprod13 > 0 & Site == 'Barro_Colorado_Island' & pooledtolerance == 'G')

# Plot predicted values under each regression

prodppfd <- predict(bci_prod_gap_ppfd)
prodnoppfd <- predict(bci_prod_gap)

with(subset(pdat, massprod13 > 0 & Site == 'Barro_Colorado_Island' & pooledtolerance == 'G'), 
     plot(log10(diameter3), log10(massprod13)))
abline(bci_prod_gap, col='red')
abline(bci_prod_gap_ppfd, col='blue')

with(subset(pdat, massprod13 > 0 & Site == 'Barro_Colorado_Island' & pooledtolerance == 'S'), 
     plot(log10(diameter3), log10(massprod13)))
abline(bci_prod_shade, col='red')
abline(bci_prod_shade_ppfd, col='blue')

confint(bci_prod_gap)
confint(bci_prod_gap_ppfd)
confint(bci_prod_shade)
confint(bci_prod_shade_ppfd)

# Attempt to do sequentially
bci_prod_gap_light <- lm(log10(massprod13) ~ log10(PPFD), data=bcigapdat)
lightresid <- bci_prod_gap_light$residuals
bci_prod_gap_both <- lm(bci_prod_gap_light$residuals ~ log10(diameter3), data=bcigapdat)

bcigapdat <- subset(pdat, massprod13 > 0 & Site == 'Barro_Colorado_Island' & pooledtolerance == 'G' & !is.na(PPFD))
bci_prod_gap_ppfd <- lm(log10(massprod13) ~ log10(diameter3) + log10(PPFD), data=bcigapdat)
bci_prod_gap <- lm(log10(massprod13) ~ log10(diameter3), data=bcigapdat)

bci_prod_gap <- broom::augment(bci_prod_gap)
bci_prod_gap_ppfd <- broom::augment(bci_prod_gap_ppfd)

gapplotdat <- data.frame(logdiameter = bci_prod_gap$log10.diameter3.,
                         logproduction = bci_prod_gap$log10.massprod13.,
                         fitted1 = bci_prod_gap$.fitted,
                         fitted2 = bci_prod_gap_ppfd$.fitted)

ggplot(gapplotdat, aes(x=logdiameter)) +
  geom_point(aes(y=logproduction)) +
  geom_line(aes(y=fitted1), color = 'red') +
  geom_line(aes(y=fitted2), color = 'blue') +
  theme_minimal()

with(bcigapdat, plot(log10(diameter3), log10(massprod13)))
abline(bci_prod_gap, col='red')
abline(bci_prod_gap_both, col='blue')

# Check order of variables
bci_prod_gap_ppfd <- lm(log10(massprod13) ~ log10(diameter3) + log10(Openness), data=bcigapdat)
bci_prod_gap_ppfd2 <- lm(log10(massprod13) ~ log10(Openness) + log10(diameter3), data=bcigapdat)

anova(bci_prod_gap_ppfd)
anova(bci_prod_gap_ppfd2)
summary(bci_prod_gap_ppfd)
summary(bci_prod_gap_ppfd2)

# Mixed model incl species as a random effect
pdat$sciname <- paste(pdat$Genus,pdat$Species)
bcigapmixed1 <- lmer(log10(massprod13) ~ log10(biomass3) + (1|sciname), data = pdat,
                    subset = massprod13 > 0 & Site == 'Barro_Colorado_Island' & tolerance == 'G' & !is.na(PPFD))
bcigapmixed2 <- lmer(log10(massprod13) ~ log10(biomass3) + log10(Openness) + (1|sciname), data = pdat,
     subset = massprod13 > 0 & Site == 'Barro_Colorado_Island' & tolerance == 'G' & !is.na(PPFD))

summary(bcigapmixed1)
confint(bcigapmixed1, method='boot', nsim=999)
summary(bcigapmixed2)
confint(bcigapmixed2, method='boot', nsim=999)

pist1 <- lm(log10(massprod13) ~ log10(biomass3), data = pdat,
                     subset = massprod13 > 0 & sciname=='Pinus strobus' & !is.na(PPFD))
pist2 <- lm(log10(massprod13) ~ log10(biomass3) + log10(Openness), data = pdat,
                     subset = massprod13 > 0 & sciname=='Pinus strobus' & !is.na(PPFD))

# t test of slope deltas by species

slope_delta <- function(x) {
  lm1 <- lm(log10(massprod13) ~ log10(diameter3), data=x)
  lm2 <- lm(log10(massprod13) ~ log10(diameter3) + log10(Openness), data=x)
  data.frame(slope1=lm1$coef[2], slope2=lm2$coef[2], delta=lm1$coef[2] - lm2$coef[2], n_indiv = nrow(x))
}

deltas <- pdat %>%
  filter(!is.na(PPFD), massprod13 > 0, !is.na(tolerance), sciname != ' ') %>%
  group_by(Site, tolerance, sciname) %>%
  do(slope_delta(.)) %>%
  filter(n_indiv >= 10)

print(deltas,n=40)

ggplot(deltas, aes(x=tolerance, y=delta)) + 
  geom_boxplot() + facet_wrap(~ Site) + 
  theme_bw() + theme(panel.grid=element_blank())

########
# 20 April. Include all trees, and include shade tolerance group as a term in the model.

lmdata <- subset(pdat, massprod13 > 0 & Site == 'Barro_Colorado_Island' & !is.na(pooledtolerance))

bci_prod_all_ppfd <- lm(log10(massprod13) ~ log10(diameter3) + tolerance : log10(PPFD) + log10(PPFD), 
                        data = pdat, subset = massprod13 > 0 & Site == 'Barro_Colorado_Island' & !is.na(tolerance))

# 21 April. Use the wright pca score as a continuous predictor.
bci_prod_all_ppfd <- lm(log10(massprod13) ~ log10(diameter3) + pca_scores : log10(PPFD), 
                        data = pdat, subset = massprod13 > 0 & Site == 'Barro_Colorado_Island' & !is.na(pca_scores))


# Show total and individual production per unit ppfd received.