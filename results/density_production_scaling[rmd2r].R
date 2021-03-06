##### NOTE: On line 26, you may need to change the file path. Otherwise it should all work.


#' ---	
#' title: "Density and production scaling relationships"	
#' author: "Quentin D. Read"	
#' date: "March 21, 2017"	

#' ## Version history	
#' 	
#' * 21 March, part 2: added density scaling by diameter in addition to by mass, to check against Enquist's numbers. Also did a separate version of the density scaling with trees that had 0 growth thrown out.	
#' * 21 March: all relevant plots are changed to histograms, and the individual production scaling is added. There is also a table at the bottom of the document.	
#' * 20 March: changed plots to PDF instead of CDF, edited axis ticks.	
#' * 17 March: changed axis labels and made y-axes logarithmic	
#' * 16 March: created document	
#'
#' # Methods	
#' 	
#' ## Preparing data	
#' 	
#' The biomass data are year 3. The production data are the increment of biomass between years 1 and 3, halved.	
#' 	
#' 	

library(dplyr)
dat <- read.csv('C:/Users/Q/Google Drive/ForestLight/data/forestlightmaster.csv', stringsAsFactors = FALSE)	
	
pdat <- dat %>% 	
  mutate(diameter3 = pmax(Year3_DGH, Year3_DBH, na.rm=TRUE),	
         biomass3 = pmax(BiomassDGH_year3_kg, BiomassDBH_year3_kg, na.rm=TRUE),	
         biomass2 = pmax(BiomassDGH_year2_kg, BiomassDBH_year2_kg, 0, na.rm=TRUE),	
         biomass1 = pmax(BiomassDGH_year1_kg, BiomassDBH_year1_kg, 0, na.rm=TRUE),	
         massprod23 = pmax(biomass3 - biomass2, 0, na.rm=TRUE),	
         massprod12 = pmax(biomass2 - biomass1, 0, na.rm=TRUE),	
         massprod13 = pmax((biomass3 - biomass1)/2, 0, na.rm=TRUE)) %>% 	
  select(Site, diameter3, biomass3, biomass2, massprod23, massprod12, massprod13, CII, tolerance)	
bcidat <- pdat %>% filter(Site == 'Barro_Colorado_Island', !is.na(CII), !is.na(biomass3))	
harvdat <- pdat %>% filter(Site == 'Harvard_Forest_LTER', !is.na(CII), !is.na(biomass3))	
#' 	
#' 	
#' ## Functions needed to run analysis	
#' 	
#' This function fits both continuous Pareto and continuous lognormal to the data (input as a single vector of values), and returns their parameters plus all data needed to plot them. This returns the pdf as well.	
#' 	
#' 	
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
  	
  return(list(plotdat = plotdat, 	
              plfit = plfit_dat, 	
              lognormfit = lognormfit_dat, 	
              plpdf = data.frame(x = plfit_dat$x, y = pl_pdf),	
              lognormpdf = data.frame(x = lognormfit_dat$x, y = lognorm_pdf),	
              xmin = xmin_pl, 	
              alpha = pars_pl$pars,	
              xmin_lognorm = xmin_lognorm,	
              pars_lognorm = pars_lognorm$pars))	
}	
#' 	
#' 	
#' This function implements Ethan White's log-binning algorithm. For our purposes, $x$ is biomass, $y$ is production, and $n$ is the desired number of bins. Modified on 20 March to also include 1-dimensional log-binning algorithm (for the density scaling bins).	
#' 	
#' 	
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
#' 	
#' 	
#' This is a verbal description of the binning algorithm given in R code above.	
#' 	
#' 1. Log transform the year 3 biomass data	
#' 2. Bin into equal width bins (on the logarithmic axis)	
#' 3. Back transform the bin edges to linear axis.	
#' 4. Get the midpoint of each bin. This is the X-value on the plot.	
#' 5. Sum the biomass production of each bin, expressed as (year 3 biomass - year 1 biomass)/2	
#' 6. Divide each bin's summed production value by the linear width of each bin (this will be the Y-value on the plot).	




	
#' This function plots the log bins as a histogram. Both x and y-axes are logarithmic. It can also plot a log-linear regression fit.	
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
                label = paste('Slope:', round(lm(I(log10(bin_value)) ~ I(log10(bin_midpoint)), data=dat)$coef[2], 2)),	
                hjust = 0, vjust = -0.25)	
  }	
  return(p)	
}	
#' 	
#' 	
#' This function, added on 20 March, plots the PDFs (not CDFs) of the Pareto and lognormal fits on a histogram of the data. It is necessary to make bins because that is the only way to plot the raw data behind the PDF.	
#' 	
#' 	
plotbinsandfits <- function(pl, bindat, xl, yl = 'log PDF', title, subtitle) {	
  library(ggplot2)	
  th <- theme_bw() + theme(panel.grid = element_blank())	
  	
  expr1 <- as.character(as.expression(substitute(	
    "Pareto:"~~alpha == a, list(a = round(pl$alpha, 2)))))	
  expr2 <- as.character(as.expression(substitute(	
    "Lognormal:"~~mu == m*","~~sigma == s, list(m = round(pl$pars_lognorm[1], 2), 	
                                                s = round(pl$pars_lognorm[2], 2)))))	
  	
  bindat <- transform(bindat, bin_value = bin_value/sum(bin_count))	
  bindat <- subset(bindat, bin_value > 0)	
  y_min <- 10^floor(log10(min(bindat$bin_value, na.rm = TRUE)))	
  y_max <- max(bindat$bin_value, na.rm = TRUE) * 1.1	
  	
  p <- ggplot(bindat) + 	
    scale_x_log10(name = xl,	
                  breaks = scales::trans_breaks("log10", function(x) 10^x),	
                  labels = scales::trans_format("log10", scales::math_format(10^.x)),	
                  expand = c(0,0)) +	
    scale_y_log10(name = yl,	
                  limits = c(y_min, y_max),	
                  breaks = scales::trans_breaks("log10", function(x) 10^x),	
                  labels = scales::trans_format("log10", scales::math_format(10^.x)),	
                  expand=c(0,0)) +	
    geom_rect(aes(xmin = bin_min, xmax = bin_max, ymin = 0, ymax = bin_value), alpha = 0.5) +	
    geom_line(data = subset(pl$plpdf, y>0), aes(x,y), color = 'forestgreen') +	
    geom_line(data = subset(pl$lognormpdf, y>0), aes(x,y), color = 'indianred') +	
    geom_text(x = -Inf, y = -Inf, label = expr1, parse = TRUE, hjust = 0, vjust = -1.5) +	
    geom_text(x = -Inf, y = -Inf, label = expr2, parse = TRUE, hjust = 0, vjust = -0.25) +	
    th +	
    labs(x = xl, y = 'log PDF') +	
    ggtitle(title, subtitle)	
  return(p)	
}	
#' 	
#' 	
#' # Results	
#' 	
#' ## Density scaling relationships	
#' 	
#' Below, I fit all the power laws at once, then plot the results. **Note: For all plots below, the green line is the Pareto fit and the red line is the lognormal fit.** Note here that the Pareto $\alpha$ parameter is equivalent to -slope, so, for example, $\alpha = 1.35$ corresponds to a slope of -1.35.	
#' 	
#' 	
numbins <- 20	
	
bci_abund_all <- powerlawfit(bcidat$biomass3)	
bci_abund_shade <- powerlawfit(subset(bcidat, tolerance == 'S')$biomass3)	
bci_abund_inter <- powerlawfit(subset(bcidat, tolerance == 'I')$biomass3)	
bci_abund_gap <- powerlawfit(subset(bcidat, tolerance == 'G')$biomass3)	
harv_abund_all <- powerlawfit(harvdat$biomass3)	
harv_abund_shade <- powerlawfit(subset(harvdat, tolerance == 'S')$biomass3)	
harv_abund_inter <- powerlawfit(subset(harvdat, tolerance == 'I')$biomass3)	
harv_abund_gap <- powerlawfit(subset(harvdat, tolerance == 'G')$biomass3)	
bci_logbin_all <- with(bcidat, logbin(x=biomass3, y=NULL, n = numbins))	
bci_logbin_shade <- with(subset(bcidat, tolerance == 'S'), logbin(x=biomass3, y=NULL, n = numbins))	
bci_logbin_inter <- with(subset(bcidat, tolerance == 'I'), logbin(x=biomass3, y=NULL, n = 15))	
bci_logbin_gap <- with(subset(bcidat, tolerance == 'G'), logbin(x=biomass3, y=NULL, n = 15))	
harv_logbin_all <- with(harvdat, logbin(x=biomass3, y=NULL, n = numbins))	
harv_logbin_shade <- with(subset(harvdat, tolerance == 'S'), logbin(x=biomass3, y=NULL, n = numbins))	
harv_logbin_inter <- with(subset(harvdat, tolerance == 'I'), logbin(x=biomass3, y=NULL, n = 10))	
harv_logbin_gap <- with(subset(harvdat, tolerance == 'G'), logbin(x=biomass3, y=NULL, n = 8))	
#' 	
#' 	
#' 	
label1 <- 'Individual mass (kg)'	
	
plotbinsandfits(bci_abund_all, bci_logbin_all, label1, 'log PDF', 'Density scaling, BCI', 'all species')	
plotbinsandfits(bci_abund_shade, bci_logbin_shade, label1, 'log PDF', 'Density scaling, BCI', 'shade-tolerant species')	
plotbinsandfits(bci_abund_inter, bci_logbin_inter, label1,'log PDF', 'Density scaling, BCI', 'intermediate species')	
plotbinsandfits(bci_abund_gap, bci_logbin_gap, label1,'log PDF', 'Density scaling, BCI', 'gap species')	
plotbinsandfits(harv_abund_all, harv_logbin_all, label1, 'log PDF', 'Density scaling, Harvard', 'all species')	
plotbinsandfits(harv_abund_shade, harv_logbin_shade, label1, 'log PDF', 'Density scaling, Harvard', 'shade-tolerant species')	
plotbinsandfits(harv_abund_inter, harv_logbin_inter, label1, 'log PDF', 'Density scaling, Harvard', 'intermediate species')	
plotbinsandfits(harv_abund_gap, harv_logbin_gap, label1, 'log PDF', 'Density scaling, Harvard', 'gap species')	
#' 	
#' 	
#' ## Density scaling relationships by diameter	
#' 	
#' Added 21 March. This is mathematically identical to the density scaling by biomass. The slope for BCI should be around -2 if it matches Enquist's method.	
#' 	
#' 	
numbins <- 20	
	
bci_abund_all <- powerlawfit(bcidat$diameter3)	
bci_abund_shade <- powerlawfit(subset(bcidat, tolerance == 'S')$diameter3)	
bci_abund_inter <- powerlawfit(subset(bcidat, tolerance == 'I')$diameter3)	
bci_abund_gap <- powerlawfit(subset(bcidat, tolerance == 'G')$diameter3)	
harv_abund_all <- powerlawfit(harvdat$diameter3)	
harv_abund_shade <- powerlawfit(subset(harvdat, tolerance == 'S')$diameter3)	
harv_abund_inter <- powerlawfit(subset(harvdat, tolerance == 'I')$diameter3)	
harv_abund_gap <- powerlawfit(subset(harvdat, tolerance == 'G')$diameter3)	
bci_logbin_all <- with(bcidat, logbin(x=diameter3, y=NULL, n = numbins))	
bci_logbin_shade <- with(subset(bcidat, tolerance == 'S'), logbin(x=diameter3, y=NULL, n = numbins))	
bci_logbin_inter <- with(subset(bcidat, tolerance == 'I'), logbin(x=diameter3, y=NULL, n = 15))	
bci_logbin_gap <- with(subset(bcidat, tolerance == 'G'), logbin(x=diameter3, y=NULL, n = 15))	
harv_logbin_all <- with(harvdat, logbin(x=diameter3, y=NULL, n = numbins))	
harv_logbin_shade <- with(subset(harvdat, tolerance == 'S'), logbin(x=diameter3, y=NULL, n = numbins))	
harv_logbin_inter <- with(subset(harvdat, tolerance == 'I'), logbin(x=diameter3, y=NULL, n = 10))	
harv_logbin_gap <- with(subset(harvdat, tolerance == 'G'), logbin(x=diameter3, y=NULL, n = 8))	
#' 	
#' 	
labeld <- 'Diameter (cm)'	
	
plotbinsandfits(bci_abund_all, bci_logbin_all, labeld, 'log PDF', 'Diameter scaling, BCI', 'all species')	
plotbinsandfits(bci_abund_shade, bci_logbin_shade, labeld, 'log PDF', 'Diameter scaling, BCI', 'shade-tolerant species')	
plotbinsandfits(bci_abund_inter, bci_logbin_inter, labeld,'log PDF', 'Diameter scaling, BCI', 'intermediate species')	
plotbinsandfits(bci_abund_gap, bci_logbin_gap, labeld,'log PDF', 'Diameter scaling, BCI', 'gap species')	
plotbinsandfits(harv_abund_all, harv_logbin_all, labeld, 'log PDF', 'Diameter scaling, Harvard', 'all species')	
plotbinsandfits(harv_abund_shade, harv_logbin_shade, labeld, 'log PDF', 'Diameter scaling, Harvard', 'shade-tolerant species')	
plotbinsandfits(harv_abund_inter, harv_logbin_inter, labeld, 'log PDF', 'Diameter scaling, Harvard', 'intermediate species')	
plotbinsandfits(harv_abund_gap, harv_logbin_gap, labeld, 'log PDF', 'Diameter scaling, Harvard', 'gap species')	
#' 	
#' 	
#' 	
#' ## Density scaling with zero-growth trees excluded	
#' 	
#' People who previously analyzed this dataset, such as Enquist, probably threw out the trees with zero biomass production recorded (whether it was measurement error or because the tree was dead/dying). Since they must be thrown out for the production scaling fits, it makes sense to also throw them out for the density scaling fits. Here are the fits from above, for both biomass and diameter, repeated with the zero-production trees thrown out before fitting. Note that the presence of zero values for production makes no difference for the *binned* production fits, because the zeroes don't add anything to the bin values. 	
#' 	
#' ### By mass	
#' 	
#' 	
numbins <- 20	
	
# Throw out trees with zeroes in the production column.	
bcisub <- subset(bcidat, massprod13 > 0)	
harvsub <- subset(harvdat, massprod13 > 0)	
	
bci_abund_all <- powerlawfit(bcisub$biomass3)	
bci_abund_shade <- powerlawfit(subset(bcisub, tolerance == 'S')$biomass3)	
bci_abund_inter <- powerlawfit(subset(bcisub, tolerance == 'I')$biomass3)	
bci_abund_gap <- powerlawfit(subset(bcisub, tolerance == 'G')$biomass3)	
harv_abund_all <- powerlawfit(harvsub$biomass3)	
harv_abund_shade <- powerlawfit(subset(harvsub, tolerance == 'S')$biomass3)	
harv_abund_inter <- powerlawfit(subset(harvsub, tolerance == 'I')$biomass3)	
harv_abund_gap <- powerlawfit(subset(harvsub, tolerance == 'G')$biomass3)	
bci_logbin_all <- with(bcisub, logbin(x=biomass3, y=NULL, n = numbins))	
bci_logbin_shade <- with(subset(bcisub, tolerance == 'S'), logbin(x=biomass3, y=NULL, n = numbins))	
bci_logbin_inter <- with(subset(bcisub, tolerance == 'I'), logbin(x=biomass3, y=NULL, n = 15))	
bci_logbin_gap <- with(subset(bcisub, tolerance == 'G'), logbin(x=biomass3, y=NULL, n = 15))	
harv_logbin_all <- with(harvsub, logbin(x=biomass3, y=NULL, n = numbins))	
harv_logbin_shade <- with(subset(harvsub, tolerance == 'S'), logbin(x=biomass3, y=NULL, n = numbins))	
harv_logbin_inter <- with(subset(harvsub, tolerance == 'I'), logbin(x=biomass3, y=NULL, n = 10))	
harv_logbin_gap <- with(subset(harvsub, tolerance == 'G'), logbin(x=biomass3, y=NULL, n = 8))	
#' 	
#' 	
#' \newpage	
#' 	
#' 	
label1 <- 'Individual mass (kg)'	
	
plotbinsandfits(bci_abund_all, bci_logbin_all, label1, 'log PDF', 'Density scaling, BCI', 'all species (no zeroes)')	
plotbinsandfits(bci_abund_shade, bci_logbin_shade, label1, 'log PDF', 'Density scaling, BCI', 'shade-tolerant species (no zeroes)')	
plotbinsandfits(bci_abund_inter, bci_logbin_inter, label1,'log PDF', 'Density scaling, BCI', 'intermediate species (no zeroes)')	
plotbinsandfits(bci_abund_gap, bci_logbin_gap, label1,'log PDF', 'Density scaling, BCI', 'gap species (no zeroes)')	
plotbinsandfits(harv_abund_all, harv_logbin_all, label1, 'log PDF', 'Density scaling, Harvard', 'all species (no zeroes)')	
plotbinsandfits(harv_abund_shade, harv_logbin_shade, label1, 'log PDF', 'Density scaling, Harvard', 'shade-tolerant species (no zeroes)')	
plotbinsandfits(harv_abund_inter, harv_logbin_inter, label1, 'log PDF', 'Density scaling, Harvard', 'intermediate species (no zeroes)')	
plotbinsandfits(harv_abund_gap, harv_logbin_gap, label1, 'log PDF', 'Density scaling, Harvard', 'gap species (no zeroes)')	
#' 	
#' 	
#' \newpage	
#' 	
#' ### By diameter	
#' 	
#' 	
numbins <- 20	
	
# Throw out trees with zeroes in the production column.	
bcisub <- subset(bcidat, massprod13 > 0)	
harvsub <- subset(harvdat, massprod13 > 0)	
	
bci_abund_all <- powerlawfit(bcisub$diameter3)	
bci_abund_shade <- powerlawfit(subset(bcisub, tolerance == 'S')$diameter3)	
bci_abund_inter <- powerlawfit(subset(bcisub, tolerance == 'I')$diameter3)	
bci_abund_gap <- powerlawfit(subset(bcisub, tolerance == 'G')$diameter3)	
harv_abund_all <- powerlawfit(harvsub$diameter3)	
harv_abund_shade <- powerlawfit(subset(harvsub, tolerance == 'S')$diameter3)	
harv_abund_inter <- powerlawfit(subset(harvsub, tolerance == 'I')$diameter3)	
harv_abund_gap <- powerlawfit(subset(harvsub, tolerance == 'G')$diameter3)	
bci_logbin_all <- with(bcisub, logbin(x=diameter3, y=NULL, n = numbins))	
bci_logbin_shade <- with(subset(bcisub, tolerance == 'S'), logbin(x=diameter3, y=NULL, n = numbins))	
bci_logbin_inter <- with(subset(bcisub, tolerance == 'I'), logbin(x=diameter3, y=NULL, n = 15))	
bci_logbin_gap <- with(subset(bcisub, tolerance == 'G'), logbin(x=diameter3, y=NULL, n = 15))	
harv_logbin_all <- with(harvsub, logbin(x=diameter3, y=NULL, n = numbins))	
harv_logbin_shade <- with(subset(harvsub, tolerance == 'S'), logbin(x=diameter3, y=NULL, n = numbins))	
harv_logbin_inter <- with(subset(harvsub, tolerance == 'I'), logbin(x=diameter3, y=NULL, n = 10))	
harv_logbin_gap <- with(subset(harvsub, tolerance == 'G'), logbin(x=diameter3, y=NULL, n = 8))	
#' 	
#' 	
#' 	
labeld <- 'Diameter (cm)'	
	
plotbinsandfits(bci_abund_all, bci_logbin_all, labeld, 'log PDF', 'Diameter scaling, BCI', 'all species (no zeroes)')	
plotbinsandfits(bci_abund_shade, bci_logbin_shade, labeld, 'log PDF', 'Diameter scaling, BCI', 'shade-tolerant species (no zeroes)')	
plotbinsandfits(bci_abund_inter, bci_logbin_inter, labeld,'log PDF', 'Diameter scaling, BCI', 'intermediate species (no zeroes)')	
plotbinsandfits(bci_abund_gap, bci_logbin_gap, labeld,'log PDF', 'Diameter scaling, BCI', 'gap species (no zeroes)')	
plotbinsandfits(harv_abund_all, harv_logbin_all, labeld, 'log PDF', 'Diameter scaling, Harvard', 'all species (no zeroes)')	
plotbinsandfits(harv_abund_shade, harv_logbin_shade, labeld, 'log PDF', 'Diameter scaling, Harvard', 'shade-tolerant species (no zeroes)')	
plotbinsandfits(harv_abund_inter, harv_logbin_inter, labeld, 'log PDF', 'Diameter scaling, Harvard', 'intermediate species (no zeroes)')	
plotbinsandfits(harv_abund_gap, harv_logbin_gap, labeld, 'log PDF', 'Diameter scaling, Harvard', 'gap species (no zeroes)')	
#' 	
#' 	
#' 	
#' ## Individual production by individual mass scaling relationships	
#' 	
#' This has been modified on 21 March to use the simpler method of fitting a regression to the log-transformed production and biomass values for each tree. However, the major caveat to this set of relationships is that the power law cannot be fit to 0 values. The 0 values probably represent trees with very low growth, but it is impossible to say. Therefore, I had to throw them out, but they comprise a pretty big portion of the trees (23% of trees in BCI and 20% of trees in HF).	
#' 	
#' 	
bci_prod_all <- lm(log10(massprod13) ~ log10(biomass3), data = bcidat, subset = massprod13 > 0)	
bci_prod_shade <- lm(log10(massprod13) ~ log10(biomass3), data = bcidat, subset = massprod13 > 0 & tolerance == 'S')	
bci_prod_inter <- lm(log10(massprod13) ~ log10(biomass3), data = bcidat, subset = massprod13 > 0 & tolerance == 'I')	
bci_prod_gap <- lm(log10(massprod13) ~ log10(biomass3), data = bcidat, subset = massprod13 > 0 & tolerance == 'G')	
harv_prod_all <- lm(log10(massprod13) ~ log10(biomass3), data = harvdat, subset = massprod13 > 0)	
harv_prod_shade <- lm(log10(massprod13) ~ log10(biomass3), data = harvdat, subset = massprod13 > 0 & tolerance == 'S')	
harv_prod_inter <- lm(log10(massprod13) ~ log10(biomass3), data = harvdat, subset = massprod13 > 0 & tolerance == 'I')	
harv_prod_gap <- lm(log10(massprod13) ~ log10(biomass3), data = harvdat, subset = massprod13 > 0 & tolerance == 'G')	
#' 	
#' 	
#' For each site, I plotted each of the regression relationships by guild. Here, I have excluded the points with no tolerance score from displaying on the plot. However, the overall regression line does include those points in its calculation.	
#' 	
#' 	
	
bcislopes <- data.frame(biomass3 = -Inf, 	
                        massprod13 = Inf, 	
                        txt = c(paste('All trees:',round(bci_prod_all$coef[2], 2)),	
                                paste('Shade-tolerant:',round(bci_prod_shade$coef[2], 2)),	
                                paste('Intermediate:',round(bci_prod_inter$coef[2], 2)),	
                                paste('Gap:',round(bci_prod_gap$coef[2], 2))	
                                ))	
harslopes <- data.frame(biomass3 = -Inf, 	
                        massprod13 = Inf, 	
                        txt = c(paste('All trees:',round(harv_prod_all$coef[2], 2)),	
                                paste('Shade-tolerant:',round(harv_prod_shade$coef[2], 2)),	
                                paste('Intermediate:',round(harv_prod_inter$coef[2], 2)),	
                                paste('Gap:',round(harv_prod_gap$coef[2], 2))	
                                ))	
	
ggplot(subset(bcidat, massprod13 > 0), aes(x = biomass3, y = massprod13)) +	
  geom_point(aes(color = tolerance), alpha = 0.5, data = subset(bcidat, massprod13 > 0 & !is.na(tolerance))) +	
  stat_smooth(aes(group = tolerance, color = tolerance), method = 'lm', se = FALSE, data = subset(bcidat, massprod13 > 0 & !is.na(tolerance))) +	
  stat_smooth(method = 'lm', se = FALSE, color = 'black') +	
  scale_x_log10(name = 'Individual mass (kg)',	
                breaks = scales::trans_breaks("log10", function(x) 10^x),	
                labels = scales::trans_format("log10", scales::math_format(10^.x)),	
                expand = c(0,0)) +	
  scale_y_log10(name = expression(paste('Individual production (kg y'^-1,')', sep = '')),	
                breaks = scales::trans_breaks("log10", function(x) 10^x),	
                labels = scales::trans_format("log10", scales::math_format(10^.x))) +	
  geom_text(x = -Inf, y = Inf, label = bcislopes$txt[1], hjust = 0, vjust = 1) +	
  geom_text(x = -Inf, y = Inf, label = bcislopes$txt[2], hjust = 0, vjust = 2) +	
  geom_text(x = -Inf, y = Inf, label = bcislopes$txt[3], hjust = 0, vjust = 3) +	
  geom_text(x = -Inf, y = Inf, label = bcislopes$txt[4], hjust = 0, vjust = 4) +	
  theme_bw() + theme(panel.grid = element_blank()) +	
  ggtitle('Individual production scaling', 'BCI')	
	
ggplot(subset(harvdat, massprod13 > 0), aes(x = biomass3, y = massprod13)) +	
  geom_point(aes(color = tolerance), alpha = 0.5, data = subset(harvdat, massprod13 > 0 & !is.na(tolerance))) +	
  stat_smooth(aes(group = tolerance, color = tolerance), method = 'lm', se = FALSE, data = subset(harvdat, massprod13 > 0 & !is.na(tolerance))) +	
  stat_smooth(method = 'lm', se = FALSE, color = 'black') +	
  scale_x_log10(name = 'Individual mass (kg)',	
                breaks = scales::trans_breaks("log10", function(x) 10^x),	
                labels = scales::trans_format("log10", scales::math_format(10^.x)),	
                expand = c(0,0)) +	
  scale_y_log10(name = expression(paste('Individual production (kg y'^-1,')', sep = '')),	
                breaks = scales::trans_breaks("log10", function(x) 10^x),	
                labels = scales::trans_format("log10", scales::math_format(10^.x))) +	
  geom_text(x = -Inf, y = Inf, label = harslopes$txt[1], hjust = 0, vjust = 1) +	
  geom_text(x = -Inf, y = Inf, label = harslopes$txt[2], hjust = 0, vjust = 2) +	
  geom_text(x = -Inf, y = Inf, label = harslopes$txt[3], hjust = 0, vjust = 3) +	
  geom_text(x = -Inf, y = Inf, label = harslopes$txt[4], hjust = 0, vjust = 4) +	
  theme_bw() + theme(panel.grid = element_blank()) +	
  ggtitle('Individual production scaling', 'Harvard Forest')	
#' 	
#' 	
#' ## Total binned production scaling relationships	
#' 	
#' For binning, I tried to create 20 bins if possible. Then, I simply fit a regression line to the relationship between $\log_{10} production$ and $\log_{10} mass$ which is the slope that we should use to evaluate the predictions of EE. Previously I tried to do a power law fit after binning, but that was dumb because power law fit is only for unbinned data!	
#' 	
#' 	
numbins <- 20	
	
bci_logbin <- with(bcidat, logbin(x=biomass3, y=massprod13, n=numbins))	
bci_logbin_shade <- with(subset(bcidat, tolerance == 'S'), logbin(x=biomass3, y=massprod13, n=numbins))	
bci_logbin_inter <- with(subset(bcidat, tolerance == 'I'), logbin(x=biomass3, y=massprod13, n=15))	
bci_logbin_gap <- with(subset(bcidat, tolerance == 'G'), logbin(x=biomass3, y=massprod13, n=15))	
harv_logbin <- with(harvdat, logbin(x=biomass3, y=massprod13, n=numbins))	
harv_logbin_shade <- with(subset(harvdat, tolerance == 'S'), logbin(x=biomass3, y=massprod13, n=numbins))	
harv_logbin_inter <- with(subset(harvdat, tolerance == 'I'), logbin(x=biomass3, y=massprod13, n=10))	
harv_logbin_gap <- with(subset(harvdat, tolerance == 'G'), logbin(x=biomass3, y=massprod13, n=8))	
	
bci_lm_all <- lm(log10(bin_value) ~ log10(bin_midpoint), data = bci_logbin)	
bci_lm_shade <- lm(log10(bin_value) ~ log10(bin_midpoint), data = bci_logbin_shade)	
bci_lm_inter <- lm(log10(bin_value) ~ log10(bin_midpoint), data = bci_logbin_inter)	
bci_lm_gap <- lm(log10(bin_value) ~ log10(bin_midpoint), data = bci_logbin_gap)	
harv_lm_all <- lm(log10(bin_value) ~ log10(bin_midpoint), data = harv_logbin)	
harv_lm_shade <- lm(log10(bin_value) ~ log10(bin_midpoint), data = harv_logbin_shade)	
harv_lm_inter <- lm(log10(bin_value) ~ log10(bin_midpoint), data = harv_logbin_inter)	
harv_lm_gap <- lm(log10(bin_value) ~ log10(bin_midpoint), data = harv_logbin_gap)	
	
#' 	
#' 	
label2 <- expression(paste('Total production (kg y'^-1,')', sep = ''))	
	
plotlogbin(bci_logbin, label1, label2, 'Production bins, BCI', 'all species', reg = TRUE)	
plotlogbin(bci_logbin_shade, label1, label2, 'Production bins, BCI', 'shade-tolerant species', reg = TRUE)	
plotlogbin(bci_logbin_inter, label1, label2, 'Production bins, BCI', 'intermediate species', reg = TRUE)	
plotlogbin(bci_logbin_gap, label1, label2, 'Production bins, BCI', 'gap species', reg = TRUE)	
plotlogbin(harv_logbin, label1, label2, 'Production bins, Harvard', 'all species', reg = TRUE)	
plotlogbin(harv_logbin_shade, label1, label2, 'Production bins, Harvard', 'shade-tolerant species', reg = TRUE)	
plotlogbin(harv_logbin_inter, label1, label2, 'Production bins, Harvard', 'intermediate species', reg = TRUE)	
plotlogbin(harv_logbin_gap, label1, label2, 'Production bins, Harvard', 'gap species', reg = TRUE)	
#' 	
#' 	
