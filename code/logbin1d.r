# Plotting empirical PDF and Pareto/lognormal fits for a vector of values.

library(poweRlaw)

dat <- bcidat$biomass3

p1 <- powerlawfit(dat)

# Empirical PDF:
df1 <- p1$plotdat
df1$cdf <- 1 - df1$y

df2 <- cbind(p1$plpdf, ylognorm = p1$lognormpdf$y)

ggplot(df2, aes(x=x,y=y)) + scale_x_log10() + scale_y_log10() + geom_line() + geom_line(aes(y=ylognorm), color = 'red')

df3 <- cbind(p1$plfit, ylognorm = p1$lognormfit$y)

ggplot(df3, aes(x=x,y=y)) + scale_x_log10() + scale_y_log10() + geom_line() + geom_line(aes(y=ylognorm), color = 'red')

# Must bin the data to superimpose on the plots of PDF fits.
# Just use a simple cutting algorithm to cut into linear bins.
n <- 100
bin_edges <- seq(min(dat), max(dat), length.out = n + 1) # get edges of bins
xbin <- rep(NA, length(dat))                           # create data structure to assign trees to bins
b <- bin_edges                                             # add a little to the biggest bin temporarily
b[length(b)] <- b[length(b)] + 1                           # (so that the biggest single tree is put in a bin)
for (i in 1:length(dat)) {
  xbin[i] <- sum(dat[i] >= b)                          # assign each tree to a bin
}
bin_midpoints <- numeric(n)
for (i in 1:n) {
  bin_midpoints[i] <- mean((bin_edges[i:(i+1)]))        # backtransform bin edges to linear, and get midpoints
}
bin_widths <- diff(bin_edges)                           # get linear width of each bin

bin_factor <- factor(xbin, levels=1:n)                  # convert bin to factor (required to deal with zeroes if present)
bin_counts <- table(bin_factor)      
bin_pdf <- bin_counts/sum(bin_counts)


ggplot(df2, aes(x=x,y=y)) + scale_x_log10() + scale_y_log10() + geom_line() + geom_line(aes(y=ylognorm), color = 'red') +
  geom_step(data=data.frame(x=bin_midpoints, y = as.numeric(bin_pdf)))


# Logbin for 1 dimensional data
logbin1d <- function(x, n) {
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
  
  return(data.frame(bin_midpoint = bin_midpoints,            # divide production by width for each bin, and return result!
                    bin_value = as.numeric(bin_counts/bin_widths), # also add bin min and max for bar plot purposes
                    bin_count = as.numeric(bin_counts),
                    bin_min = 10^bin_edges[1:n],
                    bin_max = 10^bin_edges[2:(n+1)]))
  
}

ggplot(df2, aes(x=x,y=y)) + scale_x_log10() + scale_y_log10() + geom_line() + geom_line(aes(y=ylognorm), color = 'red') +
  geom_point(data=df4, aes(x=bin_midpoint, y=bin_value/sum(bin_count)))



##################################

# Function to plot the logbin1d, along with its fits

plotbinsandfits <- function(pl, bindat, xl, title, subtitle) {
  library(ggplot2)
  th <- theme_bw() + theme(panel.grid = element_blank())
  
  expr1 <- as.character(as.expression(substitute(
    "Pareto:"~~alpha == a, list(a = round(pl$alpha, 2)))))
  expr2 <- as.character(as.expression(substitute(
    "Lognormal:"~~mu == m*","~~sigma == s, list(m = round(pl$pars_lognorm[1], 2), 
                                                s = round(pl$pars_lognorm[2], 2)))))
  
  bindat <- transform(bindat, pdf_val = bin_value/sum(bin_count))
  y_min <- 10^floor(log10(min(bindat$pdf_val, na.rm = TRUE)))
  y_max <- max(bindat$pdf_val, na.rm = TRUE) * 1.1
  
  p <- ggplot(bindat) + 
    scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                  labels = scales::trans_format("log10", scales::math_format(10^.x)),
                  expand = c(0,0)) +
    scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                  labels = scales::trans_format("log10", scales::math_format(10^.x)),
                  expand=c(0,0)) +
    geom_rect(aes(xmin = bin_min, xmax = bin_max, ymin = 0, ymax = pdf_val), alpha = 0.5) +
    geom_line(data = subset(pl$plpdf, y>0), aes(x,y), color = 'forestgreen') +
    geom_line(data = subset(pl$lognormpdf, y>0), aes(x,y), color = 'indianred') +
    geom_text(x = -Inf, y = -Inf, label = expr1, parse = TRUE, hjust = 0, vjust = -1.5) +
    geom_text(x = -Inf, y = -Inf, label = expr2, parse = TRUE, hjust = 0, vjust = -0.25) +
    th +
    labs(x = xl, y = 'log PDF') +
    ggtitle(title, subtitle)
  return(p)
}



##################################################################################

# Start from scratch.
# Follow "recipe" for 1d and 2d binning.

# 1d binning:
# 1. calculate pdf from raw data, and save parameters
# 2. use 1d logbin algorithm to make a histogram out of the raw data, dividing by the number of trees to get a density plot
# 3. plot them on the same plot

# 2d binning:
# 1. use 2d logbin algorithm to make a histogram out of the raw data, summing bin values


###fxns###

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
plotpowerlawfit <- function(pl, xl, title, subtitle) {
  library(ggplot2)
  th <- theme_bw() + theme(panel.grid = element_blank())
  
  expr1 <- as.character(as.expression(substitute(
    "Pareto:"~~alpha == a, list(a = round(pl$alpha, 2)))))
  expr2 <- as.character(as.expression(substitute(
    "Lognormal:"~~mu == m*","~~sigma == s, list(m = round(pl$pars_lognorm[1], 2), 
                                                s = round(pl$pars_lognorm[2], 2)))))
  
  p <- ggplot(pl$plotdat, aes(x = x, y = y)) +
    scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                  labels = scales::trans_format("log10", scales::math_format(10^.x))) +
    scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                  labels = scales::trans_format("log10", scales::math_format(10^.x))) +
    geom_point() +
    geom_line(data = subset(pl$plfit, y>0), color = 'forestgreen') +
    geom_line(data = subset(pl$lognormfit, y>0), color = 'indianred') +
    geom_text(x = -Inf, y = -Inf, label = expr1, parse = TRUE, hjust = 0, vjust = -1.5) +
    geom_text(x = -Inf, y = -Inf, label = expr2, parse = TRUE, hjust = 0, vjust = -0.25) +
    th +
    labs(x = xl, y = 'log CDF') +
    ggtitle(title, subtitle)
  return(p)
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

##########
### 2D ###
##########

rm(list=ls())

source('~/GitHub/forestlight/code/mergeshade.r')

pdat <- dat %>% 
  mutate(biomass3 = pmax(BiomassDGH_year3_kg, BiomassDBH_year3_kg, na.rm=TRUE),
         biomass2 = pmax(BiomassDGH_year2_kg, BiomassDBH_year2_kg, 0, na.rm=TRUE),
         biomass1 = pmax(BiomassDGH_year1_kg, BiomassDBH_year1_kg, 0, na.rm=TRUE),
         massprod23 = pmax(biomass3 - biomass2, 0, na.rm=TRUE),
         massprod12 = pmax(biomass2 - biomass1, 0, na.rm=TRUE),
         massprod13 = pmax((biomass3 - biomass1)/2, 0, na.rm=TRUE)) %>% 
  select(Site, biomass3, biomass2, massprod23, massprod12, massprod13, CII, tolerance)
bcidat <- pdat %>% filter(Site == 'Barro_Colorado_Island', !is.na(CII), !is.na(biomass3))
harvdat <- pdat %>% filter(Site == 'Harvard_Forest_LTER', !is.na(CII), !is.na(biomass3))

# Use mass production at BCI for all trees as the example for 2D binning.

x <- bcidat$biomass3
y <- bcidat$massprod13
n <- 20

bin1 <- logbin(x,y,n)
fit1 <- powerlawfit(bin1$bin_value)
plotbinsandfits(fit1, bin1, 'mass','production', 'bci','all')


# Make a pseudo-bin vector.
# Imagine that the bins are divided evenly. (this is sketchy but makes the math work)
pseudobins <- rep(bin1$bin_value, bin1$bin_count)

# Fit PDFs to the pseudobins.
fit1 <- powerlawfit(pseudobins)

# Plot the bins and the pdf fit of pseudobins on the same plot.

