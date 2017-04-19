# Sample r directly from a power law.
library(actuar)

set.seed(8720604)

x <- rpareto1(1e6, 1, 1) # one million trees, yo. 
r <- sort(x)

b0 <- 0.1
k <- 1

b <- b0 * r^2
m <- k * r^(8/3)

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

densbin_r <- logbin(r, NULL, 50)
densbin_m <- logbin(m, NULL, 50)
prodbin_r <- logbin(r, b, 50)
prodbin_m <- logbin(m, b, 50)

numbins <- c(10,50,100,1000)

densbin_m <- list()
densbin_r <- list()
for (i in 1:length(numbins)) {
  densbin_r[[i]] <- cbind(nbins = rep(numbins[i],numbins[i]), logbin(r, NULL, numbins[i]))
  densbin_m[[i]] <- cbind(nbins = rep(numbins[i],numbins[i]), logbin(m, NULL, numbins[i]))
}

densbin_r <- do.call('rbind', densbin_r)
densbin_m <- do.call('rbind', densbin_m)

library(ggplot2)

ggplot(densbin_r, aes(x=bin_midpoint, y=bin_value, group=factor(nbins), color=factor(nbins))) +
  geom_line() + theme_minimal() + scale_y_log10() + scale_x_log10(limits=c(1, 1000))

ggplot(densbin_m, aes(x=bin_midpoint, y=bin_value, group=factor(nbins), color=factor(nbins))) +
  geom_line() + theme_minimal() + scale_y_log10() + scale_x_log10(limits=c(1,1e8))
