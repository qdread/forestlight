# Overlaid histograms for relative light, as an alternative to the 3x3 plots in the Rmarkdown document.
# QDR 28 June 2017

numbins <- 10 # only for visualizing


shade_lowlight_bin <-  with(subset(shadedat, light_group == 'Low light'), logbin(x=dbh, y=NULL, n = numbins))
shade_intlight_bin <-  with(subset(shadedat, light_group == 'Intermediate light'), logbin(x=dbh, y=NULL, n = numbins))
shade_highlight_bin <-  with(subset(shadedat, light_group == 'High light'), logbin(x=dbh, y=NULL, n = numbins))

int_lowlight_bin <-  with(subset(intdat, light_group == 'Low light'), logbin(x=dbh, y=NULL, n = numbins))
int_intlight_bin <-  with(subset(intdat, light_group == 'Intermediate light'), logbin(x=dbh, y=NULL, n = numbins))
int_highlight_bin <-  with(subset(intdat, light_group == 'High light'), logbin(x=dbh, y=NULL, n = numbins))

gap_lowlight_bin <-  with(subset(gapdat, light_group == 'Low light'), logbin(x=dbh, y=NULL, n = numbins))
gap_intlight_bin <-  with(subset(gapdat, light_group == 'Intermediate light'), logbin(x=dbh, y=NULL, n = numbins))
gap_highlight_bin <-  with(subset(gapdat, light_group == 'High light'), logbin(x=dbh, y=NULL, n = numbins))

# Concatenate each light level into one df.
lowlight_bins <- rbind(transform(shade_lowlight_bin, guild = 'shade'),
                       transform(int_lowlight_bin, guild = 'intermediate'),
                       transform(gap_lowlight_bin, guild = 'gap'))
intlight_bins <- rbind(transform(shade_intlight_bin, guild = 'shade'),
                       transform(int_intlight_bin, guild = 'intermediate'),
                       transform(gap_intlight_bin, guild = 'gap'))
highlight_bins <- rbind(transform(shade_highlight_bin, guild = 'shade'),
                       transform(int_highlight_bin, guild = 'intermediate'),
                       transform(gap_highlight_bin, guild = 'gap'))

raw_numbers <- c(lowlight_bins$bin_value, intlight_bins$bin_value, highlight_bins$bin_value)/50
y_min <- 10^floor(log10(min(raw_numbers)))
y_max <- max(raw_numbers) * 1.1
y_values <- -4:3

xl2 <- 'Diameter (cm)'
yl2 <- expression(paste('Density (individuals ha'^-1,')', sep=''))

th1 <- theme(plot.title = element_text(size = 12), plot.subtitle = element_text(size = 10))

plotlogbin_overlaid <- function(dat, xl, yl, plottitle, plotsubtitle=NULL, reg = FALSE, cutoff = NA, y_min=0.1, y_max=53, x_min=1.1, x_max=141, plotarea=50, y_values=-1:3, x_values = 0:2) {
  
  dat <- transform(dat, bin_value = bin_value/plotarea) # kg per hectare.
  
  p <- ggplot(dat, aes(xmin=bin_min, xmax=bin_max, ymin=0, ymax=bin_value)) + 
    geom_rect(alpha = 1, aes(group=guild, fill=guild), color = 'black') +
    scale_x_log10(name = xl, expand = c(0,0),
                  breaks = 10^x_values,
                  labels = as.character(10^x_values), limits=c(x_min,x_max)) +
    scale_y_log10(name = yl, expand = c(0,0), limits = c(y_min, y_max),
                  breaks = 10^y_values,
                  labels = as.character(10^y_values)) +
    panel_border(colour = 'black') + 
    ggtitle(plottitle, plotsubtitle)
  if (reg) {
    p <- p +
      stat_smooth(method = 'lm', se = FALSE, color = 'forestgreen', size = 2,
                  aes(x = bin_midpoint, y = bin_value)) +
      geom_text(x = -Inf, y = -Inf, 
                label = paste('Slope without cutoff:', 
                              round(lm(I(log10(bin_value)) ~ I(log10(bin_midpoint)), data=dat)$coef[2], 2)),
                hjust = 0, vjust = -1.5)
    
    if (!is.na(cutoff)) {
      p <- p +
        stat_smooth(method = 'lm', se = FALSE, color = 'goldenrod', size = 2,
                    aes(x = bin_midpoint, y = bin_value), data = subset(dat, bin_midpoint <= cutoff)) +
        geom_text(x = -Inf, y = -Inf, 
                  label = paste('Slope with cutoff:', 
                                round(lm(I(log10(bin_value)) ~ I(log10(bin_midpoint)), data=subset(dat, bin_midpoint <= cutoff))$coef[2], 2)),
                  hjust = 0, vjust = -0.2)
    }
  }
  return(p)
}

th1 <- theme(legend.position = c(0.8, 0.8))

overlay_low <- plotlogbin_overlaid(lowlight_bins, xl2, yl2,  
                                   'Low light', NULL, reg = FALSE, cutoff = NA, y_min = y_min, y_max = y_max, x_max = 141, y_values = y_values) + 
  th1
overlay_mid <- plotlogbin_overlaid(intlight_bins, xl2, yl2,  
                                   'Intermediate light', NULL, reg = FALSE, cutoff = NA, y_min = y_min, y_max = y_max, x_max = 141, y_values = y_values) + 
  th1
overlay_high <- plotlogbin_overlaid(highlight_bins, xl2, yl2,  
                                    'High light', NULL, reg = FALSE, cutoff = NA, y_min = y_min, y_max = y_max, x_max = 141, y_values = y_values) +
  th1

pg <- plot_grid(overlay_low, overlay_mid, overlay_high, nrow = 3)
ggsave('C:/Users/Q/google_drive/ForestLight/figs/bci50hectare_overlaidhistograms.png', pg, height = 9, width = 4, dpi = 400)


####################

# Redo log bin with same width.

logbin_setwidth <- function(x, y = NULL, width) {
  logx <- log10(x)                                           # log transform x value (biomass)
  bin_edges <- seq(min(logx), max(logx) + width, by = width) # get edges of bins
  n <- length(bin_edges) - 1
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

numbins <- 10
shade_lowlight_bin <-  with(subset(shadedat, light_group == 'Low light'), logbin(x=dbh, y=NULL, n = numbins))
low_binwidth <- with(shade_lowlight_bin, log10(bin_max)-log10(bin_min))[1]
int_lowlight_binset <- with(subset(intdat, light_group == 'Low light'), logbin_setwidth(x=dbh, y=NULL, width = low_binwidth))
gap_lowlight_binset <- with(subset(gapdat, light_group == 'Low light'), logbin_setwidth(x=dbh, y=NULL, width = low_binwidth))

shade_intlight_bin <-  with(subset(shadedat, light_group == 'Intermediate light'), logbin(x=dbh, y=NULL, n = numbins))
int_binwidth <- with(shade_intlight_bin, log10(bin_max)-log10(bin_min))[1]
int_intlight_binset <- with(subset(intdat, light_group == 'Intermediate light'), logbin_setwidth(x=dbh, y=NULL, width = int_binwidth))
gap_intlight_binset <- with(subset(gapdat, light_group == 'Intermediate light'), logbin_setwidth(x=dbh, y=NULL, width = int_binwidth))

shade_highlight_bin <-  with(subset(shadedat, light_group == 'High light'), logbin(x=dbh, y=NULL, n = numbins))
high_binwidth <- with(shade_highlight_bin, log10(bin_max)-log10(bin_min))[1]
int_highlight_binset <- with(subset(intdat, light_group == 'High light'), logbin_setwidth(x=dbh, y=NULL, width = high_binwidth))
gap_highlight_binset <- with(subset(gapdat, light_group == 'High light'), logbin_setwidth(x=dbh, y=NULL, width = high_binwidth))


# Concatenate each light level into one df.
lowlight_bins <- rbind(transform(shade_lowlight_bin, guild = 'shade'),
                       transform(int_lowlight_binset, guild = 'intermediate'),
                       transform(gap_lowlight_binset, guild = 'gap'))
intlight_bins <- rbind(transform(shade_intlight_bin, guild = 'shade'),
                       transform(int_intlight_binset, guild = 'intermediate'),
                       transform(gap_intlight_binset, guild = 'gap'))
highlight_bins <- rbind(transform(shade_highlight_bin, guild = 'shade'),
                        transform(int_highlight_binset, guild = 'intermediate'),
                        transform(gap_highlight_binset, guild = 'gap'))

raw_numbers <- c(lowlight_bins$bin_value, intlight_bins$bin_value, highlight_bins$bin_value)/50
y_min <- 10^floor(log10(min(raw_numbers)))
y_max <- max(raw_numbers) * 1.1
y_values <- -4:3

xl2 <- 'Diameter (cm)'
yl2 <- expression(paste('Density (individuals ha'^-1,')', sep=''))

th1 <- theme(plot.title = element_text(size = 12), plot.subtitle = element_text(size = 10), legend.position = c(0.8, 0.8))


overlay_low <- plotlogbin_overlaid(lowlight_bins, xl2, yl2,  
                                   'Low light', NULL, reg = FALSE, cutoff = NA, y_min = y_min, y_max = y_max, x_max = 141, y_values = y_values) + 
  th1
overlay_mid <- plotlogbin_overlaid(intlight_bins, xl2, yl2,  
                                   'Intermediate light', NULL, reg = FALSE, cutoff = NA, y_min = y_min, y_max = y_max, x_max = 141, y_values = y_values) + 
  th1
overlay_high <- plotlogbin_overlaid(highlight_bins, xl2, yl2,  
                                    'High light', NULL, reg = FALSE, cutoff = NA, y_min = y_min, y_max = y_max, x_max = 141, y_values = y_values) +
  th1

pg <- plot_grid(overlay_low, overlay_mid, overlay_high, nrow = 3)
ggsave('C:/Users/Q/google_drive/ForestLight/figs/bci50hectare_overlaidhistograms.png', pg, height = 9, width = 4, dpi = 400)

