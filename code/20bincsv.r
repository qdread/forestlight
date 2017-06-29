# Create CSV with binned data for analysis.

logbin_setedges <- function(x, y = NULL, edges) {
  logx <- log10(x)                                           # log transform x value (biomass)
  bin_edges <- log10(c(edges$bin_min, edges$bin_max[length(edges$bin_max)]))
  n <- length(edges$bin_min)
  logxbin <- rep(NA, length(logx))                           # create data structure to assign trees to bins
  b <- bin_edges                                             # add a little to the biggest bin temporarily
  b[length(b)] <- b[length(b)] + 1                           # (so that the biggest single tree is put in a bin)
  for (i in 1:length(logx)) {
    logxbin[i] <- sum(logx[i] >= b)                          # assign each tree to a bin
  }
  bin_midpoints <- edges$bin_midpoint
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


# Create density scaling for all trees.
n_bins <- 20

bin_lims <- c(floor(min(alltreedat$dbh)), ceiling(max(alltreedat$dbh)))
bin_edges <- seq(log10(bin_lims[1]), log10(bin_lims[2]), length.out=n_bins + 1)

bin_dims <- data.frame(bin_min = 10^bin_edges[1:n_bins], bin_max = 10^bin_edges[2:(n_bins + 1)])
bin_dims <- transform(bin_dims, bin_midpoint = (bin_min+bin_max)/2)

unclassifieddat <- subset(alltreedat, is.na(tol_wright))

shade_bin_dens <- with(shadedat, logbin_setedges(x=dbh, y=NULL, edges=bin_dims))
shade_bin_prod <- with(shadedat, logbin_setedges(x=dbh, y=production34, edges=bin_dims))
shade_bin_light <- with(shadedat, logbin_setedges(x=dbh, y=light_received, edges=bin_dims))

shade_bin_all <- cbind(guild = 'shade', shade_bin_dens[,c('bin_midpoint','bin_min','bin_max','bin_count')],
                       density_bin_value = shade_bin_dens$bin_value/50,
                       production_bin_value = shade_bin_prod$bin_value/50, 
                       light_bin_value = shade_bin_light$bin_value/50)

int_bin_dens <- with(intdat, logbin_setedges(x=dbh, y=NULL, edges=bin_dims))
int_bin_prod <- with(intdat, logbin_setedges(x=dbh, y=production34, edges=bin_dims))
int_bin_light <- with(intdat, logbin_setedges(x=dbh, y=light_received, edges=bin_dims))

int_bin_all <- cbind(guild = 'intermediate', int_bin_dens[,c('bin_midpoint','bin_min','bin_max','bin_count')],
                       density_bin_value = int_bin_dens$bin_value/50,
                       production_bin_value = int_bin_prod$bin_value/50, 
                       light_bin_value = int_bin_light$bin_value/50)

gap_bin_dens <- with(gapdat, logbin_setedges(x=dbh, y=NULL, edges=bin_dims))
gap_bin_prod <- with(gapdat, logbin_setedges(x=dbh, y=production34, edges=bin_dims))
gap_bin_light <- with(gapdat, logbin_setedges(x=dbh, y=light_received, edges=bin_dims))

gap_bin_all <- cbind(guild = 'gap', gap_bin_dens[,c('bin_midpoint','bin_min','bin_max','bin_count')],
                       density_bin_value = gap_bin_dens$bin_value/50,
                       production_bin_value = gap_bin_prod$bin_value/50, 
                       light_bin_value = gap_bin_light$bin_value/50)

unclassified_bin_dens <- with(unclassifieddat, logbin_setedges(x=dbh, y=NULL, edges=bin_dims))
unclassified_bin_prod <- with(unclassifieddat, logbin_setedges(x=dbh, y=production34, edges=bin_dims))
unclassified_bin_light <- with(unclassifieddat, logbin_setedges(x=dbh, y=light_received, edges=bin_dims))

unclassified_bin_all <- cbind(guild = 'unclassified', unclassified_bin_dens[,c('bin_midpoint','bin_min','bin_max','bin_count')],
                       density_bin_value = unclassified_bin_dens$bin_value/50,
                       production_bin_value = unclassified_bin_prod$bin_value/50, 
                       light_bin_value = unclassified_bin_light$bin_value/50)

bin_all <- rbind(shade_bin_all, int_bin_all, gap_bin_all, unclassified_bin_all)
write.csv(bin_all, 'data/bci_binned.csv', row.names = FALSE)
