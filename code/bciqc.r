# Quality Control run on raw BCI data
# Follows Meakem et al. 2017
# QDR 11 July 2017

# Load Condit's BCI data and Nadja's light data.

fp <- 'C:/Users/Q/Dropbox/projects/forestlight/'

growth8590 <- read.delim(file.path(fp, 'BCI_light/growth_final8590.txt'), stringsAsFactors = FALSE)
growth9095 <- read.delim(file.path(fp, 'BCI_light/growth_final9095.txt'), stringsAsFactors = FALSE)

# Load census data. 1 = 1980, 2 = 1985, 3 = 1990, 4 = 1995, 5 = 2000, 6 = 2005, 7 = 2010.

for (i in 1:7) {
  load(file.path(fp, paste0('bcidata/bci.full', i, '.rdata')))
  load(file.path(fp, paste0('bcidata/bci.stem', i, '.rdata')))
}

# Get rid of tree ferns and strangler figs
# Tree ferns are already removed
species_to_remove <- c('ficubu','ficuc1','ficuc2','ficuci','ficupe','ficupo')
remove_spp <- function(dat) subset(dat, !sp %in% c(species_to_remove, toupper(species_to_remove)))

# Apply taper correction for trees not measured at dbh because of their buttress, returning an estimate of what the true dbh would be.

# Apply algorithm that tracks stems between censuses

# Remove the patch of secondary forest in the BCI plot

# Calculate the area that Nadja has light values for, so we can determine the right denominator for the per-hectare values.