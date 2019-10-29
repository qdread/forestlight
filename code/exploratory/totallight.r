# Check total insolation of bci versus total amount of par intercepted by trees in bci.
# Updated 24 Oct 2019 for new data.

# Total amount of light in W received by the 50-hectare plot at BCI

insolation <- function(lat) {
  lat <- lat * pi/180 # to radians
  y <- sin(lat)
  0.25 * 1367 * (1 - 0.482 * (3*y^2 - 1)/2)
}

# Insolation at BCI, 9.2 degrees N
(insol_bci <- insolation(9.2))

# That is in W m-2 so we multiply by the area of the 50-hectare plot, in m2 (excluding young and edge)

area_bci <- 42.84 * 100^2
(light_bci <- insol_bci * area_bci) # On average, 180 MW of light hit our study area of the BCI plot (averaging over night and all seasons). This ignores clouds and such.


# Load tree data to find total light received
load('~/google_drive/ForestLight/data/rawdataobj_alternativecluster.r')

# Total light received in 1995
(capture_bci <- sum(alltree_light_95$light_received)) # ~ 210 MW of light

# Assume 70% efficiency
capture_bci * 0.7 # ~ 150 MW of light
