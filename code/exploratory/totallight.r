# Check total insolation of bci versus total amount of par intercepted by trees in bci.

# Total amount of light in W received by the 50-hectare plot at BCI

insolation <- function(lat) {
  lat <- lat * pi/180 # to radians
  y <- sin(lat)
  0.25 * 1367 * (1 - 0.482 * (3*y^2 - 1)/2)
}

# Insolation at BCI, 9.2 degrees N
(insol_bci <- insolation(9.2))

# That is in W m-2 so we multiply by the area of the 50-hectare plot, in m2.

area_bci <- 50 * 100^2
(light_bci <- insol_bci * area_bci) # On average, 200 million watts hit BCI (averaging over night and all seasons). This ignores clouds and such.

sum(bcicensusdat$light_received, na.rm = TRUE)
