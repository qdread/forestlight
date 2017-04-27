# Calculation of light availability from tree location and neighbors

# Function to get tree height and crown dimensions from dbh
# Use same parameters for all trees, taken from Bohlman and O'Brien

tp <- function(dbh) {
  h <- exp(.438 + .595 * log(dbh))    # Height
  cd <- exp(-.157 + .702 * log(dbh))  # Crown depth
  cr <- exp(-.438 + .658 * log(dbh))  # Crown radius
  cV <- exp(-.681 + 2.02 * log(dbh))  # Crown volume
  c(h=h, cd=cd, cr=cr, cV=cV)
}

# Function to find area of "screen" above height h on a tree with given dbh.

Siplus <- function(dbh, h) {
  tps <- tp(dbh)
  
}