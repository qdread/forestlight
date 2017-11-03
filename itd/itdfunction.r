# Function to run a single timestep of the ITD model

ITD <- function(P, M, A) {
  # 1. calculate critical canopy height z*
  
  # add the crown areas of each cohort starting with the tallest
  M <- M[order(M$z, decreasing = TRUE), ]
  totalcrownA <- cumsum(M$crownA)
  # find the height of the last cohort added before plot area was reached
  # if plot area is not reached, all cohorts are in the canopy
  zstar <- if (sum(M$crownA) >= A) M$z[which(totalcrownA >= A)[1]]
  else min(M$z)
  # create logical for canopy status
  M$overstory <- M$z >= zstar
  
  # 2. implement mortality
  
  # assign canopy or understory mortality parameter to each cohort in each species
  mu <- apply(cbind(M$overstory, P$muL[M$species], P$muD[M$species]), 1, function(x) if (x[1]) x[2] else x[3])
  # get new abundance for each cohort
  M$w <- M$w * (1 - mu)
  
  # 3. implement growth
  
  # assign canopy or understory mortality parameter to each cohort in each species
  G <- apply(cbind(M$overstory, P$GL[M$species], P$GD[M$species]), 1, function(x) if (x[1]) x[2] else x[3])
  # get new diameter for each cohort
  M$dbh <- M$dbh + G
  # get new height for each cohort with height-dbh allometry
  M$z <- P$alpha[M$species] * M$dbh ^ (P$beta[M$species])
  # get new crown area for each cohort with crown area-dbh allometry
  M$crownA <- M$w * pi * (1/10000) * (P$phi[M$species] * M$dbh) ^ 2
  # sort the list of cohorts by height
  M <- M[order(M$z, decreasing = TRUE), ]
  # add the crown areas of each cohort starting with the tallest
  totalcrownA <- cumsum(M$crownA)
  # find the height of the last cohort added before plot area was reached
  # if plot area is not reached, all cohorts are in the canopy
  zstar <- ifelse(sum(M$crownA) >= A, M$z[which(totalcrownA >= A)[1]], min(M$z))
  # create logical for canopy status
  M$overstory <- M$z >= zstar
  
  # 4. implement reproduction
  
  # calculate proportion exposed crown area for each species
  ACanopy <- with(subset(M, overstory == TRUE), sum(crownA))
  AjCanopy <- rep(0, nrow(P))
  for (j in 1:nrow(P)) {
    Aj <- with(subset(M, species == j), sum(crownA))
    AjCanopy[j] <- Aj - with(subset(M, species == j & overstory == FALSE), sum(crownA))
  }
  AjCanopy[is.na(AjCanopy)] <- 0
  # multiply proportion by fecundity to get number of individuals in next cohort
  wnew <- P$F * (AjCanopy/ACanopy)
  seedlings <- data.frame(species = 1:nrow(P),
                          cohort = max(M$cohort) + 1,
                          dbh = 0,
                          w = wnew,
                          z = 0,
                          crownA = 0,
                          overstory = FALSE)
  # add new cohorts to output matrix
  M <- rbind(M, seedlings)
  return(M)
}