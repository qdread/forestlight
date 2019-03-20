# Recreate fg's 

library(dplyr)

# Load Nadja's data
fgbci <- read.table('~/google_drive/ForestLight/data/Ruger/fgroups_dynamics_new.txt', stringsAsFactors = FALSE)

# Correct functional groups so that: 1 fast, 2 pioneer, 3 slow, 4 breeder, 5 intermediate
# Old 1,2,3,4,5 --> New 2,3,1,4,5
fgbci$fg5 <- match(fgbci$fg5, c(2,3,1,4,5))

# Find percentiles divided by 2
cutoffs_1 <- with(fgbci, quantile(X1new, c(0.1,0.9))) / 2
cutoffs_2 <- with(fgbci, quantile(X2new, c(0.1,0.9))) / 2

# Separate by 45 degree angle on each PCA axis
fg_recreated <- with(fgbci, case_when(
  X1new > X2new & X1new > -X2new ~ 3,
  X1new > X2new & X1new < -X2new ~ 4,
  X1new < X2new & X1new > -X2new ~ 2,
  X1new < X2new & X1new < -X2new ~ 1,
))

# Confusion matrix to see if this worked
table(fg_recreated, fgbci$fg5) # This more or less recreates functional groups 1-4 but functional group 5 still needs to be created. Some are still in error.

# Rescale such that the cutoffs are equivalent to -1 and 1 on each axis

# Divide all species scores by absolute value 
# I was not sure what she meant by this so I interpreted it as the average of the absolute values of the two cutoffs.
X1divided <- fgbci$X1new / mean(abs(cutoffs_1))
X2divided <- fgbci$X2new / mean(abs(cutoffs_2))

fg_recreated <- ifelse(sqrt(X1divided^2 + X2divided^2) <= 1, 5, fg_recreated) # Assign fg5 as anything within a circle of radius 1.

# Final confusion matrix
table(fg_recreated, fgbci$fg5) # This confusion matrix shows that I'm close to recreating her method but not quite.

plot(fgbci$X1new, fgbci$X2new, col = as.numeric(fgbci$fg5 == fg_recreated) + 2) # Red ones are the errors.
