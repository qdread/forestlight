# Create toy dataset.
# -------------------

gdrive_path <- 'C:/Users/Q/google_drive' # Change to google drive path.
load(file.path(gdrive_path, 'ForestLight/data/data_22jan2018/rawdataobj_22jan.r'))

# Toy dataset of 1000 random dbh-production pairs from 1995.
set.seed(121)
sample_dat <- alltreedat[[3]][sample(nrow(alltreedat[[3]]), 1000), c('dbh_corr','production')]
sample_dat <- sample_dat[order(sample_dat$dbh_corr), ]
row.names(sample_dat) <- NULL

sample_dat_small <- alltreedat[[3]][sample(nrow(alltreedat[[3]]), 100), c('dbh_corr','production')]
sample_dat_small <- sample_dat_small[order(sample_dat_small$dbh_corr), ]
row.names(sample_dat_small) <- NULL
names(sample_dat_small) <- c('size','production')

# Method 1 to get total production: Directly bin total production. This method gives correct results.
# ---------------------------------------------------------------------------------------------------

# Source "logbin" function
source(file.path(gdrive_path, 'ForestLight/logbin.r'))

# Bin density
density_bin <- logbin(x = sample_dat$dbh_corr, n = 10)
# Bin total production
total_production_bin <- logbin(x = sample_dat$dbh_corr, y = sample_dat$production, n = 10)

# Confirm that method 1 gives correct results.

sum(sample_dat$production)

# Back calculate by multiplying by bin width again.
# This is the same value as the actual total production.
sum(with(total_production_bin, bin_value * (bin_max - bin_min)))


# Method 2. Get median production for each density bin and multiply by the number of individuals in each bin.
# This method gives too low of a result.
# -----------------------------------------------------------------------------------------------------------

# Get median production for each density bin
production_median_fn <- function(dat, bins) {
  bin_medians <- numeric(nrow(bins))
  for (i in 1:nrow(bins)) {
    min_i <- bins$bin_min[i]
    max_i <- ifelse(i < nrow(bins), bins$bin_max[i], bins$bin_max[i] + 1) # Make sure the largest individual is in a bin
    bin_medians[i] <- median(dat$production[dat$dbh_corr >= min_i & dat$dbh_corr < max_i])
  }
  bin_medians
}

production_medians <- production_median_fn(sample_dat, density_bin)

total_production_product <- density_bin$bin_count * production_medians / with(total_production_bin, bin_max - bin_min)

sum(total_production_product * with(total_production_bin, bin_max - bin_min)) # The answer is 7723 which is much lower than the true value

data.frame(method_1 = total_production_bin$bin_value, method_2 = total_production_product) # Method 2 always lower

# (I also checked this with geometric mean and it gives ~ the same result as the median)

production_median_fn <- function(dat, bins) {
  bin_medians <- numeric(nrow(bins))
  for (i in 1:nrow(bins)) {
    min_i <- bins$bin_min[i]
    max_i <- ifelse(i < nrow(bins), bins$bin_max[i], bins$bin_max[i] + 1) # Make sure the largest individual is in a bin
    bin_medians[i] <- mean(dat$production[dat$dbh_corr >= min_i & dat$dbh_corr < max_i])
  }
  bin_medians
}