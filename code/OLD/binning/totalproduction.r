# Correction factor for total production

prodtotals <- obs_totalprod %>%
  group_by(fg, year) %>%
  summarize(total = sum(bin_value * (bin_max-bin_min)))


dbh_pred <- exp(seq(log(1.2), log(315), length.out = 101))
mids <- function(a) a[-length(a)] + diff(a)/2
dbh_binwidth <- diff(dbh_pred)

dbh_df <- data.frame(dbh = mids(dbh_pred), binwidth = dbh_binwidth)

prodintegrals <- pred_totalprod %>%
  left_join(dbh_df) %>%
  mutate(int = q50 / binwidth) %>%
  group_by(year, dens_model, prod_model, fg) %>%
  summarize(int = sum(int, na.rm = TRUE))
  
library(pracma)  

# Trapezoidal integration is used.
prodintegrals <- pred_totalprod %>%
  group_by(year, dens_model, prod_model, fg) %>%
  summarize_at(vars(starts_with('q')), funs(trapz(x = dbh, y = .)))

plot_totalprod(year_to_plot = 1995,
               fg_names = c('all'),
               model_fit_density = 'weibull', 
               model_fit_production = 'powerlaw',
               y_limits = c(0.01, 500),
               y_breaks = c(0.1, 1, 10, 100),
               preddat = pred_totalprod %>% mutate_at(vars(starts_with('q')), funs(./0.7617)))


# Get the true numbers from the raw data.
area_core <- 42.84
library(purrr)
alltreeprod <- map_dbl(alltreedat, function(x) sum(x$production))
fgprod <- map_dfr(alltreedat, function(x) {
  x %>%
    mutate(fg = if_else(is.na(fg), 'unclassified', paste0('fg',fg))) %>%
    group_by(fg) %>%
    summarize(production = sum(production))
})

prodtotals <- rbind(data.frame(year = rep(seq(1985,2010,5), each=6), fgprod),
                    data.frame(year = seq(1985,2010,5), fg = 'alltree', production = alltreeprod)) %>%
  arrange(year, fg)

write.csv(prodtotals, file = 'C:/Users/Q/google_drive/ForestLight/data/production_total.csv', row.names = FALSE)
