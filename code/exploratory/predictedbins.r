# Function to do random draws and rebin to get a predicted value distribution for total energy

predicted_randomdraws <- function(fit, dbh_pred, dens_form, prod_form, x_min = NULL, n_indiv = 1, ll = 1.1, ul = 316, n_bin = 20, pars_to_get, delete_samples = NULL) {
  require(purrr)
  require(truncdist)
  require(extremevalues)
  powerlaw_exp_log <- function(x, a, b, c, beta0, beta1) exp(-beta0) * x^beta1 * (-a * x ^ -b + c)
  powerlaw_log <- function(x, beta0, beta1) exp(-beta0) * x^beta1
  
  pars <- as.data.frame(do.call('cbind', extract(fit, c(pars_to_get, 'sigma'))))
  
  if (!is.null(delete_samples)) {
    pars <- pars[-delete_samples,]
  }
  
  if (dens_form == 'pareto') {
    dbhs_alldraws <- pmap(pars, function(alpha, ...) rpareto(n_indiv, xm = x_min, alpha = alpha))
  }
  if (dens_form == 'weibull') {
    dbhs_alldraws <- pmap(pars, function(m, n, ...) rtrunc(n_indiv, 'weibull', a = ll, b = ul, shape = m, scale = n))
  }
  if (prod_form == 'powerlaw') {
    prods_alldraws <- lapply(1:length(dbhs_alldraws), function(i) powerlaw_log(dbhs_alldraws[[i]], beta0 = pars$beta0[i], beta1  = pars$beta1[i]))
  }
  if (prod_form == 'powerlawexp') {
    prods_alldraws <- lapply(1:length(dbhs_alldraws), function(i) powerlaw_exp_log(dbhs_alldraws[[i]], a = pars$a[i], b = pars$b[i], c = pars$c[i], beta0 = pars$beta0[i], beta1  = pars$beta1[i]))
  }
  
  fitted_bin_alldraws <- map2(dbhs_alldraws, prods_alldraws, ~ logbin_setedges(.x, .y, exp(seq(log(ll),log(ul),length.out=n_bin+1))))
  
  # Get some quantiles from each of these.
  qprobs <- c(0.025, 0.05, 0.25, 0.5, 0.75, 0.95, 0.975)
  
  binvalues <- do.call(rbind, map(fitted_bin_alldraws, 'bin_value'))
  binquantiles <- apply(binvalues, 2, quantile, probs = qprobs) %>% 
    t %>% 
    as.data.frame %>%
    setNames(c('q025', 'q05', 'q25', 'q50', 'q75', 'q95', 'q975'))
  cbind(fitted_bin_alldraws[[1]][,c('bin_midpoint','bin_min','bin_max')], binquantiles)

}

# Function to load a fit and get the predicted bins out of it.

get_predbins <- function(dens_model, prod_model, fg, year, xmin, n, total_production, use_subset = FALSE) {
  require(rstan)
  require(loo)
  
  # Load CSVs as stanfit object
  print('Loading stan fit . . .')
  fp <- '~/forestlight/old_stanoutput/jul2018'
  files <- paste0('fit_', dens_model, 'x', prod_model, '_', fg, '_', year, '_', 1:3, '.csv')
  if (use_subset) files <- paste0('ss', files) # Use the 25K subset if needed.
  fit <- read_stan_csv(file.path(fp,files))
  
  prod_model <- ifelse(prod_model == 'power', 'powerlaw', 'powerlawexp')
  
  pareto_par <- c('alpha')
  weibull_par <- c('m', 'n')
  powerlaw_par <- c('beta0', 'beta1')
  powerlawexp_par <- c('beta0', 'beta1', 'a', 'b', 'c')
  
  if (dens_model == 'pareto') get_pars <- pareto_par else get_pars <- weibull_par
  if (prod_model == 'powerlaw') get_pars <- c(get_pars, powerlaw_par) else get_pars <- c(get_pars, powerlawexp_par)
  
  
  # Get fitted values with credible intervals and prediction intervals
  print('Calculating fitted values and prediction intervals . . .')
  pred_interval <- predicted_randomdraws(fit, dbh_pred, dens_form = dens_model, prod_form = prod_model, x_min = xmin, n_indiv = n, pars_to_get = get_pars, n_bin = 20)
  pred_interval <- data.frame(year = year, dens_model = dens_model, prod_model = prod_model, fg = fg, pred_interval)

  return(pred_interval)
}

# Log bin function
logbin_setedges <- function(x, y = NULL, bin_edges) {
  logx <- log10(x)                                           # log transform x value (biomass)
  bin_edges <- log10(bin_edges)
  n <- length(bin_edges) - 1
  logxbin <- rep(NA, length(logx))                           # create data structure to assign trees to bins
  b <- bin_edges                                             # add a little to the biggest bin temporarily
  b[length(b)] <- b[length(b)] + 1                           # (so that the biggest single tree is put in a bin)
  for (i in 1:length(logx)) {
    logxbin[i] <- sum(logx[i] >= b)                          # assign each tree to a bin
  }
  bin_midpoints <- 10^(bin_edges[-1] - diff(bin_edges)/2)
  bin_widths <- diff(10^bin_edges)                              # get linear width of each bin
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
