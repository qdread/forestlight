# Compile all the data needed to recreate the plots into some easily interpretable data frames
# 05 Sep 2017

# Raw data

# All trees, all years
alltreedat # List of 6 data frames, 1985-2010 inclusive.
# Trees with light availability estimated by Nadja in 1990 and 1995
alltree_light_90
alltree_light_95

# Binned data

alltreedat_dens_logbin # and all other names from allyears_names
alltree_light_90_dens_logbin # and all other names from names1990, names1995

alltreedat_prod_logbin #and all other names
alltree_light_90_prod_logbin #and all other names

alltree_light_90_light_logbin #and all other names from 1990 and 1995
alltree_light_90_lightdens_logbin #ditto

# Coefficients with confidence intervals

allyear_dens_coeff
allyear_deltaaic
dens_light_coeff
deltaaiclight
prod_slopes
prod_light_slopes
light_diam_slopes
prod_bin_slopes
light_bin_slopes


# Save everything into an R object for exporting.
names_list <- c('alltreedat', 'alltree_light_90', 'alltree_light_95', 
                paste0(allyears_names, '_dens_logbin'), 
                paste0(allyears_names, '_prod_logbin'),
                paste0(c(names1990, names1995), '_light_logbin'),
                paste0(c(names1990, names1995), '_lightdens_logbin'),
                'allyear_dens_coeff','allyear_deltaaic','dens_light_coeff','deltaaiclight',
                'prod_slopes','prod_light_slopes','light_diam_slopes','prod_bin_slopes','light_bin_slopes')
save(list = names_list, file = 'C:/Users/Q/google_drive/ForestLight/data/data_06sep/bci_data_object.RData')
