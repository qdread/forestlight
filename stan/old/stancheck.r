
z <- expand.grid(year = c(1990, 1995, 2000, 2005, 2010),
				 fg = c('fg1', 'fg2', 'fg3', 'fg4', 'fg5', 'unclassified', 'alltree'),
				 model_name = c('ppow', 'wpow', 'pexp', 'wexp', 'pbert', 'wbert'),
				 stringsAsFactors = FALSE)

i <- 1:nrow(z)
todo <- which(!file.exists(file.path('~/forestlight/stanoutput', paste0('fit_', z$model_name[i], '_', z$fg[i], '_', z$year[i], '.r'))))
paste(todo, collapse = ',')