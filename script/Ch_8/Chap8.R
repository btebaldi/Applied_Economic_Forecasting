rm(list = ls())
setwd('E:/Book/Ch_8') # replace with own directory

library(BMR)
library(data.table)
library(ggplot2)
library(Metrics)
library(vars)

### 8.4 Simulated Data ###
sim.data <- fread('bvar_simulated_ch8_sec4.csv')
sim.data <- sim.data[101:600, .(x, y)]

sim.est <- sim.fore <- list()

# VAR(1)
sim.est$var1 <- VAR(sim.data[1:400], p = 1)
sim.fore$var1 <- predict(sim.est$var1, n.ahead = 100)
sim.fore$var1 <- data.table('x' = sim.fore$var1$fcst$x[, 1],
                            'y' = sim.fore$var1$fcst$y[, 1])

# VAR(4)
sim.est$var4 <- VAR(sim.data[1:400], p = 4)
sim.fore$var4 <- predict(sim.est$var4, n.ahead = 100)
sim.fore$var4 <- data.table('x' = sim.fore$var4$fcst$x[, 1],
                            'y' = sim.fore$var4$fcst$y[, 1])

# BVAR(1) - Minnesota prior
sim.est$bvar1 <- BVARM(sim.data[1:400], c(0, 0), p = 1,
                       HP1 = 0.2, HP2 = 0.99, HP3 = 1)
sim.fore$bvar1 <- forecast(sim.est$bvar1, periods = 100)
sim.fore$bvar1 <- data.table(sim.fore$bvar1$MeanForecast)
setnames(sim.fore$bvar1, c('x', 'y'))

# BVAR(4)
sim.est$bvar4 <- BVARM(sim.data[1:400], c(0, 0), p = 4,
                       HP1 = 0.2, HP2 = 0.99, HP3 = 1)
sim.fore$bvar4 <- forecast(sim.est$bvar4, periods = 100)
sim.fore$bvar4 <- data.table(sim.fore$bvar4$MeanForecast)
setnames(sim.fore$bvar4, c('x', 'y'))

# Evaluation
sim.error <- list('x' = data.table('Model' = c('var1', 'var4',
                                               'bvar1', 'bvar4')),
                  'y' = data.table('Model' = c('var1', 'var4',
                                               'bvar1', 'bvar4')))
for (model in c('var1', 'var4', 'bvar1', 'bvar4')) {
  sim.error$x[Model == model,
              c('RMSE', 'MAE') :=list(rmse(sim.data[401:500, x],
                                           sim.fore[[model]][, x]),
                                       mae(sim.data[401:500, x],
                                           sim.fore[[model]][, x]))]
  sim.error$y[Model == model,
              c('RMSE', 'MAE') := list(rmse(sim.data[401:500, y],
                                            sim.fore[[model]][, y]),
                                       mae(sim.data[401:500, y],
                                           sim.fore[[model]][, y]))]
}

### 8.5 Euro Area GDP Growth ###
gdp.data <- fread('dataonly_eugrowth_ch8_sec5.csv')

gdp.est <- gdp.fore <- list()

# Forecast 2003 - 2006
# BVAR(1)
gdp.est$bvar1 <- BVARM(gdp.data[1:68, .(g_fr, g_ger, g_it)],
                       c(0, 0, 0), p = 1, HP1 = 0.2, HP2 = 0.99, HP3 = 1)
gdp.fore$bvar1 <- forecast(gdp.est$bvar1, periods = 16)
gdp.fore$bvar1 <- data.table(gdp.fore$bvar1$MeanForecast)
setnames(gdp.fore$bvar1, c('g_fr', 'g_ger', 'g_it'))

# BVAR(2)
gdp.est$bvar2 <- BVARM(gdp.data[1:68, .(g_fr, g_ger, g_it)],
                       c(0, 0, 0), p = 2, HP1 = 0.2, HP2 = 0.99, HP3 = 1)
gdp.fore$bvar2 <- forecast(gdp.est$bvar2, periods = 16)
gdp.fore$bvar2 <- data.table(gdp.fore$bvar2$MeanForecast)
setnames(gdp.fore$bvar2, c('g_fr', 'g_ger', 'g_it'))

# Evaluation
gdp.error <- list('g_fr' = data.table('Model' = c('bvar1', 'bvar2')),
                  'g_ger' = data.table('Model' = c('bvar1', 'bvar2')),
                  'g_it' = data.table('Model' = c('bvar1', 'bvar2')))
for (model in c('bvar1', 'bvar2')) {
  for (var in c('g_fr', 'g_ger', 'g_it')) {
    gdp.error[[var]][Model == model,
                     c('RMSE', 'MAE') := list(rmse(gdp.data[[var]][69:84],
                                                   gdp.fore[[model]][[var]]),
                                              mae(gdp.data[[var]][69:84],
                                                  gdp.fore[[model]][[var]]))]
  }
}

# Forecast 2007 - 2012
# BVAR(1)
gdp.est$bvar1 <- BVARM(gdp.data[1:84, .(g_fr, g_ger, g_it)],
                       c(0, 0, 0), p = 1, HP1 = 0.2, HP2 = 0.99, HP3 = 1)
gdp.fore$bvar1 <- forecast(gdp.est$bvar1, periods = 24)
gdp.fore$bvar1 <- data.table(gdp.fore$bvar1$MeanForecast)
setnames(gdp.fore$bvar1, c('g_fr', 'g_ger', 'g_it'))

# BVAR(2)
gdp.est$bvar2 <- BVARM(gdp.data[1:84, .(g_fr, g_ger, g_it)],
                       c(0, 0, 0), p = 2, HP1 = 0.2, HP2 = 0.99, HP3 = 1)
gdp.fore$bvar2 <- forecast(gdp.est$bvar2, periods = 24)
gdp.fore$bvar2 <- data.table(gdp.fore$bvar2$MeanForecast)
setnames(gdp.fore$bvar2, c('g_fr', 'g_ger', 'g_it'))

# Evaluation
gdp.error <- list('g_fr' = data.table('Model' = c('bvar1', 'bvar2')),
                  'g_ger' = data.table('Model' = c('bvar1', 'bvar2')),
                  'g_it' = data.table('Model' = c('bvar1', 'bvar2')))
for (model in c('bvar1', 'bvar2')) {
  for (var in c('g_fr', 'g_ger', 'g_it')) {
    gdp.error[[var]][Model == model,
                     c('RMSE', 'MAE') := list(rmse(gdp.data[[var]][85:108],
                                                   gdp.fore[[model]][[var]]),
                                              mae(gdp.data[[var]][85:108],
                                                  gdp.fore[[model]][[var]]))]
  }
}

### 8.6 Multi-country Inflation Rates ###
inf.data <- fread('dataonly_inflation_eurostat_ch8_sec6.csv')
inf.data[, c('dateid01', 'dates') := NULL]

inf.est <- inf.fore <- list()

# VAR(4)
inf.est$var4 <- VAR(inf.data[1:115], p = 4)
inf.fore$var4 <- predict(inf.est$var4, n.ahead = 76)
inf.fore$var4 <- data.table('dp_fra' = inf.fore$var4$fcst$dp_fra[, 1],
                            'dp_ger' = inf.fore$var4$fcst$dp_ger[, 1],
                            'dp_ita' = inf.fore$var4$fcst$dp_ita[, 1],
                            'dp_spa' = inf.fore$var4$fcst$dp_spa[, 1],
                            'dp_uk' = inf.fore$var4$fcst$dp_uk[, 1],
                            'dp_usa' = inf.fore$var4$fcst$dp_usa[, 1])

# BVAR(4)
inf.est$bvar4 <- BVARM(inf.data[1:115], c(0, 0, 0, 0, 0, 0), p = 4,
                       HP1 = 0.2, HP2 = 0.99, HP3 = 1)
inf.fore$bvar4 <- forecast(inf.est$bvar4, periods = 76)
inf.fore$bvar4 <- data.table(inf.fore$bvar4$MeanForecast)
setnames(inf.fore$bvar4, names(inf.data))

# Evaluation
inf.error <- list('dp_fra' = data.table('Model' = c('var4', 'bvar4')),
                  'dp_ger' = data.table('Model' = c('var4', 'bvar4')),
                  'dp_ita' = data.table('Model' = c('var4', 'bvar4')),
                  'dp_spa' = data.table('Model' = c('var4', 'bvar4')),
                  'dp_uk' = data.table('Model' = c('var4', 'bvar4')),
                  'dp_usa' = data.table('Model' = c('var4', 'bvar4')))

for (model in c('var4', 'bvar4')) {
  for (var in c('dp_fra', 'dp_ger', 'dp_ita', 'dp_spa', 'dp_uk', 'dp_usa')) {
    inf.error[[var]][Model == model,
                     c('RMSE', 'MAE') := list(rmse(inf.data[[var]][116:191],
                                                   inf.fore[[model]][[var]]),
                                              mae(inf.data[[var]][116:191],
                                                  inf.fore[[model]][[var]]))]
  }
}