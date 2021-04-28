rm(list = ls())
setwd('E:/Book/Ch_12')

library(data.table)
library(midasr)

### 12.5 Simulated Data ###
sim.data <- fread('simulated_data_ch12_sec5.csv')
sim.data[, c('date', 'year') := list(as.Date(date, format = '%m/%d/%Y'),
                                     rep(1950:2012, each = 4))]

# Create yearly values of x and y as the sum of quarterly values
sim.data[, c('X', 'Y') := list(sum(x), sum(y)), by = year]

### 12.5.1 Bridge ###

# Low frequency Bridge equation
bridge.lf <- lm(Y ~ X, data = unique(sim.data[year < 2012, .(year, X, Y)]))

# AR(2) model for x_H & flash estimates for y_L
bridge.flash <- rep(0 ,4)
for (n in 1:3) {
  ar.fit <- ar(sim.data[1:(248 + n), x], order.max = 2, method = 'ols')
  ar.predict <- as.numeric(predict(ar.fit, n.ahead = 4 - n)$pred)
  bridge.flash[n] <- bridge.lf$coefficients[1] + bridge.lf$coefficients[2] *
    sum(sim.data[249:(249 + n - 1), x], sum(ar.predict))
}
bridge.flash[4] <- bridge.lf$coefficients[1] +
  bridge.lf$coefficients[2] * sim.data[252, X]

### 12.5.2 MIDAS ###
y.lf <- ts(unique(sim.data[year < 2012, Y]))
trend <- 1:length(y.lf)

# MIDAS & U-MIDAS
for (n in 1:4) {
  x.hf <- ts(sim.data[1:248 + n - 1, x], frequency = 4)
  
  # MIDAS
  midas.fit <- midas_r(y.lf ~ trend + mls(y.lf, 1, 1) +
                         fmls(x.hf, 3, 4, nealmon),
                       start = list(x.hf = rep(0, 4)))
  print('MIDAS')
  print(predict(midas.fit, list(x.hf = sim.data[1:248 + n, x])))
  
  # U-MIDAS
  um.fit <- midas_r(y.lf ~ trend + mls(y.lf, 1, 1) +
                         fmls(x.hf, 3, 4), start = list(x.hf = rep(0, 4)))
  print('U-MIDAS')
  print(predict(um.fit, list(x.hf = sim.data[1:248 + n, x])))
}

### 12.6 Nowcasting US GDP Growth ###

gdp.hf <- fread('us_gdp_monthly_ch12_sec6.csv',
                select = c('date', 'ci', 'e', 'ip', 'inc', 'sa'))
gdp.hf[, date := as.Date(date, format = '%m/%d/%Y')]
gdp.lf <- fread('us_gdp_quarterly_ch12_sec6.csv')

gdp.hf[, qtr := rep(1:112, each = 3)]
gdp.hf[, c('CI', 'E', 'IP', 'INC', 'SA') :=
         list(sum(ci), sum(e), sdum(ip), sum(inc), sum(sa)), by = qtr]

# Case 1: End in 2008Q3
# Case 2: End in 2011Q3
# Case 3: End in 2012Q3
for (end.date in c('2008-09-01', '2011-09-01', '2012-09-01')) {
  gdp.growth <- unique(gdp.hf[date <= end.date, .(CI, E, IP, INC, SA)])
  
  n.lf <- dim(gdp.growth)[1]
  n.hf <- n.lf * 3
  
  for (var in c('CI', 'E', 'IP', 'INC', 'SA')) {
    bridge.fit <- lm(gdp.lf[1:n.lf, y] ~ gdp.growth[[var]])
    
    # Bridge equation coefficients
    print(var)
    print(bridge.fit)
  }
  
  # MIDAS & U-MIDAS
  y.lf <- gdp.lf[1:n.lf, y]
  trend <- 1:length(y.lf)
  for (var in c('CI', 'E', 'IP', 'INC', 'SA')) {
    for (n in 1:3) {
      x.hf <- ts(gdp.hf[1:n.hf + n - 1][[var]], frequency = 3)
      
      # MIDAS
      midas.fit <- midas_r(y.lf ~ trend + mls(y.lf, 1, 1) +
                             fmls(x.hf, 2, 3, nealmon),
                           start = list(x.hf = rep(0, 3)))
      print('MIDAS')
      print(predict(midas.fit, list(x.hf = gdp.hf[1:n.hf + n][[var]])))
      
      # U-MIDAS
      um.fit <- midas_r(y.lf ~ trend + mls(y.lf, 1, 1) +
                          fmls(x.hf, 2, 3), start = list(x.hf = rep(0, 3)))
      print('U-MIDAS')
      print(predict(um.fit, list(x.hf = gdp.hf[1:n.hf + n][[var]])))
    }
  }
}