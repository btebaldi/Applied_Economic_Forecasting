rm(list = ls())
setwd('E:/Book/Ch_13')

# TPRF package source
install_github('IshmaelBelghazi/ThreePass')

library(data.table)
library(devtools)
library(forecast)
library(ggplot2)
library(Metrics)
library(psych)
library(vars)
library(T3PRF)

### 13.5 Simulated Data ###
# Quarterly target + Monthly factors
sim.dgp <- list('qtr' = fread('simulated_qtr_ch13_sec5.csv'),
                'mon' = fread('simulated_mon_ch13_sec5.csv'))

sim.dgp$y <- data.table('Date' = as.Date(sim.dgp$qtr[, date], '%Y-%m-%d'),
                        'yt' = sim.dgp$qtr[, yt],
                        'yt.new' = sim.dgp$qtr[2:160, yt])
sim.dgp$qtr[, c('date', 'yt') := NULL]
sim.dgp$mon <- fread('simulated_mon_ch13_sec5.csv')
sim.dgp$mon[, date := NULL]

## Benchmark 1 - SW
sim.dgp$pc <- principal(sim.dgp$qtr, nfactors = 3, rotate = 'varimax')$scores
sim.dgp$y[, c('PC1', 'PC2', 'PC3') := list(sim.dgp$pc[, 1], sim.dgp$pc[, 2],
                                           sim.dgp$pc[, 3])]
sim.dgp$sw <- lm(yt.new ~ PC1 + PC2 + PC3, data = sim.dgp$y[1:99])
sim.dgp$pred$sw <- predict(sim.dgp$sw, newdata = sim.dgp$y[100:159])

## Benchmark 2 - AR(1)
for (t in 1:60) {
  sim.dgp$ar <- arima(sim.dgp$y[1:(99 + t), yt], order = c(1, 0, 0))
  sim.dgp$pred$ar[t] <- as.numeric(predict(sim.dgp$ar, n.ahead = 1)$pred)
}

## Mixed-frequency TPRF
for (t in 1:60) {
  sim.dgp$ymf <- data.table('yt' = sim.dgp$y[1:(99 + t), yt.new])
  sim.dgp$ymf[, c('FS1', 'FS2', 'FS3') :=
                list(TPRF(sim.dgp$mon[seq(1, (99 + t) * 3 - 2, 3)],
                          sim.dgp$ymf, Z = sim.dgp$ymf)$factors[, 2],
                     TPRF(sim.dgp$mon[seq(2, (99 + t) * 3 - 1, 3)],
                          sim.dgp$ymf, Z = sim.dgp$ymf)$factors[, 2],
                     TPRF(sim.dgp$mon[seq(3, (99 + t) * 3, 3)],
                          sim.dgp$ymf, Z = sim.dgp$ymf)$factors[, 2])]
  sim.dgp$tprf <- lm(yt ~ FS1 + FS2 + FS3, data = sim.dgp$ymf)
  sim.dgp$pred$tprf[t] <- sim.dgp$tprf$fitted.values[99 + t]
}

## Evaluation 
model <- c('tprf', 'sw', 'ar')

# RMSE
for (m in model) {
  print(rmse(sim.dgp$y[101:160, yt], sim.dgp$pred[[m]]))
}

sim.dgp$yhat <- data.table('Y' = c(sim.dgp$y[101:160, yt], sim.dgp$pred$tprf,
                                   sim.dgp$pred$sw, sim.dgp$pred$ar),
                           'tag' = rep(c('YT', 'FORECAST3PRFS',
                                         'FORECASTPCAS', 'FORECASTARS'),
                                       each = 60),
                           'Date' = sim.dgp$y[101:160, Date])
ggplot(data = sim.dgp$yhat, aes(x = Date, y = Y, group = tag, colour = tag)) +
  theme_classic() + theme(plot.title = element_text(hjust = 0.5)) + 
  geom_line() + scale_x_date() + xlab('') + ylab('')

# Diebold-Mariano Test
for (m in 1:2) {
  for (n in (m + 1):3) {
    print(dm.test(sim.dgp$pred[[model[m]]] - sim.dgp$y[101:160, yt],
                  sim.dgp$pred[[model[n]]] - sim.dgp$y[101:160, yt]))
  }
}

### 13.6 Forecasting GDP Growth ###
g.gdp <- list()
g.gdp$us <- g.gdp$ea <- g.gdp$pol <- list()

## US - 1985 to 2013
g.gdp$us$data <- fread('usgdp_qtr_ch13_sec6.csv')
g.gdp$us$y <- g.gdp$us$data[, .(date, gdp, infl, irate)]
g.gdp$us$y[, c('gdp.new', 'infl.new', 'irate.new') :=
             list(g.gdp$us$y[2:114, gdp], g.gdp$us$y[2:114, infl],
                  g.gdp$us$y[2:114, irate])]
g.gdp$us$data[, c('date', 'gdp', 'infl', 'irate') := NULL]

# VAR - 2 lags
VARselect(g.gdp$us$y[, .(gdp, infl, irate)])
for (t in 1:51) {
  g.gdp$us$var <- VAR(g.gdp$us$y[1:(62 + t), .(gdp, infl, irate)], p = 2)
  g.gdp$us$pred$var[t] <- as.numeric(predict(g.gdp$us$var,
                                             n.ahead = 1)$fcst$gdp[, 1])
}

# SW - 4 PCs
g.gdp$us$pc <- principal(g.gdp$us$data, nfactors = 4,
                         rotate = 'varimax')$scores
g.gdp$us$y[, c('PC1', 'PC2', 'PC3', 'PC4') := list(g.gdp$us$pc[, 1],
                                                   g.gdp$us$pc[, 2],
                                                   g.gdp$us$pc[, 3],
                                                   g.gdp$us$pc[, 4])] 
g.gdp$us$sw <- lm(gdp.new ~ PC1 + PC2 + PC3 + PC4, data = g.gdp$us$y[1:62])
g.gdp$us$pred$sw <- as.numeric(predict(g.gdp$us$sw,
                                       newdata = g.gdp$us$y[63:113]))

# SWTF - 4 PCs with higher correlation factors
g.gdp$us$tf <- data.table('gdp' = g.gdp$us$y[, gdp])
for (var in names(g.gdp$us$data)) {
  if (abs(cor(g.gdp$us$data[[var]], g.gdp$us$y[, gdp])) >= 0.4) {
    g.gdp$us$tf <- cbind(g.gdp$us$tf, g.gdp$us$data[[var]])
  }
}
g.gdp$us$tf[, gdp := NULL]
g.gdp$us$tpc <- principal(g.gdp$us$tf, nfactors = 4,
                          rotate = 'varimax')$scores
g.gdp$us$y[, c('TPC1', 'TPC2', 'TPC3', 'TPC4') :=
             list(g.gdp$us$tpc[, 1], g.gdp$us$tpc[, 2],
                  g.gdp$us$tpc[, 3], g.gdp$us$tpc[, 4])]
g.gdp$us$swtf <- lm(gdp.new ~ TPC1 + TPC2 + TPC3 + TPC4,
                    data = g.gdp$us$y[1:62])
g.gdp$us$pred$swtf <- as.numeric(predict(g.gdp$us$swtf,
                                         newdata = g.gdp$us$y[63:113]))

# FAVAR - 3 PCs & 1 lag
VARselect(g.gdp$us$y[, .(gdp, infl, irate, PC1, PC2, PC3)])
for (t in 1:51) {
  g.gdp$us$favar <- VAR(g.gdp$us$y[1:(62 + t),
                                   .(gdp, infl, irate, PC1, PC2, PC3)])
  g.gdp$us$pred$favar[t] <- as.numeric(predict(g.gdp$us$favar,
                                               n.ahead = 1)$fcst$gdp[, 1])
}

# TPRF - 4 factors
for (t in 1:51) {
  g.gdp$us$rf <- data.table('gdp' = g.gdp$us$y[1:(62 + t), gdp.new])
  g.gdp$us$rf[, rf1 := TPRF(g.gdp$us$data[1:(62 + t)], gdp,
                            Z = gdp)$factors[, 2]]
  g.gdp$us$rf[, e1 := lm(gdp ~ rf1)$residuals]
  g.gdp$us$rf[, c('rf1', 'rf2') :=
                list(TPRF(g.gdp$us$data[1:(62 + t)], gdp,
                          Z = cbind(gdp, e1))$factors[, 2],
                     TPRF(g.gdp$us$data[1:(62 + t)], gdp,
                          Z = cbind(gdp, e1))$factors[, 3])]
  g.gdp$us$rf[, e2 := lm(gdp ~ rf1 + rf2)$residuals]
  g.gdp$us$rf[, c('rf1', 'rf2', 'rf3') :=
                list(TPRF(g.gdp$us$data[1:(62 + t)], gdp,
                          Z = cbind(gdp, e1, e2))$factors[, 2],
                     TPRF(g.gdp$us$data[1:(62 + t)], gdp,
                          Z = cbind(gdp, e1, e2))$factors[, 3],
                     TPRF(g.gdp$us$data[1:(62 + t)], gdp,
                          Z = cbind(gdp, e1, e2))$factors[, 4])]
  g.gdp$us$rf[, e3 := lm(gdp ~ rf1 + rf2 + rf3)$residuals]
  g.gdp$us$rf[, c('rf1', 'rf2', 'rf3', 'rf4') :=
                list(TPRF(g.gdp$us$data[1:(62 + t)], gdp,
                          Z = cbind(gdp, e1, e2, e3))$factors[, 2],
                     TPRF(g.gdp$us$data[1:(62 + t)], gdp,
                          Z = cbind(gdp, e1, e2, e3))$factors[, 3],
                     TPRF(g.gdp$us$data[1:(62 + t)], gdp,
                          Z = cbind(gdp, e1, e2, e3))$factors[, 4],
                     TPRF(g.gdp$us$data[1:(62 + t)], gdp,
                          Z = cbind(gdp, e1, e2, e3))$factors[, 5])]
  g.gdp$us$tprf <- lm(gdp ~ rf1 + rf2 + rf3 + rf4,
                      data = g.gdp$us$rf[1:(62 + t)])
  g.gdp$us$pred$tprf[t] <- g.gdp$us$tprf$fitted.values[(62 + t)]
}

## Plots
# First PC vs. First TPC
g.gdp$us$pcall <- data.table('Date' = as.Date(g.gdp$us$y[, date], '%m/%d/%Y'),
                             'PC' = c(g.gdp$us$y[, PC1], g.gdp$us$y[, TPC1]),
                             'tag' = rep(c('FIRSTPC', 'FIRSTPCTARG'),
                                         each = 114))
ggplot(data = g.gdp$us$pcall, aes(x = Date, y = PC,
                                  group = tag, colour = tag)) +
  theme_classic() + theme(plot.title = element_text(hjust = 0.5)) + 
  geom_line() + scale_x_date() + xlab('') + ylab('')

# Forecasts
g.gdp$us$yhat <- data.table('Date' = g.gdp$us$pcall[64:114, Date],
                            'Y' = c(g.gdp$us$y[64:114, gdp],
                                    g.gdp$us$pred$var, g.gdp$us$pred$sw,
                                    g.gdp$us$pred$swtf, g.gdp$us$pred$favar,
                                    g.gdp$us$pred$tprf),
                            'tag' = rep(c('Y', 'VAR', 'SW', 'SWTF',
                                          'FAVAR', 'TPRF'), each = 51))
ggplot(data = g.gdp$us$yhat, aes(x = Date, y = Y,
                                 group = tag, colour = tag)) +
  theme_classic() + theme(plot.title = element_text(hjust = 0.5)) + 
  geom_line() + scale_x_date() + xlab('') + ylab('')

## Evaluation
model <- c('var', 'sw', 'swtf', 'favar', 'tprf')

# RMSE
g.gdp$us$RMSE <- data.table()
for (m in model) {
  g.gdp$us$rmse <- data.table('Model' = toupper(m),
                              'RMSE' = rmse(g.gdp$us$y[64:114, gdp],
                                            g.gdp$us$pred[[m]]))
  g.gdp$us$RMSE <- rbindlist(list(g.gdp$us$RMSE, g.gdp$us$rmse))
}
print(g.gdp$us$RMSE)

# Diebold-Mariano Test
g.gdp$us$DM <- data.table()
for (m in 1:4) {
  for (n in (m + 1):5) {
    us.dm <- dm.test(g.gdp$us$pred[[model[m]]] - g.gdp$us$y[64:114, gdp],
                     g.gdp$us$pred[[model[n]]] - g.gdp$us$y[64:114, gdp])
    g.gdp$us$dm <- data.table('M1' = model[m], 'M2' = model[n],
                              'DM' = us.dm$p.value)
    g.gdp$us$DM <- rbindlist(list(g.gdp$us$DM, g.gdp$us$dm))
  }
}
print(g.gdp$us$DM)

## EA & Poland - 1995 to 2013
for (geo in c('ea', 'pol')) {
  g.gdp[[geo]]$data <- fread(paste0(geo, 'gdp_qtr_ch13_sec6.csv'))
  g.gdp[[geo]]$y <- g.gdp[[geo]]$data[, .(date, gdp, infl, irate)]
  g.gdp[[geo]]$y[, c('gdp.new', 'infl.new', 'irate.new') :=
                   list(g.gdp[[geo]]$y[2:74, gdp], g.gdp[[geo]]$y[2:74, infl],
                        g.gdp[[geo]]$y[2:74, irate])]
  g.gdp[[geo]]$data[, c('date', 'gdp', 'infl', 'irate') := NULL]
  
  # VAR - 1 lag
  VARselect(g.gdp[[geo]]$y[, .(gdp, infl, irate)])
  for (t in 1:35) {
    g.gdp[[geo]]$var <- VAR(g.gdp[[geo]]$y[1:(38 + t), .(gdp, infl, irate)])
    g.gdp[[geo]]$pred$var[t] <- as.numeric(predict(g.gdp[[geo]]$var,
                                                   n.ahead = 1)$fcst$gdp[, 1])
  }
  
  # SW - 4 PCs
  g.gdp[[geo]]$pc <- principal(g.gdp[[geo]]$data, nfactors = 4,
                               rotate = 'varimax')$scores
  g.gdp[[geo]]$y[, c('PC1', 'PC2', 'PC3', 'PC4') :=
                   list(g.gdp[[geo]]$pc[, 1], g.gdp[[geo]]$pc[, 2],
                        g.gdp[[geo]]$pc[, 3], g.gdp[[geo]]$pc[, 4])] 
  g.gdp[[geo]]$sw <- lm(gdp.new ~ PC1 + PC2 + PC3 + PC4,
                        data = g.gdp[[geo]]$y[1:38])
  g.gdp[[geo]]$pred$sw <- as.numeric(predict(g.gdp[[geo]]$sw,
                                         newdata = g.gdp[[geo]]$y[39:73]))
  
  # SWTF - 4 PCs with higher correlation factors
  g.gdp[[geo]]$tf <- data.table('gdp' = g.gdp[[geo]]$y[, gdp])
  for (var in names(g.gdp[[geo]]$data)) {
    if (abs(cor(g.gdp[[geo]]$data[[var]], g.gdp[[geo]]$y[, gdp])) >= 0.4) {
      g.gdp[[geo]]$tf <- cbind(g.gdp[[geo]]$tf, g.gdp[[geo]]$data[[var]])
    }
  }
  g.gdp[[geo]]$tf[, gdp := NULL]
  g.gdp[[geo]]$tpc <- principal(g.gdp[[geo]]$tf, nfactors = 4,
                            rotate = 'varimax')$scores
  g.gdp[[geo]]$y[, c('TPC1', 'TPC2', 'TPC3', 'TPC4') :=
               list(g.gdp[[geo]]$tpc[, 1], g.gdp[[geo]]$tpc[, 2],
                    g.gdp[[geo]]$tpc[, 3], g.gdp[[geo]]$tpc[, 4])]
  g.gdp[[geo]]$swtf <- lm(gdp.new ~ TPC1 + TPC2 + TPC3 + TPC4,
                      data = g.gdp[[geo]]$y[1:38])
  g.gdp[[geo]]$pred$swtf <- as.numeric(predict(g.gdp[[geo]]$swtf,
                                           newdata = g.gdp[[geo]]$y[39:73]))
  
  # FAVAR - 3 PCs & 1 lag
  VARselect(g.gdp[[geo]]$y[, .(gdp, infl, irate, PC1, PC2, PC3)])
  for (t in 1:35) {
    g.gdp[[geo]]$favar <- VAR(g.gdp[[geo]]$y[1:(38 + t),
                                     .(gdp, infl, irate, PC1, PC2, PC3)])
    g.gdp[[geo]]$pred$favar[t] <-
      as.numeric(predict(g.gdp[[geo]]$favar, n.ahead = 1)$fcst$gdp[, 1])
  }
  
  # TPRF - 4 factors
  for (t in 1:35) {
    g.gdp[[geo]]$rf <- data.table('gdp' = g.gdp[[geo]]$y[1:(38 + t), gdp.new])
    g.gdp[[geo]]$rf[, rf1 := TPRF(g.gdp[[geo]]$data[1:(38 + t)], gdp,
                              Z = gdp)$factors[, 2]]
    g.gdp[[geo]]$rf[, e1 := lm(gdp ~ rf1)$residuals]
    g.gdp[[geo]]$rf[, c('rf1', 'rf2') :=
                      list(TPRF(g.gdp[[geo]]$data[1:(38 + t)], gdp,
                                Z = cbind(gdp, e1))$factors[, 2],
                           TPRF(g.gdp[[geo]]$data[1:(38 + t)], gdp,
                                Z = cbind(gdp, e1))$factors[, 3])]
    g.gdp[[geo]]$rf[, e2 := lm(gdp ~ rf1 + rf2)$residuals]
    g.gdp[[geo]]$rf[, c('rf1', 'rf2', 'rf3') :=
                  list(TPRF(g.gdp[[geo]]$data[1:(38 + t)], gdp,
                            Z = cbind(gdp, e1, e2))$factors[, 2],
                       TPRF(g.gdp[[geo]]$data[1:(38 + t)], gdp,
                            Z = cbind(gdp, e1, e2))$factors[, 3],
                       TPRF(g.gdp[[geo]]$data[1:(38 + t)], gdp,
                            Z = cbind(gdp, e1, e2))$factors[, 4])]
    g.gdp[[geo]]$rf[, e3 := lm(gdp ~ rf1 + rf2 + rf3)$residuals]
    g.gdp[[geo]]$rf[, c('rf1', 'rf2', 'rf3', 'rf4') :=
                  list(TPRF(g.gdp[[geo]]$data[1:(38 + t)], gdp,
                            Z = cbind(gdp, e1, e2, e3))$factors[, 2],
                       TPRF(g.gdp[[geo]]$data[1:(38 + t)], gdp,
                            Z = cbind(gdp, e1, e2, e3))$factors[, 3],
                       TPRF(g.gdp[[geo]]$data[1:(38 + t)], gdp,
                            Z = cbind(gdp, e1, e2, e3))$factors[, 4],
                       TPRF(g.gdp[[geo]]$data[1:(38 + t)], gdp,
                            Z = cbind(gdp, e1, e2, e3))$factors[, 5])]
    g.gdp[[geo]]$tprf <- lm(gdp ~ rf1 + rf2 + rf3 + rf4,
                        data = g.gdp[[geo]]$rf[1:(38 + t)])
    g.gdp[[geo]]$pred$tprf[t] <- g.gdp[[geo]]$tprf$fitted.values[(38 + t)]
  }
  
  ## Plots
  # First PC vs. First TPC
  g.gdp[[geo]]$pcall <- data.table('Date' = as.Date(g.gdp[[geo]]$y[, date],
                                                    '%m/%d/%Y'),
                                   'PC' = c(g.gdp[[geo]]$y[, PC1],
                                            g.gdp[[geo]]$y[, TPC1]),
                                   'tag' = rep(c('FIRSTPC', 'FIRSTPCTARG'),
                                               each = 74))
  ggplot(data = g.gdp[[geo]]$pcall, aes(x = Date, y = PC,
                                        group = tag, colour = tag)) +
    theme_classic() + theme(plot.title = element_text(hjust = 0.5)) + 
    geom_line() + scale_x_date() + xlab('') + ylab('')
  
  # Forecasts
  g.gdp[[geo]]$yhat <- data.table('Date' = g.gdp[[geo]]$pcall[40:74, Date],
                                  'Y' = c(g.gdp[[geo]]$y[40:74, gdp],
                                          g.gdp[[geo]]$pred$var,
                                          g.gdp[[geo]]$pred$sw,
                                          g.gdp[[geo]]$pred$swtf,
                                          g.gdp[[geo]]$pred$favar,
                                          g.gdp[[geo]]$pred$tprf),
                                  'tag' = rep(c('Y', 'VAR', 'SW', 'SWTF',
                                                'FAVAR', 'TPRF'), each = 35))
  ggplot(data = g.gdp[[geo]]$yhat, aes(x = Date, y = Y,
                                       group = tag, colour = tag)) +
    theme_classic() + theme(plot.title = element_text(hjust = 0.5)) + 
    geom_line() + scale_x_date() + xlab('') + ylab('')
  
  model <- c('var', 'sw', 'swtf', 'favar', 'tprf')
  
  # RMSE
  g.gdp[[geo]]$RMSE <- data.table()
  for (m in model) {
    g.gdp[[geo]]$rmse <- data.table('Model' = toupper(m),
                                    'RMSE' = rmse(g.gdp[[geo]]$y[40:74, gdp],
                                                  g.gdp[[geo]]$pred[[m]]))
    g.gdp[[geo]]$RMSE <- rbindlist(list(g.gdp[[geo]]$RMSE, g.gdp[[geo]]$rmse))
  }
  print(g.gdp[[geo]]$RMSE)
  
  # Diebold-Mariano Test
  g.gdp[[geo]]$DM <- data.table()
  for (m in 1:4) {
    for (n in (m + 1):5) {
      dm <- dm.test(g.gdp[[geo]]$pred[[model[m]]] -
                      g.gdp[[geo]]$y[40:74, gdp],
                    g.gdp[[geo]]$pred[[model[n]]] -
                      g.gdp[[geo]]$y[40:74, gdp])
      g.gdp[[geo]]$dm <- data.table('M1' = model[m], 'M2' = model[n],
                                    'DM' = dm$p.value)
      g.gdp[[geo]]$DM <- rbindlist(list(g.gdp[[geo]]$DM, g.gdp[[geo]]$dm))
    }
  }
  print(g.gdp[[geo]]$DM)
}