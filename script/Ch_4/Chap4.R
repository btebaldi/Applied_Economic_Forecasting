rm(list = ls())
setwd('E:/Book/Ch_4') # replace with own directory

library(data.table)
library(forecast)
library(xtable)

### 4.9 Examples using simulated data ###
sim.data <- fread('simulated_datac4.csv')

sim.fit <- list()
formula.list <- c('y ~ x', 'y ~ D + x + D * x', 'y ~ D + x + D * x + y1 + x1')

yhat <- list()
yhat.rec <- list()
for (f in 1:3) {
  sim.fit[[f]] <- lm(formula.list[f], data = sim.data[102:301])
  yhat[[f]] <- predict(sim.fit[[f]], newdata = sim.data[302:501])
  
  yhat.rec[[f]] <- rep(0, 200)
  for (i in 1:200) {
    ols.rec <- lm(formula.list[f], data = sim.data[102:(300 + i)])
    yhat.rec[[f]][i] <- predict(ols.rec, newdata = sim.data[301 + i])
  }
}

# RMSE & MAE
yhat$Y <- cbind(yhat[[1]], yhat[[2]], yhat[[3]],
                yhat.rec[[1]], yhat.rec[[2]], yhat.rec[[3]])
RMSE <- sqrt(colSums((yhat$Y - sim.data[302:501, y])^2) / 200)
MAE <- colSums(abs(yhat$Y - sim.data[302:501, y])) / 200
error.table <- data.table('Forecasting Model' = c('Linear Regression',
                                                  'Dummy variable model',
                                                  'DGP'),
                          'RMSE' = RMSE[1:3], 'RMSE recursive' = RMSE[4:6],
                          'MAE' = MAE[1:3], 'MAE recursive' = MAE[4:6])
print(xtable(error.table), include.rownames = F)

## Pairwise comparison
# Diebold-Mariano test
dm.mat <- list()
model <- list(c(1, 2), c(1, 3), c(2, 3))

for (m in 1:3) {
  dm.mat[[m]] <- dm.test(sim.fit[[model[[m]][1]]]$residuals,
                         sim.fit[[model[[m]][2]]]$residuals)
}
dm.table <- data.table('Compared Models' = c('DM test statistics', 'P-Value'),
                       'M1 vs M2' = c(dm.mat[[1]]$statistic,
                                      dm.mat[[1]]$p.value),
                       'M1 vs M3' = c(dm.mat[[2]]$statistic,
                                      dm.mat[[2]]$p.value),
                       'M2 vs M3' = c(dm.mat[[3]]$statistic,
                                      dm.mat[[3]]$p.value))
print(xtable(dm.table), include.rownames = F)

# Morgan-Granger-Newbold test
mgn.table <- data.table()

for (m in 1:2) {
  for (M in (m + 1):3) {
    err.sum <- sim.fit[[m]]$residuals + sim.fit[[M]]$residuals
    err.dif <- sim.fit[[m]]$residuals - sim.fit[[M]]$residuals
    mgn <- sum(err.sum * err.dif) / sqrt(sum(err.sum ^ 2) * sum(err.dif ^ 2))
    
    # Compuate p-values for the MGN statistic
    tval <- mgn / sqrt((1 - mgn ^ 2) / 198)
    pval <- 2 * pt(abs(tval), 198, lower.tail = F)
    
    mgn.result <- data.table('t-stat' = tval, 'p-val' = pval)
    mgn.table <- rbindlist(list(mgn.table, mgn.result))
  }
}
mgn.table <- t(mgn.table)
mgn.table <- cbind(c('MGN test statistics', 'P-value'), data.table(mgn.table))
setnames(mgn.table, names(dm.table))
print(xtable(dm.table), include.colnames = T)

# Recursive
# Diebold-Mariano test
DM.mat <- list()
model <- list(c(1, 2), c(1, 3), c(2, 3))

for (m in 1:3) {
  DM.mat[[m]] <- dm.test(yhat.rec[[model[[m]][1]]] - sim.data[302:501, y],
                         yhat.rec[[model[[m]][2]]] - sim.data[302:501, y])
}
DM.table <- data.table('Compared Models' = c('DM test statistics', 'P-Value'),
                       'M1 vs M2' = c(DM.mat[[1]]$statistic,
                                      DM.mat[[1]]$p.value),
                       'M1 vs M3' = c(DM.mat[[2]]$statistic,
                                      DM.mat[[2]]$p.value),
                       'M2 vs M3' = c(DM.mat[[3]]$statistic,
                                      DM.mat[[3]]$p.value))

# Morgan-Granger-Newbold test
MGN.table <- data.table()

for (m in 1:2) {
  for (M in (m + 1):3) {
    err.sum <- yhat.rec[[m]] + yhat.rec[[M]] - 2 * sim.data[302:501, y]
    err.dif <- yhat.rec[[m]] - yhat.rec[[M]]
    mgn <- sum(err.sum * err.dif) / sqrt(sum(err.sum ^ 2) * sum(err.dif ^ 2))
    
    # Compuate p-values for the MGN statistic
    tval <- mgn / sqrt((1 - mgn ^ 2) / 198)
    pval <- 2 * pt(abs(tval), 198, lower.tail = F)
    
    mgn.result <- data.table('t-stat' = tval, 'p-val' = pval)
    MGN.table <- rbindlist(list(MGN.table, mgn.result))
  }
}
MGN.table <- t(MGN.table)
MGN.table <- cbind(c('MGN test statistics', 'P-value'), data.table(MGN.table))
setnames(MGN.table, names(dm.table))
print(xtable(rbind(DM.table, MGN.table)), include.rownames = F)


### 4.10.1 Forecasting Euro Area GDP ###
eu.gdp <- fread('ex2_forecastevaluation_gdp.csv')

gdp.fit <- list()
gdp.hat <- list()
formula.list <- c('Model1' = 'y ~ ipr + su + pr + sr',
                  'Model2' = 'y ~ ipr + su + sr',
                  'Model3' = 'y ~ ipr + su',
                  'Model4' = 'y ~ ipr + pr + sr',
                  'ARDL' = 'y ~ y1 + ipr + su1 + sr1')
for (f in names(formula.list)) {
  
  # Simple forecasts
  gdp.fit[[f]] <- lm(formula.list[[f]], data = eu.gdp[1:44])
  gdp.hat$Simple[[f]] <- predict(gdp.fit[[f]], newdata = eu.gdp[45:70])
  
  # Recursive forecasts
  gdp.hat$Rec[[f]] <- rep(0, 26)
  for (i in 1:26) {
    ols.rec <- lm(formula.list[[f]], data = eu.gdp[1:(43 + i)])
    gdp.hat$Rec[[f]][i] <- predict(ols.rec, newdata = eu.gdp[44 + i])
  }
}
summary(gdp.fit$ARDL) # ARDL Model 2

## Diebold-Mariano test
gdp.DM <- list()
for (m in c('Model1', 'Model2', 'Model3', 'Model4')) {
  gdp.DM[[m]] <- dm.test(gdp.hat$Simple[[m]] - eu.gdp[45:70, y],
                         gdp.hat$Simple$ARDL - eu.gdp[45:70, y])
}
DM.table <- data.frame('vs M1' = c(gdp.DM$Model1$statistic,
                                   gdp.DM$Model1$p.value),
                       'vs M2' = c(gdp.DM$Model2$statistic,
                                   gdp.DM$Model2$p.value),
                       'vs M3' = c(gdp.DM$Model3$statistic,
                                   gdp.DM$Model3$p.value),
                       'vs M4' = c(gdp.DM$Model4$statistic,
                                   gdp.DM$Model4$p.value))
rownames(DM.table) <- c('Test Stat.', 'p-value')
print(xtable(DM.table))

# RMSE & MAE
gdp.hat$Y <- cbind(gdp.hat$Simple$Model1, gdp.hat$Simple$Model2,
                   gdp.hat$Simple$Model3, gdp.hat$Simple$Model4,
                   gdp.hat$Simple$ARDL)
RMSE <- sqrt(colSums((gdp.hat$Y - eu.gdp[45:70, y])^2) / 26)
MAE <- colSums(abs(gdp.hat$Y - eu.gdp[45:70, y])) / 26
error.mat <- rbind(RMSE, MAE)
colnames(error.mat) <- c('Model 1', 'Model 2', 'Model 3',
                         'Model 4', 'ARDL M2')
print(xtable(error.mat))

## Morgan-Granger-Newbold test
gdp.MGN <- list('Simple' = data.table(), 'Recursive' = data.table())

for (f in c('Simple', 'Rec')) {
  for (m in c('Model1', 'Model2', 'Model3', 'Model4')) {
    err.sum <- gdp.hat[[f]][[m]] + gdp.hat[[f]]$ARDL - 2 * eu.gdp[45:70, y]
    err.dif <- gdp.hat[[f]][[m]] - gdp.hat[[f]]$ARDL
    mgn <- sum(err.sum * err.dif) / sqrt(sum(err.sum ^ 2) * sum(err.dif ^ 2))
    
    # Compute p-values for the MGN statistic
    tval <- mgn / sqrt((1 - mgn ^ 2) / 24)
    pval <- 2 * pt(abs(tval), 24, lower.tail = F)
    
    mgn.result <- data.table('t-stat' = tval, 'p-val' = pval)
    gdp.MGN[[f]] <- rbindlist(list(gdp.MGN[[f]], mgn.result))
  }
  gdp.MGN[[f]] <- t(gdp.MGN[[f]])
  rownames(gdp.MGN[[f]]) <- c('Test Stat.', 'p-value')
  colnames(gdp.MGN[[f]]) <- c('vs M1', 'vs M2', 'vs M3', 'vs M4')
  print(xtable(gdp.MGN[[f]]))
}


### 4.10.2 Forecasting US GDP ###
us.gdp <- fread('ex2_forecastevaluation_gdp_us.csv')

GDP.fit <- list()
GDP.hat <- list()

formula.list <- c('Model1' = 'y ~ ipr + su + pr + sr',
                  'Model2' = 'y ~ ipr + su + sr',
                  'Model3' = 'y ~ ipr + su',
                  'Model4' = 'y ~ ipr + pr + sr',
                  'ARDL' = 'y ~ y1 + ipr + su1 + sr1')
for (f in names(formula.list)) {
  
  # Simple forecasts
  GDP.fit[[f]] <- lm(formula.list[[f]], data = us.gdp[1:88])
  GDP.hat$Simple[[f]] <- predict(GDP.fit[[f]], newdata = us.gdp[89:114])
  
  # Recursive forecasts
  GDP.hat$Rec[[f]] <- rep(0, 26)
  for (i in 1:26) {
    ols.rec <- lm(formula.list[[f]], data = us.gdp[1:(87 + i)])
    GDP.hat$Rec[[f]][i] <- predict(ols.rec, newdata = us.gdp[88 + i])
  }
}

summary(GDP.fit$Model4)
summary(GDP.fit$ARDL)

# RMSE & MAE
GDP.hat$Y <- cbind(GDP.hat$Simple$Model1, GDP.hat$Simple$Model2,
                   GDP.hat$Simple$Model3, GDP.hat$Simple$Model4,
                   GDP.hat$Simple$ARDL)
RMSE <- sqrt(colSums((GDP.hat$Y - us.gdp[89:114, y])^2) / 26)
MAE <- colSums(abs(GDP.hat$Y - us.gdp[89:114, y])) / 26
error.mat <- rbind(RMSE, MAE)
colnames(error.mat) <- c('Model 1', 'Model 2', 'Model 3',
                         'Model 4', 'ARDL M2')
print(xtable(error.mat))

## Diebold-Mariano test
GDP.DM <- list()
DM.table <- list()

for (m in c('Model1', 'Model2', 'Model3', 'Model4')) {
  # Simple
  GDP.DM$Simple[[m]] <- dm.test(GDP.hat$Simple[[m]] - us.gdp[89:114, y],
                                GDP.hat$Simple$ARDL - us.gdp[89:114, y])
  
  # Recursive
  GDP.DM$Rec[[m]] <- dm.test(GDP.hat$Rec[[m]] - us.gdp[89:114, y],
                             GDP.hat$Rec$ARDL - us.gdp[89:114, y])
}

for (f in c('Simple', 'Rec')) {
  DM.table[[f]] <- data.frame('vs M1' = c(GDP.DM[[f]]$Model1$statistic,
                                          GDP.DM[[f]]$Model1$p.value),
                              'vs M2' = c(GDP.DM[[f]]$Model2$statistic,
                                          GDP.DM[[f]]$Model2$p.value),
                              'vs M3' = c(GDP.DM[[f]]$Model3$statistic,
                                          GDP.DM[[f]]$Model3$p.value),
                              'vs M4' = c(GDP.DM[[f]]$Model4$statistic,
                                          GDP.DM[[f]]$Model4$p.value))
  rownames(DM.table[[f]]) <- c('Test Stat.', 'p-value')
  print(xtable(DM.table[[f]]))
}

## Morgan-Granger-Newbold test
GDP.MGN <- list('Simple' = data.table(), 'Rec' = data.table())

for (f in c('Simple', 'Rec')) {
  for (m in c('Model1', 'Model2', 'Model3', 'Model4')) {
    err.sum <- GDP.hat[[f]][[m]] + GDP.hat[[f]]$ARDL- 2 * us.gdp[89:114, y]
    err.dif <- GDP.hat[[f]][[m]] - GDP.hat[[f]]$ARDL
    mgn <- sum(err.sum * err.dif) / sqrt(sum(err.sum ^ 2) * sum(err.dif ^ 2))
    
    # Compute p-values for the MGN statistic
    tval <- mgn / sqrt((1 - mgn ^ 2) / 24)
    pval <- 2 * pt(abs(tval), 24, lower.tail = F)
    
    mgn.result <- data.table('t-stat' = tval, 'p-val' = pval)
    GDP.MGN[[f]] <- rbindlist(list(GDP.MGN[[f]], mgn.result))
  }
  GDP.MGN[[f]] <- t(GDP.MGN[[f]])
  rownames(GDP.MGN[[f]]) <- c('Test Stat.', 'p-value')
  colnames(GDP.MGN[[f]]) <- c('vs M1', 'vs M2', 'vs M3', 'vs M4')
  print(xtable(GDP.MGN[[f]]))
}

### 4.10.3 Default Risk ###
oas.data <- fread('default_risk.csv')
oas.data[, Date := as.Date(Date, '%m/%d/%Y')]
oas.data[, OASnew := OAS[2:216]]
oas.data <- oas.data[1:215]

model <- c('OASnew ~ OAS + VIX', 'OASnew ~ OAS + SENT',
           'OASnew ~ OAS + PMI', 'OASnew ~ OAS + sp500')

oas.fit <- list()
BIC <- HQIC <- RMSE <- rep(0, 4)
forecast.OS <- data.table()

# Model forecasts
for (m in 1:4) {
  oas.fit[[m]] <- lm(model[m], data = oas.data[Date <= '2010-11-01'])
  print(summary(oas.fit[[m]]))
  print(dwtest(oas.fit[[m]]))
  
  BIC[m] <- BIC(oas.fit[[m]])
  HQIC[m] <- -2 * logLik(oas.fit[[m]]) + 6 * log(log(155))
  
  oas.forecast <- predict(oas.fit[[m]],
                          newdata = oas.data[Date >= '2010-12-01' &
                                               Date <= '2013-11-01'])
  RMSE[m] <- rmse(oas.data[Date >= '2011-01-01' & Date <= '2013-12-01', OAS],
                  oas.forecast)
  forecast.os <- data.table('OAS' = oas.data[Date >= '2011-01-01' &
                                               Date <= '2013-12-01', OAS],
                            'OAShat' = oas.forecast, 'model' = m)
  forecast.OS <- rbindlist(list(forecast.OS, forecast.os))
}

## Morgan-Granger-Newbold test
forecast.OS[, error := OAShat - OAS]
n.OOS <- forecast.OS[, .N] / 4
MGN.summary <- data.table()
for (m in 1:3) {
  for (M in (m + 1):4) {
    err.sum <- forecast.OS[model == m, error] + forecast.OS[model == M, error]
    err.dif <- forecast.OS[model == m, error] - forecast.OS[model == M, error]
    mgn <- sum(err.sum * err.dif) / sqrt(sum(err.sum ^ 2) * sum(err.dif ^ 2))
    
    # Compuate p-values for the MGN statistic
    tval <- mgn / sqrt((1 - mgn ^ 2) / (n.OOS - 2))
    pval <- 2 * pt(abs(tval), n.OOS - 2, lower.tail = F)
    
    # Indicator of the superior model
    indicator <- ifelse(tval > 0, 2, 1)
    
    mgn.result <- data.table('Model 1' = model[m], 'Model 2' = model[M],
                             't-stat' = tval, 'p-val' = pval,
                             'Superior' = indicator)
    MGN.summary <- rbindlist(list(MGN.summary, mgn.result))
  }
}
print(MGN.summary)