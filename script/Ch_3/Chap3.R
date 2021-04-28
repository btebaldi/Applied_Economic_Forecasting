rm(list = ls())
setwd('E:/Book/Ch_3') # replace with own directory

library(data.table)
library(ggplot2)
library(lmtest)
library(Metrics)
library(stargazer)
library(xtable)

### 3.6 Simulated Data ###
sim.data <- fread('simulated_datac3.csv')
sim.est <- sim.data[102:301]

ic.table <- data.table()
formula.list <- paste0('y0 ~ D + x0 + D * x0 ',
                       c('+ y1', '+ y1 + y2', '+ x1', '+ x1 + y1',
                         '+ x1 + y1 + y2', '+ x1 + x2', '+ x1 + x2 + y1',
                         '+ x1 + x2 + y1 + y2'))

## Identify ARDL(1, 1) as the best model
for (f in 1:8) {
  ardl.fit <- lm(formula.list[f], data = sim.est)
  ic <- data.table(cbind(AIC(ardl.fit), BIC(ardl.fit)))
  ic.table <- rbindlist(list(ic.table, ic))
}

ardl.fit <- lm(formula.list[4], data = sim.est)
summary(ardl.fit)

## Simple & recursive forecasts
yhat <- list()
yhat$y <- predict(ardl.fit, newdata = sim.data[302:501])

yhat$yrec <- rep(0, 200)
for(i in 1:200) {
  ardl.rec <- lm(formula.list[4], data = sim.data[102:(300 + i)])
  yhat$yrec[i] <- predict(ardl.rec, newdata = sim.data[301 + i])
}

# RMSE & MAE
yhat$Y <- cbind(yhat$y, yhat$yrec)
RMSE <- sqrt(colSums((yhat$Y - sim.data[302:501, y0])^2) / 200)
MAE <- colSums(abs(yhat$Y - sim.data[302:501, y0])) / 200
error.mat <- rbind(RMSE, MAE)
colnames(error.mat) <- c('Simple', 'Recursive')
print(xtable(error.mat), include.rownames = T, include.colnames = T)


### 3.7.1 Forecasting Euro Area GDP ###
eu.gdp <- fread('ex2_dynamic_gdp.csv')
gdp.ardl <- list()

# Model 1
gdp.ardl$fit1 <- lm(y ~ y1 + ipr + ipr1 + su + su1 + sr + sr1,
                    data = eu.gdp[2:60])
summary(gdp.ardl$fit1)

# Model 2
gdp.ardl$fit2 <- lm(y ~ y1 + ipr + su1 + sr1, data = eu.gdp[2:60])
summary(gdp.ardl$fit2)

eu.gdp[2:60, c('eps1', 'eps2') := list(gdp.ardl$fit1$residuals,
                                       gdp.ardl$fit2$residuals)]
eu.gdp[3:60, c('eps11', 'eps21') := list(eu.gdp[2:59, eps1],
                                         eu.gdp[2:59, eps2])]
eu.gdp[4:60, c('eps12', 'eps22') := list(eu.gdp[2:58, eps1],
                                         eu.gdp[2:58, eps2])]

# Breusch-Godfrey serial correlation LM test
bg.test <- list()
bg.test$fit1 <- lm(eps1 ~ y1 + ipr + ipr1 + su + su1 + sr + sr1 +
                     eps11 + eps12, data = eu.gdp[4:60])
summary(bg.test$fit1)

bg.test$fit2 <- lm(eps2 ~ y1 + ipr + su1 + sr1 + eps21 + eps22,
                   data = eu.gdp[4:60])
summary(bg.test$fit2)

## Forecasting with ARDL & dummies
gdp.ardl$fitD <- lm(y ~ y1 + ipr + su1 + sr1 + Dea + D2000s,
                    data = eu.gdp[2:60])
summary(gdp.ardl$fitD)

gdp.hat <- list()
gdp.fit <- lm(y ~ ipr + su + sr, data = eu.gdp[1:60]) # Model 2 - no dummy
gdp.hat$ghat <- predict(gdp.fit, newdata = eu.gdp[61:70]) 

gdp.fit <- lm(y ~ ipr + su + sr + Dea + D2000s, data = eu.gdp[1:60])
gdp.hat$ghat3 <- predict(gdp.fit, newdata = eu.gdp[61:70]) # Model 2.3 - dummy

gdp.hat$gardl <- predict(gdp.ardl$fitD, newdata = eu.gdp[61:70]) # ARDL

gdp.plot <- data.table('yhat' = rbindlist(list(data.table(eu.gdp[61:70, y]),
                                               data.table(gdp.hat$ghat),
                                               data.table(gdp.hat$ghat3),
                                               data.table(gdp.hat$gardl))),
                       'label' = rep(c('Y', 'YFOREG2_NEW',
                                       'YFOREG2_3', 'Y1STPARDL2'), each = 10))

ggplot(gdp.plot, aes(x = rep(1:10, 4), y = yhat, linetype = label)) +
  geom_line() + xlab('') + ylab('') + theme(legend.title = element_blank()) 

# RMSE & MAE
gdp.hat$Y <- cbind(gdp.hat$ghat, gdp.hat$ghat3, gdp.hat$gardl)
RMSE <- sqrt(colSums((gdp.hat$Y - eu.gdp[61:70, y])^2) / 10)
MAE <- colSums(abs(gdp.hat$Y - eu.gdp[61:70, y])) / 10
error.mat <- rbind(RMSE, MAE)
colnames(error.mat) <- c('Model 2', 'Model 2.3', 'ARDL Model 2')
print(xtable(error.mat), include.rownames = T, include.colnames = T)

## Recursive forecasts
gdp.hat$grec <- gdp.hat$grecD <- gdp.hat$grecARDL <- rep(0, 10)

formula.list <- list('grec' = 'y ~ ipr + su + sr',
                     'grecD' = 'y ~ ipr + su + sr + Dea + D2000s',
                     'grecARDL' = 'y ~ y1 + ipr + su1 + sr1 + Dea + D2000s')
for (f in names(formula.list)) {
  for (i in 1:10) {
    ardl.rec <- lm(formula.list[[f]], data = eu.gdp[1:(59 + i)])
    gdp.hat[[f]][i] <- predict(ardl.rec, newdata = eu.gdp[60 + i])
  }
}

# RMSE & MAE
gdp.hat$Yrec <- cbind(gdp.hat$grec, gdp.hat$grecD, gdp.hat$grecARDL)
RMSE <- sqrt(colSums((gdp.hat$Yrec - eu.gdp[61:70, y])^2) / 10)
MAE <- colSums(abs(gdp.hat$Yrec - eu.gdp[61:70, y])) / 10
error.mat <- rbind(RMSE, MAE)
colnames(error.mat) <- c('Model 2', 'Model 2.3', 'ARDL Model 2')
print(xtable(error.mat), include.rownames = T, include.colnames = T)

## Different forms of dynamic models
ardl.formula <- list('ardl' = paste0('y ~ y1 + ipr + ipr1 ',
                                     '+ su + su1 + sr + sr1 - 1'),
                     'static' = 'y ~ ipr + su + sr - 1',
                     'ar1' = 'y ~ y1 - 1',
                     'li' = 'y ~ ipr1 + su1 + sr1 - 1',
                     'dl' = 'y ~ ipr + ipr1 + su + su1 + sr + sr1 - 1')

## ARDL + Static + AR(1) + Leading indicator + Distributed lag
for (f in c('ardl', 'static', 'ar1', 'li', 'dl')) {
  print(summary(lm(ardl.formula[[f]], data = eu.gdp[2:70])))
}


### 3.7.2 Forecasting US GDP ###
us.gdp <- fread('ex2_dynamic_gdp_us.csv')

## ARDL with dummy variables
gdp.ARDL <- list()

# ARDL Model 1
gdp.ARDL$fit1 <- lm(y ~ y1 + ipr + ipr1 + su + su1 + sr + sr1 +
                      Dfincris + D2000s, data = us.gdp[2:104])
summary(gdp.ARDL$fit1)

# ARDL Model 2
gdp.ARDL$fit2 <- lm(y ~ y1 + ipr + su1 + sr1 + Dfincris + D2000s,
                    data = us.gdp[2:104])
summary(gdp.ARDL$fit2)

# Breusch-Godfrey serial correlation LM test
bgtest(gdp.ARDL$fit1)
bgtest(gdp.ARDL$fit2)

## Forecasts
gdp.ARDL$ghat1 <- predict(gdp.ARDL$fit1, newdata = us.gdp[105:114]) 
gdp.ARDL$ghat2 <- predict(lm(y ~ ipr + su + sr, data = us.gdp[1:104]),
                          newdata = us.gdp[105:114])
gdp.ARDL$ghat3 <- predict(lm(y ~ ipr + su + sr + Dfincris + D2000s,
                             data = us.gdp[1:104]), newdata = us.gdp[105:114])

gdp.plot <- data.table('yhat' = rbindlist(list(data.table(us.gdp[105:114, y]),
                                               data.table(gdp.ARDL$ghat1),
                                               data.table(gdp.ARDL$ghat2),
                                               data.table(gdp.ARDL$ghat3))),
                       'label' = rep(c('Y', 'YFORARDL1',
                                       'YFOREG2_NEW', 'YFOREG_3'), each = 10))

ggplot(gdp.plot, aes(x = rep(1:10, 4), y = yhat, linetype = label)) +
  geom_line() + xlab('') + ylab('') + theme(legend.title = element_blank()) 

# RMSE & MAE
gdp.ARDL$Y <- cbind(gdp.ARDL$ghat2, gdp.ARDL$ghat3, gdp.ARDL$ghat1)
RMSE <- sqrt(colSums((gdp.ARDL$Y - us.gdp[105:114, y])^2) / 10)
MAE <- colSums(abs(gdp.ARDL$Y - us.gdp[105:114, y])) / 10
error.mat <- rbind(RMSE, MAE)
colnames(error.mat) <- c('Model 2', 'Model 2.3', 'ARDL Model 1')
print(xtable(error.mat), include.rownames = T, include.colnames = T)

## Recursive forecasts
gdp.ARDL$grec <- gdp.ARDL$grecD <- gdp.ARDL$grecARDL <- rep(0, 10)
formula.list <- list('grec' = 'y ~ ipr + su + sr',
                     'grecD' = 'y ~ ipr + su + sr + Dfincris + D2000s',
                     'grecARDL' = paste0('y ~ y1 + ipr + ipr1 + su + su1 ',
                                         '+ sr + sr1 + Dfincris + D2000s'))
for (f in names(formula.list)) {
  for (i in 1:10) {
    ardl.rec <- lm(formula.list[[f]], data = us.gdp[1:(103 + i)])
    gdp.ARDL[[f]][i] <- predict(ardl.rec, newdata = us.gdp[104 + i])
  }
}

# RMSE & MAE
gdp.ARDL$Yrec <- cbind(gdp.ARDL$grecD, gdp.ARDL$grecARDL, gdp.ARDL$grec)
RMSE <- sqrt(colSums((gdp.ARDL$Yrec - us.gdp[105:114, y])^2) / 10)
MAE <- colSums(abs(gdp.ARDL$Yrec - us.gdp[105:114, y])) / 10
error.mat <- rbind(RMSE, MAE)
colnames(error.mat) <- c('Model 2', 'Model 2.3', 'ARDL Model 1')
print(xtable(error.mat), include.rownames = T, include.colnames = T)

## Different forms of dynamic models
ARDL.summary <- data.table()
ARDL.formula <- list('ardl' = paste0('y ~ y1 + ipr + ipr1 + su + su1 ',
                                     '+ sr + sr1 + Dfincris + D2000s - 1'),
                     'static' = 'y ~ ipr + su + sr + Dfincris + D2000s - 1',
                     'ar1' = 'y ~ y1 + Dfincris + D2000s - 1',
                     'li' = 'y ~ ipr1 + su1 + sr1 + Dfincris + D2000s - 1',
                     'dl' = paste0('y ~ ipr + ipr1 + su + su1 + sr + sr1 ',
                                   '+ Dfincris + D2000s - 1'))

## ARDL + Static + AR(1) + Leading indicator + Distributed lag
for (f in c('ardl', 'static', 'ar1', 'li', 'dl')) {
  gdp.ARDL[[f]] <- lm(ARDL.formula[[f]], data = us.gdp[2:104])
  ardl.summary <- summary(gdp.ARDL[[f]])
  ARDL.summary <-
    rbindlist(list(ARDL.summary,
                   data.table(t(c(ardl.summary$fstatistic[1],
                                  AIC(gdp.ARDL[[f]]), BIC(gdp.ARDL[[f]]),
                                  ardl.summary$r.squared,
                                  ardl.summary$adj.r.squared,
                                  9 - ardl.summary$fstatistic[2])))))
}
setnames(ARDL.summary, c('F-stat', 'AIC', 'Schwarz',
                         'R2', 'adj R2', '# of restr.'))
ARDL.summary <- cbind('Model' = c('ARDL', 'Static Regression',
                                  'AR(1) model', 'Leading indicator model',
                                  'Distributed lag model'), ARDL.summary)

### 3.7.3 Default Risk ###
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

# ARDL table
stargazer(oas.fit[[1]], oas.fit[[2]], oas.fit[[3]], oas.fit[[4]])  