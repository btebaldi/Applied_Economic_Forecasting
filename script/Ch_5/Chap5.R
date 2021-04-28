rm(list = ls())
setwd('E:/Book/Ch_5')

library(data.table)
library(FinTS)
library(ggplot2)
library(lmtest)
library(tseries)
library(xtable)

### 5.13 Simulated Data ###
arma.sim <- list()

# AR(2)
arma.sim$y1 <- arima.sim(model = list(ar = c(0.25, 0.25)), n = 600)

# MA(2)
arma.sim$y2 <- arima.sim(model = list(ma = c(0.25, -0.4)), n = 600)

# ARMA(2, 2)
arma.sim$y3 <- arima.sim(model = list(ar = c(-0.2, -0.1),
                                      ma = c(0.6, -0.15)), n = 600)

for (var in c('y1', 'y2', 'y3')) {
  ts.plot(arma.sim[[var]][101:600],
          gpars = list(main = toupper(var), ylab = ''))
}

# ACF & PACF
acf.sim <- pacf.sim <- list()
for (var in c('y1', 'y2', 'y3')) {
  acf.sim[[var]] <- acf(arma.sim[[var]][101:600], lag.max = 12)
  pacf.sim[[var]] <- pacf(arma.sim[[var]][101:600], lag.max = 12)
}

## ARMA selection
arma.spec <- list()
arma.spec$y1 <- arma.spec$y2 <- arma.spec$y3 <-
  list('AIC' = data.table(), 'BIC' = data.table())

# y1
for (ar.lag in 0:10) {
  arma.stat <- rep(0, 2)
  arma.fit <- arima(arma.sim$y1[101:600], order = c(ar.lag, 0, 0))
  arma.stat[1] <- arma.fit$aic
  arma.stat[2] <- -2 * arma.fit$loglik + ar.lag * log(499)
  
  arma.spec$y1$AIC <- rbindlist(list(arma.spec$y1$AIC,
                                     data.table(t(arma.stat[1]))))
  arma.spec$y1$BIC <- rbindlist(list(arma.spec$y1$BIC,
                                     data.table(t(arma.stat[2]))))
}

# y2 / y3
for (var in c('y2', 'y3')) {
  for (ar.lag in 0:3) {
    arma.stat <- rep(0, 8)
    for (ma.lag in 0:3) {
      arma.fit <- arima(arma.sim[[var]][101:600],
                        order = c(ar.lag, 0, ma.lag))
      arma.stat[ma.lag + 1] <- arma.fit$aic
      arma.stat[ma.lag + 5] <- -2 * arma.fit$loglik +
        (ar.lag + ma.lag) * log(499)
    }
    arma.spec[[var]]$AIC <- rbindlist(list(arma.spec[[var]]$AIC,
                                       data.table(t(arma.stat[1:4]))))
    arma.spec[[var]]$BIC <- rbindlist(list(arma.spec[[var]]$BIC,
                                       data.table(t(arma.stat[5:8]))))
  }
}

## Estimation
arma.fit <- list()

arma.fit$y1 <- arima(arma.sim$y1[101:600], order = c(2, 0, 0))
arma.fit$y2 <- arima(arma.sim$y2[101:600], order = c(0, 0, 2))
arma.fit$y3 <- arima(arma.sim$y3[101:600], order = c(2, 0, 2))

## Diagnostic tests
arma.test <- list()
for (var in c('y1', 'y2', 'y3')) {
  eps <- arma.fit[[var]]$residuals
  eps1 <- diff(eps, 1)
  eps2 <- diff(eps, 2)
  
  # Breusch-Godfrey serial correlation LM test
  arma.test[[var]]$bg <- lm(eps[3:500] ~ eps1[2:499] + eps2)
  
  # ARCH test
  arma.test[[var]]$arch <- ArchTest(arma.fit[[var]]$residuals, lags = 1)
  
  # White test
  arma.test[[var]]$white <- white.test(arma.fit[[var]]$residuals)
  
  # Jarque-Bera test
  arma.test[[var]]$jb <- jarque.bera.test(arma.fit[[var]]$residuals)
}

## Residuals
arma.plot <- list()
for (var in c('y1', 'y2', 'y3')) {
  arma.fit[[var]]$yhat <- arma.sim[[var]] - arma.fit[[var]]$residuals
  arma.fit[[var]]$acf <- acf(arma.fit[[var]]$residuals, lag.max = 12)
  arma.fit[[var]]$pacf <- pacf(arma.fit[[var]]$residuals, lag.max = 12)
  
  arma.plot[[var]]$yhat <-
    data.table('y' = as.numeric(cbind(arma.sim[[var]][101:600],
                                      arma.fit[[var]]$yhat,
                                      arma.fit[[var]]$residuals)),
               'label' = rep(c('Actual', 'Fitted', 'Residual'), each = 500))
}

for (c in 1:3) {
  print(ggplot(arma.plot[[c]]$yhat,
               aes(x = rep(101:600, 3), y = y, linetype = label)) +
          geom_line() + xlab('') + ylab('') +
          theme(legend.title = element_blank()))
}

## Misspecification
arma.fit$y1m <- arima(arma.sim$y1[101:600], order = c(1, 0, 0))
arma.fit$y2m <- arima(arma.sim$y2[101:600], order = c(0, 0, 1))
arma.fit$y3m <- arima(arma.sim$y3[101:600], order = c(1, 0, 1))

for (var in c('y1m', 'y2m', 'y3m')) {
  arma.fit[[var]]$acf <- acf(arma.fit[[var]]$residuals, lag.max = 12)
  arma.fit[[var]]$pacf <- pacf(arma.fit[[var]]$residuals, lag.max = 12)
}

## Subsample analysis
arma.sub <- list()
for (var in c('y1', 'y2', 'y3')) {
  for (ar.lag in 0:10) {
    arma.stat <- rep(0, 22)
    for (ma.lag in 0:10) {
      arma.fit <- arima(arma.sim[[var]][101:200],
                        order = c(ar.lag, 0, ma.lag))
      arma.stat[ma.lag + 1] <- arma.fit$aic
      arma.stat[ma.lag + 12] <- -2 * arma.fit$loglik +
        (ar.lag + ma.lag) * log(99)
    }
    arma.sub[[var]]$AIC <- rbindlist(list(arma.sub[[var]]$AIC,
                                          data.table(t(arma.stat[1:11]))))
    arma.sub[[var]]$BIC <- rbindlist(list(arma.sub[[var]]$BIC,
                                          data.table(t(arma.stat[12:22]))))
  }
}

## Forecasting
arma.forecast <- list()

# Dynamic
arma.fit$y1f <- arima(arma.sim$y1[101:500], order = c(1, 0, 0))
arma.fit$y2f <- arima(arma.sim$y2[101:500], order = c(0, 0, 1))
arma.fit$y3f <- arima(arma.sim$y3[101:500], order = c(1, 0, 1))
for (var in c('y1f', 'y2f', 'y3f')) {
  arma.forecast[[var]]$dynamic <- predict(arma.fit[[var]], n.ahead = 100)$pred
}

# Static
arma.forecast$y1f$static <- arma.forecast$y2f$static <-
  arma.forecast$y3f$static <- rep(0, 100)
for (c in 1:100) {
  arma.fit$y1f <- arima(arma.sim$y1[101:(499 + c)], order = c(1, 0, 0))
  arma.forecast$y1f$static[c] <- predict(arma.fit$y1f, n.ahead = 1)$pred
  
  arma.fit$y2f <- arima(arma.sim$y2[101:(499 + c)], order = c(0, 0, 1))
  arma.forecast$y2f$static[c] <- predict(arma.fit$y2f, n.ahead = 1)$pred
  
  arma.fit$y3f <- arima(arma.sim$y1[101:(499 + c)], order = c(1, 0, 1))
  arma.forecast$y3f$static[c] <- predict(arma.fit$y3f, n.ahead = 1)$pred
}
  
### 5.14.1 US FFR ###
arima.ffr <- fread('arma_tbill.csv')

plot.label <- c('86', '88', '90', '92', '94', '96', '98', '00',
                '02', '04', '06', '08', '10', '12')
plot(arima.ffr[, r], type = 'n', main = 'R', xlab = '', ylab = '', xaxt = 'n')
axis(1, at = (1:14) * 24, labels = plot.label)
lines(arima.ffr[, r], type = 'l')

## ACF / PACF
acf.ffr <- list()
acf.ffr$r <- list('acf' = acf(arima.ffr[1:216, r], lag.max = 24),
                  'pacf' = pacf(arima.ffr[1:216, r], lag.max = 24))
acf.ffr$dr <- list('acf' = acf(arima.ffr[2:216, Dr], lag.max = 24),
                   'pacf' = pacf(arima.ffr[2:216, Dr], lag.max = 24))

## ADF test
adf.ffr <- list('r' = adf.test(arima.ffr[1:216, r], k = 1),
                'dr' = adf.test(arima.ffr[2:216, Dr], k = 1))

## ARIMA selection
ic.ffr <- list('AIC' = data.table(), 'BIC' = data.table())
for (ar.lag in 0:12) {
  arma.stat <- rep(0, 6)
  for (ma.lag in 0:2) {
    arma.fit <- arima(arima.ffr[2:216, Dr], order = c(ar.lag, 0, ma.lag))
    arma.stat[ma.lag + 1] <- arma.fit$aic
    arma.stat[ma.lag + 4] <- -2 * arma.fit$loglik +
      (ar.lag + ma.lag) * log(214)
  }
  ic.ffr$AIC <- rbindlist(list(ic.ffr$AIC, data.table(t(arma.stat[1:3]))))
  ic.ffr$BIC <- rbindlist(list(ic.ffr$BIC, data.table(t(arma.stat[4:6]))))
}
setnames(ic.ffr$AIC, c('MA0', 'MA1', 'MA2'))
ic.ffr$AIC[, AR := (0:12)]
setnames(ic.ffr$BIC, c('MA0', 'MA1', 'MA2'))
ic.ffr$BIC[, AR := (0:12)]

arima.fit <- list()

# ARMA(5, 2) model for Dr
arima.fit$dr <- arima(arima.ffr[2:216, Dr], order = c(5, 0, 2))
print(arima.fit$dr)
acf.ffr$resid <- list('acf' = acf(arima.fit$dr$residuals, lag.max = 24),
                      'pacf' = pacf(arima.fit$dr$residuals, lag.max = 24))

## Diagnostic tests
arima.ffr[2:216, 'eps' := as.numeric(arima.fit$dr$residuals)]
arima.ffr[3:216, 'eps1' := diff(arima.ffr[2:216, eps], 1)]
arima.ffr[4:216, 'eps2' := diff(arima.ffr[3:216, eps1], 1)]

# Breusch-Godfrey serial correlation LM test
summary(lm(eps ~ eps1 + eps2, data = arima.ffr[4:216]))

# ARCH test
ArchTest(arima.fit$dr$residuals, lags = 1)

# White test
white.test(arima.fit$dr$residuals)

# AR(11) for r
arima.fit$r <- arima(arima.ffr[1:216, r], order = c(11, 0, 0))
print(arima.fit$r)
acf.ffr$resid <- list('acf' = acf(arima.fit$dr$residuals, lag.max = 24),
                      'pacf' = pacf(arima.fit$dr$residuals, lag.max = 24))

## Diagnostic tests
arima.ffr[1:216, eps := as.numeric(arima.fit$r$residuals)]
arima.ffr[2:216, 'eps1' := diff(arima.ffr[1:216, eps], 1)]
arima.ffr[3:216, 'eps2' := diff(arima.ffr[2:216, eps1], 1)]

# Breusch-Godfrey serial correlation LM test
summary(lm(eps ~ eps1 + eps2, data = arima.ffr[3:216]))

# ARCH test
ArchTest(arima.fit$r$residuals, lags = 1)

# White test
white.test(arima.fit$r$residuals)

## 1-step ahead forecasts
forecast.ffr <- list()

# Differences
# Dynamic
forecast.ffr$dr$dynamic <- as.numeric(predict(arima.fit$dr,
                                              n.ahead = 48)$pred)
# Static
forecast.ffr$dr$static <- rep(0, 48)
for (c in 1:48) {
  ffr.fit <- arima(arima.ffr[2:(215 + c), Dr], order = c(5, 0, 2))
  forecast.ffr$dr$static[c] <- predict(ffr.fit, n.ahead = 1)$pred
}
forecast.ffr$dr$dynamic <- forecast.ffr$dr$dynamic + arima.ffr[216:263, r]
forecast.ffr$dr$static <- forecast.ffr$dr$static + arima.ffr[216:263, r]

# Levels
# Dynamic
forecast.ffr$r$dynamic <- as.numeric(predict(arima.fit$r, n.ahead = 48)$pred)
# Static
forecast.ffr$r$static <- rep(0, 48)
for (c in 1:48) {
  ffr.fit <- arima(arima.ffr[1:(215 + c), r], order = c(11, 0, 0))
  forecast.ffr$r$static[c] <- predict(ffr.fit, n.ahead = 1)$pred
}

# RMSE & MAE
forecast.ffr$rhat <- cbind(forecast.ffr$r$dynamic, forecast.ffr$dr$dynamic,
                           forecast.ffr$r$static, forecast.ffr$dr$static)
RMSE <- sqrt(colSums((forecast.ffr$rhat - arima.ffr[217:264, r])^2) / 48)
MAE <- colSums(abs(forecast.ffr$rhat - arima.ffr[217:264, r])) / 48
error.mat <- rbind(RMSE, MAE)
print(xtable(error.mat))

## Crisis period
arima.fit$drc <- arima(arima.ffr[2:272, Dr], order = c(5, 0, 2))
arima.fit$rc <- arima(arima.ffr[1:272, r], order = c(11, 0, 0))

# Differences
# Dynamic
forecast.ffr$drc$dynamic <- as.numeric(predict(arima.fit$drc,
                                               n.ahead = 64)$pred)
# Static
forecast.ffr$drc$static <- rep(0, 64)
for (c in 1:64) {
  ffr.fit <- arima(arima.ffr[2:(271 + c), Dr], order = c(5, 0, 2))
  forecast.ffr$drc$static[c] <- predict(ffr.fit, n.ahead = 1)$pred
}
forecast.ffr$drc$dynamic <- forecast.ffr$drc$dynamic + arima.ffr[273:336, r]
forecast.ffr$drc$static <- forecast.ffr$drc$static + arima.ffr[273:336, r]

# Levels
# Dynamic
forecast.ffr$rc$dynamic <- as.numeric(predict(arima.fit$rc,
                                              n.ahead = 64)$pred)
# Static
forecast.ffr$rc$static <- rep(0, 64)
for (c in 1:64) {
  ffr.fit <- arima(arima.ffr[1:(271 + c), r], order = c(11, 0, 0))
  forecast.ffr$rc$static[c] <- predict(ffr.fit, n.ahead = 1)$pred
}

# RMSE & MAE
forecast.ffr$rhatc <- cbind(forecast.ffr$rc$dynamic, forecast.ffr$drc$dynamic,
                            forecast.ffr$rc$static, forecast.ffr$drc$static)
RMSE <- sqrt(colSums((forecast.ffr$rhatc - arima.ffr[273:336, r])^2) / 64)
MAE <- colSums(abs(forecast.ffr$rhatc - arima.ffr[273:336, r])) / 64
error.mat <- rbind(RMSE, MAE)
print(xtable(error.mat))

### 5.14.2 U.S. Inventories ###
arima.inven <- fread('arma_inven.csv')

plot.label <- c('86', '88', '90', '92', '94', '96', '98', '00',
                '02', '04', '06', '08', '10', '12')
plot(arima.inven[, rcpi], type = 'n', main = 'R',
     xlab = '', ylab = '', xaxt = 'n')
axis(1, at = (1:14) * 8, labels = plot.label)
lines(arima.inven[, rcpi], type = 'l')

## ACF / PACF
acf.inven <- list('acf' = acf(arima.inven[, rcpi], lag.max = 15),
                  'pacf' = pacf(arima.inven[, rcpi], lag.max = 15))

# ADF test
adf.inven <- adf.test(arima.inven[, rcpi], k = 1)
print(adf.inven)

## ARIMA selection
ic.inven <- list('AIC' = data.table(), 'BIC' = data.table())
for (ar.lag in 0:12) {
  arma.stat <- rep(0, 6)
  for (ma.lag in 0:2) {
    arma.fit <- arima(arima.inven[1:72, rcpi], order = c(ar.lag, 0, ma.lag))
    arma.stat[ma.lag + 1] <- arma.fit$aic
    arma.stat[ma.lag + 4] <- -2 * arma.fit$loglik +
      (ar.lag + ma.lag) * log(71)
  }
  ic.inven$AIC <- rbindlist(list(ic.inven$AIC, data.table(t(arma.stat[1:3]))))
  ic.inven$BIC <- rbindlist(list(ic.inven$BIC, data.table(t(arma.stat[4:6]))))
}
setnames(ic.inven$AIC, c('MA0', 'MA1', 'MA2'))
ic.inven$AIC[, AR := 0:12]
setnames(ic.inven$BIC, c('MA0', 'MA1', 'MA2'))
ic.inven$BIC[, AR := (0:12)]

# ARMA(3, 2)
arima.fit$rcpi <- arima(arima.inven[1:72, rcpi], order = c(3, 0, 2))
print(arima.fit$rcpi)
acf.inven$resid <- list('acf' = acf(arima.fit$rcpi$residuals, lag.max = 12),
                        'pacf' = pacf(arima.fit$rcpi$residuals, lag.max = 12))

## Diagnostic tests
arima.inven[1:72, eps := as.numeric(arima.fit$rcpi$residuals)]
arima.inven[2:72, 'eps1' := diff(arima.ffr[1:72, eps], 1)]
arima.inven[3:72, 'eps2' := diff(arima.ffr[2:72, eps1], 1)]

# Breusch-Godfrey serial correlation LM test
summary(lm(eps ~ eps1 + eps2, data = arima.inven[3:72]))

# ARCH test
ArchTest(arima.fit$rcpi$residuals, lags = 1)

# White test
white.test(arima.fit$rcpi$residuals)

## Forecasts
forecast.inven <- list()

# 1-step ahead
forecast.inven$rcpi$dynamic <- as.numeric(predict(arima.fit$rcpi,
                                                  n.ahead = 16)$pred)
forecast.inven$rcpi$static <- rep(0, 16)
for (c in 1:16) {
  inven.fit <- arima(arima.inven[1:(71 + c), rcpi], order = c(3, 0, 2))
  forecast.inven$rcpi$static[c] <- predict(inven.fit, n.ahead = 1)$pred
}

# 2-step ahead
forecast.inven$rcpi$multi <- rep(0, 16)
for (c in 1:16) {
  inven.fit <- arima(arima.inven[1:(70 + c), rcpi], order = c(3, 0, 2))
  forecast.inven$rcpi$multi[c] <- predict(inven.fit, n.ahead = 2)$pred[2]
}

# RMSE & MAE
forecast.inven$rcpihat <- cbind(forecast.inven$rcpi$dynamic,
                                forecast.inven$rcpi$static,
                                forecast.inven$rcpi$multi)
RMSE <- sqrt(colSums((forecast.inven$rcpihat -
                        arima.inven[73:88, rcpi])^2) / 16)
MAE <- colSums(abs(forecast.inven$rcpihat - arima.inven[73:88, rcpi])) / 16
error.mat <- rbind(RMSE, MAE)
print(xtable(error.mat))

## Crisis period
arima.fit$rcpic <- arima(arima.inven[1:100, rcpi], order = c(7, 0, 0))

forecast.inven$rcpic$dynamic <- as.numeric(predict(arima.fit$rcpic,
                                                   n.ahead = 12)$pred)
forecast.inven$rcpic$static <- rep(0, 12)
for (c in 1:12) {
  inven.fit <- arima(arima.inven[1:(100 + c), rcpi], order = c(7, 0, 0))
  forecast.inven$rcpic$static[c] <- predict(inven.fit, n.ahead = 1)$pred
}

# 2-step ahead
forecast.inven$rcpic$multi <- rep(0, 12)
for (c in 1:12) {
  inven.fit <- arima(arima.inven[1:(100 + c), rcpi], order = c(7, 0, 0))
  forecast.inven$rcpic$multi[c] <- predict(inven.fit, n.ahead = 2)$pred[2]
}

# RMSE & MAE
forecast.inven$rcpichat <- cbind(forecast.inven$rcpic$dynamic,
                                forecast.inven$rcpic$static,
                                forecast.inven$rcpic$multi)
RMSE <- sqrt(colSums((forecast.inven$rcpichat -
                        arima.inven[101:112, rcpi])^2) / 12)
MAE <- colSums(abs(forecast.inven$rcpichat - arima.inven[101:112, rcpi])) / 12
error.mat <- rbind(RMSE, MAE)
print(xtable(error.mat))