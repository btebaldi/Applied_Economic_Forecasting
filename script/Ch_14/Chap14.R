rm(list = ls())
setwd('E:/Book/Ch_14')

library(data.table)
library(forecast)
library(ggplot2)
library(midasr)
library(rugarch)

daily.sp <- list()

### 14.6.1 ARCH-type Models ###
# S&P 500 Index
daily.sp$data <- fread('OxfordManRealizedVolatility.csv',
                       select = c('Date', 'Return'))
daily.sp$data[, Date := as.Date(as.character(Date), '%Y%m%d')]
ggplot(daily.sp$data, aes(Date)) + geom_line(aes(y = Return)) +
  theme_classic() + theme(plot.title = element_text(hjust = 0.5)) +
  xlab('') + ylab('') + ggtitle('S&P 500 Daily Returns')
ggplot(daily.sp$data, aes(Date)) + geom_line(aes(y = RV)) +
  theme_classic() + theme(plot.title = element_text(hjust = 0.5)) +
  xlab('') + ylab('') + ggtitle('S&P 500 Index Realized Volatility (5-min)')

T <- dim(daily.sp$data)[1]
tt <- dim(daily.rv$data[Date <= '2006-12-31'])[1]
model <- c('sGARCH', 'eGARCH', 'gjrGARCH')
freq <- c('Day', 'Week', 'Month', 'SixMonth')

## GARCH / EGARCH / GJR-GARCH
# Jan. 2000 to Dec. 2006
for (g in model) {
  garch.spec <- ugarchspec(variance.model = list(model = g,
                                                 garchOrder = c(1, 1)),
                           mean.model = list(armaOrder = c(0, 0),
                                             include.mean = F))
  daily.sp[[g]] <- ugarchfit(data = daily.sp$data[Date <= '2006-12-31',
                                                  Return], garch.spec)
}

# Jan. 2000 to Oct. 2008
for (g in model) {
  garch.spec <- ugarchspec(variance.model = list(model = g,
                                                 garchOrder = c(1, 1)),
                           mean.model = list(armaOrder = c(0, 0),
                                             include.mean = F))
  daily.sp[[g]] <- ugarchfit(data = daily.sp$data[Date <= '2008-10-31',
                                                  Return], garch.spec)
}

# Jan. 2000 to Dec. 2009
for (g in model) {
  garch.spec <- ugarchspec(variance.model = list(model = g,
                                                 garchOrder = c(1, 1)),
                           mean.model = list(armaOrder = c(0, 0),
                                             include.mean = F))
  daily.sp[[g]] <- ugarchfit(data = daily.sp$data[Date <= '2009-12-31',
                                                  Return], garch.spec)
}

# Forecast evaluation
garch.pred <- list('Day' = data.table(), 'Week' = data.table(),
                   'Month' = data.table(), 'SixMonth' = data.table())

ptm <- proc.time()
for (g in model) {
  g.spec <- ugarchspec(variance.model = list(model = g, garchOrder = c(1, 1)),
                       mean.model = list(armaOrder = c(0, 0),
                                         include.mean = F))
  
  # daily rolling window
  g.pred <- list('Day' = rep(0, T - tt), 'Week' = rep(0, T - tt),
                 'Month' = rep(0, T - tt), 'SixMonth' = rep(0, T - tt))
 
  for (t in 1:(T - tt)) {
    g.fit <- ugarchfit(data = daily.sp$data[t:(t + tt - 1), Return],
                       garch.spec)
    g.pred$Day[t] <-
      as.numeric(ugarchforecast(g.fit, n.ahead = 1)@forecast$sigmaFor)
   
    g.fit <- ugarchfit(data = daily.sp$data[t:(t + tt - 5), Return],
                       garch.spec)
    g.pred$Week[t] <-
      as.numeric(ugarchforecast(g.fit, n.ahead = 5)@forecast$sigmaFor)[5]
    
    garch.fit <- ugarchfit(data = daily.sp$data[t:(t + tt - 20), Return],
                           garch.spec)
    g.pred$Month[t] <-
      as.numeric(ugarchforecast(garch.fit,
                                n.ahead = 20)@forecast$sigmaFor)[20]
    
    garch.fit <- ugarchfit(data = daily.sp$data[t:(t + tt - 120), Return],
                           garch.spec)
    g.pred$SixMonth[t] <-
      as.numeric(ugarchforecast(garch.fit,
                                n.ahead = 120)@forecast$sigmaFor)[120]
  }
  
  for (h in freq) {
    g.pred[[h]] <- data.table('sigma.hat' = g.pred[[h]], 'GARCH' = g)
    garch.pred[[h]] <- rbindlist(list(garch.pred[[h]], g.pred[[h]]))
  }
}
proc.time() - ptm

sigma.pred <- rbindlist(list(garch.pred$Day, garch.pred$Week,
                             garch.pred$Month, garch.pred$SixMonth))
sigma.pred[, Frequency := rep(freq, each = t * 3)]
for (g in model) {
  for (h in freq) {
    sigma.mse <- (sigma.pred[GARCH == g & Frequency == h, sigma.hat]^2 -
                    daily.sp$data[(tt + 1):T, Return]^2)^2
    sigma.qlike <- log(sigma.pred[GARCH == g & Frequency == h, sigma.hat]^2) +
      daily.sp$data[(tt + 1):T, Return]^2 /
      sigma.pred[GARCH == g & Frequency == h, sigma.hat]^2
    sigma.pred[GARCH == g & Frequency == h,
               c('MSE', 'QLIKE') := list(sigma.mse, sigma.qlike)]
  }
}
write.csv(sigma.pred, 'sigma_hat_garch.csv', row.names = F)

# DM test
for (h in freq) {
  sigma.freq <- sigma.pred[Frequency == h]
  for (loss in c('MSE', 'QLIKE')) {
    print(dm.test(sigma.freq[GARCH == 'sGARCH'][[loss]],
                  sigma.freq[GARCH == 'eGARCH'][[loss]]))
    print(dm.test(sigma.freq[GARCH == 'sGARCH'][[loss]],
                  sigma.freq[GARCH == 'gjrGARCH'][[loss]]))
    print(dm.test(sigma.freq[GARCH == 'eGARCH'][[loss]],
                  sigma.freq[GARCH == 'gjrGARCH'][[loss]]))
  }
}

### 14.6.2 Realized Volatility ###
daily.rv <- list()
daily.rv$data <- fread('OxfordManRealizedVolatility.csv',
                       select = c('Date', 'RV'))
daily.rv$data[, c('Date', 'RV.D', 'RV.W', 'RV.M', 'RV.S') :=
                list(as.Date(as.character(Date), '%Y%m%d'),
                     c(0, RV[1:(T - 1)]),
                     c(rep(0, 5), rollapply(RV[1:(T - 1)], 5, sum)),
                     c(rep(0, 22), rollapply(RV[1:(T - 1)], 22, sum)),
                     c(rep(0, 125), rollapply(RV[1:(T - 1)], 125, sum)))]

# Pre-crisis sample
tt <- dim(daily.rv$data[Date <= '2006-12-31'])[1]

# Crisis sample
tt <- dim(daily.rv$data[Date <= '2008-10-31'])[1]

# Post-crisis sample
tt <- dim(daily.rv$data[Date <= '2009-12-31'])[1]


# MIDAS-RV
daily.rv$midas <- midas_r(RV ~ fmls(RV.D, 3, 1, nealmon),
                          start = list(RV.D = c(0, 0, 0)),
                          data = daily.rv$data[2:tt])

# HAR
daily.rv$har <- lm(RV ~ RV.D + RV.W + RV.M, data = daily.rv$data[22:tt])

summary(daily.rv$midas)
summary(daily.rv$har)

## Multiple horizons
daily.rv$MIDAS <- daily.rv$HAR <- list()

# One-day-ahead
daily.rv$MIDAS$Day <- forecast(daily.rv$midas,
                               newdata = daily.rv$data[tt:(T - 1)])$mean
daily.rv$HAR$Day <- as.numeric(predict(daily.rv$har,
                                       newdata = daily.rv$data[tt:(T - 1)]))

# One-week-ahead
daily.rv$midas <- midas_r(RV.W ~ fmls(RV.D, 3, 1, nealmon),
                          start = list(RV.D = rep(0, 3)),
                          data = daily.rv$data[6:tt])
daily.rv$MIDAS$Week <-
  forecast(daily.rv$midas, newdata = daily.rv$data[(tt - 4):(T - 5)])$mean
daily.rv$HAR$Week <-
  as.numeric(predict(daily.rv$har, newdata = daily.rv$data[(tt - 4):(T - 5)]))

# One-month-ahead
daily.rv$midas <- midas_r(RV.M ~ fmls(RV.D, 3, 1, nealmon),
                          start = list(RV.D = rep(0, 3)),
                          data = daily.rv$data[23:tt])
daily.rv$MIDAS$Month <-
  forecast(daily.rv$midas, newdata = daily.rv$data[(tt - 21):(T - 22)])$mean
daily.rv$HAR$Month <-
  as.numeric(predict(daily.rv$har,
                     newdata = daily.rv$data[(tt - 21):(T - 22)]))

# Six-month-ahead
daily.rv$midas <- midas_r(RV.S ~ fmls(RV.D, 3, 1, nealmon),
                          start = list(RV.D = rep(0, 3)),
                          data = daily.rv$data[126:tt])
daily.rv$MIDAS$SixMonth <-
  forecast(daily.rv$midas, newdata = daily.rv$data[(tt - 124):(T - 125)])$mean
daily.rv$HAR$SixMonth <-
  as.numeric(predict(daily.rv$har,
                     newdata = daily.rv$data[(tt - 124):(T - 125)]))

## DM test
daily.rv$oos <- daily.rv$data[(tt + 1):T, RV]
for (h in freq) {
  daily.rv$MSE[[h]] <- data.table('MIDAS' = (daily.rv$MIDAS[[h]] -
                                               daily.rv$oos)^2,
                                  'HAR' = (daily.rv$HAR[[h]] -
                                             daily.rv$oos)^2)
  daily.rv$QLIKE[[h]] <- data.table('MIDAS' = log(daily.rv$MIDAS[[h]]) +
                                      daily.rv$oos / daily.rv$MIDAS[[h]],
                                    'HAR' = log(daily.rv$HAR[[h]]) +
                                      daily.rv$oos / daily.rv$HAR[[h]])
  print(dm.test(daily.rv$MSE[[h]]$MIDAS, daily.rv$MSE[[h]]$HAR))
  print(dm.test(daily.rv$QLIKE[[h]]$MIDAS, daily.rv$QLIKE[[h]]$HAR))
}