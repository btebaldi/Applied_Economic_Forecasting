rm(list = ls())
setwd('E:/Book/Ch_6')

library(data.table)
library(ggplot2)
library(vars)

### 6.9 Simulated Data ###
sim.data <- fread('var_simulated_ch6_sec9.csv')

# Time series plot
ggplot(sim.data[101:500], aes(x = 101:500, y = y)) + geom_line() +
  theme_bw() + xlab('') + ylab('') + ggtitle('Y')
ggplot(sim.data[101:500], aes(x = 101:500, y = x)) + geom_line() +
  theme_bw() + xlab('') + ylab('') + ggtitle('X')

# VAR lag order
VARselect(sim.data, lag.max = 12)

# VAR(1)
summary(VAR(sim.data, p = 1))

# VAR(4)
summary(VAR(sim.data, p = 4))

# VAR(1) diagnostic tests
var.sim <- VAR(sim.data, p = 1)

var.test <- list()

# VAR(1) residual normality test
var.test$norm <- normality.test(var.sim)

# VAR(1) residual heteroskedasticity test
var.test$arch <- arch.test(var.sim)

# VAR(1) residual serial correlation test
var.test$serial <- serial.test(var.sim, lags.pt = 2)

# Static and dynamic forecasts
var.pred <- predict(var.sim, n.ahead = 100)
var.hat <- cbind(var.pred$fcst$y[, 1], var.pred$fcst$x[, 1])
RMSE <- sqrt(colSums((var.hat - sim.data[501:600])^2) / 100)
MAE <- colSums(abs(var.hat - sim.data[501:600])) / 100 
error.mat <- rbind(RMSE, MAE)

# ARMA fit
arma.sim <- arma.pred <- list()
arma.sim$y <- arima(sim.data[101:500, y], order = c(1, 0, 0))
arma.sim$x <- arima(sim.data[101:500, x], order = c(1, 0, 3))

arma.pred$y <- predict(arma.sim$y, n.ahead = 100)
arma.pred$x <- predict(arma.sim$x, n.ahead = 100)

arma.hat <- cbind(as.numeric(arma.pred$y$pred), as.numeric(arma.pred$x$pred))
RMSE <- sqrt(colSums((arma.hat - sim.data[501:600])^2) / 100)
MAE <- colSums(abs(arma.hat - sim.data[501:600])) / 100 
error.mat <- rbind(RMSE, MAE)

# IRF
irf.sim <- irf(var.sim, impulse = c('y', 'x'), response = c('y', 'x'))
plot(irf.sim)

### 6.10 GDP Growth in the Euro Area ###
eu.data <- fread('eu_growth_ch6_sec10.csv')
eu.data[, date := as.Date(date, format = '%Y-%m-%d')]

# Endogenous variables
eu.endo <- eu.data[, .(g_fr, g_ger, g_it)]

# Exogenous variables
eu.exo <- eu.data[, .(g_us, d_1988q1, d_1993q1, d_2000q1)]

# Time-series plot
eu.plot <- melt(eu.data[, .(date, g_fr, g_ger, g_it)], id = 'date')
ggplot(data = eu.plot, aes(x = date, y = value, linetype = variable)) +
  geom_line() + theme_bw() + theme(legend.title = element_blank()) +
  scale_x_date() + xlab('') + ylab('i (%)')  

# VAR lag order
VARselect(eu.endo[1:72, ], exogen = eu.exo[1:72, ], lag.max = 4)

eu.var <- list()
# VAR(1)
eu.var$l1 <- VAR(eu.endo[1:72, ], exogen = eu.exo[1:72, ], lag = 1)

eu.test <- list()
eu.test$norm <- normality.test(eu.var$l1) # normality test
eu.test$arch <- arch.test(eu.var$l1) # heteroskedasticity test

eu.irf <- irf(eu.var$l1, impulse = c('g_fr', 'g_ger', 'g_it'),
              response = c('g_fr', 'g_ger', 'g_it'))
plot(eu.irf)

# VAR(2)
eu.var$l2 <- VAR(eu.endo[1:72, ], exogen = eu.exo[1:72, ], lag = 2)
causality(eu.var$l2) # Granger causality test

# Static and dynamic forecasts
eu.pred <- list()
eu.pred$l1 <- predict(eu.var$l1, n.ahead = 12, dumvar = eu.exo[73:84, ])
eu.pred$l2 <- predict(eu.var$l2, n.ahead = 12, dumvar = eu.exo[73:84, ])

### 6.11 Monetary Transmission Mechanism ###
us.data <- fread('usmonetary_ch6_sec11.csv')
us.endo <- us.data[, .(ly, lp, i)]
us.exo <- us.data[, .(d0109, d0110)]

# VAR lag order
VARselect(us.endo[13:270], lag.max = 12)

us.var <- us.test <- us.irf <- list()

us.var$l2 <- VAR(us.endo[13:270], lag = 2) # VAR(2)
summary(us.var$l2)
us.irf$l2 <- irf(us.var$l2, impulse = c('i', 'lp', 'ly'),
                 response = c('i', 'lp', 'ly'))
plot(us.irf$l2)

us.var$l2D <- VAR(us.endo[13:270], exogen = us.exo[13:270], lag = 2) 
summary(us.var$l2D) # Dummy variables for 9/11

us.test$norm <- normality.test(us.var$l2D) # normality test
us.test$arch <- arch.test(us.var$l2D) # heteroskedasticity test

# Pre-crisis
us.endo <- us.data[, .(ly, lp, i, lpcm, lsmtr, lsmnbr)]
us.exo <- us.data[, .(d0109, d0110, d9101, d0111)]

us.var$l1 <- VAR(us.endo[13:270, ], exogen = us.exo[13:270, ], lag = 1)
summary(us.var$l1)
us.irf$l1 <- irf(us.var$l1, impulse = 'i',
                 response = c('lpcm', 'lp', 'ly', 'i', 'lsmtr', 'lsmnbr'))
plot(us.irf$l1)