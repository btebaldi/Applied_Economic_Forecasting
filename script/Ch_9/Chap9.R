rm(list = ls())
setwd('E:/Book/Ch_9')

library(data.table)
library(forecast)
library(ggplot2)
library(Metrics)
library(tsDyn)

### 9.7 Simulated Data ###

y <- list('TAR' = rep(0, 601), 'STAR' = rep(0, 601))
e <- rnorm(601)

# y0 = 1.5 or y0 = 0.5
y$TAR[1] <- y$STAR[1] <- 1.5

for (t in 2:601) {
  # TAR DGP
  y$TAR[t] <- (0.5 - 0.3 * y$TAR[t - 1]) * (1 - (y$TAR[t - 1] > 1)) +
    (0.1 + 0.8 * y$TAR[t - 1]) * (y$TAR[t - 1] > 1)
  
  # STAR DGP
  y$STAR[t] <- 2 + 0.9 * y$STAR[t - 1] +
    (3 - 0.8 * y$STAR[t - 1]) / (1 + exp(-5 * (y$STAR[t - 1] - 1)))
}
y$TAR <- y$TAR[102:601] + e[102:601]
y$STAR <- y$STAR[102:601] + e[102:601]

# Plots
plot(y$TAR, type = 'l', xlab = '', ylab = '')
plot(y$STAR, type = 'l', xlab = '', ylab = '')

# Linearity tests
delta.lin.test(y$TAR)
delta.lin.test(y$STAR)

## Estimation
mod.test <- list()

# TAR
mod.test$tar <- setar(y$TAR[1:400], m = 1, d = 1)

# STAR
mod.test$star <- star(y$STAR[1:400], m = 1)

## Forecasting
frc.test <- lapply(mod.test, predict, n.ahead = 100)

### 9.8 US Industrial Production Growth ###
ind.prod <- fread('indprod_ch9.csv')
prod <- list()
prod$est <- ind.prod[date >= '1930-01-01' & date <= '2000-01-01', indpro]
prod$frc <- ind.prod[date >= '2000-04-01' & date <= '2014-10-01', indpro]

# Parameter estimates
mod.prod <- list()
mod.prod$ar <- arima(prod$est, order = c(2, 0, 0))
mod.prod$tar <- setar(prod$est, m = 2, d = 1)
mod.prod$star <- star(prod$est, m = 2)

# Recursive 1-step ahead
frc.prod <- list('ar' = rep(0, 59), 'tar' = rep(0, 59), 'star' = rep(0, 59))

# AR(2)
prod$est <- ind.prod[date >= '1930-01-01' & date <= '2000-01-01', indpro]
for (t in 1:59) {
  mod.prod$ar <- arima(prod$est, order = c(2, 0, 0))
  frc.prod$ar[t] <- as.numeric(predict(mod.prod$ar, n.ahead = 1)$pred)
  prod$est <- c(prod$est, frc.prod$ar[t])
}
frc.prod$ar <- predict(mod.prod$ar, n.ahead = 59)
  
# TAR
prod$est <- ind.prod[date >= '1930-01-01' & date <= '2000-01-01', indpro]
for (t in 1:59) {
  mod.prod$tar <- setar(prod$est, m = 2, d = 1, thDelay = 1)
  frc.prod$tar[t] <- as.numeric(predict(mod.prod$tar, n.ahead = 1))
  prod$est <- c(prod$est, frc.prod$tar[t])
}
frc.prod$tar <- predict(mod.prod$tar, n.ahead = 59)
  
# STAR
prod$est <- ind.prod[date >= '1930-01-01' & date <= '2000-01-01', indpro]
for (t in 1:59) {
  mod.prod$star <- star(prod$est, m = 2, thDelay = 1)
  frc.prod$star[t] <- as.numeric(predict(mod.prod$star, n.ahead = 1))
  prod$est <- c(prod$est, frc.prod$star[t])
}
frc.prod$star <- predict(mod.prod$star, n.ahead = 59)

# RMSE
rmse.prod <- lapply(frc.prod, rmse, actual = prod$frc)

# Diebold-Mariano test
error.prod <- lapply(frc.prod, `-`, prod$frc)
dm.test(error.prod$ar, error.prod$tar) # AR vs TAR
dm.test(error.prod$star, error.prod$ar) # STAR vs AR
dm.test(error.prod$star, error.prod$tar) # STAR vs TAR

# Plot
frc.plot <- data.table('Y' = prod$frc, 'tag' = 'Y')
for (model in c('ar', 'tar', 'star')) {
  frc.plot <- rbindlist(list(frc.plot,
                             data.table('Y' = frc.prod[[model]],
                                        'tag' = paste0('Y_', toupper(model)))))
}
frc.plot[tag == 'Y_AR', tag := 'Y_LINEAR']

ggplot(frc.plot, aes(x = rep(1:59, 4), y = Y, linetype = tag)) +
  geom_line() + xlab('') + ylab('') + theme(legend.title = element_blank())