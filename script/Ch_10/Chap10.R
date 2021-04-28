rm(list = ls())
setwd('E:/Book/Ch_10')

library(data.table)
library(ggplot2)
library(MSwM)

### 10.5 Simulated Data ###
sim.data <- fread('simulated_ch10_sec5.csv')
sim.data[, x := 1:600]
ggplot(data = sim.data, aes(x = x, y = y)) + theme_classic() +
  geom_line() + xlab('') + ylab('')

mod.sim <- list()

# Baseline - AR(2)
mod.sim$ar <- arima(sim.data[101:600, y], order = c(2, 0, 0))

## Markov regime switching
mod.sim$lm <- lm(y ~ y.L1 + y.L2,
                 data = data.table('y' = sim.data[103:600, y],
                                   'y.L1' = sim.data[102:599, y],
                                   'y.L2' = sim.data[101:598, y]))

# Regime change in the intercept and the first AR lag
mod.sim$mswm <- msmFit(mod.sim$lm, k = 2, sw = c(T, T, F, T),
                       control = list(parallel = F))

## Forecasting
mod.sim$mswm@transMat %*% as.matrix(mod.sim$mswm@Coef[, 1:2]) %*%
  matrix(c(1, sim.data[599, y]), nrow = 2) +
  mod.sim$mswm@Coef[1, 3] * sim.data[598, y]

### 10.6 US Industrial Production ###
ind.prod <- fread('indprod_ch9.csv')
prod <- list()
prod$est <- ind.prod[date >= '1930-01-01' & date <= '2000-01-01', indpro]
prod$frc <- ind.prod[date >= '2000-04-01' & date <= '2014-10-01', indpro]

mod.prod <- list()

# Baseline - AR(2)
mod.prod$ar <- arima(prod$est, order = c(2, 0, 0))

## Markov regime switching
prod$ms <- data.table('y' = prod$est[3:281],
                      'y.L1' = prod$est[2:280], 'y.L2' = prod$est[1:279])
mod.prod$lm <- lm(y ~ y.L1 + y.L2, data = prod$ms)

# Regime change in the intercept and the first AR lag
mod.prod$mswm <- msmFit(mod.prod$lm, k = 2,
                        sw = c(T, T, F, T), control = list(parallel = F))

## Forecasting
mod.prod$mswm@transMat %*% as.matrix(mod.prod$mswm@Coef[, 1:2]) %*%
  matrix(c(1, prod$ms[279, y.L1]), nrow = 2) +
  mod.prod$mswm@Coef[1, 3] * prod$ms[279, y.L2]