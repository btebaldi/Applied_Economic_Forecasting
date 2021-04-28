rm(list = ls())
setwd('E:/Book/Ch_7') # replace with own directory

library(data.table)
library(ggplot2)
library(vars)

### 7.9 Simulated Data ###
sim.data <- fread('simulated_cointegration.csv')
sim.data <- sim.data[, .(x, y)]

## Plots
ggplot(sim.data, aes(x = 101:600, y = x)) + geom_line() + theme_bw() +
  xlab('') + ylab('') + ggtitle('X')
ggplot(sim.data, aes(x = 101:600, y = y)) + geom_line() + theme_bw() +
  xlab('') + ylab('') + ggtitle('X')

## Correlogram
sim.acf <- sim.pacf <- list()
sim.acf$x <- acf(sim.data[, x], lag.max = 20)
sim.acf$y <- acf(sim.data[, y], lag.max = 20)
sim.pacf$x <- pacf(sim.data[, x], lag.max = 20)
sim.pacf$y <- pacf(sim.data[, y], lag.max = 20)

## ADF Unit Root Test
sim.df <- list()
sim.df$y <- ur.df(sim.data[, y], type = 'trend', lags = 1)
sim.df$dy <- ur.df(diff(sim.data[, y]), lags = 1)
summary(sim.df$y)
summary(sim.df$dy)

## Forecasting
sim.est <- sim.fore <- list()

# VAR(2) in levels
sim.est$var <- VAR(sim.data[1:400], p = 2, type = 'both')
sim.fore$var <- predict(sim.est$var, n.ahead = 100)
sim.fore$var <- data.table('x' = sim.fore$var$fcst$x[, 1],
                           'y' = sim.fore$var$fcst$y[, 1])

# VAR(1) in differences
sim.est$dvar <- VAR(apply(sim.data[1:400], 2, diff), p = 1)
sim.fore$dvar <- predict(sim.est$dvar, n.ahead = 100)
sim.fore$dvar <- data.table('x' = sim.fore$dvar$fcst$x[, 1],
                            'y' = sim.fore$dvar$fcst$y[, 1])
sim.fore$dvar <- cumsum(rbind(sim.data[400], sim.fore$dvar))
sim.fore$dvar <- sim.fore$dvar[2:101]

# VECM
sim.est$vec <- ca.jo(sim.data[1:400], type = 'eigen', ecdet = 'const',
                     K = 2, spec = 'longrun')
sim.fore$vec <- predict(cajools(sim.est$vec), n.ahead = 100)
sim.fore$vec <- sim.data[401:500] + cumsum(sim.fore$vec)

sim.plot <- data.table('X' = sim.data[401:500, x],
                       'X_VAR_LEVEL' = sim.fore$var[, x],
                       'X_VAR_DIFF' = sim.fore$var[, x],
                       'X_VECM' = sim.fore$vec[, x],
                       'Y' = sim.data[401:500, y],
                       'Y_VAR_LEVEL' = sim.fore$var[, y],
                       'Y_VAR_DIFF' = sim.fore$dvar[, y],
                       'Y_VECM' = sim.fore$vec[, y])

## X
# VAR in levels
xhat.plot <- data.table(rbindlist(list(data.table(sim.plot[, X]),
                                       data.table(sim.plot[, X_VAR_LEVEL]))),
                        rep(c('X', 'X_VAR_LEVEL'), each = 100))

# VAR in differences
xhat.plot <- data.table(rbindlist(list(data.table(sim.plot[, X]),
                                       data.table(sim.plot[, X_VAR_DIFF]))),
                        rep(c('X', 'X_VAR_DIFF'), each = 100))

# VECM
xhat.plot <- data.table(rbindlist(list(data.table(sim.plot[, X]),
                                       data.table(sim.plot[, X_VECM]))),
                        rep(c('X', 'X_VAR_VECM'), each = 100))

setnames(xhat.plot, c('xhat', 'label'))
ggplot(xhat.plot, aes(x = rep(501:600, 2), y = xhat, linetype = label)) +
  geom_line() + xlab('') + ylab('') + theme(legend.title = element_blank())

## Y
# VAR in levels
yhat.plot <- data.table(rbindlist(list(data.table(sim.plot[, Y]),
                                       data.table(sim.plot[, Y_VAR_LEVEL]))),
                        rep(c('Y', 'Y_VAR_LEVEL'), each = 100))

# VAR in differences
yhat.plot <- data.table(rbindlist(list(data.table(sim.plot[, Y]),
                                       data.table(sim.plot[, Y_VAR_DIFF]))),
                        rep(c('Y', 'Y_VAR_DIFF'), each = 100))

# VECM
yhat.plot <- data.table(rbindlist(list(data.table(sim.plot[, Y]),
                                       data.table(sim.plot[, Y_VECM]))),
                        rep(c('Y', 'Y_VAR_VECM'), each = 100))

setnames(yhat.plot, c('yhat', 'label'))
ggplot(yhat.plot, aes(x = rep(501:600, 2), y = yhat, linetype = label)) +
  geom_line() + xlab('') + ylab('') + theme(legend.title = element_blank())

### 7.10 UK Term Structure ###
rate.data <- fread('dataonly_uktermstructure.csv')

plot.label <- c('1985', '1990', '1995', '2000', '2005', '2010')
for (r in c('r1', 'r3', 'r6', 'r24')) {
  plot(rate.data[[r]], type = 'n', xlab = '', ylab = '',
       xaxt = 'n', ylim = c(0, 16))
  axis(1, at = c(0:5) * 66, labels = plot.label)
  lines(rate.data[[r]], lty = 1)
}

rate.acf <- list()
for (r in c('r1', 'r3', 'r6', 'r24')) {
  rate.acf[[r]] <- list('ACF' = acf(rate.data[[r]], lag.max = 20),
                        'PACF' = pacf(rate.data[[r]], lag.max = 20))
}

VARselect(rate.data[1:216, .(r1, r3, r6, r24)], lag.max = 10)

# VAR(4)
rate.est <- list()
rate.est$var <- VAR(rate.data[1:216, .(r1, r3, r6, r24)], p = 4)
rate.est$resid <- data.table(resid(rate.est$var))

for (r in c('r1', 'r3', 'r6', 'r24')) {
  plot(rate.est$resid[[r]], type = 'n', xlab = '', ylab = '', xaxt = 'n')
  axis(1, at = c(0:5) * 66, labels = plot.label)
  lines(rate.est$resid[[r]], lty = 1)
  title(paste0(toupper(r), ' Residuals'))
}

# VAR(3) in differences
rate.est$dvar <- VAR(apply(rate.data[1:216, .(r1, r3, r6, r24)], 2, diff),
                     p = 3)

# VECM
rate.est$vec <- ca.jo(rate.data[1:216, .(r1, r3, r6, r24)], type = 'eigen',
                      ecdet = 'const', K = 3, spec = 'longrun')

## Forecasting
rate.fore <- list()

# VAR(4) in levels
rate.fore$var <- predict(rate.est$var, n.ahead = 120)
rate.fore$var <- data.table('r1' = rate.fore$var$fcst$r1[, 1],
                            'r3' = rate.fore$var$fcst$r3[, 1],
                            'r6' = rate.fore$var$fcst$r6[, 1],
                            'r24' = rate.fore$var$fcst$r24[, 1])

# VAR(3) in differences
rate.fore$dvar <- predict(rate.est$dvar, n.ahead = 120)
rate.fore$dvar <- data.table('r1' = rate.fore$dvar$fcst$r1[, 1],
                             'r3' = rate.fore$dvar$fcst$r3[, 1],
                             'r6' = rate.fore$dvar$fcst$r6[, 1],
                             'r24' = rate.fore$dvar$fcst$r24[, 1])
rate.fore$dvar <- cumsum(rbind(rate.data[216, .(r1, r3, r6, r24)],
                               rate.fore$dvar))
rate.fore$dvar <- rate.fore$dvar[2:121]

# VECM
rate.fore$vec <- predict(cajools(rate.est$vec), n.ahead = 120)
rate.fore$vec <- rate.data[217:336, .(r1, r3, r6, r24)] +
  cumsum(rate.fore$vec)

# Plots
for (r in c('r1', 'r3', 'r6', 'r24')) {
  rhat.plot <- data.table(rbindlist(list(data.table(rate.data[[r]][217:336]),
                                         data.table(rate.fore$var[[r]]),
                                         data.table(rate.fore$dvar[[r]]),
                                         data.table(rate.fore$vec[[r]]))))
  rhat.plot[, label := rep(c('r', 'r_VAR_level', 'r_VAR_diff', 'r_VECM'),
                           each = 120)]
  setnames(rhat.plot, c('rhat', 'label'))
  
  ggplot(rhat.plot[label == 'r' | label == 'r_VAR_level'],
         aes(x = rep(1:120, 2), y = rhat, linetype = label)) +
    geom_line() + xlab('') + ylab('') + theme(legend.title = element_blank())
  
  ggplot(rhat.plot[label == 'r' | label == 'r_VAR_diff'],
         aes(x = rep(1:120, 2), y = rhat, linetype = label)) +
    geom_line() + xlab('') + ylab('') + theme(legend.title = element_blank())
  
  ggplot(rhat.plot[label == 'r' | label == 'r_VECM'],
         aes(x = rep(1:120, 2), y = rhat, linetype = label)) +
    geom_line() + xlab('') + ylab('') + theme(legend.title = element_blank())
}

### 7.11 CCI & CLI ###
index.data <- fread('dataonly_uslei.csv')

plot.label <- c('86', '88', '90', '92', '94', '96', '98', '00', '02')
for (i in c('cci', 'cli')) {
  plot(index.data[[i]][313:540], type = 'n', xlab = '', ylab = '',
       xaxt = 'n', ylim = c(50, 100))
  axis(1, at = (0:8) * 24 + 18, labels = plot.label)
  lines(index.data[[i]][313:540], lty = 1)
  title(paste0(toupper(i), ' (2004 = 100)'))
}

# VAR(4)
index.est <- list()
index.est$var <- VAR(index.data[313:540, .(cci, cli)], p = 4)

# VAR(3) in differences
index.est$dvar <- VAR(apply(index.data[313:540, .(cci, cli)], 2, diff), p = 3)

# VECM
index.est$vec <- ca.jo(index.data[313:540, .(cci, cli)], type = 'eigen',
                       ecdet = 'const', K = 2, spec = 'longrun')

## Forecasting
index.fore <- list()

# VAR(4) in levels
index.fore$var <- predict(index.est$var, n.ahead = 36)
index.fore$var <- data.table(cbind(index.fore$var$fcst$cci[, 1],
                                   index.fore$var$fcst$cli[, 1]))
setnames(index.fore$var, c('cci', 'cli'))

# VAR(3) in differences
index.fore$dvar <- predict(index.est$dvar, n.ahead = 36)
index.fore$dvar <- data.table(cbind(index.fore$dvar$fcst$cci[, 1],
                                    index.fore$dvar$fcst$cli[, 1]))
setnames(index.fore$dvar, c('cci', 'cli'))
index.fore$dvar <- cumsum(rbind(index.data[540, .(cci, cli)],
                                index.fore$dvar))
index.fore$dvar <- index.fore$dvar[2:37]

# VECM
index.fore$vec <- predict(cajools(index.est$vec), n.ahead = 36)
index.fore$vec <- index.data[541:576, .(cci, cli)] + cumsum(index.fore$vec)

# Plots
for (i in c('cci', 'cli')) {
  ihat.plot <- data.table(rbindlist(list(data.table(index.data[[i]][541:576]),
                                         data.table(index.fore$var[[i]]),
                                         data.table(index.fore$dvar[[i]]),
                                         data.table(index.fore$vec[[i]]))))
  ihat.plot[, label := rep(c('i', 'i_VAR_level', 'i_VAR_diff', 'i_VECM'),
                           each = 36)]
  setnames(ihat.plot, c('ihat', 'label'))
  
  ggplot(ihat.plot[label == 'i' | label == 'i_VAR_level'],
         aes(x = rep(1:36, 2), y = ihat, linetype = label)) +
    geom_line() + xlab('') + ylab('') + theme(legend.title = element_blank())
  
  ggplot(ihat.plot[label == 'i' | label == 'i_VAR_diff'],
         aes(x = rep(1:36, 2), y = ihat, linetype = label)) +
    geom_line() + xlab('') + ylab('') + theme(legend.title = element_blank())
  
  ggplot(ihat.plot[label == 'i' | label == 'i_VECM'],
         aes(x = rep(1:36, 2), y = ihat, linetype = label)) +
    geom_line() + xlab('') + ylab('') + theme(legend.title = element_blank())
}


## Rolling window
index.roll <- list('var' = data.table(),
                   'dvar' = data.table(),
                   'vec' = data.table())
for (j in 313:516) {
  pred.var <- predict(VAR(index.data[j:(j + 59), .(cci, cli)],
                          p = 4), n.ahead = 1)
  pred.var <- data.table('cci' = pred.var$fcst$cci[, 1],
                         'cli' = pred.var$fcst$cli[, 1])
  index.roll$var <- rbindlist(list(index.roll$var, pred.var))
  
  pred.dvar <- predict(VAR(apply(index.data[j:(j + 59), .(cci, cli)],
                                 2, diff), p = 3), n.ahead = 1)
  pred.dvar <- data.table('cci' = pred.dvar$fcst$cci[, 1],
                          'cli' = pred.dvar$fcst$cli[, 1])
  index.roll$dvar <- rbindlist(list(index.roll$dvar, pred.dvar))
  
  pred.vec <- ca.jo(index.data[j:(j + 59), .(cci, cli)], type = 'eigen',
                    ecdet = 'const', K = 2, spec = 'longrun')
  pred.vec <- predict(cajools(pred.vec), n.ahead = 1)
  pred.vec <- data.table(t(pred.vec[58, ]))
  setnames(pred.vec, c('cci', 'cli'))
  index.roll$vec <- rbindlist(list(index.roll$vec, pred.vec))
}

## Crisis period
# VAR(3)
index.est <- list()
index.est$var <- VAR(index.data[313:576, .(cci, cli)], p = 3)

# VAR(2) in differences
index.est$dvar <- VAR(apply(index.data[313:576, .(cci, cli)], 2, diff), p = 2)

# VECM
index.est$vec <- ca.jo(index.data[313:576, .(cci, cli)], type = 'eigen',
                       ecdet = 'const', K = 2, spec = 'longrun')

## Forecasting
index.fore <- list()

# VAR(3) in levels
index.fore$var <- predict(index.est$var, n.ahead = 76)
index.fore$var <- data.table(cbind(index.fore$var$fcst$cci[, 1],
                                   index.fore$var$fcst$cli[, 1]))
setnames(index.fore$var, c('cci', 'cli'))

# VAR(2) in differences
index.fore$dvar <- predict(index.est$dvar, n.ahead = 76)
index.fore$dvar <- data.table(cbind(index.fore$dvar$fcst$cci[, 1],
                                    index.fore$dvar$fcst$cli[, 1]))
setnames(index.fore$dvar, c('cci', 'cli'))
index.fore$dvar <- cumsum(rbind(index.data[576, .(cci, cli)],
                                index.fore$dvar))
index.fore$dvar <- index.fore$dvar[2:77]

# VECM
index.fore$vec <- predict(cajools(index.est$vec), n.ahead = 76)
index.fore$vec <- index.data[577:652, .(cci, cli)] + cumsum(index.fore$vec)

# Plots
for (i in c('cci', 'cli')) {
  ihat.plot <- data.table(rbindlist(list(data.table(index.data[[i]][577:652]),
                                         data.table(index.fore$var[[i]]),
                                         data.table(index.fore$dvar[[i]]),
                                         data.table(index.fore$vec[[i]]))))
  ihat.plot[, label := rep(c('i', 'i_VAR_level', 'i_VAR_diff', 'i_VECM'),
                           each = 76)]
  setnames(ihat.plot, c('ihat', 'label'))
  
  ggplot(ihat.plot[label == 'i' | label == 'i_VAR_level'],
         aes(x = rep(1:76, 2), y = ihat, linetype = label)) +
    geom_line() + xlab('') + ylab('') + theme(legend.title = element_blank())
  
  ggplot(ihat.plot[label == 'i' | label == 'i_VAR_diff'],
         aes(x = rep(1:76, 2), y = ihat, linetype = label)) +
    geom_line() + xlab('') + ylab('') + theme(legend.title = element_blank())
  
  ggplot(ihat.plot[label == 'i' | label == 'i_VECM'],
         aes(x = rep(1:76, 2), y = ihat, linetype = label)) +
    geom_line() + xlab('') + ylab('') + theme(legend.title = element_blank())
}


## Rolling window
index.roll <- list('var' = data.table(),
                   'dvar' = data.table(),
                   'vec' = data.table())
for (j in 313:592) {
  pred.var <- predict(VAR(index.data[j:(j + 59), .(cci, cli)],
                          p = 3), n.ahead = 1)
  pred.var <- data.table('cci' = pred.var$fcst$cci[, 1],
                         'cli' = pred.var$fcst$cli[, 1])
  index.roll$var <- rbindlist(list(index.roll$var, pred.var))
  
  pred.dvar <- predict(VAR(apply(index.data[j:(j + 59), .(cci, cli)],
                                 2, diff), p = 2), n.ahead = 1)
  pred.dvar <- data.table('cci' = pred.dvar$fcst$cci[, 1],
                          'cli' = pred.dvar$fcst$cli[, 1])
  index.roll$dvar <- rbindlist(list(index.roll$dvar, pred.dvar))
  
  pred.vec <- ca.jo(index.data[j:(j + 59), .(cci, cli)], type = 'eigen',
                    ecdet = 'const', K = 2, spec = 'longrun')
  pred.vec <- predict(cajools(pred.vec), n.ahead = 1)
  pred.vec <- data.table(t(pred.vec[58, ]))
  setnames(pred.vec, c('cci', 'cli'))
  index.roll$vec <- rbindlist(list(index.roll$vec, pred.vec))
}