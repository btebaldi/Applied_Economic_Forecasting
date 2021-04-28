rm(list = ls())
setwd('E:/Book/Ch_11')

library(data.table)
library(dlm)
library(ggplot2)
library(Metrics)
library(psych)

### 11.4 TVP Regression ###
tvp.sim <- list()

tvp.sim$data <- fread('tvp_simulated_ch11_sec4.csv')
tvp.sim$data[, date := as.Date(date, '%Y-%m-%d')]

ggplot(data = tvp.sim$data, aes(x = date, y = x1)) + theme_classic() +
  theme(plot.title = element_text(hjust = 0.5)) + geom_line() +
  scale_x_date() + xlab('') + ylab('') + ggtitle('X1')

ggplot(data = tvp.sim$data, aes(x = date, y = b1_0)) + theme_classic() +
  theme(plot.title = element_text(hjust = 0.5)) + geom_line() +
  scale_x_date() + xlab('') + ylab('') + ggtitle('B1')

ggplot(data = tvp.sim$data, aes(x = date, y = y_all)) + theme_classic() +
  theme(plot.title = element_text(hjust = 0.5)) + geom_line() +
  scale_x_date() + xlab('') + ylab('') + ggtitle('Y')

tvp.sim$mod <- dlmModReg(tvp.sim$data[, .(x1, x2, x3, x4)])
tvp.sim$filt <- dlmFilter(tvp.sim$data[, y_all], mod = tvp.sim$mod)
tvp.sim$rmse <- list()
for (t in 1:357) {
  tvp.sim$rmse[[t]] <- sqrt(diag(tvp.sim$filt$U.C[[t]] %*%
    diag(tvp.sim$filt$D.C[t, ]^2) %*% t(tvp.sim$filt$U.C[[t]])))
}
tvp.sim$rmse <- matrix(unlist(tvp.sim$rmse), ncol = 5, byrow = T)
tvp.sim$beta <- data.table(tvp.sim$filt$m)
setnames(tvp.sim$beta, c('B0', 'B1', 'B2', 'B3', 'B4'))

for (b in 1:4) {
  plot(tvp.sim$beta[[b + 1]], type = 'l', xlab = '', ylab = '',
       main = paste0('Smoothed SV', b, ' State Estimate'))
}

lm(tvp.sim$data[, b1_0] ~ tvp.sim$beta[2:357, B1])

## OLS vs. TVP
tvp.sim$ols <- predict(lm(y_all ~ x1 + x2 + x3 + x4,
                          data = tvp.sim$data[1:327]),
                       newdata = tvp.sim$data[328:356, .(x1, x2, x3, x4)])

tvp.sim$pred <- data.table('Date' = tvp.sim$data[328:356, date],
                           'Y' = c(tvp.sim$data[328:356, y_all],
                                   tvp.sim$ols, tvp.sim$filt$f[328:356]),
                           'Tag' = c(rep('Y', 29), rep('YF_OLS_CHW', 29),
                                     rep('YF_TVP_CHW', 29)))
ggplot(data = tvp.sim$pred, aes(x = Date, y = Y, group = Tag, colour = Tag)) +
  theme_classic() + theme(plot.title = element_text(hjust = 0.5)) + 
  geom_line() + scale_x_date() + xlab('') + ylab('')

### 11.5.1 US GDP Growth ###
us.gdp <- list()

us.gdp$data <- fread('gdp_us_ch11_sec5.csv')
us.gdp$ols <- arima(us.gdp$data[, g_us], order = c(2, 0, 0))

us.gdp$Y <- data.table('Date' = as.Date(us.gdp$data[3:108, date],
                                        '%m/%d/%Y'),
                       'y' = us.gdp$data[3:108, g_us],
                       'y.L1' = us.gdp$data[2:107, g_us],
                       'y.L2' = us.gdp$data[1:106, g_us])
us.gdp$mod <- dlmModReg(us.gdp$Y[, .(y.L1, y.L2)])
us.gdp$filt <- dlmFilter(us.gdp$Y[, y], mod = us.gdp$mod)

us.gdp$pred <- data.table('yhat' = rep(0, 192))

## Recursive + Rolling / OLS + TVP
for (t in 1:48) {
  # Recursive OLS
  us.gdp$ols <- arima(us.gdp$data[1:(59 + t), g_us], order = c(2, 0, 0))
  us.gdp$pred[t, yhat := as.numeric(predict(us.gdp$ols)$pred)]
  
  # Recursive TVP
  us.gdp$mod <- dlmModReg(us.gdp$Y[1:(57 + t), .(y.L1, y.L2)])
  us.gdp$filt <- dlmFilter(us.gdp$Y[1:(57 + t), y], mod = us.gdp$mod)
  us.gdp$pred[t + 48, yhat := us.gdp$filt$f[57 + t]]
  
  # Rolling OLS
  us.gdp$ols <- arima(us.gdp$data[t:(59 + t), g_us], order = c(2, 0, 0))
  us.gdp$pred[t + 96, yhat := as.numeric(predict(us.gdp$ols)$pred)]
  
  # Rolling TVP
  us.gdp$mod <- dlmModReg(us.gdp$Y[t:(57 + t), .(y.L1, y.L2)])
  us.gdp$filt <- dlmFilter(us.gdp$Y[t:(57 + t), y], mod = us.gdp$mod)
  us.gdp$pred[t + 144, yhat := us.gdp$filt$f[58]]
}
us.gdp$pred[, tag := c(rep('YF_OLS_CHW', 48), rep('YF_TVP_CHW', 48),
                       rep('YF_OLS_ROLW', 48), rep('YF_TVP_ROLW', 48))]
us.gdp$pred <- rbindlist(list(us.gdp$pred,
                              data.table('yhat' = us.gdp$data[61:108, g_us],
                                         'tag' = 'Y')))
us.gdp$pred[, Date := us.gdp$Y[59:106, Date]]

for (label in c('YF_OLS_CHW', 'YF_TVP_CHW', 'YF_OLS_ROLW', 'YF_TVP_ROLW')) {
  ggplot(data = us.gdp$pred[tag == 'Y' | tag == label],
         aes(x = Date, y = yhat, group = tag, colour = tag)) +
    theme_classic() + theme(plot.title = element_text(hjust = 0.5)) + 
    geom_line() + scale_x_date() + xlab('') + ylab('')
  print(rmse(us.gdp$pred[tag == 'Y', yhat], us.gdp$pred[tag == label, yhat]))
}

### 11.5.2 Dynamic Factor Models ###
g.gdp <- list()

g.gdp$data <- fread('gdp_ch11_sec5.csv')

# One factor
g.gdp$pc <- principal(g.gdp$data[, .(g_fr, g_ger, g_it, g_us)],
                      nfactors = 1, rotate = 'varimax')
g.gdp$data[, pc1 := g.gdp$pc$scores]

# Two factors
g.gdp$pc <- principal(g.gdp$data[, .(g_fr, g_ger, g_it, g_us)],
                      nfactors = 2, rotate = 'varimax')
g.gdp$data[, c('PC1', 'PC2', 'PC1.L1', 'PC2.L1') :=
             list(g.gdp$pc$scores[, 1], g.gdp$pc$scores[, 2],
                  c(0, g.gdp$pc$scores[1:107, 1]),
                  c(0, g.gdp$pc$scores[1:107, 2]))]

g.gdp$mod <- dlmModReg(g.gdp$data[, .(g_us, g_fr, g_ger, g_it, pc1)])
g.gdp$of <- dlmFilter(g.gdp$data[, g_us], mod = g.gdp$mod)

g.gdp$mod <- dlmModReg(g.gdp$data[2:108, .(g_us, g_fr, g_ger, g_it,
                                           PC1, PC2, PC1.L1, PC2.L1)])
g.gdp$tf <- dlmFilter(g.gdp$data[, g_us], mod = g.gdp$mod)

g.gdp$factors <- data.table('Factor' = c(g.gdp$data[, PC1],
                                         g.gdp$data[, PC2],
                                         g.gdp$tf$m[2:109, 2],
                                         g.gdp$tf$m[2:109, 4]),
                            'tag' = c(rep('TWOFACT_1', 108),
                                      rep('TWOFACT_2', 108),
                                      rep('F_SV1_TWOF', 108),
                                      rep('F_SV3_TWOF', 108)),
                            'Date' = as.Date(g.gdp$data[, date], '%Y-%m-%d'))
ggplot(data = g.gdp$factors, aes(x = Date, y = Factor,
                                 group = tag, colour = tag)) +
  theme_classic() + theme(plot.title = element_text(hjust = 0.5)) + 
  geom_line() + scale_x_date() + xlab('') + ylab('')

summary(lm(g.gdp$factors[tag == 'TWOFACT_1', Factor] ~
             g.gdp$factors[tag == 'F_SV1_TWOF', Factor] +
             g.gdp$factors[tag == 'F_SV3_TWOF', Factor]))

g.gdp$ols <- lm(g_us ~ PC1 + PC2, data = g.gdp$data[1:96])
g.gdp$pred <- data.table('yhat' = c(g.gdp$data[97:108, g_us],
                                    predict(g.gdp$ols,
                                            newdata = g.gdp$data[97:108]),
                                    g.gdp$tf$f[97:108]),
                         'tag' = c(rep('Y', 12), rep('G_US_OLS2F_CHW', 12),
                                   rep('G_US_KF2F_CHW', 12)),
                         'Date' = g.gdp$factors[97:108, Date])
ggplot(data = g.gdp$pred, aes(x = Date, y = yhat,
                              group = tag, colour = tag)) +
  theme_classic() + theme(plot.title = element_text(hjust = 0.5)) + 
  geom_line() + scale_x_date() + xlab('') + ylab('')

## RMSE
for (g in c('g_us', 'g_fr', 'g_ger', 'g_it')) {
  # KF two factors
  g.gdp$tf <- dlmFilter(g.gdp$data[[g]], mod = g.gdp$mod)
  
  # OLS + PC
  g.gdp$ols <- lm(paste0(g, ' ~ PC1 + PC2'), data = g.gdp$data[1:96])
  g.gdp$ols <- predict(g.gdp$ols, newdata = g.gdp$data[97:108])
  
  print(c(rmse(g.gdp$data[97:108][[g]], g.gdp$tf$f[97:108]),
          rmse(g.gdp$data[97:108][[g]], g.gdp$ols)))
}