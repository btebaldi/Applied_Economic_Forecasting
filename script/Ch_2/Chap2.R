
# setup -------------------------------------------------------------------
rm(list = ls())


library(readxl)
# library(data.table)
library(dplyr)
library(ggplot2)
library(lmtest)
library(strucchange)
library(tseries)
library(xtable)

source("./lib/helper.r")

# 2.8 Simulated Data ------------------------------------------------------
sim.data <- list()
sim.data$full <- read_excel("./database/simulated_model_ch2_sec8.xlsx")
sim.data$est <- sim.data$full[102:301,] # estimation sample
sim.data$fore <- sim.data$full[302:401,] # forecasting sample


# Table 2.7.1 - Breusch-Pagan-Godfrey homoskedasticity test ---------------
# OLS model fit
ols.fit <- lm(y ~ x, data = sim.data$est)
summary(ols.fit)

# By "hand" calculations
tbl <- sim.data$est 
tbl$res <- ols.fit$residuals
tbl$res2 <- tbl$res^2

BPG.mdl <- lm(res2 ~ x, data = tbl)
summary(BPG.mdl)

BPG.mdl.sta <- summary(BPG.mdl)$r.squared * length(BPG.mdl$fitted.values)
BPG.mdl.df <- (length(BPG.mdl$coefficients)-1)
BPG.mdl.pvalue <- pchisq(BPG.mdl.sta,
       df = BPG.mdl.df,
       lower.tail = FALSE)

cat(sprintf("Breusch-Pagan test\n\nBP = %6.4f, df = %d, p-value = %8.6f\n",
            BPG.mdl.sta,BPG.mdl.df,BPG.mdl.pvalue))

# Using lmtest function
lmtest::bptest(ols.fit)

rm(list = c("BPG.mdl.sta", "BPG.mdl.df", "BPG.mdl.pvalue", "tbl"))


# Figure 2.8.1 - Normality test & Histogram -------------------------------

# Normality test & Histogram
sim.data$est$eps = ols.fit$residuals
hist(sim.data$est$eps, xlab = '', ylab = '', main = '') 
sim.data$est %>% 
  ggplot() +
  geom_histogram(aes(x=eps), binwidth = 2, fill="#FFFFFF00", colour="black")

jarque.bera.test(sim.data$est$eps) 



# Table 2.8.1: White heteroskedasticity test ------------------------------
sim.data$est <- sim.data$est %>% mutate(x2 = x^2, eps2 = eps^2)
white.fit <- lm(eps2 ~ x2 + x, data = sim.data$est) 
summary(white.fit)
 
white.fit.sta <- summary(white.fit)$r.squared * length(white.fit$fitted.values)
white.fit.df <- (length(white.fit$coefficients)-1)
white.fit.pvalue <- pchisq(white.fit.sta,
                         df = white.fit.df,
                         lower.tail = FALSE)


cat(sprintf("White heteroskedasticity test\n\nWH = %6.4f, df = %d, p-value = %8.6f\n",
            white.fit.sta, white.fit.df, white.fit.pvalue))

rm(list = c("white.fit.sta", "white.fit.df", "white.fit.pvalue", "white.fit"))
# Table 2.8.2 Breusch-Godfrey LM test -------------------------------------

# Durbin-Watson test
dwtest(ols.fit) 

CH.mdl1 <- lm(y~x,data = sim.data$est)
CH.mdl1.T1 <- lm(y~x,data = sim.data$est[1:49, ])
CH.mdl1.T2 <- lm(y~x,data = sim.data$est[50:200, ])
k <- length(CH.mdl1$coefficients)
CH.df1 <- k
CH.df2 <- length(CH.mdl1$residuals) - 2*k 
CH <- (RSS(CH.mdl1) - RSS(CH.mdl1.T1) - RSS(CH.mdl1.T2))/(RSS(CH.mdl1.T1) + RSS(CH.mdl1.T2)) * (CH.df2)/(CH.df1)
CH.df2.pvalue <- pf(CH,
   df1 = CH.df1,
   df2 = CH.df2,
   lower.tail = FALSE)

# logLikelihood test
tbl <- sim.data$est %>% mutate(D = as.integer(row_number()>49))

mdl.1 <- lm(y ~ x*D, data = tbl)
LR.stat <- -2*(logLik(CH.mdl1) - logLik(mdl.1))
LR.df <- length(mdl.1$coefficients) - length(CH.mdl1$coefficients)
LR.pvalue <- pchisq(LR.stat, df = LR.df, lower.tail = FALSE)

# Wald Test

W <- Wald.test(mdl = CH.mdl1, mdl.free = mdl.1)
cat(sprintf("Chow Breakpoint Test: %d\n%20s\t%6.3f\tProb. F(%d,%d)\t%5.3f\n%20s\t%6.3f\tProb. ChiSq(%d)\t%5.3f\n%20s\t%6.3f\tProb. ChiSq(%d)\t%5.3f\n", 102+49,
            "F-Statistic", CH, CH.df1, CH.df2, CH.df2.pvalue,
            "Log likelihood ratio", LR.stat, LR.df, LR.pvalue,
            "Wald Statistic", W$Statistic, W$q, W$pValue))

CH.mdl1.T1 <- lm(y~x,data = sim.data$est[1:99, ])
CH.mdl1.T2 <- lm(y~x,data = sim.data$est[100:200, ])
k <- length(CH.mdl1$coefficients) 
CH.df1 <- k
CH.df2 <- length(CH.mdl1$residuals) - 2*k 
CH <- (RSS(CH.mdl1) - RSS(CH.mdl1.T1) - RSS(CH.mdl1.T2))/(RSS(CH.mdl1.T1) + RSS(CH.mdl1.T2)) * (CH.df2)/(CH.df1)
CH.df2.pvalue <- pf(CH,
                    df1 = CH.df1,
                    df2 = CH.df2,
                    lower.tail = FALSE)

# logLikelihood test
tbl <- sim.data$est %>% mutate(D = as.integer(row_number()>99))

mdl.1 <- lm(y ~ x*D, data = tbl)
LR.stat <- -2*(logLik(CH.mdl1) - logLik(mdl.1))
LR.df <- length(mdl.1$coefficients) - length(CH.mdl1$coefficients)
LR.pvalue <- pchisq(LR.stat, df = LR.df, lower.tail = FALSE)

W <- Wald.test(mdl = CH.mdl1, mdl.free = mdl.1)
cat(sprintf("Chow Breakpoint Test: %d\n%20s\t%6.3f\tProb. F(%d,%d)\t%5.3f\n%20s\t%6.3f\tProb. ChiSq(%d)\t%5.3f\n%20s\t%6.3f\tProb. ChiSq(%d)\t%5.3f\n", 102+49,
            "F-Statistic", CH, CH.df1, CH.df2, CH.df2.pvalue,
            "Log likelihood ratio", LR.stat, LR.df, LR.pvalue,
            "Wald Statistic", W$Statistic, W$q, W$pValue))



CH.mdl1.T1 <- lm(y~x,data = sim.data$est[1:149, ])
CH.mdl1.T2 <- lm(y~x,data = sim.data$est[150:200, ])
k <- length(CH.mdl1$coefficients) 
CH.df1 <- k
CH.df2 <- length(CH.mdl1$residuals) - 2*k 
CH <- (RSS(CH.mdl1) - RSS(CH.mdl1.T1) - RSS(CH.mdl1.T2))/(RSS(CH.mdl1.T1) + RSS(CH.mdl1.T2)) * (CH.df2)/(CH.df1)
CH.df2.pvalue <- pf(CH,
                    df1 = CH.df1,
                    df2 = CH.df2,
                    lower.tail = FALSE)

# logLikelihood test
tbl <- sim.data$est %>% mutate(D = as.integer(row_number()>149))

mdl.1 <- lm(y ~ x*D, data = tbl)
LR.stat <- -2*(logLik(CH.mdl1) - logLik(mdl.1))
LR.df <- length(mdl.1$coefficients) - length(CH.mdl1$coefficients)
LR.pvalue <- pchisq(LR.stat, df = LR.df, lower.tail = FALSE)

W <- Wald.test(mdl = CH.mdl1, mdl.free = mdl.1)
cat(sprintf("Chow Breakpoint Test: %d\n%20s\t%6.3f\tProb. F(%d,%d)\t%5.3f\n%20s\t%6.3f\tProb. ChiSq(%d)\t%5.3f\n%20s\t%6.3f\tProb. ChiSq(%d)\t%5.3f\n", 102+49,
            "F-Statistic", CH, CH.df1, CH.df2, CH.df2.pvalue,
            "Log likelihood ratio", LR.stat, LR.df, LR.pvalue,
            "Wald Statistic", W$Statistic, W$q, W$pValue))


rm(list = c("CH", "CH.df1", "CH.df2", "CH.df2.pvalue", "mdl.1", "LR.stat", "LR.df", "LR.pvalue"))
rm(list = c("tbl", "k", "CH.mdl1", "CH.mdl1.T1", "CH.mdl1.T2", "BPG.mdl", "W"))

# Chow break point test
# sctest(y ~ x, data = sim.data$est, type = 'Chow', from = 0.1, to = 0.255)
# sctest(y ~ x, data = sim.data$est, type = 'Chow', from = 0.25, to = 0.5)
# sctest(y ~ x, data = sim.data$est, type = 'Chow', from = 0.5, to = 0.75)
# 
# fs <- Fstats(y ~ x, data = sim.data$est, from = 0.25, to = NULL)
# fs <- Fstats(y ~ x, data = sim.data$est, from = 0.5, to = NULL)
# fs <- Fstats(y ~ x, data = sim.data$est, from = 0.75, to = NULL)

tbl <- tibble(res1 = sim.data$full$y - (sim.data$full$x * ols.fit$coefficients["x"] + ols.fit$coefficients["(Intercept)"]),
              x = sim.data$full$x)
 
tbl <- tbl %>% mutate(res1_1 = lag(res1,1),
               res1_2 = lag(res1,2),
               res1_3 = lag(res1,3),
               res1_4 = lag(res1,4)) %>% 
  filter(row_number() >= 102 & row_number() <= 301)

tbl <- tibble(res1 = ols.fit$residuals, x = ols.fit$model$x) %>% mutate(res1_1 = lag(res1,1),
                                                                        res1_2 = lag(res1,2),
                                                                        res1_3 = lag(res1,3),
                                                                        res1_4 = lag(res1,4))

tbl[is.na(tbl)] = 10

BG.mdl <- lm(res1 ~ 1 + x + res1_1 + res1_2 + res1_3 + res1_4, data=tbl)
names(BG.mdl$coefficients) = c("(Intercept)", "x", "Res_1", "Res_2", "Res_3", "Res_4")
BG.mdl.stat <- summary(BG.mdl)$r.squared * length(BG.mdl$residuals)
BG.mdl.stat
lmtest::bgtest(y ~ x, order = 4, data = sim.data$est, fill = 10)


# Bai-Perron test
breakpoints(y ~ x, data = sim.data$est, h = 0.15)

# Dummy variable
olsD.fit <- lm(y ~ D + x + D * x, data = sim.data$est)
summary(olsD.fit)

# Simple forecasts - with / no dummy
yhat <- list()
yhat$y <- predict(ols.fit, newdata = sim.data$fore) # no dummy
yhat$yD <- predict(olsD.fit, newdata = sim.data$fore) # with dummy

yhat$yD.se <- sqrt(sum(olsD.fit$residuals^2) / 198)
yhat$yD.up <- yhat$yD + 1.96 * yhat$yD.se
yhat$yD.low <- yhat$yD - 1.96 * yhat$yD.se

# Plot - y / yhat1 / yhat2
yhat.plot <- data.table('yhat' = rbindlist(list(data.table(sim.data$fore[, y]),
                                                data.table(yhat$y),
                                                data.table(yhat$yD))),
                        'label' = rep(c('Y', 'YHAT1', 'YHAT2'), each = 100))

ggplot(yhat.plot, aes(x = rep(302:401, 3), y = yhat, linetype = label)) +
  geom_line() + xlab('') + ylab('') + theme(legend.title = element_blank()) 

# Plot - yhat2
yhat.plot <- data.table('yhat' = rbindlist(list(data.table(yhat$yD),
                                                data.table(yhat$yD.up),
                                                data.table(yhat$yD.low))),
                        'label' = rep(c('YHAT2', 'YHAT2_UP', 'YHAT2_LOW'),
                                      each = 100))

ggplot(yhat.plot, aes(x = rep(302:401, 3), y = yhat, linetype = label)) +
  geom_line() + xlab('') + ylab('') + theme(legend.title = element_blank()) 

# Recursive
yhat$y.rec <- yhat$yD.rec <- yhat$yD.recse <- rep(0, 100)

for (i in 1:100) {
  ols.rec <- lm(y ~ x, data = sim.data$full[102:(300 + i)])
  yhat$y.rec[i] <- predict(ols.rec, newdata = sim.data$full[301 + i])
  
  olsD.rec <- lm(y ~ D + x + D * x, data = sim.data$full[102:(300 + i)])
  yhat$yD.rec[i] <- predict(olsD.rec, newdata = sim.data$full[301 + i])
  yhat$yD.recse[i] <- sqrt(sum(olsD.rec$residuals^2) / (197 + i))
}

# Plot - recursive forecasts
yrec.plot <- data.table('yhat' = rbindlist(list(data.table(sim.data$fore[, y]),
                                                data.table(yhat$y.rec),
                                                data.table(yhat$yD.rec))),
                        'label' = rep(c('Y', 'YHAT1_REC', 'YHAT2_REC'),
                                      each = 100))

# Plot - recursive forecasts with dummy
ggplot(yrec.plot, aes(x = rep(302:401, 3), y = yhat, linetype = label)) +
  geom_line() + xlab('') + ylab('') + theme(legend.title = element_blank()) 

yrec.plot <- data.table('yhat' =
                          rbindlist(list(data.table(yhat$yD.rec),
                                         data.table(yhat$yD.rec +
                                                      1.96 * yhat$yD.recse),
                                         data.table(yhat$yD.rec -
                                                      1.96 * yhat$yD.recse))),
                        'label' =
                          rep(c('YHAT2_REC', 'YHAT2_REC_UP', 'YHAT2_REC_LOW'),
                              each = 100))

ggplot(yrec.plot, aes(x = rep(302:401, 3), y = yhat, linetype = label)) +
  geom_line() + xlab('') + ylab('') + theme(legend.title = element_blank()) 

# RMSE & MAE
yhat$Y <- cbind(yhat$y, yhat$yD.rec)
RMSE <- sqrt(colSums((yhat$Y - sim.data$fore[, y])^2) / 100)
MAE <- colSums(abs(yhat$Y - sim.data$fore[, y])) / 100
error.mat <- rbind(RMSE, MAE)
colnames(error.mat) <- c('Simple', 'Recursive')
print(xtable(error.mat), include.rownames = T, include.colnames = T)

### 2.9.1 Forecasting Euro Area GDP ###
eu.gdp <- fread('ex2_misspecification_gdp.csv')
gdp.fit <- lm(y ~ ipr + su + sr, data = eu.gdp)

# Breusch Pagan Test
bptest(gdp.fit)

# White test
eu.gdp[, c('eps', 'eps2', 'ipr2', 'su2', 'sr2') :=
         list(gdp.fit$residuals, gdp.fit$residuals^2, ipr^2, su^2, sr^2)]
white.fit <- lm(eps2 ~ ipr + ipr2 + ipr * su + ipr * sr +
                  su + su2 + su * sr + sr + sr2, data = eu.gdp)
summary(white.fit)

# Durbin-Watson test
dwtest(gdp.fit)  

# Breusch-Godfrey test
bgtest(y ~ ipr + su + sr, data = eu.gdp, order = 2, type = c('Chisq', 'F'))
bgtest(y ~ ipr + su + sr, data = eu.gdp, order = 3, type = c('Chisq', 'F'))

# Normality test & Histogram
hist(eu.gdp[, eps], xlab = '', ylab = '', main = '') 
jarque.bera.test(eu.gdp[, eps]) 

# Chow break point test
sctest(y ~ ipr + su + sr, data = eu.gdp, type = 'Chow', from = 0.1, to = 0.3)
sctest(y ~ ipr + su + sr, data = eu.gdp, type = 'Chow', from = 0.3, to = 0.7)
sctest(y ~ ipr + su + sr, data = eu.gdp, type = 'Chow', from = 0.7, to = 1)

fs <- Fstats(y ~ x, data = eu.gdp, from = 0.3, to = 1)
fs <- Fstats(y ~ x, data = eu.gdp, from = 0.7, to = 1)

# Bai-Perron test
breakpoints(y ~ ipr + su + sr, data = eu.gdp, h = 0.15)

# Recursive estimation
gdp.rr <- recresid(y ~ ipr + su + sr, data = eu.gdp$full)
plot(gdp.rr, type = 'l')

# Dummy variable
gdpD.fit <- list()
gdpD.formula <- c('y ~ ipr + su + sr + Dea', 'y ~ ipr + su + sr + D2000s',
                  'y ~ ipr + su + sr + Dea + D2000s')
for (model in 1:3) {
  gdpD.fit[[model]] <- lm(gdpD.formula[model], data = eu.gdp)
  print(summary(gdpD.fit[[model]]))
}

# Forecasting
gdp.hat <- list()
gdp.fit <- lm(y ~ ipr + su + sr, data = eu.gdp[1:60]) # Model 2 - no dummy
gdp.hat$ghat <- predict(gdp.fit, newdata = eu.gdp[61:70]) 

gdp.fit <- lm(y ~ ipr + su + sr + Dea + D2000s, data = eu.gdp[1:60])
gdp.hat$ghat3 <- predict(gdp.fit, newdata = eu.gdp[61:70]) # Model 2.3 - dummy

gdp.plot <- data.table('yhat' = rbindlist(list(data.table(eu.gdp[61:70, y]),
                                               data.table(gdp.hat$ghat),
                                               data.table(gdp.hat$ghat3))),
                       'label' = rep(c('Y', 'YFOREG2_NEW', 'YFOREG2_3'),
                                     each = 10))

ggplot(gdp.plot, aes(x = rep(1:10, 3), y = yhat, linetype = label)) +
  geom_line() + xlab('') + ylab('') + theme(legend.title = element_blank()) 

# RMSE & MAE
gdp.hat$Y <- cbind(gdp.hat$ghat, gdp.hat$ghat3)
RMSE <- sqrt(colSums((gdp.hat$Y - eu.gdp[61:70, y])^2) / 10)
MAE <- colSums(abs(gdp.hat$Y - eu.gdp[61:70, y])) / 10
error.mat <- rbind(RMSE, MAE)
colnames(error.mat) <- c('Model 2', 'Model 2.3')
print(xtable(error.mat), include.rownames = T, include.colnames = T)

### 2.9.2 Forecasting US GDP ###
us.gdp <- fread('ex2_misspecification_gdp_us.csv')
gdp.fit <- lm(y ~ ipr + su + sr, data = us.gdp)

# Breusch Pagan Test
bptest(gdp.fit)

# White test
us.gdp[, c('eps', 'eps2', 'ipr2', 'su2', 'sr2') :=
         list(gdp.fit$residuals, gdp.fit$residuals^2, ipr^2, su^2, sr^2)]
white.fit <- lm(eps2 ~ ipr + ipr2 + ipr * su + ipr * sr +
                  su + su2 + su * sr + sr + sr2, data = us.gdp)
summary(white.fit)

# Durbin-Watson test
dwtest(gdp.fit)  

# Normality test & Histogram
hist(us.gdp[, eps], xlab = '', ylab = '', main = '') 
jarque.bera.test(us.gdp[, eps]) 

# Chow break point test
sctest(y ~ ipr + su + sr, data = us.gdp, type = 'Chow', from = 0.55, to = 0.6)
sctest(y ~ ipr + su + sr, data = us.gdp, type = 'Chow', from = 0.75, to = 0.8)

fs <- Fstats(y ~ x, data = us.gdp, from = 0.55, to = 1)
fs <- Fstats(y ~ x, data = us.gdp, from = 0.75, to = 1)

# Bai-Perron test
breakpoints(y ~ ipr + su + sr, data = eu.gdp, h = 0.15)

# Recursive estimation
gdp.rr <- recresid(y ~ ipr + su + sr, data = us.gdp)
plot(gdp.rr, type = 'l')

# Dummy variable
gdpD.fit <- list()
gdpD.formula <- c('y ~ ipr + su + sr + Dfincris', 'y ~ ipr + su + sr + D2000s',
                  'y ~ ipr + su + sr + Dfincris + D2000s')
for (model in 1:3) {
  gdpD.fit[[model]] <- lm(gdpD.formula[model], data = us.gdp)
}
summary(gdpD.fit[[1]])
summary(gdpD.fit[[2]])
summary(gdpD.fit[[3]])

## Forecasting
gdp.hat <- list()

# Model 2 - no dummy
gdp.fit <- lm(y ~ ipr + su + sr, data = us.gdp[1:104]) 
gdp.hat$ghat <- predict(gdp.fit, newdata = us.gdp[105:114]) 

# Model 2.3 - dummy
gdp.fit <- lm(y ~ ipr + su + sr + Dfincris + D2000s, data = us.gdp[1:104])
gdp.hat$ghat3 <- predict(gdp.fit, newdata = us.gdp[105:114]) 

gdp.plot <- data.table('yhat' = rbindlist(list(data.table(us.gdp[105:114, y]),
                                               data.table(gdp.hat$ghat),
                                               data.table(gdp.hat$ghat3))),
                       'label' = rep(c('Y', 'YFOREG2_3', 'YRFOREG2_3'),
                                     each = 10))

ggplot(gdp.plot, aes(x = rep(1:10, 3), y = yhat, linetype = label)) +
  geom_line() + xlab('') + ylab('') + theme(legend.title = element_blank()) 

# Recursive
for (i in 1:10) {
  ols.rec <- lm(y ~ ipr + su + sr, data = us.gdp[1:(103 + i)])
  gdp.hat$rec[i] <- predict(ols.rec, newdata = us.gdp[104 + i])
  
  olsD.rec <- lm(y ~ ipr + su + sr + Dfincris + D2000s,
                 data = us.gdp[1:(103 + i)])
  gdp.hat$rec3[i] <- predict(olsD.rec, newdata = us.gdp[104 + i])
}

# RMSE & MAE
gdp.hat$Y <- cbind(gdp.hat$ghat, gdp.hat$ghat3) # simple
RMSE <- sqrt(colSums((gdp.hat$Y - us.gdp[105:114, y])^2) / 10)
MAE <- colSums(abs(gdp.hat$Y - us.gdp[105:114, y])) / 10
error.mat <- rbind(RMSE, MAE)

# Recursive RMSE & MAE
gdp.hat$Yrec <- cbind(gdp.hat$rec, gdp.hat$rec3) # recursive
RMSE <- sqrt(colSums((gdp.hat$Yrec - us.gdp[105:114, y])^2) / 10)
MAE <- colSums(abs(gdp.hat$Yrec - us.gdp[105:114, y])) / 10
error.mat <- rbind(error.mat, RMSE, MAE)
rownames(error.mat) <- c('Simple RMSE', 'Simple MAE',
                         'Recursive RMSE', 'Recursive MAE')
colnames(error.mat) <- c('Model 2', 'Model 2.3')
print(xtable(error.mat), include.rownames = T, include.colnames = T)

### 2.9.3 Default Risk ###
default.risk <- fread('default_risk.csv')
default.risk[, Date := as.Date(Date, format = '%m/%d/%Y')] 
default.risk[, OAS := OAS[2:216]]
default.risk <- default.risk[1:215]

# Dummy and interaction term
default.risk[Date >= '2008-01-01' & Date < '2010-01-01', D := 1]
default.risk[Date < '2008-01-01' | Date >= '2010-01-01', D:= 0]
default.risk[, c('VIX.D', 'SENT.D', 'PMI.D', 'sp500.D') :=
               list(VIX * D, SENT * D, PMI * D, sp500 * D)]

# Dummy
default.D <- list('M1' = 'OAS ~ VIX + D', 'M2' = 'OAS ~ SENT + D',
                  'M3' = 'OAS ~ PMI + D', 'M4' = 'OAS ~ sp500 + D')

# Interaction
default.I <- list('M1' = 'OAS ~ VIX + D + VIX.D',
                  'M2' = 'OAS ~ SENT + D + SENT.D',
                  'M3' = 'OAS ~ PMI + D + PMI.D',
                  'M4' = 'OAS ~ sp500 + D + sp500.D')

for (m in c('M1', 'M2', 'M3', 'M4')) {
  fit.D <- lm(default.D[[m]], data = default.risk)
  print(summary(fit.D))
  print(coeftest(fit.D, vcov = NeweyWest(fit.D, lag = 12)))
}

for (m in c('M1', 'M2', 'M3', 'M4')) {
  fit.I <- lm(default.I[[m]], data = default.risk)
  print(summary(fit.I))
  print(coeftest(fit.I, vcov = NeweyWest(fit.I, lag = 12)))
}