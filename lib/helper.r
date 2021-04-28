RSS <- function(lm){
  return(sum(resid(lm)^2))
}

Wald.test <- function(mdl, mdl.free){
  
  ret1 <- (mdl.free$df.residual) * (RSS(mdl) - RSS(mdl.free))/RSS(mdl.free)
  q <- length(mdl.free$coefficients) - length(mdl$coefficients) 
  ret2 <- pchisq(ret1, df = q, lower.tail = FALSE)
  ret <- list("Statistic" = ret1, "pValue" = ret2, q = q)
}
