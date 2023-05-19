library(MASS)
library(survival)

source("./function/estBeta.r")

## simulation for the number of covariates > 10

generateData.p <- function(n = 5000, p = 15){ # censor = 50%
  
  # strata: age
  eta0j <- seq(2.40, 4.08, 0.04)  # corresponding to eta0j for ages 18-60
  
  # covariates
  age <- 17 + sample.int(43, n, replace = TRUE)
  Z <- mvrnorm(n, mu = rep(0, ceiling(p/2)), Sigma = matrix(diag(rep(c(1, 0.8, 0.7, 0.8, 1), length.out = ceiling(p/2))), ncol = ceiling(p/2)))
  for (q in 1:(p-ceiling(p/2))){
    Z <- cbind.data.frame(Z, rbinom(n, 1, 0.5))
  }
  colnames(Z) <- paste("Z", 1:p, sep = "")
  
  # survival times
  eta0 <- eta0j[age-17]
  alpha <- rep(c(0.1, 0.16, 0.12, 0.16, 0.1), length.out = p)
  # D ~ Exp(lambda); 1/lambda = eta0+alpha1*Z1+alpha2*Z2
  D <- rexp(n, rate = 1/(eta0 + as.matrix(Z)%*%alpha))
  summary(D)
  
  # censoring times
  lambda0j <- seq(0.1, 0.31, 0.005)
  theta <- rep(c(0.07, 0.1, 0.12, 0.08, 0.1), length.out = p)
  lambda = lambda0j[age-17] * exp(as.matrix(Z)%*%theta)
  C <- rexp(n, rate = lambda)
  # calculate the censoring rate
  # sum(C<D)/n
  
  # combine data
  dat <- cbind.data.frame("ID" = 1:n, "age" = age, Z, 
                          "X" = apply(cbind(D,C), 1, min, na.rm = TRUE), # Xi = Di^Ci
                          "deltaX" = ifelse(is.na(ifelse(D <= C, 1, 0)), 0, ifelse(D <= C, 1, 0))) # delta_i = I(Di <= Ci)
  return(dat)
}

trueBeta.p <- function(L, p){ # L = 2nd and 3rd quantile
  n = 10000000
  
  # strata: age
  eta0j <- seq(2.40, 4.08, 0.04) # corresponding to eta0j for ages 18-60
  
  # covariates
  age <- 17 + sample.int(43, n, replace = TRUE)
  Z <- mvrnorm(n, mu = rep(0, ceiling(p/2)), Sigma = matrix(diag(rep(c(1, 0.8, 0.7, 0.8, 1), length.out = ceiling(p/2))), ncol = ceiling(p/2)))
  for (q in 1:(p-ceiling(p/2))){
    Z <- cbind.data.frame(Z, rbinom(n, 1, 0.5))
  }
  colnames(Z) <- paste("Z", 1:p, sep = "")
  
  # survival times
  eta0 <- eta0j[age-17]
  alpha <- rep(c(0.1, 0.16, 0.12, 0.16, 0.1), length.out = p)
  D <- rexp(n, rate = 1/(eta0 + as.matrix(Z)%*%alpha))
  DL <- unlist(lapply(D, function(x){min(x, L)}))
  dat <- cbind.data.frame(age, Z, DL)
  mod <- lm(as.formula(paste("DL ~ factor(age) + ", paste0(colnames(Z), collapse = " + "), " - 1")), dat = dat)
  beta <- coef(mod)[colnames(Z)]
  
  return(list("beta" = beta, "mu0" = coef(mod)[1:43]))
}

sim.p <- function(L, p, n.sim = 1000){
  beta <- trueBeta.p(L, p)$beta
  
  est <- data.frame(matrix(ncol = p, nrow = n.sim))
  se <- data.frame(matrix(ncol = p, nrow = n.sim))
  cover <- data.frame(matrix(ncol = p, nrow = n.sim))
  for (iter in 1:n.sim){
    dat <- generateData.p(n = 5000, p = p)
    res.temp <- estBeta(dat, Xname = "X", deltaXname = "deltaX", 
                        Znames = paste("Z", 1:p, sep = ""), ZCnames = paste("Z", 1:p, sep = ""), 
                        strname = "age", L = L)
    betahat <- res.temp$betahat
    var <- res.temp$var
    ci <- ifelse(beta > betahat-1.96*sqrt(var) & beta < betahat+1.96*sqrt(var), 1, 0)
    est[iter, ] <- betahat
    se[iter, ] <- sqrt(var)
    cover[iter, ] <- ci
  }
  
  ## metrics
  # bias
  bias <- apply(est, 2, mean) - beta
  # empirical standard deviation (ESD)
  esd <- sqrt(apply(apply(est, 1, function(x){(x - beta)^2}), 1, mean))
  # average asymptotic standard error (ASE)
  ase <- apply(se, 2, mean)
  # empirical coverage probabilities (CP): use ASE to create a confidence interval, and test if the true value lies within the interval
  cp <- apply(cover, 2, mean)
  
  res <- cbind.data.frame("L" = L, "Parameter" = round(beta, digits = 3), 
                          "Bias" = round(bias, digits = 3), 
                          "ESD" = round(esd, digits = 3), 
                          "ASE" = round(ase, digits = 3), 
                          "CP" = round(cp, digits = 3))
  
  return(res)
}

set.seed(0408)
res.final <- data.frame()
Ls = c(2.5,5)
for (k in 1:length(Ls)){
  res.final <- rbind.data.frame(res.final, sim.p(L = Ls[k], p = 15, n.sim = 1000))
}
write.csv(res.final, "./results/TableB1.csv", row.names = FALSE)
