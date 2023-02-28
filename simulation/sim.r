library(survival)
library(MASS)

generateData <- function(n, 
                         censor){ # censor can take values "15%", "30%", "45%", or "60%"
  
  # strata: age
  eta0j <- seq(2.40, 4.08, 0.04)  # corresponding to eta0j for ages 18-60
  
  # covariates
  age <- 17 + sample.int(43, n, replace = TRUE)
  Z1 <- rbinom(n, size = 1, prob = 0.5)
  Z2 <- rnorm(n, mean = 0, sd = 1)
  
  # survival times
  eta0 <- eta0j[age-17]
  alpha1 = 0.9
  alpha2 = 0.4
  # D ~ Exp(lambda); 1/lambda = eta0+alpha1*Z1+alpha2*Z2
  D <- rexp(n, rate = 1/(eta0+alpha1*Z1+alpha2*Z2))
  
  # censoring times
  if(censor == "15%"){
    lambda0j <- seq(0.0184, 0.07, 0.0008)
    theta1 = 0.4
    theta2 = 0.1
    lambda = lambda0j[age-17] * exp(theta1*Z1+theta2*Z2)
    C <- rexp(n, rate = lambda)
  }
  else if (censor == "30%"){
    lambda0j <- seq(0.062, 0.160, 0.002)
    theta1 = 0.2
    theta2 = -0.2
    lambda = lambda0j[age-17] * exp(theta1*Z1+theta2*Z2)
    C <- rexp(n, rate = lambda)
  }
  else if (censor == "45%"){
    lambda0j <- seq(0.1, 0.31, 0.005)
    theta1 = 0.3
    theta2 = 0.2
    lambda = lambda0j[age-17] * exp(theta1*Z1+theta2*Z2)
    C <- rexp(n, rate = lambda)
  }
  else if (censor == "60%"){
    lambda0j <- seq(0.125, 0.65, 0.0125)
    theta1 = 0.4
    theta2 = -0.4
    lambda = lambda0j[age-17] * exp(theta1*Z1+theta2*Z2)
    C <- rexp(n, rate = lambda)
  }
  # calculate the censoring rate
  # sum(C<D)/n
  
  # combine data
  dat <- cbind.data.frame("ID" = 1:n, "age" = age, "Z1" = Z1, "Z2" = Z2, 
                          "X" = apply(cbind(D,C), 1, min, na.rm = TRUE), # Xi = Di^Ci
                          "deltaX" = ifelse(is.na(ifelse(D <= C, 1, 0)), 0, ifelse(D <= C, 1, 0))) # delta_i = I(Di <= Ci)
  return(dat)
}

trueBeta <- function(L){ # L = 2nd and 3rd quantile
  n = 10000000
  
  # strata: age
  eta0j <- seq(2.40, 4.08, 0.04) # corresponding to eta0j for ages 18-60
  
  # covariates
  age <- 17 + sample.int(43, n, replace = TRUE)
  Z1 <- rbinom(n, size = 1, prob = 0.5)
  Z2 <- rnorm(n, mean = 0, sd = 1)
  
  # survival times
  eta0 <- eta0j[age-17]
  alpha1 = 0.9
  alpha2 = 0.4
  
  # D ~ Exp(lambda); 1/lambda = eta0+alpha1*Z1+alpha2*Z2
  D <- rexp(n, rate = 1/(eta0+alpha1*Z1+alpha2*Z2))
  DL <- unlist(lapply(D, function(x){min(x, L)}))
  mod <- lm(DL ~ factor(age) + Z1 + Z2 - 1)
  beta <- coef(mod)[c("Z1","Z2")]
  
  return(list("beta" = beta, "mu0" = coef(mod)[1:43]))
}

estBeta <- function(dat, Xname = "X", deltaXname = "deltaX", Znames = c("Z1", "Z2"), strname = "age", L){
  dat$Y <- apply(cbind(dat[, Xname], L), 1, min) # Yi = Di^Ci^L
  dat$deltaY <- ifelse(dat[, deltaXname] == 1, 1, ifelse(L < dat[, Xname], 1, 0)) # delta_iki = I(Di^L <= Ci) = delta_i
  modC <- coxph(as.formula(paste("Surv(", Xname, ", 1-", deltaXname, ") ~ ", paste0(Znames, collapse = " + "), " + strata(", strname, ")")), 
                data = dat, ties = "breslow")
  dat.new <- dat
  dat.new[, Xname] <- dat.new$Y
  dat$W <- dat$deltaY/predict(modC, newdata = dat.new, type = "survival")
  n = nrow(dat)
  p = length(Znames)
  J = length(unique(dat[, strname]))
  strvalue <- sort(unique(dat[, strname]))
  Zbar <-  matrix(ncol = p, nrow = J)
  for (j in 1:J){
    Zbar[j, ] <- colSums(as.data.frame(dat[dat[, strname] == strvalue[j], Znames]) * dat$W[dat[, strname] == strvalue[j]]) / sum(dat$W[dat[, strname] == strvalue[j]])
  }
  Z <- dat[, Znames]
  Zres <- Z - Zbar[match(dat[, strname], strvalue), ]
  B <- apply(Zres * dat$W * dat$Y, 2, sum)
  A <- matrix(0, ncol = p, nrow = p)
  for (i in 1:n){
    A = A + t(as.matrix(Zres[i, ])) %*% as.matrix(Zres[i, ]) * dat$W[i]
  }
  betahat <- solve(A) %*% B
  mu0 <- rep(0, J)
  for(j in 1:J){
    dat.temp <- dat[dat[, strname] == strvalue[j], ]
    mu0[j] <- sum(dat.temp$W * (dat.temp$Y - as.matrix(dat.temp[, Znames]) %*% betahat)) / sum(dat.temp$W)
  }
  Avar <- matrix(0, ncol = p, nrow = p)
  for (i in 1:n){
    Avar = Avar + t(as.matrix(Zres[i, ])) %*% as.matrix(Zres[i, ]) * dat$W[i] / n
  }
  Bvar <- matrix(0, ncol = p, nrow = p)
  for (i in 1:nrow(Zres)){
    temp <- Zres[i, ] * dat$W[i] * (dat$Y[i] - mu0[match(dat[i, strname], strvalue)] - as.matrix(Z[i, ]) %*% betahat)
    Bvar = Bvar + t(as.matrix(temp)) %*% as.matrix(temp) / n
  }
  var <- diag(t(solve(Avar)) %*% Bvar %*% solve(Avar)) / n
  
  return(list("betahat" = betahat, "var" = var))
}

sim <- function(censor, L, ns = c(1250, 2500, 5000, 10000), n.sim = 500){
  beta <- trueBeta(L)$beta
  p <- length(beta)
  res <- data.frame(matrix(ncol = 8, nrow = length(ns)*p))
  colnames(res) <- c("L", "Censoring", "Parameter", "n", "Bias", "ESD", "ASE", "CP")
  res$L <- L
  res$Censoring <- censor
  res$Parameter <- rep(round(beta, digits = 3), each = length(ns))
  res$n <- rep(ns, 2)
  
  for (s in 1:length(ns)){
    est <- data.frame(matrix(ncol = p, nrow = n.sim))
    se <- data.frame(matrix(ncol = p, nrow = n.sim))
    cover <- data.frame(matrix(ncol = p, nrow = n.sim))
    for (iter in 1:n.sim){
      dat <- generateData(n = ns[s], censor = censor)
      res.temp <- estBeta(dat, Xname = "X", deltaXname = "deltaX", Znames = c("Z1", "Z2"), strname = "age", L = L)
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
    
    res$Bias[seq(from = s, to = s + length(ns) * (p-1), by = length(ns))] <- round(bias, digits = 3)
    res$ESD[seq(from = s, to = s + length(ns) * (p-1), by = length(ns))] <- round(esd, digits = 3)
    res$ASE[seq(from = s, to = s + length(ns) * (p-1), by = length(ns))] <- round(ase, digits = 3)
    res$CP[seq(from = s, to = s + length(ns) * (p-1), by = length(ns))] <- round(cp, digits = 3)
  }
  
  return(res)
}

# Table 1
set.seed(0227)
res.final <- data.frame()
ns = c(1250,2500,5000,10000)
Ls = c(2.5,5)
for (k in 1:length(Ls)){
  res.final <- rbind.data.frame(res.final, sim(censor = "15%", L = Ls[k], ns = ns, n.sim = 500))
  res.final <- rbind.data.frame(res.final, sim(censor = "30%", L = Ls[k], ns = ns, n.sim = 500))
}
write.csv(res.final, "res_low.csv", row.names = FALSE)

# Table 2
set.seed(0315)
res.final <- data.frame()
ns = c(2500,5000,10000,20000)
Ls = c(2.5,5)
for (k in 1:length(Ls)){
  res.final <- rbind.data.frame(res.final, sim(censor = "45%", L = Ls[k], ns = ns, n.sim = 500))
  res.final <- rbind.data.frame(res.final, sim(censor = "60%", L = Ls[k], ns = ns, n.sim = 500))
}
write.csv(res.final, "res_high.csv", row.names = FALSE)
