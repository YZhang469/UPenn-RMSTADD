source("../function/estBeta.r")

generateData <- function(n, 
                         censor){ # censor can take values "25%", or "50%"
  
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
  if(censor == "25%"){
    lambda0j <- seq(0.054, 0.096, 0.001)
    theta1 = 0.4
    theta2 = 0.1
    lambda = lambda0j[age-17] * exp(theta1*Z1+theta2*Z2)
    C <- rexp(n, rate = lambda)
  }
  else if (censor == "50%"){
    lambda0j <- seq(0.15, 0.36, 0.005)
    theta1 = 0.2
    theta2 = -0.2
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

sim <- function(censor, L, ns = c(1250, 2500, 5000, 10000), n.sim = 1000){
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
      res.temp <- estBeta(dat, Xname = "X", deltaXname = "deltaX", Znames = c("Z1", "Z2"), ZCnames = c("Z1", "Z2"), strname = "age", L = L)
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

res.final <- data.frame()
ns = c(1250,2500,5000,10000)
Ls = c(2.5,5)
for (k in 1:length(Ls)){
  res.final <- rbind.data.frame(res.final, sim(censor = "25%", L = Ls[k], ns = ns, n.sim = 1000))
  res.final <- rbind.data.frame(res.final, sim(censor = "50%", L = Ls[k], ns = ns, n.sim = 1000))
}
write.csv(res.final, "res_main.csv", row.names = FALSE)
