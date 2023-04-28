library(haven)
library(survival)
library(ggplot2)
library(dplyr)
library(survRM2)

## data cleaning
dat <- read_sas("postlt_yuan_15jul2022.sas7bdat")
dat <- dat[dat$status_1 == 0 & dat$AGE_DON >= 12, ]
dat <- dat[!is.na(dat$don_bili), ]

## functions
estStrAddModel <- function(dat, Znames, modC){
  dat.new <- dat
  dat.new$X <- dat.new$Y
  dat$W <- dat$deltaY/predict(modC, newdata = dat.new, type = "survival")
  
  Z <- dat[, Znames]
  Z <- apply(Z, 2, as.numeric)
  Zbar <-  matrix(ncol = ncol(Z), nrow = length(unique(dat$stratum)))
  for (j in 1:length(unique(dat$stratum))){
    Zbar[j, ] <- apply(Z, 2, function(x){sum(dat$W[dat$stratum == unique(dat$stratum)[j]] * x[dat$stratum == unique(dat$stratum)[j]])}) / sum(dat$W[dat$stratum == unique(dat$stratum)[j]])
  }
  Zres <- Z - Zbar[match(dat$stratum, unique(dat$stratum)), ]
  B <- apply(Zres * dat$W * dat$Y, 2, sum)
  A <- matrix(0, ncol = ncol(Z), nrow = ncol(Z))
  for (i in 1:nrow(dat)){
    A = A + as.matrix(Zres[i, ]) %*% t(as.matrix(Zres[i, ])) * dat$W[i]
  }
  betahat <- solve(A) %*% B
  
  mu0 <- rep(0, length(unique(dat$stratum)))
  for(j in 1:length(mu0)){
    dat.temp <- dat[dat$stratum == unique(dat$stratum)[j], ]
    Z.temp <- apply(dat.temp[, Znames], 2, as.numeric)
    mu0[j] <- sum(dat.temp$W * (dat.temp$Y - Z.temp %*% betahat)) / sum(dat.temp$W)
  }
  mu0_df <- cbind.data.frame("stratum" = unique(dat$stratum), "mu0" = mu0)
  
  return(list("betahat" = betahat, "mu0hat" = mu0_df))
}

estStrMultModel <- function(dat, Znames, modC){
  dat.new <- dat
  dat.new$X <- dat.new$Y
  dat.new$W <- 1/predict(modC, newdata = dat.new, type = "survival")
  dat.new$W.new <- dat.new$deltaY * dat.new$Y * dat.new$W
  dat.new$delta.new <- 1
  dat.new$X.new <- 1
  dat.new$Off <- -log(dat.new$Y)
  mod.new <- coxph(as.formula(paste("Surv(X.new, delta.new) ~ ", paste(Znames, collapse = " + "), " + offset(-log(Y)) + strata(stratum)")), data = dat.new)
  betahat <- coefficients(mod.new)
  num <- tapply(dat.new$deltaY * dat.new$W * dat.new$Y, as.factor(dat.new$stratum), sum)
  denom <- tapply(dat.new$deltaY * dat.new$W * exp(as.matrix(dat.new[, Znames]) %*% betahat), as.factor(dat.new$stratum), sum)
  mu0 <- num/denom
  mu0_df <- cbind.data.frame("stratum" = unique(dat.new$stratum), "mu0" = mu0[!is.na(mu0)])
  
  return(list("betahat" = betahat, "mu0hat" = mu0_df))
}

IOC <- function(predicted, time, event, dat, modC, method){ # method = "Harrell" or "Uno" (with IPTW)
  if (method == "Uno"){
    weights <- 1/predict(modC, newdata = dat, type = "survival")
  }
  num = 0
  denom = 0
  for (i in 1:nrow(dat)){
    if (method == "Harrell"){
      num <- num + sum(event[i] * ifelse(time[i] < time, 1, 0) * ifelse(predicted[i] < predicted, 1, 0))
      denom <- denom + sum(event[i] * ifelse(time[i] < time, 1, 0))
    }
    else if (method == "Uno"){
      dat.temp <- dat
      dat.temp$X <- time[i]
      weights.temp <- 1/predict(modC, newdata = dat.temp, type = "survival")
      num <- num + sum(event[i] * weights[i] * weights.temp * ifelse(time[i] < time, 1, 0) * ifelse(predicted[i] < predicted, 1, 0))
      denom <- denom + sum(event[i] * weights[i] * weights.temp * ifelse(time[i] < time, 1, 0))
    }
  }
  return(num/denom)
}

brier_rmst <- function(predicted, dat, modC, L){
  dat.new <- dat
  dat.new$X <- dat.new$Y
  weights <- dat$deltaY/predict(modC, newdata = dat.new, type = "survival")
  weights[weights > 50] = 50
  B <- mean((dat$Y - predicted)^2 * weights)
  return(B)
}

## k-fold cross-validation
CV <- function(dat, Znames, ZCnames, L = 5*365, K = 2){
  
  # c-index
  c <- data.frame(matrix(nrow = K, ncol = 4))
  # Brier score
  b <- data.frame(matrix(nrow = K, ncol = 4))
  colnames(c) = colnames(b) = c("StrAdd", "Add", "StrMulti", "Cox")
  
  # data
  dat$Y <- apply(cbind(dat$X, L), 1, min) # Yi = Di^Ci^L
  dat$deltaY <- ifelse(dat$GF == 1, 1, ifelse(L < dat$X, 1, 0)) # deltaY_i = I(Di^L <= Ci)
  dat$stratum <- factor(paste(dat$CTR_CODE, ceiling(dat$AGE/10)*10, sep = "/"))
  dat <- dat[!dat$stratum %in% names(which(tapply(dat$deltaY, dat$stratum, sum) < K)), ]
  
  dat$creat1_dialysis <- dat$creat1 * dat$dialysis
  dat$yr_lt <- dat$yr_lt - 2010
  # donor age group, 1: 12-20, 2: 21-40, 3: 41-60, 4: 61-92
  dat$AGE_DON00 <- ifelse(dat$AGE_DON <= 20, 1, 0)
  dat$AGE_DON20 <- ifelse(dat$AGE_DON > 20 & dat$AGE_DON <= 40, 1, 0) # reference
  dat$AGE_DON40 <- ifelse(dat$AGE_DON > 40 & dat$AGE_DON <= 60, 1, 0)
  dat$AGE_DON60 <- ifelse(dat$AGE_DON > 60, 1, 0)
  dat$don_hgt1 <- (dat$don_hgt-170)/10 # cm
  dat$don_wgt1 <- (dat$don_wgt-80)/10 # kg
  
  # sample within stratum*deltaY (make sure the training set has an event for each stratum; try to use stratified block randomization)
  dat$group <- factor(paste(dat$stratum, dat$deltaY, sep = "/"))
  dat <- dat[sample(nrow(dat)), ]
  dat <- arrange(dat, group)
  
  index <- sapply(table(dat$group), function(x){rep(1:K, length.out = x)})
  dat$index <- unlist(index)
  for (k in 1:K){
    train <- dat[dat$index == k, ]
    test <- dat[dat$index != k, ]
    Z.test <- test[, Znames]
    
    # stratified additive model
    modC <- coxph(as.formula(paste("Surv(X, 1-GF) ~ ", paste0(ZCnames, collapse = " + "), " + strata(stratum)", sep = "")), 
                  data = dat, ties = "breslow")
    est <- estStrAddModel(train, Znames, modC)
    betahat <- est$betahat
    mu0_df <- est$mu0hat
    Yhat <- mu0_df$mu0[match(test$stratum, mu0_df$stratum)] + as.matrix(Z.test) %*% betahat
    c[k,1] <- IOC(Yhat, test$Y, test$deltaY, test, modC, "Harrell")
    b[k,1] <- brier_rmst(Yhat, test, modC, L)/(365^2)
    
    # additive model without stratification (weighted least squares regression)
    modC.nonstrat <- coxph(as.formula(paste("Surv(X, 1-GF) ~ ", paste0(ZCnames, collapse = " + "), sep = "")), 
                           data = train, ties = "breslow")
    weights.nonstrat <- train$deltaY/predict(modC.nonstrat, newdata = train, type = "survival")
    mod.nonstrat <- lm(as.formula(paste("Y ~ ", paste(Znames, collapse = " + "))), data = train, weights = weights.nonstrat)
    Yhat.nonstrat <- predict(mod.nonstrat, newdata = test)
    c[k,2] <- IOC(Yhat.nonstrat, test$Y, test$deltaY, test, modC.nonstrat, "Harrell")
    b[k,2] <- brier_rmst(Yhat.nonstrat, test, modC.nonstrat, L)/(365^2)
    
    # multiplicative model with stratification (Wang et al., 2019)
    est.mult <- estStrMultModel(train, Znames, modC)
    betahat.mult <- est.mult$betahat
    mu0_df.mult <- est.mult$mu0hat
    Yhat.mult <- mu0_df.mult$mu0[match(test$stratum, mu0_df$stratum)] + as.matrix(Z.test) %*% betahat.mult
    c[k,3] <- IOC(Yhat.mult, test$Y, test$deltaY, test, modC.mult, "Harrell")
    b[k,3] <- brier_rmst(Yhat.mult, test, modC, L)/(365^2)
    
    # cox model: calculate RMST or assess hazards directly
    mod.coxph <- coxph(as.formula(paste("Surv(X, GF) ~ ", paste(Znames, collapse = " + "))), data = train)
    c[k,4] <- concordance(mod.coxph, newdata = test, ymin = 0, ymax = L)$concordance
    # Brier score for survival probability: Brier(object = Surv(test$X, test$GF), pre_sp = predict(mod.coxph, newdata = test), t_star = L)
    # to be comparable to the model of RMST, convert to the RMST version of Brier score
    rmst.cox <- survival:::survmean(survfit(mod.coxph, newdata = test), rmean = L)$matrix[, 5]
    b[k,4] <- brier_rmst(rmst.cox, test, modC.nonstrat, L)/(365^2)
  }
  
  return(list("c" = apply(c, 2, mean), "b" = apply(b, 2, mean)))
}

## results
Znames <- c("female", "dialysis", "creat1", "creat1_dialysis", "diabetes", "albumin3", "working_lt", # recipient covariate (age as stratum)
            "diag_HCV", "diag_ahn", "diag_chol_cirr", "diag_mal_neo", "diag_met_dis", # reference: diag_nonchol_cirr
            "don_female", "don_Black", "don_hisp", "don_Asian", # donor covariate
            "AGE_DON00", "AGE_DON40", "AGE_DON60", 
            "don_cod_anoxia", "don_cod_cva", "don_cod_other", 
            "don_DCD", "don_hgt1", "don_wgt1", "don_smoke", "don_coke", 
            "don_creat", "don_partsplit", 
            "abo_mat2", # perfect match between donor and recipient
            "cmv_DposRneg", "cmv_DposRpos", "cmv_DnegRpos") # hypothesis on the effect of DposRneg
ZCnames <- c("female", "dialysis", "creat1", "creat1_dialysis", "diabetes", "albumin3", "working_lt", # recipient covariate (age as stratum)
            "diag_HCV", "diag_ahn", "diag_chol_cirr", "diag_mal_neo", "diag_met_dis", # reference: diag_nonchol_cirr
            "yr_lt", 
            "don_female", "don_Black", "don_hisp", "don_Asian", # donor covariate
            "AGE_DON00", "AGE_DON40", "AGE_DON60", 
            "don_cod_anoxia", "don_cod_cva", "don_cod_other", 
            "don_DCD", "don_hgt1", "don_wgt1", "don_smoke", "don_coke", 
            "don_creat", "don_partsplit", 
            "abo_mat2", # perfect match between donor and recipient
            "cmv_DposRneg", "cmv_DposRpos", "cmv_DnegRpos") # hypothesis on the effect of DposRneg

set.seed(0715)
score_5 <- CV(dat, Znames, L = 5*365, K = 2)
score_3 <- CV(dat, Znames, L = 3*365, K = 2)
score_1 <- CV(dat, Znames, L = 365, K = 2)
diag <- cbind.data.frame("Metric" = rep(c("C", "B"), each = 3), "L" = rep(c(5,3,1), 2), 
                         rbind(score_5$c, score_3$c, score_1$c, 
                               score_5$b, score_3$b, score_1$b))
write.csv(diag, "./results/Table4.csv", row.names = FALSE)
