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
