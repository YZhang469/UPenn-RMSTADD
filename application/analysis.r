library(haven)
library(survival)
library(ggplot2)

source("../function/estBeta.r")

## data cleaning
dat <- read_sas("postlt_yuan_15jul2022.sas7bdat")
dat <- dat[dat$status_1 == 0 & dat$AGE_DON >= 12, ]
dat <- dat[!is.na(dat$don_bili), ]
# define stratum = age * center
dat$stratum <- factor(paste(dat$CTR_CODE, ceiling(dat$AGE/10)*10, sep = "/"))
# recode covariates
dat$creat1_dialysis <- dat$creat1 * dat$dialysis
dat$AGE_DON00 <- ifelse(dat$AGE_DON <= 20, 1, 0)
dat$AGE_DON20 <- ifelse(dat$AGE_DON > 20 & dat$AGE_DON <= 40, 1, 0) # reference
dat$AGE_DON40 <- ifelse(dat$AGE_DON > 40 & dat$AGE_DON <= 60, 1, 0)
dat$AGE_DON60 <- ifelse(dat$AGE_DON > 60, 1, 0)
dat$don_hgt1 <- (dat$don_hgt-170)/10 # cm
dat$don_wgt1 <- (dat$don_wgt-80)/10 # kg

Znames <- c("female", "dialysis", "creat1", "creat1_dialysis", "diabetes", "albumin3", "working_lt", # recipient covariate (age as stratum)
            "diag_HCV", "diag_ahn", "diag_chol_cirr", "diag_mal_neo", "diag_met_dis", # reference: diag_nonchol_cirr
            "don_female", "don_Black", "don_hisp", "don_Asian", # donor covariate
            "AGE_DON00", "AGE_DON40", "AGE_DON60", 
            "don_cod_anoxia", "don_cod_cva", "don_cod_other", 
            "don_DCD", "don_hgt1", "don_wgt1", "don_smoke", "don_coke", 
            "don_creat", "don_partsplit", 
            "abo_mat2", # perfect match between donor and recipient
            "cmv_DposRneg", "cmv_DposRpos", "cmv_DnegRpos") # hypothesis on the effect of DposRneg

est <- estBeta(dat = dat, Xname = "X", deltaXname = "GF", Znames = Znames, strname = "stratum", L = 5*365)
betahat <- est$betahat
se <- sqrt(est$var)
res <- cbind.data.frame("est" = betahat, "se" = se, "p" = pnorm(-abs(betahat/se), mean = 0, sd = 1, lower.tail = TRUE) * 2)
write.csv(res, "res.csv", row.names = TRUE)

## histogram of mu0
png(file = "hist_mu0.png", width = 12, height = 10, units = "cm", res = 600)
hist_mu0 <- hist(mu0/365, breaks = seq(floor(min(mu0/365)), ceiling(max(mu0/365)), 0.2))
hist_mu0$density <- hist_mu0$counts/sum(hist_mu0$counts)*100
plot(hist_mu0, freq = FALSE, ylab = "Proportion (%)", xlab = expression(paste(hat(mu)[0], " (years)")), main = "", 
     xlim = c(0, 5), col = "lightblue")
dev.off()
mu0_df$age <- sub(".*/", "", mu0_df$stratum)
mu0_df$age <- paste(as.numeric(mu0_df$age) - 9, mu0_df$age, sep = " - ")
mu0_df$center <- sub("/.*", "", mu0_df$stratum)
# by age
png(file = "hist_mu0_by_age.png", width = 12, height = 18, units = "cm", res = 600)
ggplot(mu0_df[!mu0_df$age %in% c("11 - 20", "81 - 90"), ], aes(x = mu0/365)) + 
  geom_histogram(color = "black", fill = "lightblue") + 
  facet_wrap(~ age, ncol = 2) + 
  labs(y = "Count", x = expression(paste(hat(mu)[0], " (years)")), fill = "Age")
dev.off()
