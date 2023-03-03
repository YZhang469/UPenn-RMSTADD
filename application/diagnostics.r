library(haven)
library(survival)
library(ggplot2)
library(dplyr)
library(survRM2)

source("../function/estBeta.r")

## data cleaning
dat <- read_sas("postlt_yuan_15jul2022.sas7bdat")
dat <- dat[dat$status_1 == 0 & dat$AGE_DON >= 12, ]
dat <- dat[!is.na(dat$don_bili), ]
