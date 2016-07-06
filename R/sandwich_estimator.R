#Trying out the sandwich estimator (vs OLS)
library(rms)

mouse.data <- readRDS("~/Documents/research/causal_challenge_repo/data/mouse_data.RDS")
var.names <- readRDS("~/Documents/research/causal_challenge_repo/data/variable_names.RDS")

#Find variables with missing data  ####
vars.missing = sapply(mouse.data, function(x){any(is.na(x))})

# Regressing 1st phenotype on genotype ####
ols1 = lm(IMPC_HEM_001_001~as.factor(geno), data=mouse.data)
#Trying with sandwich SEs to take litter into account
mouse.data$geno = factor(mouse.data$geno)
rms_ols1 = ols(IMPC_HEM_001_001~geno, data=mouse.data, x=TRUE)
rc1 = robcov(rms_ols1, cluster=mouse.data$litter)









