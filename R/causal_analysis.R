library(rms)
library(mice)
set.seed(213)
# causal interpretation genotype-phenotype

mouse.data <- readRDS("mouse_data.RDS")
var.names <- readRDS("variable_names.RDS")
mouse.data <- mouse.data[,-1]
colnames(mouse.data)[4:25] <- var.names$nam

## create a dataframe with genotype-phenotype association for complete variables 
#     13 columns, one for each knockout conditions (genotype)
#     22 rows, one for each phenotypic measurement
#     each element report the p-value for association if pval <=0.05/22, reports 1 otherwise

full_mouse <- mouse.data[,-c(13,14,16,19,23)]
geno <- full_mouse$geno
litter <- full_mouse$litter

geno_pheno1 <- matrix(NA, ncol=13, nrow=17)
rownames(geno_pheno1) <- colnames(full_mouse)[4:20] 

for(i in 4:20)
{
  modi <- ols(full_mouse[,i] ~ geno, x=TRUE)
  adj <- robcov(modi, cluster=litter)
  pval <- 2*(1-pt(abs(modi$coefficients/sqrt(diag(adj$var))), df=613))[-1]
  geno_pheno1[i-3,] <- ifelse(pval<=0.05/22, pval, 1)
}

## get pooled coefficients from MICE
imp <- mouse.data

# id for variables with missing values
id <- c(13, 14, 16, 19, 23)
method <- rep("", 25)
method[id] <- "norm"

# pred matrix
pred <- matrix(0, nrow = ncol(imp), ncol = ncol(imp))
pred[id,] <- rep(1, ncol(imp))
diag(pred) <- rep(0, ncol(imp))
pred[,3] <- 0 # dont use litter

# apply MICE
imp_res <- mice(imp, m=30, method=method, predictorMatrix=pred, printFlag=FALSE)

# create a list, each element = 1 imputed dataset

imputed_data <- vector(mode="list", length=30)

id_13 <- as.numeric(rownames(imp)[which(is.na(imp$`Red blood cell distribution width`)==TRUE)]) # id with missing values for each 5 variables
id_14 <- as.numeric(rownames(imp)[which(is.na(imp$`Neutrophil differential count`)==TRUE)])
id_16 <- as.numeric(rownames(imp)[which(is.na(imp$`Lymphocyte differential count`)==TRUE)])
id_19 <- as.numeric(rownames(imp)[which(is.na(imp$`Monocyte cell count`)==TRUE)])
id_23 <- as.numeric(rownames(imp)[which(is.na(imp$`Basophil differential count`)==TRUE)])

for(i in 1:30)
{
  impute <- imp
  impute[id_13,13] <- as.vector(unlist(imp_res$imp$`Red blood cell distribution width`[i]))
  impute[id_14,14] <- as.vector(unlist(imp_res$imp$`Neutrophil differential count`[i]))
  impute[id_16,16] <- as.vector(unlist(imp_res$imp$`Lymphocyte differential count`[i]))
  impute[id_19,19] <- as.vector(unlist(imp_res$imp$`Monocyte cell count`[i]))
  impute[id_23,23] <- as.vector(unlist(imp_res$imp$`Basophil differential count`[i]))
  imputed_data[[i]] <- impute
}

# run 5 ols, one phenotypic measurement at the time, only for phenotypic measurement with missing values!!
geno <- imp$geno
litter <- imp$litter
geno_pheno2 <- matrix(NA, ncol=13, nrow=5)
rownames(geno_pheno2) <- colnames(mouse.data)[c(13,14,16,19,23)] 

# id 13 "Red blood cell distribution width" 
W <- matrix(NA, ncol=30, nrow=14)
theta <- matrix(NA, nrow=14, ncol=30)

for(i in 1:30) # loop over imputed dataset
{
  olsi <- ols(as.vector(imputed_data[[i]][,13]) ~ geno, x=TRUE)
  theta[,i] <- olsi$coefficients
  W[,i] <- diag(vcov(robcov(olsi, cluster=litter)))
}

thetaO <- apply(theta, 1, mean) # estimate itself doesnt change a lot, SEs change more
WO <- apply(W, 1, mean)
B <- sum(apply(theta, 2, function(x) (x-thetaO)^2))/29
var_thetaO <- WO + (1+1/29)*B
pval <- 2*(1-pt(abs(thetaO/sqrt(var_thetaO)), df=613))[-1]
geno_pheno2[1,] <- ifelse(pval <= 0.05/22, pval, 1)


# id 14 "Neutrophil differential count" 
W <- matrix(NA, ncol=30, nrow=14)
theta <- matrix(NA, nrow=14, ncol=30)

for(i in 1:30) # loop over imputed dataset
{
  olsi <- ols(as.matrix(imputed_data[[i]][,14]) ~ geno, x=TRUE)
  theta[,i] <- olsi$coefficients
  W[,i] <- diag(vcov(robcov(olsi, cluster=litter)))
}

thetaO <- apply(theta, 1, mean) # estimate itself doesnt change a lot, SEs change more
WO <- apply(W, 1, mean)
B <- sum(apply(theta, 2, function(x) (x-thetaO)^2))/29
var_thetaO <- WO + (1+1/29)*B
pval <- 2*(1-pt(abs(thetaO/sqrt(var_thetaO)), df=613))[-1]
geno_pheno2[2,] <- ifelse(pval <= 0.05/22, pval, 1)



# id 16 "Lymphocyte differential count" 
W <- matrix(NA, ncol=30, nrow=14)
theta <- matrix(NA, nrow=14, ncol=30)

for(i in 1:30) # loop over imputed dataset
{
  olsi <- ols(as.matrix(imputed_data[[i]][,16]) ~ geno, x=TRUE)
  theta[,i] <- olsi$coefficients
  W[,i] <- diag(vcov(robcov(olsi, cluster=litter)))
}

thetaO <- apply(theta, 1, mean) # estimate itself doesnt change a lot, SEs change more
WO <- apply(W, 1, mean)
B <- sum(apply(theta, 2, function(x) (x-thetaO)^2))/29
pval <- 2*(1-pt(abs(thetaO/sqrt(var_thetaO)), df=613))[-1]
geno_pheno2[3,] <- ifelse(pval <= 0.05/22, pval, 1)




# id 19 "Monocyte cell count" 
W <- matrix(NA, ncol=30, nrow=14)
theta <- matrix(NA, nrow=14, ncol=30)

for(i in 1:30) # loop over imputed dataset
{
  olsi <- ols(as.matrix(imputed_data[[i]][,19]) ~ geno, x=TRUE)
  theta[,i] <- olsi$coefficients
  W[,i] <- diag(vcov(robcov(olsi, cluster=litter)))
}

thetaO <- apply(theta, 1, mean) # estimate itself doesnt change a lot, SEs change more
WO <- apply(W, 1, mean)
B <- sum(apply(theta, 2, function(x) (x-thetaO)^2))/29
var_thetaO <- WO + (1+1/29)*B
pval <- 2*(1-pt(abs(thetaO/sqrt(var_thetaO)), df=613))[-1]
geno_pheno2[4,] <- ifelse(pval <= 0.05/22, pval, 1)




# id 23 "Basophil differential count" 
W <- matrix(NA, ncol=30, nrow=14)
theta <- matrix(NA, nrow=14, ncol=30)

for(i in 1:30) # loop over imputed dataset
{
  olsi <- ols(as.matrix(imputed_data[[i]][,23]) ~ geno, x=TRUE)
  theta[,i] <- olsi$coefficients
  W[,i] <- diag(vcov(robcov(olsi, cluster=litter)))
}

thetaO <- apply(theta, 1, mean) # estimate itself doesnt change a lot, SEs change more
WO <- apply(W, 1, mean)
B <- sum(apply(theta, 2, function(x) (x-thetaO)^2))/29
var_thetaO <- WO + (1+1/29)*B
pval <- 2*(1-pt(abs(thetaO/sqrt(var_thetaO)), df=613))[-1]
geno_pheno2[5,] <- ifelse(pval <= 0.05/22, pval, 1)


## create list with phenotypic measurement associated to each genotype -> quick view
#     each element of the list represents one genotype
#     lists the associated phenotypes

geno_pheno <- rbind(geno_pheno1, geno_pheno2)
colnames(geno_pheno) <- c("geno=1550_1", "geno=1796_1", "geno=1797_1", "geno=1798_1", "geno=1799_1", "geno=3157_1", "geno=3621_1", "geno=3803_1", "geno=3805_1", "geno=3887_1", "geno=4045_1", "geno=4047_1", "geno=727_1")

quick <- vector(mode="list", length=13) # change 8 to 13 when we have predictions
names(quick) <- colnames(geno_pheno)

for(i in 1:13) # change 8 to 13
{
  pvalues <- geno_pheno[,i]
  quick[[i]] <- rownames(geno_pheno)[which(pvalues != 1)]
}

