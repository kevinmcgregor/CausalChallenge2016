#Simulation as before, except missing data from all 5 phenotypes at once

library(mice)
library(Metrics)

# read data
mouse.data <- readRDS("mouse_data.rds")

set.seed(30062016)
Nsimul <- 100

# empty matrix to save result
res_norm <- matrix(0, Nsimul, 5)
res_pmm <- matrix(0, Nsimul, 5)
colnames(res_norm) = colnames(res_pmm) = c("IMPC_HEM_027_001", "IMPC_HEM_029_001", "IMPC_HEM_031_001", "IMPC_HEM_034_001", "IMPC_HEM_038_001")

# consider only complete observations: delete genotypes with missing values
full_mouse <- mouse.data[-which(mouse.data$geno %in% c("1796_1", "3621_1", "4045_1", "3803_1", "3887_1")),]

for(i in 1:Nsimul)
{
  # temporary dataframe
  temp <- full_mouse
  
  # empty method vector
  method_norm <- c("", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "")
  method_pmm <- c("", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "")
  
  # randomly select phenotypes to be deleted
  NAgeno <- sample(unique(temp$geno)[2:9], 5)
  NAvar <- c("IMPC_HEM_027_001", "IMPC_HEM_029_001", "IMPC_HEM_031_001", "IMPC_HEM_034_001", "IMPC_HEM_038_001") 
  
  # save true values + delete in temp
  true = list(IMPC_HEM_027_001=temp[temp$geno==NAgeno[1],NAvar[1]],
              IMPC_HEM_029_001=temp[temp$geno==NAgeno[2],NAvar[2]],
              IMPC_HEM_031_001=temp[temp$geno==NAgeno[3],NAvar[3]],
              IMPC_HEM_034_001=temp[temp$geno==NAgeno[4],NAvar[4]],
              IMPC_HEM_038_001=temp[temp$geno==NAgeno[5],NAvar[5]])
  for (j in 1:5) {
    temp[temp$geno == NAgeno[j],  colnames(temp) == NAvar[j]] <- NA
  }
  
  # create data frame with all variables to impute + predictors
  impMICE <- as.data.frame(cbind(temp[,2:3], temp[,5:26]))
  
  # column ids of variables with missing values
  id <- which(colnames(impMICE) %in% NAvar)
  
  # method norm + pmm
  method_norm[id] <- "norm"
  method_pmm[id] <- "pmm"
  
  # matrix to specify which preditors to use 
  pred <- matrix(0, nrow = ncol(impMICE), ncol = ncol(impMICE))
  pred[id,] <- rep(1, ncol(impMICE))
  diag(pred) <- rep(0, ncol(impMICE))
  
  # apply MICE
  MICE_norm <- mice(impMICE, m = 30, method = method_norm, predictorMatrix = pred, printFlag = FALSE)
  MICE_pmm <- mice(impMICE, m = 30, method = method_pmm, predictorMatrix = pred, printFlag = FALSE)
  
  # get predictions and compute mse
  predic_norm = predic_pmm = list(IMPC_HEM_027_001=NULL,IMPC_HEM_029_001=NULL,IMPC_HEM_031_001=NULL,
                                  IMPC_HEM_034_001=NULL,IMPC_HEM_038_001=NULL)
  mse_norm = rep(0,5)
  mse_pmm = rep(0,5)
  for (j in 1:5) {
    predic_norm[[j]] = apply(MICE_norm$imp[id[j]][[1]], 1, mean)
    predic_pmm[[j]] = apply(MICE_pmm$imp[id[j]][[1]], 1, mean)

    mse_norm[j] <- mse(actual = true[[j]], predicted = predic_norm[[j]])
    mse_pmm[j] <- mse(actual = true[[j]], predicted = predic_pmm[[j]])
  }
  
  # save results
  res_norm[i, ] <- mse_norm
  res_pmm[i, ] <- mse_pmm
  
  cat(i, "\n")
}

save(res_norm, res_pmm, file = "simulation_all.RData")

par(mfrow = c(1, 2))
boxplot(res_norm, main = "Imputation via OLS")
boxplot(res_pmm, main = "Imputation via PMM", ylim=c(0,2))



#Trying again with mixed combination of OLS and PMM
set.seed(30062016)

res_mixed <- matrix(0, Nsimul, 5)
colnames(res_mixed) = c("IMPC_HEM_027_001", "IMPC_HEM_029_001", "IMPC_HEM_031_001", "IMPC_HEM_034_001", "IMPC_HEM_038_001")
sim_geno = matrix("", NSimul, 5)

for(i in 1:Nsimul)
{
  # temporary dataframe
  temp <- full_mouse
  
  # empty method vector
  method_mixed <- c("", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "")
  
  # randomly select phenotypes to be deleted
  NAgeno <- sample(unique(temp$geno)[2:9], 5)
  NAvar <- c("IMPC_HEM_027_001", "IMPC_HEM_029_001", "IMPC_HEM_031_001", "IMPC_HEM_034_001", "IMPC_HEM_038_001") 
  
  sim_geno[i,] = NAgeno
  
  # save true values + delete in temp
  true = list(IMPC_HEM_027_001=temp[temp$geno==NAgeno[1],NAvar[1]],
              IMPC_HEM_029_001=temp[temp$geno==NAgeno[2],NAvar[2]],
              IMPC_HEM_031_001=temp[temp$geno==NAgeno[3],NAvar[3]],
              IMPC_HEM_034_001=temp[temp$geno==NAgeno[4],NAvar[4]],
              IMPC_HEM_038_001=temp[temp$geno==NAgeno[5],NAvar[5]])
  for (j in 1:5) {
    temp[temp$geno == NAgeno[j],  colnames(temp) == NAvar[j]] <- NA
  }
  
  # create data frame with all variables to impute + predictors
  impMICE <- as.data.frame(cbind(temp[,2:3], temp[,5:26]))
  
  # column ids of variables with missing values
  id <- which(colnames(impMICE) %in% NAvar)
  
  # method norm + pmm
  method_mixed[id[c(1:3,5)]] <- "norm"
  method_mixed[id[4]] <- "pmm"
  
  # matrix to specify which preditors to use 
  pred <- matrix(0, nrow = ncol(impMICE), ncol = ncol(impMICE))
  pred[id,] <- rep(1, ncol(impMICE))
  diag(pred) <- rep(0, ncol(impMICE))
  
  # apply MICE
  MICE_mixed <- mice(impMICE, m = 30, method = method_mixed, predictorMatrix = pred, printFlag = FALSE)
  
  # get predictions and compute mse
  predic_mixed = list(IMPC_HEM_027_001=NULL,IMPC_HEM_029_001=NULL,IMPC_HEM_031_001=NULL,
                                  IMPC_HEM_034_001=NULL,IMPC_HEM_038_001=NULL)
  mse_mixed = rep(0,5)
  for (j in 1:5) {
    predic_mixed[[j]] = apply(MICE_mixed$imp[id[j]][[1]], 1, mean)
    
    mse_mixed[j] <- mse(actual = true[[j]], predicted = predic_mixed[[j]])
  }
  
  # save results
  res_mixed[i, ] <- mse_mixed
  
  cat(i, "\n")
}


save(res_mixed, sim_geno, file = "simulation_mixed.RData")

boxplot(res_mixed, main = "Imputation via OLS")

#Trying out splines
library(splines)

#Using splines on important variables for IMPC_HEM_029_001 and comparing to OLS
#find important vars
summary(lm(IMPC_HEM_029_001~.-id-geno-litter, data=full_mouse))
plot(full_mouse$IMPC_HEM_031_001, full_mouse$IMPC_HEM_029_001)
plot(full_mouse$IMPC_HEM_033_001, full_mouse$IMPC_HEM_029_001)
plot(full_mouse$IMPC_HEM_035_001, full_mouse$IMPC_HEM_029_001)
plot(full_mouse$IMPC_HEM_038_001, full_mouse$IMPC_HEM_029_001)
plot(full_mouse$IMPC_HEM_040_001, full_mouse$IMPC_HEM_029_001)

mod_31 = lm(IMPC_HEM_029_001~ns(IMPC_HEM_031_001,df=3)+
                             ns(IMPC_HEM_033_001,df=3)+
                             ns(IMPC_HEM_035_001,df=3)+
                             ns(IMPC_HEM_038_001,df=3)+
                             ns(IMPC_HEM_040_001,df=3), data=full_mouse)
pred_31 = predict(mod_31, type="response")
mse_31 = mse(pred_31, full_mouse$IMPC_HEM_029_001)

mod_31_nosplines = lm(IMPC_HEM_029_001~.-id-litter-geno, data=full_mouse)
pred_31_nosplines = predict(mod_31_nosplines, type="response")
mse_31_nosplines = mse(pred_31_nosplines, full_mouse$IMPC_HEM_029_001)


#Simulation with splines for IMPC_HEM_029_001
Nsimul=500
set.seed(30062016)

res_splines <- matrix(0, Nsimul, 1)
res_nosplines <- matrix(0, Nsimul, 1)
colnames(res_splines) = colnames(res_nosplines) = "IMPC_HEM_029_001"
#New data frame with splines basis vectors
splines_mouse = data.frame(full_mouse, ns(full_mouse$IMPC_HEM_031_001,df=3),
                                       ns(full_mouse$IMPC_HEM_033_001,df=3),
                                       ns(full_mouse$IMPC_HEM_035_001,df=3),
                                       ns(full_mouse$IMPC_HEM_038_001,df=3),
                                       ns(full_mouse$IMPC_HEM_040_001,df=3))
colnames(splines_mouse)[27:41] = paste0("S",1:15)
for(i in 1:Nsimul)
{
  # temporary dataframe
  temp <- splines_mouse
  
  # empty method vector
  method_splines <- rep("", NCOL(temp)-2)
  
  # randomly select phenotypes to be deleted
  NAgeno <- sample(unique(temp$geno)[2:9], 1)
  NAvar <- c("IMPC_HEM_029_001") 
  
  # save true values + delete in temp
  true = list(IMPC_HEM_029_001=temp[temp$geno==NAgeno,NAvar])
  
  temp[temp$geno == NAgeno,  colnames(temp) == NAvar] <- NA
  
  #Data frames with just splines
  splines_alone = temp[,27:41]
  splines_alone = cbind(splines_alone, IMPC_HEM_029_001=temp$IMPC_HEM_029_001)
  
  nosplines_data = temp[,-(27:41)]
  
  # column ids of variables with missing values
  id <- which(colnames(temp) %in% NAvar)
  
  #Linear models
  lm_splines = lm(IMPC_HEM_029_001~., data=splines_alone)
  lm_nosplines = lm(temp$IMPC_HEM_029_001~., data=nosplines_data[,-c(1,3,4)])
  
  pred_splines = predict(lm_splines, newdata=splines_alone[is.na(temp[,id]),])
  pred_nosplines = predict(lm_nosplines, newdata=temp[is.na(nosplines_data[,id]),])
  
  # save results
  res_splines[i, ] <- mse(true[[1]], pred_splines)
  res_nosplines[i, ] <- mse(true[[1]], pred_nosplines)
  
  cat(i, "\n")
}
























