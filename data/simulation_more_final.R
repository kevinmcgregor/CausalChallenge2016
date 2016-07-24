#Another 'final' simulation with some of the predictions coming from directly applying formulas

library(mice)
library(Metrics)

# read data
mouse.data <- readRDS("mouse_data.rds")
var.names = readRDS("variable_names.RDS")

# consider only complete observations: delete genotypes with missing values
full_mouse <- mouse.data[-which(mouse.data$geno %in% c("1796_1", "3621_1", "4045_1", "3803_1", "3887_1")),]

Nsimul = 100
res_direct <- matrix(0, Nsimul, 5)
colnames(res_direct) =  c("IMPC_HEM_027_001", "IMPC_HEM_029_001", "IMPC_HEM_031_001", "IMPC_HEM_034_001", "IMPC_HEM_038_001")
sim_geno = matrix("", Nsimul, 5)
wh_var_27 = which(colnames(full_mouse)=="IMPC_HEM_027_001")
wh_var_29 = which(colnames(full_mouse)=="IMPC_HEM_029_001")
wh_var_31 =  which(colnames(full_mouse)=="IMPC_HEM_031_001")

# consider only complete observations: delete genotypes with missing values
full_mouse <- mouse.data[-which(mouse.data$geno %in% c("1796_1", "3621_1", "4045_1", "3803_1", "3887_1")),]

for(i in 1:Nsimul)
{
  # temporary dataframe
  temp <- full_mouse
  
  # empty method vector
  method <- c("", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "")
  
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
  
  sim_geno[i,] = NAgeno
  
  #Predict lympho diff before MICE
  predic = predic_no29 = list(IMPC_HEM_027_001=NULL,IMPC_HEM_029_001=NULL,IMPC_HEM_031_001=NULL,
                              IMPC_HEM_034_001=NULL,IMPC_HEM_038_001=NULL)
  temp$IMPC_HEM_031_001[temp$geno==NAgeno[3]] = 
            temp$IMPC_HEM_032_001[temp$geno==NAgeno[3]]/temp$IMPC_HEM_001_001[temp$geno==NAgeno[3]]*100
  predic[[3]] = temp$IMPC_HEM_031_001[temp$geno==NAgeno[3]]
  
  #Predict RBCDW before MICE
  sdcv = NA
  sdcv = temp$IMPC_HEM_005_001*temp$IMPC_HEM_027_001/100
  lm_sdcv = lm(log(sdcv)~.-IMPC_HEM_005_001-IMPC_HEM_027_001, data=temp[,-c(1,3,4)])
  pred_sdcv = exp(predict(lm_sdcv, newdata = temp[temp$geno==NAgeno[1],]))
  temp$IMPC_HEM_027_001[temp$geno==NAgeno[1]] = pred_sdcv/temp$IMPC_HEM_005_001[temp$geno==NAgeno[1]] * 100
  predic[[1]] = temp$IMPC_HEM_027_001[temp$geno==NAgeno[1]]
  
  # create data frame with all variables to impute + predictors
  impMICE <- as.data.frame(cbind(temp[,2:3], temp[,5:26]))
  
  # column ids of variables with missing values
  id <- which(colnames(impMICE) %in% NAvar[-c(1,3)])
  
  # method norm + pmm
  method[id] <- "norm"
  
  # matrix to specify which preditors to use 
  pred <- matrix(0, nrow = ncol(impMICE), ncol = ncol(impMICE))
  pred[id,] <- rep(1, ncol(impMICE))
  diag(pred) <- rep(0, ncol(impMICE))
  
  # apply MICE
  MICE <- mice(impMICE, m = 30, method = method, predictorMatrix = pred, printFlag = FALSE)
  
  # get predictions and compute mse
  mse = rep(0,5)
  mse_no29 = rep(0,5)
  imputed_mat = impMICE
  k=1
  for (j in 1:5) {
    if (j %in% c(2,4,5)) {
      predic[[j]] = apply(MICE$imp[id[k]][[1]], 1, mean)
      #Put imputed values back into data frame
      imputed_mat[imputed_mat$geno==NAgeno[j],id[k]] = predic[[j]]
      k=k+1
    }
    mse[j] <- mse(actual = true[[j]], predicted = predic[[j]])
  }
  
  #Now, given the predictions for the other variables, predict var 29
  lm_mod = lm(IMPC_HEM_029_001~., data=imputed_mat[!imputed_mat$geno %in% NAgeno,-2])
  predic[[2]] = predict(lm_mod, newdata=imputed_mat[imputed_mat$geno==NAgeno[2],])
  
  mse[2] = mse(actual = true[[2]], predicted = predic[[2]])
  
  # save results
  res_direct[i, ] <- mse
  
  cat(i, "\n")
}

#Comparing to "final" simulation
load("simulation_final.RData")

#By variable
colMeans(res_replace29)
colMeans(res_direct)
#Overall
mean(res_replace29)
mean(res_direct)







