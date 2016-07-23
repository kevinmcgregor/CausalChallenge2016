#Final simulation (I hope)
#Trying regular MICE algorithm, except without variable 029. This variable will be predicted by itself
#afterwards with OLS

library(mice)
library(Metrics)

# read data
mouse.data <- readRDS("mouse_data.rds")
var.names = readRDS("variable_names.RDS")

# consider only complete observations: delete genotypes with missing values
full_mouse <- mouse.data[-which(mouse.data$geno %in% c("1796_1", "3621_1", "4045_1", "3803_1", "3887_1")),]

Nsimul = 100
res <- matrix(0, Nsimul, 5)
res_no29 <- matrix(0, Nsimul, 5)
colnames(res) = colnames(res_no29) = c("IMPC_HEM_027_001", "IMPC_HEM_029_001", "IMPC_HEM_031_001", "IMPC_HEM_034_001", "IMPC_HEM_038_001")
sim_geno = matrix("", Nsimul, 5)
wh_var_29 = which(colnames(full_mouse)=="IMPC_HEM_029_001")

# consider only complete observations: delete genotypes with missing values
full_mouse <- mouse.data[-which(mouse.data$geno %in% c("1796_1", "3621_1", "4045_1", "3803_1", "3887_1")),]

for(i in 1:Nsimul)
{
  # temporary dataframe
  temp <- full_mouse
  
  # empty method vector
  method <- c("", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "")
  method_no29 <- c("", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "")
  
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
  
  # create data frame with all variables to impute + predictors
  impMICE <- as.data.frame(cbind(temp[,2:3], temp[,5:26]))
  impMICE_no29 <- as.data.frame(cbind(temp[,2:3], temp[,c(5:(wh_var_29-1),(wh_var_29+1):26)]))
  
  # column ids of variables with missing values
  id <- which(colnames(impMICE) %in% NAvar)
  id_no29 = which(colnames(impMICE_no29) %in% NAvar)
  
  # method norm + pmm
  method[id] <- "norm"
  method_no29[id_no29] <- "norm"
  
  # matrix to specify which preditors to use 
  pred <- matrix(0, nrow = ncol(impMICE), ncol = ncol(impMICE))
  pred[id,] <- rep(1, ncol(impMICE))
  diag(pred) <- rep(0, ncol(impMICE))
  
  pred_no29 = pred[-wh_var_29,-wh_var_29]
  
  # apply MICE
  MICE <- mice(impMICE, m = 30, method = method, predictorMatrix = pred, printFlag = FALSE)
  MICE_no29 <- mice(impMICE_no29, m = 30, method = method_no29, predictorMatrix = pred_no29, printFlag = FALSE)
  
  # get predictions and compute mse
  predic = predic_no29 = list(IMPC_HEM_027_001=NULL,IMPC_HEM_029_001=NULL,IMPC_HEM_031_001=NULL,
                                  IMPC_HEM_034_001=NULL,IMPC_HEM_038_001=NULL)
  mse = rep(0,5)
  mse_no29 = rep(0,5)
  imputed_mat = impMICE
  for (j in 1:5) {
    predic[[j]] = apply(MICE$imp[id[j]][[1]], 1, mean)
    
    #Only have predictions for 4 vars for now
    if (j!=2) {
      k = ifelse(j==1,j,j-1)
      predic_no29[[j]] = apply(MICE_no29$imp[id_no29[k]][[1]], 1, mean)
      
      #Put imputed values back into data frame
      imputed_mat[imputed_mat$geno==NAgeno[j],id[j]] = predic_no29[[j]]
    }
    
    mse[j] <- mse(actual = true[[j]], predicted = predic[[j]])
  }
  
  #Now, given the predictions for the other variables, predict var 29
  lm_mod = lm(IMPC_HEM_029_001~., data=imputed_mat[!imputed_mat$geno %in% NAgeno,-2])
  predic_no29[[2]] = predict(lm_mod, newdata=imputed_mat[imputed_mat$geno==NAgeno[2],])
  
  for (j in 1:5) {
    mse_no29[j] = mse(actual = true[[j]], predicted = predic_no29[[j]])
  }
  
  # save results
  res[i, ] <- mse
  res_no29[i, ] <- mse_no29
  
  cat(i, "\n")
}

res_replace29 = cbind(res[,1],res_no29[,2],res[,3:5])

save(res, res_no29, res_replace29, sim_geno, file = "simulation_final.RData")

boxplot(res, main = "Imputation via MICE")
boxplot(res_no29, main = "Imputation via MICE , var29 alone")
boxplot(res_replace29, main = "Imputation via mice, \n var29 replaced in original")






