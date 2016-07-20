#Get the predictions for best MSE
# Use mice on all vars, then use it on only 4 vars (i.e. without var 29). Then pick and choose based on simulations.

library(mice)
library(Metrics)

# read data
mouse.data <- readRDS("mouse_data.rds")

#Missing variables and genotypes for which these variables are missing ####
NAvar = c("IMPC_HEM_027_001", "IMPC_HEM_029_001", "IMPC_HEM_031_001", "IMPC_HEM_034_001", "IMPC_HEM_038_001")
NAgeno = c("3621_1", "4045_1", "1796_1", "3887_1", "3803_1")

#Making sure these two vectors are in the correct order. (i.e. checking if all NA in each pair)
allNA = rep(FALSE, 5)
for (i in 1:5) {
  allNA[i] = all(is.na(mouse.data[mouse.data$geno==NAgeno[i],NAvar[i]]))
}
all(allNA) #SHOULD BE TRUE, OTHERWISE SOMETHING IS WRONG

#Setting up MICE algorithm ####
# create data frame with all variables to impute + predictors
wh_var_29 = which(colnames(mouse.data)=="IMPC_HEM_029_001")
impMICE <- as.data.frame(cbind(mouse.data[,2:3], mouse.data[,5:26]))
impMICE_no29 <- as.data.frame(cbind(mouse.data[,2:3], mouse.data[,c(5:(wh_var_29-1),(wh_var_29+1):26)]))

id = which(colnames(impMICE) %in% NAvar)
method = rep("", NCOL(impMICE))
method[id] = "norm"

id_no29 = which(colnames(impMICE_no29) %in% NAvar)
method_no29 = rep("", NCOL(impMICE_no29))
method_no29[id_no29] <- "norm"

pred <- matrix(0, nrow = ncol(impMICE), ncol = ncol(impMICE))
pred[id,] <- rep(1, ncol(impMICE))
diag(pred) <- rep(0, ncol(impMICE))

pred_no29 = pred[-wh_var_29,-wh_var_29]

#Running MICE ####
MICE <- mice(impMICE, m = 30, method = method, predictorMatrix = pred, printFlag = FALSE)
MICE_no29 <- mice(impMICE_no29, m = 30, method = method_no29, predictorMatrix = pred_no29, printFlag = FALSE)

#Get predictions ####
predic = predic_no29 = list(IMPC_HEM_027_001=NULL,IMPC_HEM_029_001=NULL,IMPC_HEM_031_001=NULL,
                            IMPC_HEM_034_001=NULL,IMPC_HEM_038_001=NULL)
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
}

#Now, given the predictions for the other variables, predict var 29
lm_mod = lm(IMPC_HEM_029_001~., data=imputed_mat[!imputed_mat$geno %in% NAgeno,-2])
predic_no29[[2]] = predict(lm_mod, newdata=imputed_mat[imputed_mat$geno==NAgeno[2],])

#Picking and choosing which predictions come from which prediciton method based on simulation results
predictions = list(IMPC_HEM_027_001=predic$IMPC_HEM_027_001,
                   IMPC_HEM_029_001=predic_no29$IMPC_HEM_029_001,
                   IMPC_HEM_031_001=predic$IMPC_HEM_031_001,
                   IMPC_HEM_034_001=predic_no29$IMPC_HEM_034_001,
                   IMPC_HEM_038_001=predic_no29$IMPC_HEM_038_001)

#Final imputed dataset
mouse.data.imp = mouse.data
for (i in 1:5) {
  mouse.data.imp[mouse.data.imp$geno==NAgeno[i],NAvar[i]] = predictions[[i]]
}
any(is.na(mouse.data.imp)) #SHOULD BE FALSE (i.e. no NAs present)

save(predictions, mouse.data.imp, file="~/Documents/research/causal_challenge_repo/Results/FINAL_PREDICTIONS.RData")






