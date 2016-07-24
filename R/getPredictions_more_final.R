#Get (more) final predictions (better than the previous "final" predictions)

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

#Predictions!
mouse.data.imp = mouse.data

#Predict lympho diff before MICE
predic = predic_no29 = list(IMPC_HEM_027_001=NULL,IMPC_HEM_029_001=NULL,IMPC_HEM_031_001=NULL,
                            IMPC_HEM_034_001=NULL,IMPC_HEM_038_001=NULL)
mouse.data.imp$IMPC_HEM_031_001[mouse.data$geno==NAgeno[3]] = 
      mouse.data$IMPC_HEM_032_001[mouse.data$geno==NAgeno[3]]/mouse.data$IMPC_HEM_001_001[mouse.data$geno==NAgeno[3]]*100
predic[[3]] = mouse.data.imp$IMPC_HEM_031_001[mouse.data$geno==NAgeno[3]]

#Predict RBCDW before MICE
sdcv = NA
sdcv = mouse.data$IMPC_HEM_005_001*mouse.data$IMPC_HEM_027_001/100
lm_sdcv = lm(log(sdcv)~.-IMPC_HEM_005_001-IMPC_HEM_027_001, data=mouse.data[,-c(1,3,4)])
pred_sdcv = exp(predict(lm_sdcv, newdata = mouse.data[mouse.data$geno==NAgeno[1],]))
mouse.data.imp$IMPC_HEM_027_001[mouse.data.imp$geno==NAgeno[1]] = 
                  pred_sdcv/mouse.data$IMPC_HEM_005_001[mouse.data$geno==NAgeno[1]] * 100
predic[[1]] = mouse.data.imp$IMPC_HEM_027_001[mouse.data.imp$geno==NAgeno[1]]


impMICE <- as.data.frame(cbind(mouse.data[,2:3], mouse.data[,5:26]))

# column ids of variables with missing values
id <- which(colnames(impMICE) %in% NAvar[-c(1,3)])
id_imp <- which(colnames(mouse.data.imp) %in% NAvar[-c(1,3)])

# method norm + pmm
method <- c("", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "")
method[id] <- "norm"

# matrix to specify which preditors to use 
pred <- matrix(0, nrow = ncol(impMICE), ncol = ncol(impMICE))
pred[id,] <- rep(1, ncol(impMICE))
diag(pred) <- rep(0, ncol(impMICE))

# apply MICE
MICE <- mice(impMICE, m = 30, method = method, predictorMatrix = pred, printFlag = FALSE)

k=1
for (j in c(2,4,5)) {
  predic[[j]] = apply(MICE$imp[id[k]][[1]], 1, mean)
  #Put imputed values back into data frame
  mouse.data.imp[mouse.data.imp$geno==NAgeno[j],id_imp[k]] = predic[[j]]
  k=k+1
}

#Putting in var 29 prediction from genotypes that never had missing values
lm_mod = lm(IMPC_HEM_029_001~., data=impMICE[!impMICE$geno %in% NAgeno,-2])
predic[[2]] = predict(lm_mod, newdata=mouse.data.imp[mouse.data.imp$geno==NAgeno[2],])
mouse.data.imp[mouse.data.imp$geno==NAgeno[2],id_imp[2]] = predic[[2]]

#Renaming
predictions = predic 

save(predictions, mouse.data.imp, file="~/Documents/research/causal_challenge_repo/Results/MORE_FINAL_PREDICTIONS.RData")












