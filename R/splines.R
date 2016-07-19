#Finding variables on which to invoke splines in order to better predict IMPC_HEM_029_001

library(mice)
library(Metrics)
library(splines)

# read data
setwd("~/Documents/research/causal_challenge_repo/data")
mouse.data <- readRDS("mouse_data.rds")

# consider only complete observations: delete genotypes with missing values
full_mouse <- mouse.data[-which(mouse.data$geno %in% c("1796_1", "3621_1", "4045_1", "3803_1", "3887_1")),]

var_summaries = list()
var29 = which(colnames(full_mouse)=="IMPC_HEM_029_001")
count=1
for (i in c(5:(var29-1), (var29+1):26)) {
  mod = lm(full_mouse[,var29]~ns(full_mouse[,i],df=3))
  var_summaries[[count]] = summary(mod)$coef
  names(var_summaries)[count] = colnames(full_mouse)[i]
  count = count+1
}

#########################
#Trying out splines on every variable and comparing to OLS
n_df = 3
set.seed(30062016)
Nsimul <- 500

# empty matrix to save result
res_ols <- matrix(0, Nsimul, 1)
res_splines <- matrix(0, Nsimul, 1)

for (i in 1:Nsimul) {
  temp <- full_mouse
  
  NAgeno <- sample(unique(temp$geno)[2:9], 1)
  NAvar <- c("IMPC_HEM_029_001") 
  
  # save true values + delete in temp
  true = list(IMPC_HEM_029_001=temp[temp$geno==NAgeno,NAvar])
  temp[temp$geno == NAgeno,  colnames(temp) == NAvar] <- NA
  
  splines_matrix = temp[,colnames(temp)==NAvar]
  for (j in c(5:(var29-1), (var29+1):26)) {
    splines_matrix = cbind(splines_matrix, ns(temp[,j],df=n_df))
  }
  colnames(splines_matrix) = c(NAvar,paste0("S",1:63))
  
  # column ids of variables with missing values
  id <- which(colnames(temp) %in% NAvar)
  
  #Linear models
  lm_ols = lm(temp$IMPC_HEM_029_001~., data=temp[,-c(1,3,4)])
  lm_splines = lm(IMPC_HEM_029_001~., data=as.data.frame(splines_matrix))
  
  pred_ols = predict(lm_ols, newdata=temp[is.na(temp[,id]),])  
  pred_splines = predict(lm_splines, newdata=as.data.frame(splines_matrix[is.na(temp[,id]),]))

  
  # save results
  res_ols[i] <- mse(true[[1]], pred_ols)
  res_splines[i] <- mse(true[[1]], pred_splines)

  cat(i, "\n")
}



#Trying out with variable 27

#Trying out splines on every variable and comparing to OLS
n_df = 3
set.seed(30062016)
Nsimul <- 500

# empty matrix to save result
res_ols <- matrix(0, Nsimul, 1)
res_splines <- matrix(0, Nsimul, 1)
var27 = which(colnames(full_mouse)=="IMPC_HEM_027_001")

for (i in 1:Nsimul) {
  temp <- full_mouse
  
  NAgeno <- sample(unique(temp$geno)[2:9], 1)
  NAvar <- c("IMPC_HEM_027_001") 
  
  # save true values + delete in temp
  true = list(IMPC_HEM_027_001=temp[temp$geno==NAgeno,NAvar])
  temp[temp$geno == NAgeno,  colnames(temp) == NAvar] <- NA
  
  splines_matrix = temp[,colnames(temp)==NAvar]
  for (j in c(5:(var27-1), (var27+1):26)) {
    splines_matrix = cbind(splines_matrix, ns(temp[,j],df=n_df))
  }
  colnames(splines_matrix) = c(NAvar,paste0("S",1:63))
  
  # column ids of variables with missing values
  id <- which(colnames(temp) %in% NAvar)
  
  #Linear models
  lm_ols = lm(temp$IMPC_HEM_027_001~., data=temp[,-c(1,3,4)])
  lm_splines = lm(IMPC_HEM_027_001~., data=as.data.frame(splines_matrix))
  
  pred_ols = predict(lm_ols, newdata=temp[is.na(temp[,id]),])  
  pred_splines = predict(lm_splines, newdata=as.data.frame(splines_matrix[is.na(temp[,id]),]))
  
  
  # save results
  res_ols[i] <- mse(true[[1]], pred_ols)
  res_splines[i] <- mse(true[[1]], pred_splines)
  
  cat(i, "\n")
}
















