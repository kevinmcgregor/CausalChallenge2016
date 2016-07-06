#Simulation for PC K-means clustering

library(mice)
library(Metrics)

# read data
mouse.data <- readRDS("mouse_data.rds")

set.seed(30062016)
Nsimul <- 500
Nclust <- 10

# empty matrix to save result
res <- matrix(NA, ncol = 4, nrow = Nsimul)

# consider only complete observations: delete genotypes with missing values
full_mouse <- mouse.data[-which(mouse.data$geno %in% c("1796_1", "3621_1", "4045_1", "3803_1", "3887_1")),]

for(i in 1:Nsimul)
{
  # temporary dataframe
  temp <- full_mouse
  
  #Calculating clusters based on first two PCs
  pheno.pc <- princomp(temp[,5:26], cor=TRUE)
  temp$kclust = factor(kmeans(pheno.pc$scores[,1:2], Nclust, nstart=25)$clust)
  
  # empty method vector
  #method_norm <- c("", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "")
  method_norm = rep("", NCOL(temp)-2)
  #method_pmm <- c("", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "")
  method_pmm = rep("", NCOL(temp)-2)
  
  # randomly select one variable among the 5 variables with missing values + one knock-out condition
  NAgeno <- sample(unique(temp$geno)[2:9], 1)
  NAvar <- sample(c("IMPC_HEM_027_001", "IMPC_HEM_029_001", "IMPC_HEM_031_001", "IMPC_HEM_034_001", "IMPC_HEM_038_001"), 1) 
  
  # save true value + delete in temp
  true <- temp[temp$geno == NAgeno, colnames(temp) == NAvar]
  temp[temp$geno == NAgeno,  colnames(temp) == NAvar] <- NA
  
  # create data frame with all variables to impute + predictors
  impMICE <- as.data.frame(cbind(temp[,2:3], temp[,5:27]))
  
  # column id of variable with missing values
  id <- which(colnames(impMICE) == NAvar)
  
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
  
  # get predictions
  predic_norm <- apply(MICE_norm$imp[id][[1]], 1, mean)
  predic_pmm <- apply(MICE_pmm$imp[id][[1]], 1, mean)
  
  # compute mse
  mse_norm <- mse(actual = true, predicted = predic_norm)
  mse_pmm <- mse(actual = true, predicted = predic_pmm)
  
  # save results
  res[i, ] <- c(NAvar, NAgeno, mse_norm, mse_pmm)
  
  print(i)
#   print("TRUE: ")
#   print(true)
#   print("PRED_NORM: ")
#   print(predic_norm)
#   print("PRED_PMM: ")
#   print(predic_pmm)
}
res[,3] <- as.numeric(res[,3])
res[,4] <- as.numeric(res[,4])
res <- as.data.frame(res)
colnames(res) <- c("Variable", "Genotype", "MSE_norm", "MSE_pmm")
write.csv(res, file = "simulation_pc_k10.csv")

results_pc <- read.csv("simulation_pc_k10.csv")

par(mfrow = c(1, 2))
boxplot(MSE_norm ~ Variable, data = results_pc, main = "Imputation via OLS")
boxplot(MSE_pmm ~ Variable, data = results_pc, main = "Imputation via PMM")

boxplot(MSE_norm ~ Variable, data = results_pc, main = "Imputation via OLS", ylim = c(0, 0.4))
boxplot(MSE_pmm ~ Variable, data = results_pc, main = "Imputation via PMM", ylim = c(0, 0.4))

#Read in simulation not using PCs to compare
results <- read.csv("simulation_mice.csv")

mean(results$MSE_norm)
mean(results_pc$MSE_norm)

mean(results$MSE_pmm)
mean(results_pc$MSE_pmm)





