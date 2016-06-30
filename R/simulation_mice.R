#################################
# simulation MICE               #
#################################

# For discussion:
#   - simulations: (1) delete one (randomly selected) variable for (randomly selected) knock-out condition 
#                      and repetively apply MICE (e.g. 1000 times) with different imputation scenarios 
#                  OR  
#                  (2) for each "simulated dataset" (e.g. 1000), delete one (randomly selected) variable for (randomly selected) knock-out condition 
#                      and repetively apply (naive) MICE i.e. always use the same imputation model
#   - measure in each simulation: MSE
#   - parameters to consider: validity of the imputation model (variable selection and form of the model)
#                             this we really don't know, we want to see whether MICE performs well in predicting missing values with a naive imputation model

install.packages("mice")
install.packages("Metrics")
library(mice)
library(Metrics)

# read data
mouse.data <- readRDS("mouse_data.rds")

set.seed(30062016)
Nsimul <- 500

# empty matrix to save result
res <- matrix(NA, ncol = 4, nrow = Nsimul)

# consider only complete observations: delete genotypes with missing values
full_mouse <- mouse.data[-which(mouse.data$geno %in% c("1796_1", "3621_1", "4045_1", "3803_1", "3887_1")),]

for(i in 1:Nsimul)
{
  # temporary dataframe
  temp <- full_mouse
  
  # empty method vector
  method_norm <- c("", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "")
  method_pmm <- c("", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "")
  
  # randomly select one variable + one knock-out condition
  NAgeno <- sample(unique(temp$geno)[2:9], 1) # 727_1
  NAvar <- sample(colnames(temp)[5:26], 1)   # IMPC_HEM_005_001
  
  # save true value + delete in temp
  true <- temp[temp$geno == NAgeno, colnames(temp) == NAvar]
  temp[temp$geno == NAgeno,  colnames(temp) == NAvar] <- NA
  
  # create data frame with all variables to impute + predictors
  impMICE <- as.data.frame(cbind(temp[,2:3], temp[,5:26]))
  
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
}
res[,3] <- as.numeric(res[,3])
res[,4] <- as.numeric(res[,4])
res <- as.data.frame(res)
colnames(res) <- c("Variable", "Genotype", "MSE_norm", "MSE_pmm")
write.csv(res, file = "simulation_mice.csv")















