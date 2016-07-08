#Same MICE simulation that Gabrielle wrote, except using quickpred to choose subset of
#variables for prediction of missing values


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
  method_norm <- c("", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "")
  method_pmm <- c("", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "")
  
  # randomly select one variable among the 5 variables with missing values + one knock-out condition
  NAgeno <- sample(unique(temp$geno)[2:9], 1)
  NAvar <- sample(c("IMPC_HEM_027_001", "IMPC_HEM_029_001", "IMPC_HEM_031_001", "IMPC_HEM_034_001", "IMPC_HEM_038_001"), 1) 
  
  # save true value + delete in temp
  true <- temp[temp$geno == NAgeno, colnames(temp) == NAvar]
  temp[temp$geno == NAgeno,  colnames(temp) == NAvar] <- NA
  
  # create data frame with all variables to impute + predictors
  #LEAVING GENOTYPE OUT FOR NOW
  impMICE <- as.data.frame(temp[,c(2,5:26)])
  
  # column id of variable with missing values
  id <- which(colnames(impMICE) == NAvar)
  
  # method norm + pmm
  method_norm[id] <- "norm"
  method_pmm[id] <- "pmm"
  
  # matrix to specify which preditors to use 
#   pred <- matrix(0, nrow = ncol(impMICE), ncol = ncol(impMICE))
#   pred[id,] <- rep(1, ncol(impMICE))
#   diag(pred) <- rep(0, ncol(impMICE))
  
  #Choose predictors for each missing variable
  pred <- quickpred(impMICE)
  pred <- cbind(0,)

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
  
  cat(i,"\n")
}
res[,3] <- as.numeric(res[,3])
res[,4] <- as.numeric(res[,4])
res <- as.data.frame(res)
colnames(res) <- c("Variable", "Genotype", "MSE_norm", "MSE_pmm")


write.csv(res, file = "simulation_quickpred.csv")

# use this command when running on gate.mcgill 
quit(save = "no")

# I initially ran the simulations for all phenotypic measurements i.e. at each iteration I would sample any of the 22 
# phenotypic measurements and assess the performance of MICE for predicting that measurement
# However, we are really only interested in the performance of MICE for the 5 variables we will have to eventually impute
# So instead of sampling from all phenotypic measurements, I sampled only from the 5 variables of interest.

####################################
# analyze results simulations MICE #
####################################

results <- read.csv("simulation_quickpred.csv")

# View distribution of MSE with norm vs pmm for the 5 variables to impute
par(mfrow = c(1, 2))
boxplot(MSE_norm ~ Variable, data = results, main = "Imputation via OLS")
boxplot(MSE_pmm ~ Variable, data = results, main = "Imputation via PMM")

boxplot(MSE_norm ~ Variable, data = results, main = "Imputation via OLS", ylim = c(0, 0.4))
boxplot(MSE_pmm ~ Variable, data = results, main = "Imputation via PMM", ylim = c(0, 0.4))

# test difference mse_norm vs mse_pmm for each variable
t.test(x = results$MSE_norm[results$Variable == "IMPC_HEM_027_001"], y = results$MSE_pmm[results$Variable == "IMPC_HEM_027_001"]) # no
t.test(x = results$MSE_norm[results$Variable == "IMPC_HEM_029_001"], y = results$MSE_pmm[results$Variable == "IMPC_HEM_029_001"]) # yes - norm is better
t.test(x = results$MSE_norm[results$Variable == "IMPC_HEM_031_001"], y = results$MSE_pmm[results$Variable == "IMPC_HEM_031_001"]) # yes - norm is better
t.test(x = results$MSE_norm[results$Variable == "IMPC_HEM_034_001"], y = results$MSE_pmm[results$Variable == "IMPC_HEM_034_001"]) # no - both very good, pmm is better
t.test(x = results$MSE_norm[results$Variable == "IMPC_HEM_038_001"], y = results$MSE_pmm[results$Variable == "IMPC_HEM_038_001"]) # yes - norm is better

#Loading original results to compare
orig_results = read.csv("simulation_mice.csv")

mean(orig_results$MSE_norm)
mean(results$MSE_norm)

mean(orig_results$MSE_pmm)
mean(results$MSE_pmm)


#Trying higher order polynomials for variable 27 prediction
init_lm = lm(IMPC_HEM_027_001~.-litter-id, data=full_mouse)
init_pred = predict(init_lm, type="response")
mse(full_mouse$IMPC_HEM_027_001, init_pred)

quad_lm = lm(IMPC_HEM_027_001~.-litter-id+I(IMPC_HEM_002_001^2)+I(IMPC_HEM_004_001^2), data=full_mouse)
quad_pred = predict(quad_lm, type="response")
mse(full_mouse$IMPC_HEM_027_001, quad_pred)

sp = smooth.spline(full_mouse$IMPC_HEM_004_001,full_mouse$IMPC_HEM_027_001, nknots = 10)
plot(mouse.data$IMPC_HEM_004_001, mouse.data$IMPC_HEM_027_001)
lines(sp, col="red")

#Trying generalized additive model
library(mgcv)
gam = gam(IMPC_HEM_027_001~s(IMPC_HEM_002_001,IMPC_HEM_004_001, IMPC_HEM_005_001, IMPC_HEM_008_001,
                             IMPC_HEM_019_001,fx=FALSE,bs="tp"), data=full_mouse)

phen.matrix = as.matrix(full_mouse[,5:26])
phen.matrix = phen.matrix[,colnames(phen.matrix)!="IMPC_HEM_027_001"]
gam2 = gam(IMPC_HEM_027_001~s(phen.matrix,fx=FALSE,bs="tp"), data=full_mouse)




