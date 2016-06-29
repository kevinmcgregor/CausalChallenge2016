#################################
# Read in data                  #
#################################
mouse.data <- readRDS("mouse_data.RDS")
var.names <- readRDS("variable_names.RDS")
names(mouse.data)

#################################
# Descriptive Data              #
#################################

table(mouse.data$geno) 
mean(table(mouse.data$geno)[2:14]) # ≈13.38 mice / knock-out conditions
sum(is.na(mouse.data)) # 67 missing data overall
length(unique(mouse.data$litter)) # 190 different litters
mean(table(mouse.data$litter)) # ≈3.23 mice / litter
range(table(mouse.data$litter))

## missing data
# knock-out 1796_1, missing IMPC_HEM_031_001 = Lymphocyte differential count
# knock-out 3621_1, missing IMPC_HEM_027_001 = Red blood cell distribution width
# knock-out 4045_1, missing IMPC_HEM_029_001 = Neutrophil differential count
# knock-out 3803_1, missing IMPC_HEM_038_001 = Basophil differential count
# knock-out 3887_1, missing IMPC_HEM_034_001 = Monocyte cell count

# correlation between phenotypic measurements within genotype
wild <- mouse.data[mouse.data$geno=="0",]
cor(wild[,5:26])
geno_1 <- mouse.data[mouse.data$geno=="1796_1",]
cor(geno_1[,5:26], use = "pairwise.complete.obs")
boxplot(IMPC_HEM_027_001 ~ geno, data = mouse.data)
boxplot(IMPC_HEM_029_001 ~ geno, data = mouse.data)
boxplot(IMPC_HEM_031_001 ~ geno, data = mouse.data)
boxplot(IMPC_HEM_034_001 ~ geno, data = mouse.data)
boxplot(IMPC_HEM_038_001 ~ geno, data = mouse.data)

#################################
# MICE naive                    #
#################################
library(mice)

# look at the distribution of the missing variables

hist(mouse.data$IMPC_HEM_031_001) # normal, some small extreme values
hist(mouse.data$IMPC_HEM_027_001) # right skewed
hist(mouse.data$IMPC_HEM_029_001) # normal, some large extreme values
hist(mouse.data$IMPC_HEM_038_001) # right skewed
hist(mouse.data$IMPC_HEM_034_001) # slightly right skewed
# extreme values for 031 and 029 -> related to litter effect??
# skewed variables: transformation or predictive mean matching

# consider interactions ???

# create data frame with all variables to impute + predictors
impMICE <- as.data.frame(cbind(mouse.data[,2:3], mouse.data[,5:26]))
colnames(impMICE)

# imputation method: pmm for skewed variables, linear model for normally distributed variables
method <- c("", "", "", "", "", "", "", "", "", "", "", "pmm", "norm", "", "norm", "", "", "pmm", "", "", "", "pmm", "", "")

# matrix to specify which preditors to use for each imputed variables: rows 12, 13, 15, 18, 22
pred <- matrix(0, nrow = ncol(impMICE), ncol = ncol(impMICE))
pred[12,] <- rep(1, ncol(impMICE))
pred[13,] <- rep(1, ncol(impMICE))
pred[15,] <- rep(1, ncol(impMICE))
pred[18,] <- rep(1, ncol(impMICE))
pred[22,] <- rep(1, ncol(impMICE))
diag(pred) <- rep(0, ncol(impMICE))

# generate m=30 imputed dataset
mod <- mice(impMICE, m = 30, method = method, predictorMatrix = pred)

# view imputations
mod$imp$IMPC_HEM_027_001

# derive predictions as the average imputed values
pred027 <- apply(mod$imp$IMPC_HEM_027_001, 1, mean)
pred029 <- apply(mod$imp$IMPC_HEM_029_001, 1, mean)
pred031 <- apply(mod$imp$IMPC_HEM_031_001, 1, mean)
pred034 <- apply(mod$imp$IMPC_HEM_034_001, 1, mean)
pred038 <- apply(mod$imp$IMPC_HEM_038_001, 1, mean)

# create complete_mouse dataset with predictions
mouse.data$id_row <- seq(1:nrow(mouse.data))
complete_mouse <- mouse.data
complete_mouse$IMPC_HEM_027_001[which(mouse.data$id_row %in% names(pred027))] <- pred027
complete_mouse$IMPC_HEM_029_001[which(mouse.data$id_row %in% names(pred029))] <- pred029
complete_mouse$IMPC_HEM_031_001[which(mouse.data$id_row %in% names(pred031))] <- pred031
complete_mouse$IMPC_HEM_034_001[which(mouse.data$id_row %in% names(pred034))] <- pred034
complete_mouse$IMPC_HEM_038_001[which(mouse.data$id_row %in% names(pred038))] <- pred038

# look at the distribution of predictions
boxplot(IMPC_HEM_027_001 ~ geno, data = complete_mouse)
boxplot(IMPC_HEM_029_001 ~ geno, data = complete_mouse)
boxplot(IMPC_HEM_031_001 ~ geno, data = complete_mouse)
boxplot(IMPC_HEM_034_001 ~ geno, data = complete_mouse)
boxplot(IMPC_HEM_038_001 ~ geno, data = complete_mouse)

#################################
# simple simulation MICE        #
#################################
library(Metrics) # mse

# delete one (randomly selected) variable for (randomly selected) knock-out condition 
#   and assess predictive accuracy via MSE

# consider only complete observations
full_mouse <- mouse.data[-which(mouse.data$geno %in% c("1796_1", "3621_1", "4045_1", "3803_1", "3887_1")),]

# randomly select one variable + one knock-out condition
set.seed(29062016)
sample(unique(full_mouse$geno)[2:9], 1) # 727_1
sample(colnames(full_mouse)[5:26], 1)   # IMPC_HEM_005_001

# save true value + delete in full_mouse
true <- full_mouse$IMPC_HEM_005_001[full_mouse$geno == "727_1"]
full_mouse$IMPC_HEM_005_001[full_mouse$geno == "727_1"] <- NA

# distribution of the missing variable
hist(full_mouse$IMPC_HEM_005_001) # normal

# apply MICE
impMICE <- as.data.frame(cbind(full_mouse[,2:3], full_mouse[,5:26]))
colnames(impMICE)
method <- c("", "", "", "", "", "", "norm", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "")
pred <- matrix(0, nrow = ncol(impMICE), ncol = ncol(impMICE))
pred[7,] <- rep(1, ncol(impMICE))
diag(pred) <- rep(0, ncol(impMICE))
mod <- mice(impMICE, m = 30, method = method, predictorMatrix = pred)

# extract predictions
pred005 <- apply(mod$imp$IMPC_HEM_005_001, 1, mean)
mse(actual = true, predicted = pred005) # not bad


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
















