################################################
## Data for the CRM Causal Inference Challenge
## To register, please visit-
## http://www.crm.umontreal.ca/2016/Genetics16/competition_e.php

#################################
# Read in data                  #
#################################
mouse.data <- readRDS("mouse_data.RDS")
var.names <- readRDS("variable_names.RDS")

names(mouse.data)
# id: is a unique identifier for each mouse
# sex:
# geno: indicates the experimental condition.
#       0 indicates wildtype
# litter: identifier for the litter of the mouse
# IMPC_HEM_XXX_XXX: Phenotype data. Scientific name given in var.names

#################################
# Descriptive Data              #
#################################

table(mouse.data$geno) 
mean(table(mouse.data$geno)[2:14]) # ≈13.38 mice / knock-out conditions
sum(is.na(mouse.data)) # 67 missing data overall
length(unique(mouse.data$litter)) # 190 different litters
mean(table(mouse.data$litter)) # ≈3.23 mice / litter
range(table(mouse.data$litter))

# link between phenotypic measurement, following Claudia
mouse <- mouse.data[mouse.data$geno==0,]

colnames(mouse)[5:26] <- var.names$nam
cor(mouse[,which(colnames(mouse) %in% c("Red blood cell count", "Hemoglobin", "Hematocrit"))])
# hematocrit correlated with red blood cell count and hemoglobin
# causal direction: rbcc -> hema, hemo -> hema ?
# might want to delete red blood cell count?

cor(mouse[,which(colnames(mouse) %in% c("Red blood cell count", "Hemoglobin", "Hematocrit", "Mean corpuscular hemoglobin", "Mean cell hemoglobin concentration"))])

cor(mouse[,which(colnames(mouse) %in% c("White blood cell count", "Neutrophil cell count", "Lymphocyte cell count", "Eosinophil cell count", "Monocyte cell count", "Basophil cell count"))])
# white blood cell count correlated with all separate cell counts: might want to delete it?
  
cor(mouse[,which(colnames(mouse) %in% c("Large Unstained Cell (LUC) count", "Large Unstained Cell (LUC) differential count", "Lymphocyte differential count", "Lymphocyte cell count"))])
# lymphocyte and LUC cell counts correlated, but not differential counts

cor(mouse[,which(colnames(mouse) %in% c("Eosinophil differential count", "Eosinophil cell count", "Basophil cell count", "Basophil differential count"))])

cor(mouse[,which(colnames(mouse) %in% c("Eosinophil differential count", "Basophil differential count", "Lymphocyte differential count", "Monocyte differential count", "Neutrophil differential count"))])

# look at all cell counts
cor(mouse[,which(colnames(mouse) %in% c("Lymphocyte cell count", "White blood cell count"))])

# differential counts sum up to 100%
mouse$`Neutrophil differential count`+mouse$`Lymphocyte differential count`+mouse$`Monocyte differential count`+mouse$`Eosinophil differential count`+mouse$`Basophil differential count`

#################################
# Pairwise Plots                #
#################################
pairs(mouse[, 5:26], pch=19, cex=0.5)


#################################
# Estimate Graph                #
#################################
library(pcalg)
library(Rgraphviz)

# LiNGAM 
A <- LINGAM(mouse.data[mouse.data$geno == 0, 8:26])

# PC Alg
pc.fit <- pc(suffStat = list(C = cor(mouse.data[c(5:26)]), n = nrow(mouse.data)),
             indepTest = gaussCItest, p=22, alpha=0.01, verbose = TRUE)
plot(pc.fit, labels=seq(1:22), main = "Estimated CPDAG")

