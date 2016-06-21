################################################
## Data for the CRM Causal Inference Challenge
## To register, please visit-
## http://www.crm.umontreal.ca/2016/Genetics16/competition_e.php

#################################
# Read in data                  #
#################################
mouse.data <- readRDS("~/Documents/research/causal_challenge_2016/mouse_data.RDS")
var.names <- readRDS("~/Documents/research/causal_challenge_2016/variable_names.RDS")

names(mouse.data)
# id: is a unique identifier for each mouse
# sex:
# geno: indicates the experimental condition.
#       0 indicates wildtype
# litter: identifier for the litter of the mouse
# IMPC_HEM_XXX_XXX: Phenotype data. Scientific name given in var.names


#################################
# Pairwise Plots                #
#################################
pairs(mouse.data[, 5:15])
pairs(mouse.data[, 16:26])


#################################
# Estimate Graph                #
#################################
library('pcalg')
library('Rgraphviz')
# LiNGAM 
A <- LINGAM(mouse.data[mouse.data$geno == 0, 8:26])

# PC Alg
pc.fit <- pc(suffStat = list(C = cor(mouse.data[which(mouse.data$geno == 0), 5:26]), n = sum(mouse.data$geno == 0)),
             indepTest = gaussCItest, alpha=0.01, labels = names(mouse.data)[5:26], verbose = TRUE)
plot(pc.fit, main = "Estimated CPDAG", cex=1.5)
