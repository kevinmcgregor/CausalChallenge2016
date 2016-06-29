# Trying PCA to see if we can identify clusters based on litter

mouse.data <- readRDS("~/Documents/research/causal_challenge_repo/data/mouse_data.RDS")
var.names <- readRDS("~/Documents/research/causal_challenge_repo/data/variable_names.RDS")

#Find variables with missing data  ####
vars.missing = sapply(mouse.data, function(x){any(is.na(x))})

#Reduced dataset with only complete variables and only phenotypes 
vars.to.remove = c(1:4,which(vars.missing))
pheno.nonmissing = as.matrix(mouse.data[,-vars.to.remove])

#Constructing principal components of complete phenotype vars ####
pheno.pc = princomp(pheno.nonmissing, cor=TRUE)

#Include scores in mouse data frame
mouse.data = data.frame(mouse.data, pheno.pc$scores[,1:2])

#Plotting by sex
plot(pheno.pc$scores[mouse.data$sex=="1",1], pheno.pc$scores[mouse.data$sex=="1",2], xlim=c(-6,7),ylim=c(-5,8),col='red')
points(pheno.pc$scores[mouse.data$sex=="0",1], pheno.pc$scores[mouse.data$sex=="0",2],col='green')

#Boxplots of PCs for different litters
bymedian <- with(mouse.data, reorder(litter, Comp.1, median))
boxplot(Comp.1~bymedian, data=mouse.data)

mean_by_litter = aggregate(Comp.1~litter, data=mouse.data, FUN=mean)
sorted_mean_by_litter = mean_by_litter[order(mean_by_litter$Comp.1),]
sorted_mean_by_litter$litter = factor(sorted_mean_by_litter$litter, levels = sorted_mean_by_litter$litter)
boxplot(Comp.1~as.factor(litter), data=sorted_mean_by_litter)

#Trying again with sum of first two PCs
mouse.data$sum.pc = mouse.data$Comp.1 + mouse.data$Comp.2
bymedian_sum = with(mouse.data, reorder(litter, sum.pc, median))
boxplot(sum.pc~bymedian_sum, data=mouse.data)

#Running predictive models with litter as predictor and comparing to same model with 2 PCs as predictors
#to see if the predictive power is the same...
pred1_mod = lm(IMPC_HEM_001_001~as.factor(litter), data=mouse.data)
pred1 = predict(pred1_mod, type="response")
mse1 = mean((pred1-mouse.data$IMPC_HEM_001_001)^2)

predpc1_mod = lm(IMPC_HEM_001_001~Comp.1+Comp.2, data=mouse.data)
predpc1 = predict(predpc1_mod, type="response")
msepc1 = mean((predpc1-mouse.data$IMPC_HEM_001_001)^2)

#Comparing litter predictive power to 2 PCs predictive power for each (non-missing) phenotype
getMSE = function(phen, litter, pcs) {
  pred = predict(lm(phen~as.factor(litter)), type="response")
  pred_pc = predict(lm(phen~pcs), type="response")
  return( c(mean((pred-phen)^2), mean((pred_pc-phen)^2)) )
}

pcs = cbind(mouse.data$Comp.1, mouse.data$Comp.2)
cols_to_include = 5:26

allphen_mse = apply(as.matrix(mouse.data[,cols_to_include]), 2, getMSE, litter=mouse.data$litter, pcs=pcs)


###################################
#Plotting different phenotypes with respect to litter with and without PCs 
#(to see if PCs make litters more homogenous)

boxplot(mouse.data$IMPC_HEM_001_001~mouse.data$litter)
res1 = residuals(lm(IMPC_HEM_001_001~Comp.1+Comp.2, data=mouse.data))
boxplot(res1~mouse.data$litter)





