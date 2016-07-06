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

res1 = residuals(lm(IMPC_HEM_001_001~Comp.1+Comp.2, data=mouse.data))
par(mfrow=c(2,1), mar=c(4,2,2,1))
boxplot(mouse.data$IMPC_HEM_001_001~mouse.data$litter)
boxplot(res1~mouse.data$litter)

#Calculate quintiles of means of 1st PC
mean_by_litter$comp1_level = cut(mean_by_litter$Comp.1, quantile(mean_by_litter$Comp.1, seq(0,1,1/20)), 
                                 labels=1:20, include.lowest=TRUE)
mouse.data$comp1_level = rep(0,NROW(mouse.data))
for (l in mean_by_litter$litter) {
  mouse.data[mouse.data$litter==l,]$comp1_level = mean_by_litter[mean_by_litter==l,]$comp1_level
}

#Trying K-means clustering based on PC1 and PC2
clust = kmeans(pheno.pc$scores[,1:3], 15, nstart=20 )
library(RColorBrewer)
library(scales)
palette(alpha(brewer.pal(9,'Set1'), 0.5))
plot(pheno.pc$scores[,1:2], col=clust$clust, pch=16)

mouse.data$kclust = factor(clust$clust)


######################################
#Deleting missing rows and testing predictive model on variables that had missing data before the deletion ####
rows.missing = apply(as.matrix(mouse.data), 1, function(x){any(is.na(x))})
mouse.nonmissing.rows = mouse.data[!rows.missing,]
#Randomly delete observations from one variable for all observations of a given genotype
set.seed(3487)
rand.deleted.geno = sample(mouse.nonmissing.rows$geno[mouse.nonmissing.rows$geno!=0],1)

#Prediction without PC category
p1 = lm(IMPC_HEM_027_001~.-Comp.1-Comp.2-sum.pc-kclust, data=mouse.nonmissing.rows[,-c(1,3,4,30)],
        subset = (mouse.nonmissing.rows$geno!=rand.deleted.geno))

p1_predicted = predict(p1,
                       newdata=mouse.nonmissing.rows[mouse.nonmissing.rows$geno==rand.deleted.geno,], 
                       type="response")
mean((p1_predicted-mouse.nonmissing.rows[mouse.nonmissing.rows$geno==rand.deleted.geno,]$IMPC_HEM_027_001)^2)

#Prediction with PC category
p1_cat = lm(IMPC_HEM_027_001~.-Comp.1-Comp.2-sum.pc-kclust, data=mouse.nonmissing.rows[,-c(1,3,4)],
        subset = (mouse.nonmissing.rows$geno!=rand.deleted.geno))

p1_predicted_cat = predict(p1_cat,
                       newdata=mouse.nonmissing.rows[mouse.nonmissing.rows$geno==rand.deleted.geno,], 
                       type="response")

mean((p1_predicted_cat-mouse.nonmissing.rows[mouse.nonmissing.rows$geno==rand.deleted.geno,]$IMPC_HEM_027_001)^2)


#Prediction using k-means clusters
p1_clust = lm(IMPC_HEM_027_001~.-Comp.1-Comp.2-sum.pc-comp1_level, data=mouse.nonmissing.rows[,-c(1,3,4)],
            subset = (mouse.nonmissing.rows$geno!=rand.deleted.geno))

p1_predicted_clust = predict(p1_clust,
                           newdata=mouse.nonmissing.rows[mouse.nonmissing.rows$geno==rand.deleted.geno,], 
                           type="response")

mean((p1_predicted_clust-mouse.nonmissing.rows[mouse.nonmissing.rows$geno==rand.deleted.geno,]$IMPC_HEM_027_001)^2)












