# following example in course/2012/ch25 - Discovering Causal Structure from Observations 

library(pcalg)
library(Rgraphviz)
library(SMPracticals)

## example data
data(mathmarks)
suffStat <- list(C=cor(mathmarks), n=nrow(mathmarks))
pc.fit <- pc(suffStat, indepTest=gaussCItest, p=ncol(mathmarks), alpha=0.005, verbose=TRUE)
plot(pc.fit, labels=colnames(mathmarks), main="Inferred DAG for mathmarks")

## apply to mouse data
mouse.data <- readRDS("mouse_data.RDS") 
var.names <- readRDS("variable_names.RDS")
mouse <- mouse.data[mouse.data$geno==0,5:26]
colnames(mouse) <- var.names$nam

cor(mouse)
suffStat <- list(C=cor(mouse), n=nrow(mouse))
pc.fit <- pc(suffStat, indepTest=gaussCItest, p=ncol(mouse), alpha=0.005, verbose=TRUE, conservative = TRUE)
plot(pc.fit, labels=seq(1:22), main="Inferred DAG for 22 phenotypic measurements", cex.lab = 2)
# hard to interpret...

## try deleting red and white blood cell count
mouse1 <- mouse[,-c(1,2)]
suffStat <- list(C=cor(mouse1), n=nrow(mouse1))
pc.fit <- pc(suffStat, indepTest=gaussCItest, p=ncol(mouse1), alpha=0.005, verbose=TRUE)
plot(pc.fit, labels=seq(3,22), main="Inferred DAG for 20 phenotypic measurements", cex.lab = 2)

## try with only differential counts
mouse2 <- mouse[,which(colnames(mouse) %in% c("Neutrophil differential count","Lymphocyte differential count", "Monocyte differential count", "Eosinophil differential count", "Basophil differential count"))]
suffStat <- list(C=cor(mouse2), n=nrow(mouse2))
pc.fit <- pc(suffStat, indepTest=gaussCItest, p=ncol(mouse2), alpha=0.005, verbose=TRUE)
plot(pc.fit, labels=colnames(mouse2), main="Inferred DAG for 5 phenotypic measurements", cex.lab = 2)

## try including genotype
mouse3 <- mouse.data[-which(mouse.data$geno %in% c("1796_1", "3621_1", "4045_1", "3803_1", "3887_1")),]
colnames(mouse3)[5:26] <- var.names$nam
mouse3 <- mouse3[,-c(1,4)]
mouse3$sex <- as.numeric(mouse3$sex)
mouse3$geno <- as.numeric(mouse3$geno)

suffStat <- list(C=cor(mouse3), n=nrow(mouse3))
pc.fit <- pc(suffStat, indepTest=gaussCItest, p=ncol(mouse3), alpha=0.005, verbose=TRUE)
plot(pc.fit, labels=c("sex", "geno", seq(1:22)), main="Inferred DAG for 22 phenotypic measurements", cex.lab = 2)

## with only differential counts
mouse4 <- mouse3[,which(colnames(mouse3) %in% c("sex", "geno", "Neutrophil differential count","Lymphocyte differential count", "Monocyte differential count", "Eosinophil differential count", "Basophil differential count"))]
suffStat <- list(C=cor(mouse4), n=nrow(mouse4))
pc.fit <- pc(suffStat, indepTest=gaussCItest, p=ncol(mouse4), alpha=0.005, verbose=TRUE)
plot(pc.fit, labels=colnames(mouse4), main="Inferred DAG for 5 phenotypic measurements", cex.lab = 2)




