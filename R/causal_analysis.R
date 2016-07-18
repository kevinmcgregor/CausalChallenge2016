# causal interpretation genotype-phenotype

mouse.data <- readRDS("mouse_data.RDS")
var.names <- readRDS("variable_names.RDS")
colnames(mouse.data)[5:26] <- var.names$nam
attach(mouse.data)

# univariate linear regressions
all.lm <- lm(as.matrix(mouse.data[,5:26]) ~ geno)
summary(all.lm)


## create a dataframe with genotype-phenotype association
#     13 columns, one for each knockout conditions (genotype)
#     22 rows, one for each phenotypic measurement
#     each element report the p-value for association if pval <=0.05, reports 0 otherwise

extract.coefficients <- function(var) # to extract coefficient matrix of each model
{
  return(as.matrix(var$coefficients))
}
extract.pvalue <- function(var) # to extract p-value for each genotype, each variable
{
  variable <- ifelse(var[,4]<=0.05, var[,4], 0)
  return(variable)
}

coeff <- lapply(summary(all.lm), extract.coefficients)
pval <- lapply(coeff, extract.pvalue)
geno_pheno <- t(as.data.frame(pval))
geno_pheno <- geno_pheno[,-1] # remove column intercept


## create list with phenotypic measurement associated to each genotype -> quick view
#     each element of the list represents one genotype
#     lists the associated phenotypes

quick <- vector(mode = "list", length = 8) # change 8 to 13 when we have predictions
names(quick) <- colnames(geno_pheno)

for(i in 1:8) # change 8 to 13
{
  pvalues <- geno_pheno[,i]
  quick[[i]] <- rownames(geno_pheno)[which(pvalues != 0)]
}



