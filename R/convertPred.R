#Convert predictions to csv for submission

load("~/Documents/research/causal_challenge_repo/Results/MORE_FINAL_PREDICTIONS.RData")

p = matrix(NA,14,5)
p[,1] = predictions[[1]] 
p[,2] = predictions[[2]]
p[,3] = c(predictions[[3]],NA)
p[,4] = predictions[[4]]
p[,5] = c(predictions[[5]],NA,NA)
colnames(p) = names(predictions)

write.csv(p, row.names=FALSE, quote=FALSE, file="~/Documents/research/causal_challenge_repo/Results/PREDICTIONS_TO_SUMIT_KM.csv")

