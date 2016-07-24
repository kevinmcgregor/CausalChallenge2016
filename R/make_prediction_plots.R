#Make nice plots for poster (prediction stuff)
library(ggplot2)

load("~/Documents/research/causal_challenge_repo/data/simulation_more_final.RData")

plot_data = data.frame(log(c(res_direct[,1],res_direct[,2],res_direct[,3],res_direct[,4],res_direct[,5])),
                       rep(1:5,each = 100))
colnames(plot_data) = c("MSE", "Variable")

#Plotting simulation results
pdf("~/Documents/research/causal_challenge_repo/plots/sim_plot.pdf",width=8,height=5)
gplot = ggplot(plot_data, aes(factor(Variable), MSE))
gplot+geom_boxplot(fill="darkgreen", color="red") +
  scale_x_discrete(labels=c("RBC dist width","Neutro diff count",
                            "Lymph diff count",
                            "Mono cell count",
                            "Baso diff count")) + theme_bw() +
  theme(axis.text.x=element_text(angle=25,size=11,hjust = 1),
        axis.title.y=element_text(size=14), axis.title.x=element_blank()) + ylab("log(MSE)") +
  ggtitle("Mean square error of simulation predictions")
dev.off()

#Plotting final imputed values for var 29 among the true values in the other genotypes
load("~/Documents/research/causal_challenge_repo/Results/MORE_FINAL_PREDICTIONS.RData")

pdf("~/Documents/research/causal_challenge_repo/plots/var29_plot.pdf",width=9,height=5)
gplot2 = ggplot(mouse.data.imp, aes(factor(geno), IMPC_HEM_029_001))
gplot2+geom_boxplot(fill="darkgreen", color="red") + 
  geom_boxplot(data=mouse.data.imp[mouse.data.imp$geno=="3621_1",],aes(factor(geno), IMPC_HEM_029_001), fill="blue") +
  theme_bw() +
  theme(axis.text.x=element_text(angle=25,size=11,hjust = 1),
        axis.title.y=element_text(size=14)) + ylab("Neutrophil differential count") + xlab("Genotype") +
  ggtitle("Predicted values for Neutrophil differential count")
dev.off()










