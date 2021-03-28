###########################################
# Comparing intervention choice simulations 
# 4/14/19
###########################################

library(ggplot2)
 # folderName <- c("AdaptiveChild","AdaptiveDesc", "ParentSet","ChildSet",
 #                 "DescSet","Random")
 # Scen <- c(10,10,10,10,10,10)
folderName <- c("ChildSet")
Scen <- c(52)
 nNodes <- 8
#folderName <- c("ChildSet","ChildSet","ChildSet")
#folderName <- c("Random","DescSet")
# Scen <- c(14, 'Scen14_100obs', 'Scen14_500obs')
# nNodes <- 9
# set directory to save files
#outPath <- setwd("C:/Users/Michele/Dropbox/Jeff Miller/Presentations/Student Seminar April 2019/Images/")
outPath <- setwd("C:/Users/Michele/Dropbox/Biostat Other/Committee Meetings/May 2019/P2Images/")
# set label for saving graphs
label <- "Scen10_DSupdate"

# initialize vectors that will be columns of data frame
method <- c()
experiment <- c()

# posterior entropy columns
avg.ent <- c(); var.ent <- c(); med.ent <- c(); halfCI.ent <- c()
# accuracy of network inference
#prop.match <- c(); halfCI.match <- c()
avg.tpr <-c(); halfCI.tpr <- c()

# Hamming distance
avg.hamMean <- c(); med.hamMean <- c();halfCI.hamMean <- c()
for(i in 1:length(folderName)){
  
  # entropy
  ent <- read.csv(paste("C:/Users/Michele/Dropbox/Jeff Miller/BDAGL/BDAGLMZ/Sim Scenarios/compStructures/",folderName[i],"/Scen",Scen[i],"/postEntropy_Results.csv",sep=""), na.strings = "NaN", header = F)
  nSIM <- nrow(ent)
  nExp <- ncol(ent)
  method <- c(method, rep(folderName[i], times = nExp))
  experiment <- c(experiment, 1:nExp)
  
  # functions related to posterior entropy
  avg.ent <- c(avg.ent, apply(ent, MARGIN = 2, function(x) mean(x,na.rm = T)))
  var.ent <- c(var.ent, apply(ent, MARGIN = 2, function(x) var(x,na.rm = T)))
  med.ent <- c(med.ent, apply(ent, MARGIN = 2, function(x) median(x,na.rm = T)))
  se.ent  <- apply(ent, MARGIN = 2, sd) / sqrt(nSIM)
  halfCI.ent <- c(halfCI.ent, (qt(0.975,nSIM-1) * se.ent))
  
  
  # accuracy of network inference
  # match <- read.csv(paste("C:/Users/Michele/Dropbox/Jeff Miller/BDAGL/BDAGLMZ/Sim Scenarios/compStructures/",folderName[i],"/Scen",Scen[i],"/matchTD_Results.csv",sep=""), na.strings = "NaN", header = F)
  #pm <- apply(match, MARGIN = 2, function(x) sum(x == 1)/length(x))
  # hc <- sapply(pm, function(pr) {1.96*sqrt((pr*(1-pr))/ nSIM)} )
  # prop.match <- c(prop.match, pm)
  # halfCI.match <- c(halfCI.match, hc)
  tpr <- read.csv(paste("C:/Users/Michele/Dropbox/Jeff Miller/BDAGL/BDAGLMZ/Sim Scenarios/compStructures/",folderName[i],"/Scen",Scen[i],"/tpr_Results.csv",sep=""), na.strings = "NaN", header = F)
  avg.tpr <- c(avg.tpr, apply(tpr, MARGIN = 2, function(x) mean(x,na.rm = T)))
  se.tpr  <- apply(tpr, MARGIN = 2, sd) / sqrt(nSIM)
  halfCI.tpr <- c(halfCI.tpr, (qt(0.975,nSIM-1) * se.tpr))
  
  # Hamming distance
  hamMean<- read.csv(paste("C:/Users/Michele/Dropbox/Jeff Miller/BDAGL/BDAGLMZ/Sim Scenarios/compStructures/", folderName[i],"/Scen",Scen[i],"/hammingMean_Results.csv",sep=""), na.strings = "NaN", header = F)
  avg.hamMean <- c(avg.hamMean, apply(hamMean, MARGIN = 2, function(x) mean(x,na.rm = T)))
  med.hamMean <- c(med.hamMean, apply(hamMean, MARGIN = 2, function(x) median(x,na.rm = T)))
  se.ham <- apply(hamMean, MARGIN = 2, sd) / sqrt(nSIM)
  halfCI.hamMean <- c(halfCI.hamMean, (qt(0.975,nSIM-1) * se.ham))
  
  
}
# combine columns into dataframe
df <- data.frame(method, experiment, 
                 avg.ent, var.ent, med.ent, se.ent, halfCI.ent,
                 #prop.match, halfCI.match,
                 avg.tpr, halfCI.tpr,
                 avg.hamMean, med.hamMean, se.ham, halfCI.hamMean)
library(plyr)

df$method <- mapvalues(df$method, from = c("ChildSet","ChildSet"), to = c("AllowRepeat","NoRepeat"))
levels(df$method)

# change names of method for plotting putposes
library(plyr)
levels(df$method)
df$method <- mapvalues(df$method, from = c("AdaptiveChild"), to = c("ChildPW"))
levels(df$method)
# plot posterior entropy across experiments
pMeanEntropy <- ggplot(data = df, aes(x=experiment, y = avg.ent, group = method))+
  geom_line(aes(color = method), size = 1)+
  geom_point(aes(color = method),size = 2)+
  geom_errorbar(aes(ymin=avg.ent-halfCI.ent, ymax=avg.ent+halfCI.ent, color = method), width=.1)+
  labs(title = "Mean Posterior Entropy",x="Experiment Number",y="Mean Posterior Entropy")+
  ylim(0,max(avg.ent))+
  scale_x_continuous(breaks= seq(from = 2, to = nExp, by = 2))+
  theme(plot.title = element_text(size = rel(2))) +
  theme(plot.title = element_text(hjust = 0.5));pMeanEntropy
ggsave(filename=paste("meanEntropy_",label,".pdf", sep=""), device = "pdf")

# plot true positive edge rate
pTPR <- ggplot(data=df, aes(x=experiment, y = avg.tpr, group = method))+
  geom_line(aes(color = method), size = 1)+
  geom_point(aes(color = method),size = 2)+
  geom_errorbar(aes(ymin=avg.tpr-halfCI.tpr, ymax=avg.tpr+halfCI.tpr, color = method), width=.1)+
  labs(title = "True Positive Rate",x="Experiment Number",y="True Positive Edge Rate")+
  ylim(0,max(avg.tpr))+
  scale_x_continuous(breaks= seq(from = 2, to = nExp, by = 2))+
  theme(plot.title = element_text(size = rel(2))) +
  theme(plot.title = element_text(hjust = 0.5));pTPR
ggsave(filename=paste("tpr_",label,".pdf", sep=""), device = "pdf")

# plot accuracy of network inference
# pAccuracy <- ggplot(data=df, aes(x=experiment, y = prop.match, group = method))+
#   geom_line(aes(color = method), size = 1)+
#   geom_point(aes(color = method),size = 2)+
#   geom_errorbar(aes(ymin=prop.match-halfCI.match, ymax=prop.match+halfCI.match, color = method), width=.1)+
#   labs(title = "Accuracy of Network Inference",x="Experiment Number",y="Prop. of Simulations Matching True Network")+
#   ylim(0,1)+
#   scale_x_continuous(breaks= seq(from = 2, to = nExp, by = 2))+
#   theme(plot.title = element_text(size = rel(2))) +
#   theme(plot.title = element_text(hjust = 0.5));pAccuracy
# ggsave(filename=paste("accuracy_",label,".pdf", sep=""), device = "pdf")

# plot median posterior entropy across experiments
pMedianEntropy <- ggplot(data = df, aes(x=experiment, y = med.ent, group = method))+
  geom_line(aes(color = method), size = 1)+
  geom_point(aes(color = method),size = 2)+
 # geom_errorbar(aes(ymin=avg.ent-halfCI.ent, ymax=avg.ent+halfCI.ent, color = method), width=.1)+
  labs(title = "Median Posterior Entropy",x="Experiment Number",y="Median Posterior Entropy")+
  ylim(0, max(med.ent))+
  scale_x_continuous(breaks= seq(from = 2, to = nExp, by = 2))+
  theme(plot.title = element_text(size = rel(2))) +
  theme(plot.title = element_text(hjust = 0.5));pMedianEntropy
ggsave(filename=paste("medianEntropy_",label,".pdf", sep=""), device = "pdf")

# plot variance of posterior entropy across experiments
pVarEntropy <- ggplot(data = df, aes(x=experiment, y = var.ent, group = method))+
  geom_line(aes(color = method), size = 1)+
  geom_point(aes(color = method),size = 2)+
  # geom_errorbar(aes(ymin=avg.ent-halfCI.ent, ymax=avg.ent+halfCI.ent, color = method), width=.1)+
  ylim(0,max(var.ent))+
  labs(title = "Variance of Posterior Entropy",x="Experiment Number",y="Variance of Posterior Entropy")+
  scale_x_continuous(breaks= seq(from = 2, to = nExp, by = 2))+
  theme(plot.title = element_text(size = rel(2))) +
  theme(plot.title = element_text(hjust = 0.5));pVarEntropy
ggsave(filename=paste("varEntropy_",label,".pdf", sep=""), device = "pdf")

# plot mean Hamming distance across experiments
pMeanHamming <- ggplot(data = df, aes(x=experiment, y = avg.hamMean, group = method))+
  geom_line(aes(color = method), size = 1)+
  geom_point(aes(color = method),size = 2)+
  geom_errorbar(aes(ymin=avg.hamMean-halfCI.hamMean, ymax=avg.hamMean+halfCI.hamMean, color = method), width=.1)+
  ylim(0,max(avg.hamMean))+
  labs(title = "Mean Hamming Distance",x="Experiment Number",y="Mean Hamming Distance")+
  scale_x_continuous(breaks= seq(from = 2, to = nExp, by = 2))+
  theme(plot.title = element_text(size = rel(2))) +
  theme(plot.title = element_text(hjust = 0.5));pMeanHamming
ggsave(filename=paste("meanHamming_",label,".pdf", sep=""), device = "pdf")

########################################
# looking at sequences of interventions
########################################
library("plotrix")
library(reshape2)

for (i in 1:length(folderName)){
  intv <- read.csv(paste("C:/Users/Michele/Dropbox/Jeff Miller/BDAGL/BDAGLMZ/Sim Scenarios/compStructures/",folderName[i],"/Scen",Scen[i],"/intvSeq_Results.csv",sep=""), na.strings = "NaN", header = F)
  perc <- t(apply(intv, MARGIN = 2, function(x) table(factor(x, levels =0:nNodes))/nSIM))
  row.names(perc) <- c(1:nExp)
  melted.perc <- melt(perc)
  colnames(melted.perc) <- c("expNum", "node", "prop")
  pIntvMat <- ggplot(data = melted.perc, aes(x = node, y=expNum, fill=prop))+
    geom_tile(color="white")+
    #scale_fill_gradient(low="#132B43", high = "#56B1F7")+
    labs(x = "Node Intervened On", y = "Experiment Number",
         title = paste(folderName[i]," Intervention Sequence",sep=""))+
    theme(plot.title = element_text(hjust = 0.5))+
    scale_y_continuous(breaks= seq(from = 1, to = nExp, by = 1))+
    scale_x_continuous(breaks= seq(from = 0, to = nNodes, by = 1))+
    geom_text(aes(node, expNum, label=prop),color = "white"); pIntvMat
  ggsave(filename=paste(folderName[i],"_IntvMat_",label,".pdf",sep=""), device = "pdf")
}

# old way using matplot
color2D.matplot(perc, main = paste(folderName[i]," Intervention Sequence",sep=""),
                xlab = paste("Node Intervened On (0 -> ",nNodes,")", sep = ""), 
                ylab = paste("Experiment Number (",nExp," <- 1)", sep=""), 
                show.values = 2, yrev= T, axes = F )

                                                                                             