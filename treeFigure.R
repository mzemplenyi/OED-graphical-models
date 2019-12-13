#############
# 11/2/19
# Script for making 8-node tree figures for Paper 2 
#############

library(reshape2)
library(ggplot2)

# folderName <- c("Random","Bninfo","Leong", "ParentSet","AdaptiveChild") # need to add Cpdag to this
# Scen <- c("129","129","129","129","129","129")

folderName <- c("Leong","AdaptiveChild","Cpdag")
Scen <- c("136","136","136")
nNodes <-8

# set directory to save files
outPath <- "C:/Users/Michele/Documents/Jeff Miller/Paper Draft/Images/Simulations"
setwd(outPath); getwd()

# set label for saving graphs
label <- 'tree8'

# initialize vectors that will be columns of data frame
  method <- c()
  scenario <- c()
  experiment <- c()
  
  # accuracy of network inference
  avg.tpr <-c(); halfCI.tpr <- c()
  avg.fpr <-c(); halfCI.fpr <- c()
  avg.tnr <-c(); halfCI.tnr <- c()
  avg.fnr <-c(); halfCI.fnr <- c()
  # Hamming distance
  avg.hamMean <- c(); med.hamMean <- c();halfCI.hamMean <- c()
  for(i in 1:length(folderName)){
    print(folderName[i])
    # Hamming distance
    hamMean<- read.csv(paste("C:/Users/Michele/Documents/Jeff Miller/BDAGL/BDAGLMZ/Sim Scenarios/compStructures/", folderName[i],"/Scen",Scen[i],"/hammingMean_Results.csv",sep=""), na.strings = "NaN", header = F)
    nSIM <- nrow(hamMean)
    nExp <- ncol(hamMean)
    method <- c(method, rep(folderName[i], times = nExp))
    scenario <- c(scenario, rep(Scen[i], times = nExp))
    experiment <- c(experiment, 1:nExp)
    
    
    avg.hamMean <- c(avg.hamMean, apply(hamMean, MARGIN = 2, function(x) mean(x,na.rm = T)))
    med.hamMean <- c(med.hamMean, apply(hamMean, MARGIN = 2, function(x) median(x,na.rm = T)))
    se.ham <- apply(hamMean, MARGIN = 2, sd) / sqrt(nSIM)
    #halfCI.hamMean <- c(halfCI.hamMean, (qt(0.975,nSIM-1) * se.ham)) # 95% CI
    halfCI.hamMean <- c(halfCI.hamMean, se.ham) # one standard error 
    
    # accuracy of network inference
    tpr <- read.csv(paste("C:/Users/Michele/Documents/Jeff Miller/BDAGL/BDAGLMZ/Sim Scenarios/compStructures/",folderName[i],"/Scen",Scen[i],"/tpr_Results.csv",sep=""), na.strings = "NaN", header = F)
    avg.tpr <- c(avg.tpr, apply(tpr, MARGIN = 2, function(x) mean(x,na.rm = T)))
    se.tpr  <- apply(tpr, MARGIN = 2, sd) / sqrt(nSIM)
    #halfCI.tpr <- c(halfCI.tpr, (qt(0.975,nSIM-1) * se.tpr)) # 95% CI
    halfCI.tpr <- c(halfCI.tpr, se.tpr) # one standard error 
    
    fpr <- read.csv(paste("C:/Users/Michele/Documents/Jeff Miller/BDAGL/BDAGLMZ/Sim Scenarios/compStructures/",folderName[i],"/Scen",Scen[i],"/fpr_Results.csv",sep=""), na.strings = "NaN", header = F)
    avg.fpr <- c(avg.fpr, apply(fpr, MARGIN = 2, function(x) mean(x,na.rm = T)))
    se.fpr  <- apply(fpr, MARGIN = 2, sd) / sqrt(nSIM)
    halfCI.fpr <- c(halfCI.fpr, (qt(0.975,nSIM-1) * se.fpr))
    
    tnr <- read.csv(paste("C:/Users/Michele/Documents/Jeff Miller/BDAGL/BDAGLMZ/Sim Scenarios/compStructures/",folderName[i],"/Scen",Scen[i],"/tnr_Results.csv",sep=""), na.strings = "NaN", header = F)
    avg.tnr <- c(avg.tnr, apply(tnr, MARGIN = 2, function(x) mean(x,na.rm = T)))
    se.tnr  <- apply(tnr, MARGIN = 2, sd) / sqrt(nSIM)
    halfCI.tnr <- c(halfCI.tnr, (qt(0.975,nSIM-1) * se.tnr))
    
    fnr <- read.csv(paste("C:/Users/Michele/Documents/Jeff Miller/BDAGL/BDAGLMZ/Sim Scenarios/compStructures/",folderName[i],"/Scen",Scen[i],"/fnr_Results.csv",sep=""), na.strings = "NaN", header = F)
    avg.fnr <- c(avg.fnr, apply(fnr, MARGIN = 2, function(x) mean(x,na.rm = T)))
    se.fnr  <- apply(fnr, MARGIN = 2, sd) / sqrt(nSIM)
    halfCI.fnr <- c(halfCI.fnr, (qt(0.975,nSIM-1) * se.fnr))
  }    
    
    # combine columns into dataframe
    df <- data.frame(method, scenario, experiment, 
                     avg.tpr, halfCI.tpr,
                     avg.fpr, halfCI.fpr,
                     avg.tnr, halfCI.tnr,
                     avg.fnr, halfCI.fnr,
                     avg.hamMean, med.hamMean, halfCI.hamMean)


# map any values for methods that need to be renamed
library(plyr)
df$method <- mapvalues(df$method, from = c("Cpdag"), to = c("TSE"))
df$method <- mapvalues(df$method, from = c("AdaptiveChild"), to = c("PW.Child"))
# remove extra bninfo lines
df <- df[-which(df$method == "Bninfo" & df$experiment == 8),]


## HAMMING DISTANCE / L1 ERROR
# plot mean Hamming distance across experiments with both method and scenario interaction
pMeanHamming <- ggplot(data = df, aes(x=experiment, y = avg.hamMean, group = interaction(method,scenario)))+
  geom_line(aes(color = method, linetype = scenario), size = 1)+
  geom_point(aes(color = method),size = 1)+
  geom_errorbar(aes(ymin=avg.hamMean-halfCI.hamMean, ymax=avg.hamMean+halfCI.hamMean, color = method), width=.1)+
  #scale_linetype_manual(values = c("solid","dashed"))+
  #  ylim(-0.5,,max(avg.hamMean+halfCI.hamMean))+
  ylim(0,8) +
  labs(title = "Mean Hamming Distance",x="Experiment Number",y="Mean Hamming Distance")+
  scale_x_continuous(breaks= seq(from = 1, to = 7, by = 1))+
#  scale_y_continuous(breaks= seq(from = 0, to = 8, by = 1))+
  theme(plot.title = element_text(size = rel(2))) +
  theme(plot.title = element_text(hjust = 0.5));pMeanHamming
ggsave(filename=paste("meanHamming_",label,".pdf", sep=""), device = "pdf")

pTPR <- ggplot(data=df, aes(x=experiment, y = avg.tpr, group = interaction(method,scenario)))+
  geom_line(aes(color = method, linetype = scenario), size = 1)+
  geom_point(aes(color = method),size = 2)+
  geom_errorbar(aes(ymin=avg.tpr-halfCI.tpr, ymax=avg.tpr+halfCI.tpr, color = method), width=.1)+
  labs(title = "True Positive Rate",x="Experiment Number",y="True Positive Edge Rate")+
  #ylim(0.5,1)+
  scale_x_continuous(breaks= seq(from = , to = nExp, by = 1))+
  scale_y_continuous(breaks= seq(from = 0.5, to = 1, by = 0.1))+
  theme(plot.title = element_text(size = rel(2))) +
  theme(plot.title = element_text(hjust = 0.5));pTPR
ggsave(filename=paste("tpr_",label,".pdf", sep=""), device = "pdf")

# plot LOG mean Hamming distance across experiments
pLogMeanHamming <- ggplot(data = df, aes(x=experiment, y = avg.hamMean, group = interaction(method,scenario)))+
  geom_line(aes(color = method, linetype = scenario), size = 1)+
  geom_point(aes(color = method),size = 2)+
  geom_errorbar(aes(ymin=avg.hamMean-halfCI.hamMean, ymax=avg.hamMean+halfCI.hamMean, color = method), width=.1)+
#  ylim(0,max(avg.hamMean+halfCI.hamMean))+
  #  ylim(0,2) +
  labs(title = "Log10 Mean Hamming Distance",x="Experiment Number",y="Log10 Mean Hamming Distance")+
  scale_x_continuous(breaks= seq(from = 2, to = nExp, by = 2))+
  scale_y_continuous(trans = "log10")+
  theme(plot.title = element_text(size = rel(2))) +
  theme(plot.title = element_text(hjust = 0.5));pLogMeanHamming
ggsave(filename=paste("logmeanHamming_",label,".pdf", sep=""), device = "pdf")

# plot median Hamming distance across experiments
pMedianHamming <- ggplot(data = df, aes(x=experiment, y = med.hamMean, group = interaction(method,scenario)))+
  geom_line(aes(color = method, linetype = scenario), size = 1)+
  geom_point(aes(color = method),size = 2)+
  #geom_errorbar(aes(ymin=avg.hamMean-halfCI.hamMean, ymax=avg.hamMean+halfCI.hamMean, color = scenario), width=.1)+
  ylim(0,max(med.hamMean))+
  labs(title = "Median Hamming Distance",x="Experiment Number",y="Median Hamming Distance")+
  scale_x_continuous(breaks= seq(from = 2, to = nExp, by = 2))+
  theme(plot.title = element_text(size = rel(2))) +
  theme(plot.title = element_text(hjust = 0.5));pMedianHamming
ggsave(filename=paste("medianHamming_",label,".pdf", sep=""), device = "pdf")



pFPR <- ggplot(data=df, aes(x=experiment, y = avg.fpr, group = interaction(method,scenario)))+
  geom_line(aes(color = method, linetype = scenario), size = 1)+
  geom_point(aes(color = method),size = 2)+
  geom_errorbar(aes(ymin=avg.fpr-halfCI.fpr, ymax=avg.fpr+halfCI.fpr, color = method), width=.1)+
  labs(title = "False Positive Rate",x="Experiment Number",y="False Positive Edge Rate")+
#  ylim(0,0.1)+
  scale_x_continuous(breaks= seq(from = 2, to = nExp, by = 2))+
  theme(plot.title = element_text(size = rel(2))) +
  theme(plot.title = element_text(hjust = 0.5));pFPR
ggsave(filename=paste("fpr_",label,".pdf", sep=""), device = "pdf")

pFNR <- ggplot(data=df, aes(x=experiment, y = avg.fnr, group = interaction(method,scenario)))+
  geom_line(aes(color = method, linetype = scenario), size = 1)+
  geom_point(aes(color = method),size = 2)+
  geom_errorbar(aes(ymin=avg.fnr-halfCI.fnr, ymax=avg.fnr+halfCI.fnr, color = method), width=.1)+
  labs(title = "False Negative Rate",x="Experiment Number",y="False Negative Edge Rate")+
#  ylim(0,1)+
  scale_x_continuous(breaks= seq(from = 2, to = nExp, by = 2))+
  theme(plot.title = element_text(size = rel(2))) +
  theme(plot.title = element_text(hjust = 0.5));pFNR
ggsave(filename=paste("fnr_",label,".pdf", sep=""), device = "pdf")

pTNR <- ggplot(data=df, aes(x=experiment, y = avg.tnr, group = interaction(method,scenario)))+
  geom_line(aes(color = method, linetype = scenario), size = 1)+
  geom_point(aes(color = method),size = 2)+
  geom_errorbar(aes(ymin=avg.tnr-halfCI.tnr, ymax=avg.tnr+halfCI.tnr, color = method), width=.1)+
  labs(title = "True Negative Rate",x="Experiment Number",y="True Negative Edge Rate")+
#  ylim(0.9,1)+
  scale_x_continuous(breaks= seq(from = 2, to = nExp, by = 2))+
  theme(plot.title = element_text(size = rel(2))) +
  theme(plot.title = element_text(hjust = 0.5));pTNR
ggsave(filename=paste("tnr_",label,".pdf", sep=""), device = "pdf")

########################################
# looking at sequences of interventions
########################################
library("plotrix")
#outpath = sprintf('C:/Users/Michele/Documents/Jeff Miller/BDAGL/BDAGLMZ/Sim Scenarios/compStructures/%s/Scen%d/',method,Scen);

for (i in 1:length(folderName)){
  intv <- read.csv(paste("C:/Users/Michele/Documents/Jeff Miller/BDAGL/BDAGLMZ/Sim Scenarios/compStructures/",folderName[i],"/Scen",Scen[i],"/intvSeq_Results.csv",sep=""), na.strings = "NaN", header = F)
  perc <- t(apply(intv, MARGIN = 2, function(x) table(factor(x, levels =0:nNodes))/nSIM))
  perc <- round(perc, digits = 2)
  row.names(perc) <- c(1:nrow(perc))
  melted.perc <- melt(perc)
  colnames(melted.perc) <- c("expNum", "node", "prop")
  pIntvMat <- ggplot(data = melted.perc, aes(x = node, y=expNum, fill=prop))+
    geom_tile(color="white")+
    scale_fill_gradientn(colours = c("white", "pink", "red"), values = c(0,0.05,1)) +
    #scale_fill_gradientn(colours = c("white", "blue", "red"), values = c(0,0.05,1)) +
    #scale_fill_gradient(low="#132B43", high = "#56B1F7")+
    labs(x = "Node Intervened On", y = "Experiment Number",
         title = paste(folderName[i]," Intervention Sequence",sep=""))+
    theme(plot.title = element_text(hjust = 0.5))+
    scale_y_continuous(breaks= seq(from = 1, to = nExp, by = 1))+
    scale_x_continuous(breaks= seq(from = 0, to = nNodes, by = 1))+
    geom_text(aes(node, expNum, label=prop),color = "white"); pIntvMat
  ggsave(filename=paste(folderName[i],"_IntvMat_",label,".pdf",sep=""), device = "pdf")
}
