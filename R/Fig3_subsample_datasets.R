
#################
#This script provides the code to generate the subsample analysis 
#of four Populus RNA-seq data sets used in Figure 3

#Matt Zinkgraf
#US Forest Servce and UC Davis Computer Science

#################

setwd(".")

library(edgeR)
library(WGCNA)
require(plyr)
library(RColorBrewer)
library(preprocessCore)
library(fields)
source("R/networkFunctions-extras-05.R")
source("R/functions.R")

allowWGCNAThreads(n=2);
options(stringsAsFactors = FALSE);

#load expression and modules
load("Data/results/Grand_analysis_expressions_mods.rdata")

lib_names<-names(data)

#build multiset
setLabels=c("TW","PT","mansfield","tsai")
multiExpr =list(TW=list(data=t(data[,1:56])), PT =list(data=t(data[,57:71])), man =list(data=t(data[,72:91])),tsai=list(data=t(data[,92:127])));

multiColor =list(TW= TW_mods,PT= PT_mods,man= man_mods,tsai=tsai_mods);

#Basic sizes of data
size = checkSets(multiExpr);
nSamples = size$nSamples;
nGenes = size$nGenes;
nSets = size$nSets;

############################
#calculate consensus network
############################

#soft thresholds were determined in ConsensusCoExpression.R
STPowers =c(8,10,14,12);

TOMinfo_vas = blockwiseIndividualTOMs(multiExpr,
                                      maxBlockSize = 40000,
                                      power= STPowers);

#This time we will save the adjacency matrix
print(system.time( {
  mods_vas = blockwiseConsensusModules(
    multiExpr,
    individualTOMInfo = TOMinfo_vas,
    power=8,
    consensusQuantile = 0,
    deepSplit = 2,
    reassignThresholdPS = 0,
    pamRespectsDendro = FALSE,
    detectCutHeight = 0.99, minModuleSize = 300,
    mergeCutHeight = 0.25, numericLabels = TRUE,
    saveIndividualTOMs=TRUE,
    individualTOMFileNames = "Data/results/individualTOM-Set%s-Block%b.RData",
    saveConsensusTOMs = FALSE,
    getTOMScalingSamples = TRUE,
    verbose = 3, indent = 2);
} ) );

collectGarbage();

#create list of all possible combinations of 4 data sets
cy<-list(c(1,2),c(1,3),c(1,4),c(2,3),c(2,4),c(3,4),c(1,2,3),c(1,3,4),c(2,3,4))

#create empty list to hold the results
sample_out<-list()


#need to do this manually and load adjcencies as needed because of memory allocation
#http://labs.genetics.ucla.edu/horvath/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/Consensus-NetworkConstruction-man.pdf

for(q in 1:9)
{
  TOMs<-list()
  if(cy[[q]][1]==1) 
  {
    load("Data/results/individualTOM-Set1-Block1.RData"); TOMs[[1]]<-tomDS; rm(tomDS);
  } else if(cy[[q]][1]==2) 
  {
    load("Data/results/individualTOM-Set2-Block1.RData"); TOMs[[1]]<-tomDS; rm(tomDS);
  }else if(cy[[q]][1]==3) 
  {
    load("Data/results/individualTOM-Set3-Block1.RData"); TOMs[[1]]<-tomDS; rm(tomDS);
  }
  
  if(cy[[q]][2]==2) 
  {
    load("Data/results/individualTOM-Set2-Block1.RData"); TOMs[[2]]<-tomDS; rm(tomDS);
  }else if(cy[[q]][2]==3) 
  {
    load("Data/results/individualTOM-Set3-Block1.RData"); TOMs[[2]]<-tomDS; rm(tomDS);
  }else if(cy[[q]][2]==4) 
  {
    load("Data/results/individualTOM-Set4-Block1.RData"); TOMs[[2]]<-tomDS; rm(tomDS);
  }
  
  if(!is.na(cy[[q]][3]))
  {
    if(cy[[q]][3]==3) 
    {
      load("Data/results/individualTOM-Set3-Block1.RData"); TOMs[[3]]<-tomDS; rm(tomDS);
    }else if(cy[[q]][3]==4) 
    {
      load("Data/results/individualTOM-Set4-Block1.RData"); TOMs[[3]]<-tomDS; rm(tomDS);
    }
  }
  
  scaleP=0.95
  set.seed(12345)
  nSets<-length(cy[[q]])
  nSamples=as.integer(1/(1-scaleP)*1000)
  
  scaleSample= sample(nGenes*(nGenes-1)/2,size=nSamples)
  
  TOMScalingSamples = list();
  
  scaleQuant =rep(1, nSets)
  
  scalePowers =rep(1, nSets)
  
  # Loop over sets
  for (set in 1:nSets)
  {
    # Select the sampled TOM entries
    TOMScalingSamples[[set]] = as.dist(TOMs[[set]])[scaleSample]
    # Calculate the 95th percentile
    scaleQuant[set] = quantile(TOMScalingSamples[[set]],
                               probs = scaleP, type = 8);
    
    # Scale the male TOM
    if (set>1)
    {
      scalePowers[set] = log(scaleQuant[1])/log(scaleQuant[set]);
      TOMs[[set]] = as.dist(TOMs[[set]])^scalePowers[set];
      collectGarbage()
    }
  }
  
  rm(TOMScalingSamples); rm(scalePowers); rm(scaleQuant);
  
  if(length(TOMs)==length(cy[[q]]) & length(TOMs)==2)
  {  
    consensusTOM = pmin(TOMs[[1]],TOMs[[2]]);
    rm(TOMs); 
    collectGarbage()
    # Clustering
    consTree = hclust(as.dist(1-consensusTOM), method = "average");
    # We like large modules, so we set the minimum module size relatively high:
    minModuleSize = 300;
    # Module identification using dynamic tree cut:
    unmergedLabels = cutreeDynamic(dendro = consTree, distM = lowerTri2matrix(1-consensusTOM),
                                   deepSplit = 2, cutHeight = 0.99,
                                   minClusterSize = minModuleSize,
                                   pamRespectsDendro = FALSE );
    unmergedColors = labels2colors(unmergedLabels)
    merge= mergeCloseModules(multiExpr[cy[[q]]], unmergedColors, cutHeight = 0.4, getNewUnassdME = TRUE, relabel = TRUE, consensusQuantile = 0.25);
    sample_out[[q]]<-table(merge$colors)
    rm(consensusTOM);
    collectGarbage()
  }else if(length(TOMs)==length(cy[[q]]) & length(TOMs)==3)
  {  
    consensusTOM = pmin(TOMs[[1]],TOMs[[2]],TOMs[[3]]);
    rm(TOMs);
    collectGarbage()
    # Clustering
    consTree = hclust(as.dist(1-consensusTOM), method = "average");
    # We like large modules, so we set the minimum module size relatively high:
    minModuleSize = 300;
    # Module identification using dynamic tree cut:
    unmergedLabels = cutreeDynamic(dendro = consTree, distM = lowerTri2matrix(1-consensusTOM),
                                   deepSplit = 2, cutHeight = 0.99,
                                   minClusterSize = minModuleSize,
                                   pamRespectsDendro = FALSE );
    unmergedColors = labels2colors(unmergedLabels)
    merge= mergeCloseModules(multiExpr[cy[[q]]], unmergedColors, cutHeight = 0.4, getNewUnassdME = TRUE, relabel = TRUE, consensusQuantile = 0.25);
    sample_out[[q]]<-table(merge$colors)
    rm(consensusTOM);
    collectGarbage()
  }
}

#save(sample_out, file="Data/results/subsamplingVascularDatasets.rdata")
load("Data/results/subsamplingVascularDatasets.rdata")

output<-data.frame(matrix(vector(), 0, 4, dimnames=list(c(), c("nSets","nModules","avgModGenes","totModGene"))), stringsAsFactors=F)

######
#summarize the sampling
for(st in 1:length(sample_out))
{
  ind<-which(names(unlist(sample_out[st]))!="grey")
  output[st,1]<-length(cy[[st]])
  output[st,2]<-length(ind)
  output[st,3]<-mean(sample_out[[st]][ind])
  output[st,4]<-sum(sample_out[[st]][ind])
}

ct=9
for(st in 1:length(multiColor))
{
  ind<-which(names(unlist(table(multiColor[st])))!="grey")
  output[(ct+st),1]<-1
  output[(ct+st),2]<-length(ind)
  output[(ct+st),3]<-mean(table(multiColor[st])[ind])
  output[(ct+st),4]<-sum(table(multiColor[st])[ind])
}

#get the consensus module assignments
load("Data/results/merge-labels.RData")

ind<-which(names(unlist(table(labels2colors(labels))))!="grey")
output[14,1]<-4
output[14,2]<-length(ind)
output[14,3]<-mean(table(labels2colors(labels))[ind])
output[14,4]<-sum(table(labels2colors(labels))[ind])  

####plot
require(plyr); require(ggplot2)

soutput<-ddply(output, .(nSets), summarise, nAvg.avg=mean(avgModGenes), nAvg.se=sd(avgModGenes)/sqrt(length(avgModGenes)), nTot.avg=mean(totModGene), nTot.se=sd(totModGene)/sqrt(length(totModGene)), nN.avg=mean(nModules), nN.se=sd(nModules)/sqrt(length(nModules)))

#rearrange soutput
soutput[5:12,]<-NA
soutput[5:8,1:3]<-soutput[1:4,c(1,4,5)]
soutput<-soutput[,-c(4,5)]
soutput[9:12,1:3]<-soutput[1:4,c(1,4,5)]
soutput<-soutput[,-c(4,5)]
soutput[1:4,4]<-"nAvgModGene"
soutput[5:8,4]<-"nTotModGene"
soutput[9:12,4]<-"nModules"

ggplot(soutput,aes(x=nSets,y=nAvg.avg)) + geom_line(stat = "identity", position="dodge") + geom_errorbar(aes(ymax=nAvg.avg+nAvg.se, ymin=nAvg.avg-nAvg.se), position=position_dodge(width=0.9),width=0.25,size=0.5)+theme_bw()+ facet_grid(V4~.,scales="free_y")

ggsave(file="Data/results/Sample_datasets.pdf",w=4,h=6)


########
# Significance of number of data sets on netowrk structure
summary(lm(log10(totModGene)~nSets, data=output))

summary(lm(log10(nModules)~nSets, data=output))

summary(lm(log10(avgModGenes)~nSets, data=output))

######
#predict the log network size size
model<-lm(log10(totModGene)~nSets, data=output)

ggplot(output,aes(x=nSets,y=log10(totModGene))) + geom_point() + geom_smooth(method = "lm", se = FALSE) + ylim(0,5) + xlim(1,20) + geom_abline(intercept = model$coefficients[1], slope = model$coefficients[2])
