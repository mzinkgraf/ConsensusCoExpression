
#################
#This script provides the code to generate Consensus CoExpression 
#analysis of Populus RNA-seq

#Matt Zinkgraf
#US Forest Service and UC Davis Computer Science

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

# #gravitropism experiment
# Data1<-read.table("Data/Gravitropism_normALLrpkm.txt",sep="\t",stringsAsFactors=F)
# gsg1<-goodSamplesGenes(t(Data1), verbose=3)
# Data1<-Data1[gsg1$goodGenes,]
# 
# 
# #PT Woody Tissue named PT_datExpr0
# load(file='Data/PT_WGCNA_data.rdata')
# Data2<-data.frame(t(PT_datExpr0))
# rm(PT_datExpr0); rm(PT_moduleColors)
# gsg2<-goodSamplesGenes(t(Data2), verbose=3)
# Data2<-Data2[gsg2$goodGenes,]
# 
# #provenance named datExpr
# load("Data/provenance_xylem.rdata")
# Data3<-data.frame(datExp0)
# rm(datExp0)
# gsg3<-goodSamplesGenes(t(Data3), verbose=3)
# Data3<-Data3[gsg3$goodGenes,]
# 
# #tsai named datExpr
# load("Data/Tsai_vascular_rpkm.rdata")
# Data4<-data.frame(datExp0)
# rm(datExp0)
# gsg4<-goodSamplesGenes(t(Data4), verbose=3)
# Data4<-Data4[gsg4$goodGenes,]
#  
# # #combine datasets
# data<-merge(Data1,Data2,by.x="row.names",by.y="row.names")
# data<-merge(data,Data3,by.x="Row.names",by.y="row.names")
# data<-merge(data,Data4,by.x="Row.names",by.y="row.names")
# row.names(data)<-data[,1]
# data<-data[,-1]
# 
# 
# #load module information from individual analyses
# 
# TW_output<-read.table("Data/tensionwood_WGCNAallGENES_moduleSize500.txt",sep="\t",header=T)
# row.names(TW_output)<-TW_output$gene
# TW_mods<-TW_output[row.names(data),2]
# 
# pt_output<-read.table("Data/PTRNAseq_WGCNAallGENES_modules.txt",sep="\t",header=T)
# row.names(pt_output)<-pt_output$gene
# PT_mods<-pt_output[row.names(data),2]
# 
# man_out<-read.table("Data/provenance_xylem_WGCNAallGENES_modules.txt",sep="\t",header=T)
# row.names(man_out)<-man_out$gene
# man_mods<-man_out[row.names(data),2]
# 
# tsai_out<-read.table("Data/Tsai_vascular_drought_module_output.txt",sep="\t",header=T)
# row.names(tsai_out)<-tsai_out$gene
# tsai_mods<-tsai_out[row.names(data),2]
# 
# save(data,TW_mods,PT_mods,man_mods,tsai_mods,file="Data/results/Grand_analysis_expressions_mods.rdata")

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

#Get soft threshold for individual experiments

powers = c(c(1:10), seq(from = 12, to=20, by=2))
# powerTables = vector(mode= "list" ,length = nSet);
# for(set in 1:nSets)
#   powerTables[[set]] =  list(data = pickSoftThreshold(multiExpr[[set]]$data, powerVector=powers,verbose = 2 )[[2]]);
# 
# # Save the results
# save(powerTables,file="Data/results/scaleFreeAnalysis-powerTables.RData");
load("Data/results/scaleFreeAnalysis-powerTables.RData")
collectGarbage();

#Re-format results for plotting
meanK = modelFit = matrix(0,length(powers), nSets);
for (set in 1:nSets)
{
  modelFit[, set] = -sign(powerTables[[set]]$data[,3])*powerTables[[set]]$data[,2];
  meanK[,set] = powerTables[[set]]$data[,5];
}
# Plot scatterplots of topology indices vs. soft-thresholding power
colors = c("black", "red", "blue", "green");


#################
#
#Plot Figure S1: Soft threshold analysis of individual data sets used in the consensus network
#
#################

pdf(file = "Data/results/MultiExp_scaleFreeTopologyAnalysis.pdf", wi = 8, h=6);
par(mfrow = c(1,2));
plot(powers, modelFit[, 1],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",
     pch = 21, col = 1, bg = 1,main = "Scale independence",ylim = range(modelFit));
addGrid();

# this line corresponds to using an R^2 cut-off of h
abline(h=0.85,col="red")
for(set in 2:nSets)
  points(powers, modelFit[, set], pch = 21, col = colors[set], bg = colors[set]);
legendClean("bottomright", legend = setLabels, pch = 21, col = colors);

# Plot of mean connectivity
plot(powers, meanK[, 1],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity",
     main = "Mean connectivity", pch = 21, col = 1, bg = 1);
addGrid();

for(set in 2:nSets)
  points(powers, meanK[, set], pch = 21, col = colors[set], bg = colors[set]);
legendClean("topright", legend = setLabels, pch = 21, col = colors);

# If plotting into a file, close it
dev.off();

############################
#calculate consensus network
############################

STPowers =c(8,10,14,12);

TOMinfo_vas = blockwiseIndividualTOMs(multiExpr,
                                  maxBlockSize = 40000,
                                  power= STPowers);

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
    saveTOMs = FALSE,
    saveConsensusTOMs = FALSE,
    getTOMScalingSamples = TRUE,
    verbose = 3, indent = 2);
} ) );

collectGarbage();

#Save the results for future use
#save(mods_vas,TOMinfo_vas, file="Data/results/Grand_mods.RData");
load("Data/results/Grand_mods.RData")

#merge similar modules
mergeCut = 0.4

merge= mergeCloseModules(multiExpr, mods_vas$unmergedColors, cutHeight = mergeCut, consensusQuantile = 0.25, getNewUnassdME = TRUE, relabel = TRUE);
labels=merge$colors;

#Save the resulting labels for future use
#save(merge, labels, file="Data/results/merge-labels.RData");
load("Data/results/merge-labels.RData")

###############
#
#Plot Figure 2A: Dendrogram and Consensus module identification
#
################

pdf(file="Data/results/consensusDendro.pdf",wi=10,5);
plotDendroAndColors(mods_vas$dendrograms[[1]], labels2colors(labels[mods_vas$goodGenes]),
                    "Consensus",
                    main ="Consensus modules",
                    dendroLabels = FALSE,
                    addGuide = TRUE, hang = 0.01, abHeight = 0.99,
                    guideHang = 0.03, marAll =c(1,5,1,1));
dev.off();



#Create Eigengene network for consensus modules
MEs0 = multiSetMEs(multiExpr, universalColors = labels2colors(labels))

MEs<-orderMEs(rbind(MEs0[[1]]$data,MEs0[[2]]$data,MEs0[[3]]$data,MEs0[[4]]$data))
row.names(MEs)<-lib_names

MEDiss = 1-cor(MEs);
# Cluster module eigengenes
METree = hclust(as.dist(MEDiss), method = "average");

###############
#
#Plot Figure 2B: Consensus and Individual Eigenegene dendrograms
#
################

pdf(file="Data/results/Experiments_ME_Dendro.pdf",w=10,h=2)
par(mfrow = c(1, 5))

plot(as.dendrogram(METree), horiz=FALSE, main = "Consensus", xlab = "", sub = "",cex = 0.8)

dataset_labels<-c("Gravitropism","Woody Tissues","Provenance","Drought")
for(r in 1:4)
{
  tmpMEDiss = 1-cor(orderMEs(MEs0[[r]]$data));
  # Cluster module eigengenes
  tmpMETree = hclust(as.dist(tmpMEDiss), method = "average");
  plot(as.dendrogram(tmpMETree), horiz=FALSE, main = dataset_labels[r], xlab = "", sub = "",cex = 0.8)
}
dev.off()


###############
#
#Plot Figure 2C: ME correlations for consensus and individual experiments
#
################

oMEs<-orderMEs(MEs0)
pdf(file="Data/results/Experiment_MEcor.pdf",w=10,h=2)
PlotExpPCsCor(oMEs,dataset_labels,IncludeSign=TRUE,setMargins=TRUE,plotConsensus=TRUE)
dev.off()


###############
#
#Plot Figure 2D: Identify global and experiment specific co-expression relationships
#
################

pdf(file="Data/results/Consensus_module_ExperMod_heatmap.pdf",w=10,h=3)
par(mfrow=c(1,length(multiColor)));
par(cex = 0.8);
for(l in 1:length(multiColor))
{
  # Convert the numeric module labels to color labels
  TWColors = multiColor[[l]][mods_vas$goodGenes]
  TWModules<-unique(TWColors)
  consColors= labels2colors(labels[mods_vas$goodGenes])
  consModules<-unique(consColors)
  
  # Numbers of individual and consensus modules
  nTWMods = length(TWModules)
  nConsMods = length(consModules)
  
  # Initialize tables of p-values and of the corresponding counts
  pTable = matrix(0, nrow = nTWMods, ncol = nConsMods);
  CountTbl = matrix(0, nrow = nTWMods, ncol = nConsMods);
  # Execute all pairwaise comparisons
  for (fmod in 1:nTWMods)
  {
    for (cmod in 1:nConsMods)
    {
      TWMembers = (TWColors == TWModules[fmod]);
      consMembers = (consColors == consModules[cmod]);
      pTable[fmod, cmod] = -log10(fisher.test(TWMembers, consMembers, alternative = "greater")$p.value);
      CountTbl[fmod, cmod] = sum(TWColors == TWModules[fmod] & consColors ==
                                   consModules[cmod])
    }
  }
  
  # Truncate small p values
  pTable[is.infinite(pTable)] = 1.3*max(pTable[is.finite(pTable)]);
  pTable[pTable>100 ] = 100 ;
  # Marginal counts (really module sizes)
  TWModTotals = apply(CountTbl, 1, sum)
  consModTotals = apply(CountTbl, 2, sum)
  
  # Use function labeledHeatmap to produce the color-coded table with all the trimmings
  gry<-which(consModules=="grey")
  labeledHeatmap(Matrix = t(pTable[,-gry]),
                 yLabels = paste(" ", consModules[-gry]),
                 xLabels = paste(" ", TWModules),
                 colorLabels = TRUE,
                 ySymbols = consModules[-gry],
                 xSymbols = TWModules,
                 colors = greenWhiteRed(100)[50:100],
                 main = dataset_labels[l],
                 cex.text = 0.8, cex.lab = 0.8, setStdMargins = FALSE);
  
  
}
dev.off()



###############
#
#Plot Figure 4: Plot the consensus eigengene network against treatments for each experiment
#
################

#######
require(ggplot2)
require(plyr)

#Gravitropism 
samples<-read.table("Data/GA_RNAandLibraries.txt",sep="\t",header=T,stringsAsFactors=F)
samples<-samples[-28,]
exper1<-read.table("Data/sample_and_library_names.csv",sep=",", row.names=1,stringsAsFactors=F,header=T)
TW_traits<-data.frame(cbind(c(exper1$Genotype,samples$Genotype),c(exper1$Sample1,samples$Sample1),c(exper1$sample2,samples$Treatment),c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2)))
names(TW_traits)<-c("genotype","woodtype","GA","year")
row.names(TW_traits)<-c(row.names(exper1),samples$Name)

TW_MEs<-merge(TW_traits,MEs,by.x="row.names",by.y="row.names")
row.names(TW_MEs)<-TW_MEs[,1]
TW_MEs<-TW_MEs[,-1]

allME<-NULL
for(j in 5:ncol(TW_MEs))
{
  EG<-TW_MEs[,c(1:4,j)]
  names(EG)[5]<-"EG"
  tmp.avg<-ddply(EG, .(genotype,woodtype,GA), summarise, tmp.avg=mean(EG), tmp.se=sd(EG)/sqrt(length(EG)))
  tmp.avg$genotype <- factor(tmp.avg$genotype, levels = c("miRNA2","717","35SARK2"))
  levels(tmp.avg$genotype)<-c("miRNA ARK2", "wild type", "OE ARK2")
  
  allME<-rbind(allME,cbind(names(TW_MEs)[j],tmp.avg))
}

#plot allME
names(allME)[1]<-"module"
allME$module<-factor(allME$module, levels = c("MEturquoise","MEbrown","MEyellow","MEblue","MEgrey"))

plot1<-ggplot(allME,aes(x=woodtype,y=tmp.avg, group=GA, colour=GA)) + geom_line() + geom_errorbar(aes(ymax=tmp.avg+tmp.se, ymin=tmp.avg-tmp.se), width=0.25,size=0.5)+theme_bw() + ylab("Gene expression (rpkm)") +xlab("wood type") + theme(text = element_text(size=8)) + theme(plot.background = element_blank() ,panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border = element_blank()) + theme(axis.line = element_line(color = 'black',size=0.5))+ scale_fill_grey( start=0.8,end=0.3) + scale_y_continuous(limits=c(-0.45,0.85)) + facet_grid(module ~ genotype)

#print plot
plot1  
ggsave(file="Data/results/GraviallMEline.pdf",plot=plot1,width=4,height=6)

#get significance of Gravitropism 
TW_ME_cor<-data.frame(matrix(vector(), 0, 4, dimnames=list(c(), c("genotype","woodtype","GA","year"))), stringsAsFactors=F)
for (e in 1:4)
{
  count=1
  for(k in 5:9)
  {
    TW_ME_cor[count,1:4]<-anova(lm(TW_MEs[,k]~genotype+woodtype+GA+year,data=TW_MEs))[1:4,5]
    row.names(TW_ME_cor)[count]<-names(TW_MEs)[k]
    count=count+1
  }
}


#############
#drought
sample_tsai<-read.table("Data/Tsai_SUT4_RNAseq.txt",sep="\t",header=T)

tsai_MEs<-merge(sample_tsai,MEs,by.x="Run_s",by.y="row.names")
row.names(tsai_MEs)<-tsai_MEs[,1]

#relevel genotype
tsai_MEs$genotype_s<-factor(tsai_MEs$genotype_s, levels = c("wild type","SUT4-RNAi transgenic"))

#relevel genotype
tsai_MEs$tissue_s<-factor(tsai_MEs$tissue_s, levels = c("bark","xylem"))

allME<-NULL
for(j in 33:ncol(tsai_MEs))
{
  EG<-tsai_MEs[,c(9,11,12,j)]
  names(EG)[4]<-"EG"
  tmp.avg<-ddply(EG, .(genotype_s,tissue_s,treatment_s), summarise, tmp.avg=mean(EG), tmp.se=sd(EG)/sqrt(length(EG)))
  
  tmp.avg$genotype_s<-factor(tmp.avg$genotype_s, levels = c("wild type","SUT4-RNAi transgenic"))
  
  #relevel genotype
  tmp.avg$tissue_s<-factor(tmp.avg$tissue_s, levels = c("bark","xylem"))
  
  
  allME<-rbind(allME,cbind(names(tsai_MEs)[j],tmp.avg))
}

#plot allME
names(allME)[1]<-"module"
allME$module<-factor(allME$module, levels = c("MEturquoise","MEbrown","MEyellow","MEblue","MEgrey"))

plot2<-ggplot(allME,aes(x=tissue_s,y=tmp.avg, group=treatment_s, colour=treatment_s)) + geom_line() + geom_errorbar(aes(ymax=tmp.avg+tmp.se, ymin=tmp.avg-tmp.se), width=0.25,size=0.5)+theme_bw() + ylab("Gene expression (rpkm)") +xlab("wood type") + theme(text = element_text(size=8)) + theme(plot.background = element_blank() ,panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border = element_blank()) + theme(axis.line = element_line(color = 'black',size=0.5))+ scale_fill_grey( start=0.8,end=0.3) + scale_y_continuous(limits=c(-0.45,0.85)) + facet_grid(module ~ genotype_s)

#print plot
plot2
ggsave(file="Data/results/Drought_allMEline.pdf",plot=plot2,width=4.5,height=6)

# get significance of Drought
tsai_ME_cor<-data.frame(matrix(vector(), 0, 3, dimnames=list(c(), c("genotype","tissue","treatment"))), stringsAsFactors=F)
for (e in 1:3)
{
  count=1
  for(k in 33:ncol(tsai_MEs))
  {
    tsai_ME_cor[count,1:3]<-anova(lm(tsai_MEs[,k]~genotype_s+tissue_s+treatment_s,data=tsai_MEs))[1:3,5]
    row.names(tsai_ME_cor)[count]<-names(tsai_MEs)[k]
    count=count+1
  }
}


#############
#PT Woody Tissues
tis<-c("phloem","phloem","xylem","phloem","xylem","phloem","xylem","phloem","xylem","phloem","xylem","phloem","xylem","phloem","xylem")
tree<-c("t0","t0","t0","t1","t1","t2","t2","t3","t3","t4","t4","t5","t5","t6","t6")
pt_lib<-names(data)[57:71]
PT_traits<-data.frame(cbind(pt_lib,tis,tree))

PT_MEs<-merge(PT_traits,MEs,by.x="pt_lib",by.y="row.names")
row.names(PT_MEs)<-PT_MEs[,1]
PT_MEs<-PT_MEs[,-1]

allME<-NULL
for(j in 3:ncol(PT_MEs))
{
  EG<-PT_MEs[,c(1,2,j)]
  names(EG)[3]<-"EG"
  tmp.avg<-ddply(EG, .(tis), summarise, tmp.avg=mean(EG), tmp.se=sd(EG)/sqrt(length(EG)))
  
  
  allME<-rbind(allME,cbind(names(PT_MEs)[j],tmp.avg))
}

#plot allME
names(allME)[1]<-"module"
allME$module<-factor(allME$module, levels = c("MEturquoise","MEbrown","MEyellow","MEblue","MEgrey"))

plot3<-ggplot(allME,aes(x=tis,y=tmp.avg,group=1)) + geom_line() + geom_errorbar(aes(ymax=tmp.avg+tmp.se, ymin=tmp.avg-tmp.se), width=0.25,size=0.5)+theme_bw() + ylab("Gene expression (rpkm)") +xlab("P. trichocarpa tissue") + theme(text = element_text(size=8)) + theme(plot.background = element_blank() ,panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border = element_blank()) + theme(axis.line = element_line(color = 'black',size=0.5))+ scale_fill_grey( start=0.8,end=0.3) + scale_y_continuous(limits=c(-0.45,0.85)) + facet_grid(module ~ .)

#print plot
plot3
ggsave(file="Data/results/WoodyTissue_allMEline.pdf",plot=plot3,width=2.5,height=6)

#get significance of woody tissues
PT_ME_cor<-data.frame(matrix(vector(), 0, 2, dimnames=list(c(), c("tissue","tree"))), stringsAsFactors=F)
for (e in 1:2)
{
  count=1
  for(k in 3:ncol(PT_MEs))
  {
    PT_ME_cor[count,1:2]<-anova(lm(PT_MEs[,k]~tis+tree,data=PT_MEs))[1:2,5]
    row.names(PT_ME_cor)[count]<-names(PT_MEs)[k]
    count=count+1
  }
}


#############
#Provenance

sample_loc<-read.table("Data/provenance_sample_locations.txt",sep="\t",header=T)

man_MEs<-merge(sample_loc,MEs,by.x="row.names",by.y="row.names")
row.names(man_MEs)<-man_MEs[,1]
man_MEs<-man_MEs[,-c(1,2)]

allME<-NULL
for(j in 5:ncol(man_MEs))
{
  
  allME<-rbind(allME,cbind(names(man_MEs)[j],man_MEs[,j],"Longitude",man_MEs$Longitude))
  allME<-rbind(allME,cbind(names(man_MEs)[j],man_MEs[,j],"Latitude",man_MEs$Latitude))
}
allME<-as.data.frame(allME)

#plot allME
names(allME)<-c("module","EG","meassure", "coordinate")

allME$EG<-as.numeric(allME$EG)
allME$coordinate<-as.numeric(allME$coordinate)
allME$module<-factor(allME$module, levels = c("MEturquoise","MEbrown","MEyellow","MEblue","MEgrey"))

plot4<-ggplot(allME,aes(x=coordinate,y=EG,group=meassure,colour=meassure)) + geom_point(size=1) +theme_bw() + stat_smooth(method="glm", se=TRUE, colour="black",size=0.5)+ ylab("Gene expression (rpkm)") +xlab("P. trichocarpa Provenance") + theme(text = element_text(size=8)) + theme(plot.background = element_blank() ,panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border = element_blank()) + theme(axis.line = element_line(color = 'black',size=0.5))+ scale_fill_grey( start=0.8,end=0.3) + facet_grid(module ~ meassure,scale="free_x")+ scale_y_continuous(limits=c(-0.45,0.85))

#print plot
plot4
ggsave(file="Data/results/Provenance_allMEline.pdf",plot4,width=4,height=6)

#get significance of provenance
man_ME_cor<-data.frame(matrix(vector(), 0, 3, dimnames=list(c(), c("longitude","latitude","year"))), stringsAsFactors=F)
for (e in 1:2)
{
  count=1
  for(k in 5:ncol(man_MEs))
  {
    man_ME_cor[count,1:3]<-anova(lm(man_MEs[,k]~Longitude + Latitude + year,data=man_MEs))[1:3,5]
    row.names(man_ME_cor)[count]<-names(man_MEs)[k]
    count=count+1
  }
}

###############################
#
# GO analysis of Populus genes using AT best blast hit
# Potri to AT is obtained from the phytozome annotation Ptrichocarpa_210
#
###############################

source("R/GO_populus.R")
GOtems<-read.csv("Data/GO_terms.txt",sep="\t",header=F)
pt<-read.table("Data/Ptrichocarpa_210_annotation_primary.txt",sep="\t")
consColors = labels2colors(labels[mods_vas$goodGenes])

out.anno<-data.frame(cbind(row.names(data),consColors))
out.anno2<-merge(out.anno,pt,by.x="V1",by.y="V2")
 
# Iblue<-which(out.anno2$consColors=="blue")
# GOblue<-atGOanalysis(out.anno2$V10[Iblue])
# #View(summary(GOblue$BP))
# 
# Ibrown<-which(out.anno2$consColors=="brown")
# GObrown<-atGOanalysis(out.anno2$V10[Ibrown])
# #View(summary(GObrown$BP))
# 
# Iturquoise<-which(out.anno2$consColors=="turquoise")
# GOturquoise<-atGOanalysis(out.anno2$V10[Iturquoise])
# 
# Iyellow<-which(out.anno2$consColors=="yellow")
# GOyellow<-atGOanalysis(out.anno2$V10[Iyellow])
# 
# Igrey<-which(out.anno2$consColors=="grey")
# GOgrey<-atGOanalysis(out.anno2$V10[Igrey])
#  
# save(GOblue,GObrown,GOturquoise,GOyellow,GOgrey,file="Data/results/Module_GO.rdata")

load("Data/results/Module_GO.rdata")

require(xlsx)
write.xlsx(summary(GOblue$BP),file="Data/results/Supplementary_Table_2.xlsx",sheetName = "GOblue_BP")
write.xlsx(summary(GOblue$MF),file="Data/results/Supplementary_Table_2.xlsx",sheetName = "GOblue_MF",append = T)
write.xlsx(summary(GOblue$CC),file="Data/results/Supplementary_Table_2.xlsx",sheetName = "GOblue_CC",append = T)

write.xlsx(summary(GObrown$BP),file="Data/results/Supplementary_Table_2.xlsx",sheetName = "GObrown_BP",append=TRUE)
write.xlsx(summary(GObrown$MF),file="Data/results/Supplementary_Table_2.xlsx",sheetName = "GObrown_MF",append=TRUE)
write.xlsx(summary(GObrown$CC),file="Data/results/Supplementary_Table_2.xlsx",sheetName = "GObrown_CC",append=TRUE)

write.xlsx(summary(GOturquoise$BP),file="Data/results/Supplementary_Table_2.xlsx",sheetName = "GOturquoise_BP",append=TRUE)
write.xlsx(summary(GOturquoise$MF),file="Data/results/Supplementary_Table_2.xlsx",sheetName = "GOturquoise_MF",append=TRUE)
write.xlsx(summary(GOturquoise$CC),file="Data/results/Supplementary_Table_2.xlsx",sheetName = "GOturquoise_CC",append=TRUE)

write.xlsx(summary(GOyellow$BP),file="Data/results/Supplementary_Table_2.xlsx",sheetName = "GOyellow_BP",append=TRUE)
write.xlsx(summary(GOyellow$MF),file="Data/results/Supplementary_Table_2.xlsx",sheetName = "GOyellow_MF",append=TRUE)
write.xlsx(summary(GOyellow$CC),file="Data/results/Supplementary_Table_2.xlsx",sheetName = "GOyellow_CC",append=TRUE)

###############
#
#Plot Figure 5: Functional enrichment of consensus module using GO analysis
#buid heatmap showing differences in GO enrichment for auxin, cell-wall, hormone and meristem
#
################


key<-"auxin|hormone|gibberelli|brassino|cytokinin|cell wall|meristem|protein localization|methylation|histone|chromatin|cellulose|lignin|xylan|xylose"

blueGO<-summary(GOblue$BP)[grep(key,summary(GOblue$BP)[,7],perl=TRUE),1:2]
names(blueGO)[2]<-"blue"

brownGO<-summary(GObrown$BP)[grep(key,summary(GObrown$BP)[,7],perl=TRUE),1:2]
names(brownGO)[2]<-"brown"

turquoiseGO<-summary(GOturquoise$BP)[grep(key,summary(GOturquoise$BP)[,7],perl=TRUE),1:2]
names(turquoiseGO)[2]<-"turquoise"

yellowGO<-summary(GOyellow$BP)[grep(key,summary(GOyellow$BP)[,7],perl=TRUE),1:2]
names(yellowGO)[2]<-"yellow"

greyGO<-summary(GOgrey$BP)[grep(key,summary(GOgrey$BP)[,7],perl=TRUE),1:2]
names(greyGO)[2]<-"grey"

GOtable<-merge(blueGO,brownGO,by.x="GOBPID",by.y="GOBPID",all=T)
GOtable<-merge(GOtable,turquoiseGO,by.x="GOBPID",by.y="GOBPID",all=T)
GOtable<-merge(GOtable,yellowGO,by.x="GOBPID",by.y="GOBPID",all=T)
GOtable<-merge(GOtable,greyGO,by.x="GOBPID",by.y="GOBPID",all=T)
GOtable<-merge(GOtable,GOtems,by.x="GOBPID",by.y="V1",all.x=T)

#generate order

o3<-grep("hormone|gibberelli|brassino|auxin|cytokinin",GOtable$V3)
o4<-grep("cell wall",GOtable$V3)
o5<-grep("meristem",GOtable$V3)
o2<-grep("protein localization",GOtable$V3)
o7<-grep("methylation|histone|chromatin",GOtable$V3)
o6<-grep("cellulose|lignin|xylan|xylose|glucomannan",GOtable$V3)
o<-c(o3,o4,o6,o5,o2,o7)

#convert to -log10(pvalue)
GOtable[is.na(GOtable)] <- 1
GOtable[,2:6]<--log10(GOtable[,2:6])


results<-data.frame(t(GOtable[o,2:6]))
names(results)<-GOtable[o,1]

#reorder modules
results<-results[c("turquoise","brown","yellow","blue","grey"),]

my_palette <- colorRampPalette(c("white", "red"))(n = 8)

vlines<-cumsum(c(length(o3),length(o4),length(o6),length(o5),length(o2),length(o7)))
#make list of group names
GOgroups<-rep(NA,ncol(results))
GOgroups[round(length(o3)/2)]<-"hormone"
GOgroups[round(length(o4)/2)+vlines[1]]<-"cell wall"
GOgroups[round(length(o6)/2)+vlines[2]]<-"cellulose | lignin | xylose"
GOgroups[round(length(o5)/2)+vlines[3]]<-"meristem"
GOgroups[round(length(o2)/2)+vlines[4]]<-"protein localization"
GOgroups[round(length(o7)/2)+vlines[5]]<-"chromatin | histone | methylation"

pdf(file="Data/results/Module_Go_enrichment_wGery.pdf",width=8,height=4)
par(mar = c(9, 8, 2, 2));
labeledHeatmap(Matrix = results,
               xLabels = GOgroups,
               yLabels = paste("ME",row.names(results),sep=""),
               yColorLabels = TRUE,
               yColorWidth = 0.05,
               ySymbols = row.names(results),
               colors = my_palette,
               #textMatrix = txt,
               setStdMargins = FALSE,
               cex.text = 0.5,
               cex.lab.y = 0.8,
               cex.lab.x = 1,
               xLabelsAngle = 45,
               zlim = c(0,15),
               main = paste("GO Enrichment of Co-expression Modules"),
               verticalSeparator.x = vlines
)

dev.off()



###############
#
#Plot Figure 6: Enrichment of DNase-seq footprints and binding of transcription factors within consensus modules. 
#
################



########################################
#chip enrichment of modules
########################################

tf_data<-read.table("Data/TF_binding.txt",sep="\t",header=T)
row.names(tf_data)<-tf_data[,1]
tf<-tf_data[row.names(data),]
tf<-tf[,c("gene","ark1B","ark2B","blrB","pcnB","revB","polB")]

#create empty dataframe
tf_binding<-data.frame(matrix(vector(), 0, 7, dimnames=list(c(), c("Module","ark1","ark2","blr","pcn","rev","pol2"))), stringsAsFactors=F)

for (j in 2:ncol(tf)){
  
  ugene<-row.names(tf)[which(tf[,j]!=0)]
  
  n_chip<-length(ugene)
  n_no_chip<-nrow(tf)-n_chip
  count=1
  for (i in unique(out.anno[,2]))
  {
    n1<-length(unique(out.anno[which(out.anno$consColors==i),1]))
    n2<-length(intersect(unique(out.anno[which(out.anno$consColors==i),1]), ugene))
    p<-phyper(n2,n_chip,n_no_chip,n1,lower.tail=FALSE)
    #print(n2/n1)
    #print(p)
    tf_binding[count,1]<-i
    tf_binding[count,j]<-p
    count<-count+1
  }
}

row.names(tf_binding)<-paste("ME",tf_binding[,1],sep="")
tf_binding<-tf_binding[c("MEturquoise","MEbrown","MEyellow","MEblue","MEgrey"),]

pdf(file="Data/results/consensus_modules_TFbinding_heatmap.pdf",w=8,h=6)
par(mar=c(5,10,5,1))
labeledHeatmap(Matrix = -log10(tf_binding[,2:ncol(tf_binding)]),
               xLabels = names(tf_binding)[2:ncol(tf_binding)],
               yLabels = paste("ME",tf_binding[,1],sep=""),
               ySymbols = paste("ME",tf_binding[,1],sep=""),
               colorLabels = TRUE,
               colors = colorRampPalette(c("white", "red"))(n = 10),
               #textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(0,15),
               main = paste("Module-binding relationships"))
dev.off()


########################################
#DNase enrichment of modules
########################################

FPG<-read.table("Data/DNase_FP_TargetGenes.txt",sep="\t",header=T)

FPB<-ddply(FPG,.(V1),summarise, n=length(V1), b=if(length(V1>0)) 1 else 0, OverlapStart=if(length(which(V8=="OverlapStart"))>0) 1 else 0, upstream=if(length(which(V8=="upstream"))>0) 1 else 0, downstream=if(length(which(V8=="downstream"))>0) 1 else 0, OverlapEnd=if(length(which(V8=="OverlapEnd"))>0) 1 else 0,inside=if(length(which(V8=="inside"))>0) 1 else 0, OverlapAll=if(length(which(V8=="OverlapAll"))>0) 1 else 0 )

names(FPB)<-c("V1", "fpN", "fpB", "fpOS", "fpuS", "fpdS", "fpOE", "fpin", "fpOA")

#add remaining genes that did not have a target
pt_genes<-read.table("Data/pt210_gene_features.txt",sep="\t",header=T)
tmp<-merge(FPB,pt_genes,by.x="V1",by.y="geneID",all.y=T)

tmp[which(is.na(tmp[,2])),2:9]=0

row.names(tmp)<-tmp[,1]
fp<-tmp[out.anno$V1,]
fp<-fp[,-c(2,3,10,11,12,13)]

#create empty dataframe
fp_binding<-data.frame(matrix(vector(), 0, 7, dimnames=list(c(), c("Module","fpOS", "fpuS", "fpdS", "fpOE", "fpin", "fpOA"))), stringsAsFactors=F)

for (j in 2:ncol(fp)){
  
  ugene<-row.names(fp)[which(fp[,j]!=0)]
  
  n_chip<-length(ugene)
  n_no_chip<-nrow(fp)-n_chip
  count=1
  for (i in unique(out.anno[,2]))
  {
    n1<-length(unique(out.anno[which(out.anno$consColors==i),1]))
    n2<-length(intersect(unique(out.anno[which(out.anno$consColors==i),1]), ugene))
    p<-phyper(n2,n_chip,n_no_chip,n1,lower.tail=FALSE)
    
    fp_binding[count,1]<-i
    fp_binding[count,j]<-p
    count<-count+1
  }
}

row.names(fp_binding)<-paste("ME",fp_binding[,1],sep="")
fp_binding<-fp_binding[c("MEturquoise","MEbrown","MEyellow","MEblue","MEgrey"),c("Module","fpuS","fpOS","fpin","fpOE","fpdS","fpOA")]

pdf(file="Data/results/consensus_modules_DNase_enrich_heatmap.pdf",w=8,h=6)
par(mar=c(5,10,5,1))
labeledHeatmap(Matrix = -log10(fp_binding[,2:ncol(fp_binding)]),
               xLabels = names(fp_binding)[2:ncol(fp_binding)],
               yLabels = paste("ME",fp_binding[,1],sep=""),
               ySymbols = paste("ME",fp_binding[,1],sep=""),
               colorLabels = TRUE,
               colors = colorRampPalette(c("white", "red"))(n = 10),
               #textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(0,15),
               main = paste("Module-DNase relationships"))
dev.off()


###############
#
#Plot Figure 7: Frequency of co-occurrence of paralogous genes within the same versus different gene co-expression modules.
#
################

row.names(out.anno)<-out.anno[,1]
para<-read.table("Data/version3_duplications.out",sep=",",header=T)

para[,9]<-out.anno[para[,1],2]
para[,10]<-out.anno[para[,2],2]
para[is.na(para)]<-"grey"

tmp<-chisq.test(para[,9],para[,10])

ord<-c("turquoise","brown","yellow","blue","grey")
resid<-as.matrix(tmp$stdres)[ord,ord]

pdf(file="Data/results/Paralog_module.pdf",w=8,h=6)
labeledHeatmap(Matrix = resid,
               xLabels =dimnames(resid)[[1]],
               yLabels =dimnames(resid)[[2]],
               ySymbols =dimnames(resid)[[2]],
               colorLabels = FALSE,
               colors= greenWhiteRed(50),
               setStdMargins = FALSE,
               cex.text= 0.5,
               zlim =c(-10,10))
dev.off()


#######################
#
#Plot Figure 8: Enrichment of censensus modules for SNPs from two GWAS studies in Populus (Porth et al., 2013; McKowm et al., 2014)
#
#######################

load("Data/GWAS_associations.rdata")
x<-unique(associations[,1])
y<-unique(associations[,2])

associations<-merge(associations,out.anno,by.x="V1",by.y="V1",all.x=T)
associations[is.na(associations[,4]),4]<-"grey"

load("Data/34KarrayGenes.rdata")
colors<-labels2colors(labels)
names(colors)<-row.names(data)

#reorder porth traits
pr<-y[1:16]
pr<-pr[order(tolower(pr))]
#reorder Mcknown trait
mn<-y[17:36]
mn<-mn[order(mn)]
trs<-c(pr,mn)

Gt<-length(colors[UarrayGenes])
Mt<-table(colors[UarrayGenes])

GWAS_enrich<-data.frame(matrix(vector(), length(trs), 5, dimnames=list(trs, c("blue","brown","turquoise","yellow","grey"))), stringsAsFactors=F)

for(i in 1:length(trs))
{
  n1<-which(associations$Trait==trs[i])
  tb<-table(associations[n1,4])
  
  if("blue" %in% names(tb)) GWAS_enrich[i,1]<-phyper(tb[["blue"]], Mt[["blue"]], (Gt-Mt[["blue"]]), length(n1),lower.tail=F) else GWAS_enrich[i,1]<-NA
  
  if("brown" %in% names(tb)) GWAS_enrich[i,2]<-phyper(tb[["brown"]], Mt[["brown"]], (Gt-Mt[["brown"]]), length(n1),lower.tail=F) else GWAS_enrich[i,2]<-NA
  
  if("turquoise" %in% names(tb)) GWAS_enrich[i,3]<-phyper(tb[["turquoise"]], Mt[["turquoise"]], (Gt-Mt[["turquoise"]]), length(n1),lower.tail=F) else GWAS_enrich[i,3]<-NA
  
  if("yellow" %in% names(tb)) GWAS_enrich[i,4]<-phyper(tb[["yellow"]], Mt[["yellow"]], (Gt-Mt[["yellow"]]), length(n1),lower.tail=F) else GWAS_enrich[i,4]<-NA
  
  if("grey" %in% names(tb)) GWAS_enrich[i,5]<-phyper(tb[["grey"]], Mt[["grey"]], (Gt-Mt[["grey"]]), length(n1),lower.tail=F) else GWAS_enrich[i,5]<-NA
  
}

GWAS_enrich[is.na(GWAS_enrich)]<-1
GWAS_enrich<-GWAS_enrich[,c("turquoise","brown","yellow","blue","grey")]

pdf(file="Data/results/arrayGenesGWAS_modules_enrichment.pdf",width=8,height=4)
par(mar = c(6, 8.5, 3, 3));
labeledHeatmap(Matrix = -log10(t(GWAS_enrich[trs,])),
               xLabels = row.names(GWAS_enrich),
               yLabels = paste("ME",names(GWAS_enrich),sep=""),
               yColorLabels = TRUE,
               yColorWidth = 0.05,
               ySymbols = paste("ME",names(GWAS_enrich),sep=""),
               colors = colorRampPalette(c("white", "blue"))(n = 10),
               #textMatrix = txt,
               setStdMargins = FALSE,
               cex.text = 0.5,
               cex.lab.y = 0.8,
               cex.lab.x = 0.5,
               xLabelsAngle = 90,
               zlim = c(0,15),
               main = paste("GWAS_enrichment"),
               verticalSeparator.x = length(pr)
)

dev.off()


########################
#
#Plot Figure S2: Network connectivity of consensus modules
#
########################

require(reshape)
kme<-consensusKME(multiExpr,labels,multiEigengenes = MEs0,consensusQuantile = 0.25,signed=F,excludeGrey=F)
conKME<-cbind(labels2colors(labels),kme)

tmp<-apply(conKME, 1, function(x) x[ which(names(conKME)==paste("consensus.kME",x[1],sep=""))] )

conKMEcombined<-as.data.frame(cbind(labels2colors(labels),tmp))
names(conKMEcombined)<-c("module","kME")

conKMEcombined[,1]<-factor(conKMEcombined[,1],levels=c("turquoise","brown","yellow","blue","grey"))

pdf(file="Data/results/Consensus_module_kME.pdf",h=4,w=6)
boxplot(abs(as.numeric(conKMEcombined[,2]))~conKMEcombined[,1],col=levels(conKMEcombined[,1]), xlab="Consensus Modules", ylab="Gene Connectivity")
dev.off()


################################
# generate supplementary table 1
################################

out<-cbind(row.names(data),multiColor$TW,multiColor$tsai,multiColor$man,multiColor$PT,consColors,conKMEcombined$kME)
out<-merge(out,pt,by.x="V1",by.y="V2")

out<-out[,-c(8,9,14)]
names(out)<-c("Gene", "Gravitropism_Modules",	"Drought_Modules",	"Provenance_Modules",	"WoodyTissue_Modules",	"Consensus_Modules",	"Consensus_kME",	"pFAM",	"Panther",	"KOG",	"EC",	"PT_GO",	"ATG",	"Alias",	"Description")

mapping<-read.table("Data/MAPPINGatGO2PotriV3.txt",sep="\t",header=T)
atGO<-data.frame(ddply(mapping, .(V2), summarise, atGOterms=list(frame.go_id)))
row.names(atGO)<-atGO[,1]

out[,16]<-apply(out,1,function(x) paste(atGO[sub("(AT\\d+G\\d+)\\.\\d+","\\1",x[13],perl=TRUE),2][[1]], collapse = "; "))
names(out)[16]<-"AT_GO"

write.table(out[,c(1:12,16,13:15)],file="Data/results/Supplementary_Table_1.txt",sep="\t")
