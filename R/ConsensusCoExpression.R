#################
#Consensus CoExpression analysis of Populus RNA-seq 
library(edgeR)
library(WGCNA) 
library(RColorBrewer)
library(preprocessCore)
library(fields)
source("R/networkFunctions-extras-05.R")
source("R/functions.R")
#source("R/NetworkFunctions-Mouse.R") #additional functions can be found here but they are not used in this script

allowWGCNAThreads(n=2);
options(stringsAsFactors = FALSE);

# #gravitropism experiment
# Data1<-read.table("Data/Gravitropism_normALLrpkm.txt",sep="\t",stringsAsFactors=F)
# 
# #PT Woody Tissue named PT_datExpr0
# load(file='Data/PT_WGCNA_data.rdata')
# Data2<-data.frame(t(PT_datExpr0))
# rm(PT_datExpr0); rm(PT_moduleColors)
# 
# #provenance named datExpr
# load("Data/provenance_xylem.rdata")
# Data3<-data.frame(datExp0)
# rm(datExp0)
# 
# #tsai named datExpr
# load("Data/Tsai_rpkm.rdata")
# Data4<-data.frame(datExp0)
# rm(datExp0)
#  
# # #combine datasets
# data<-merge(Data1,Data2,by.x="row.names",by.y="row.names")
# data<-merge(data,Data3,by.x="Row.names",by.y="row.names")
# data<-merge(data,Data4,by.x="Row.names",by.y="row.names")
# row.names(data)<-data[,1]
# data<-data[,-1]
# 
# #check genes
# gsg<-goodSamplesGenes(t(data), verbose=3)
# gsg$allOK
# 
# data<-data.frame(data[gsg$goodGenes,]) 
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
# tsai_out<-read.table("Data/Tsai_drought_module_output.txt",sep="\t",header=T)
# row.names(tsai_out)<-tsai_out$gene
# tsai_mods<-tsai_out[row.names(data),2]
# 
# save(data,TW_mods,PT_mods,man_mods,tsai_mods,file="Data/results/Grand_analysis_expressions_mods.rdata")

#load expression and modules
load("Data/results/Grand_analysis_expressions_mods.rdata")

lib_names<-names(data)

#build multiset
setLabels=c("TW","PT","mansfield","tsai")
multiExpr =list(TW=list(data=t(data[,1:56])), PT =list(data=t(data[,57:71])), man =list(data=t(data[,72:91])),tsai=list(data=t(data[,92:163])));

multiColor =list(TW= TW_mods,PT= PT_mods,man= man_mods,tsai=tsai_mods);


#Basic sizes of data
size = checkSets(multiExpr);
nSamples = size$nSamples;
nGenes = size$nGenes;
nSets = size$nSets;

powers =c(seq(1,12, by=1));
powerTables = vector(mode= "list" ,length = nSets);
for(set in 1:nSets)
  powerTables[[set]] =  list(data = pickSoftThreshold(multiExpr[[set]]$data, powerVector=powers,verbose = 2 )[[2]]);

# Save the results
#save(powerTables,file="Data/results/scaleFreeAnalysis-powerTables.RData");
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
colors = c("black", "red", "blue");
sizeGrWindow(10, 8);

pdf(file = "Data/results/scaleFreeTopologyAnalysis.pdf", wi = 10, h=8);
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

STPowers =c(8,10,12,10);

TOMinfo = blockwiseIndividualTOMs(multiExpr,
                                  maxBlockSize = 40000,
                                  power= STPowers);

print(system.time( {
  mods = blockwiseConsensusModules(
    multiExpr,
    individualTOMInfo = TOMinfo,
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
#save(mods,TOMinfo, file="Data/results/Grand_mods.RData");
load("Data/results/Grand_mods.RData")

mergeCut = 0.4
#This function is not working properly
merge= mergeCloseModules(multiExpr, mods$unmergedColors, cutHeight = mergeCut, consensusQuantile = 0.25, getNewUnassdME = TRUE, relabel = TRUE);
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
plotDendroAndColors(mods$dendrograms[[1]], labels2colors(labels[mods$goodGenes]),
                    "Consensus",
                    main ="Consensus modules",
                    dendroLabels = FALSE,
                    addGuide = TRUE, hang = 0.01, abHeight = 0.99,
                    guideHang = 0.03, marAll =c(1,5,1,1));
dev.off();



#Create Eigengene network for Consensus modules
MEs0 = multiSetMEs(multiExpr, universalColors = labels2colors(labels))

MEs<-orderMEs(rbind(MEs0[[1]]$data,MEs0[[2]]$data,MEs0[[3]]$data,MEs0[[4]]$data))
row.names(MEs)<-lib_names

MEDiss = 1-cor(MEs);
# Cluster module eigengenes
METree = hclust(as.dist(MEDiss), method = "average");

###############
#
#Plot Figure 2B: Consensus and Individual Eigenegene dendrogram
#
################

pdf(file="Data/results/Experiments_ME_Dendro.pdf",w=10,h=2)
par(mfrow = c(1, 5))

plot(as.dendrogram(METree), horiz=FALSE, main = "Consensus", xlab = "", sub = "",cex = 0.8)

dataset_labels<-c("Gravitropism","Drought","Provenance","Woody Tissues")
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
#Plot Figure 2C: Ordered ME correlations for Consensus and individual experiments
#
################

oMEs<-orderMEs(MEs0)
pdf(file="Data/results/Experiment_MEcor.pdf",w=10,h=2)
PlotExpPCsCor(oMEs,dataset_labels,IncludeSign=TRUE,setMargins=TRUE,plotConsensus=TRUE)
dev.off()


###############
#
#Plot Figure 2D: Identify Global and Experiment specific co-expression relationships
#
################

pdf(file="Data/results/Consensus_module_ExperMod_heatmap.pdf",w=10,h=3)
par(mfrow=c(1,length(multiColor)));
par(cex = 0.8);
for(l in 1:length(multiColor))
{
  # Convert the numeric module labels to color labels
  TWColors = multiColor[[l]][mods$goodGenes]
  TWModules<-unique(TWColors)
  consColors= labels2colors(labels[mods$goodGenes])
  consModules<-unique(consColors)
  
  # Numbers of TW and consensus modules
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
  
  # Truncate p values smaller than 10^{-50} to 10^{-50} 
  pTable[is.infinite(pTable)] = 1.3*max(pTable[is.finite(pTable)]);
  pTable[pTable>100 ] = 100 ;
  # Marginal counts (really module sizes)
  TWModTotals = apply(CountTbl, 1, sum)
  consModTotals = apply(CountTbl, 2, sum)
  
  #par(mar=c(20, 15, 2.7, 1)+0.3);
  # Use function labeledHeatmap to produce the color-coded table with all the trimmings
  gry<-which(consModules=="grey")
  #par(mar=c(10,15,5,1))
  labeledHeatmap(Matrix = t(pTable[,-gry]),
                 yLabels = paste(" ", consModules[-gry]),
                 xLabels = paste(" ", TWModules),
                 colorLabels = TRUE,
                 #ySymbols = paste("Cons ", consModules, ": ", consModTotals, sep=""),
                 #xSymbols = paste(setLabels[l]," ", TWModules, ": ", TWModTotals, sep=""),
                 #textMatrix = CountTbl,
                 ySymbols = consModules[-gry],
                 xSymbols = TWModules,
                 colors = greenWhiteRed(100)[50:100],
                 main = dataset_labels[l],
                 cex.text = 0.8, cex.lab = 0.8, setStdMargins = FALSE);
  
  
}
dev.off()



###############
#
#Plot Figure 4: Plot the Consensus eigengene network against treatments for each experiment
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
allME$module<-factor(allME$module, levels = c("MEyellow","MEblue","MEturquoise","MEgrey","MEbrown"))

plot1<-ggplot(allME,aes(x=woodtype,y=tmp.avg, group=GA, colour=GA)) + geom_line() + geom_errorbar(aes(ymax=tmp.avg+tmp.se, ymin=tmp.avg-tmp.se), width=0.25,size=0.5)+theme_bw() + ylab("Gene expression (rpkm)") +xlab("wood type") + theme(text = element_text(size=8)) + theme(plot.background = element_blank() ,panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border = element_blank()) + theme(axis.line = element_line(color = 'black',size=0.5))+ scale_fill_grey( start=0.8,end=0.3) + scale_y_continuous(limits=c(-0.45,0.85)) + facet_grid(module ~ genotype)

#print plot
plot1  
ggsave(file="Data/results/GraviallMEline.pdf",plot=plot1,width=4,height=6)

#############
#drought
sample_tsai<-read.table("Data/Tsai_SUT4_RNAseq.txt",sep="\t",header=T)

tsai_MEs<-merge(sample_tsai,MEs,by.x="Run_s",by.y="row.names")
row.names(tsai_MEs)<-tsai_MEs[,1]

#relevel genotype
tsai_MEs$genotype_s<-factor(tsai_MEs$genotype_s, levels = c("wild type","SUT4-RNAi transgenic"))

#relevel genotype
tsai_MEs$tissue_s<-factor(tsai_MEs$tissue_s, levels = c("leaf LPI-15","root","bark","xylem"))

allME<-NULL
for(j in 33:ncol(tsai_MEs))
{
  EG<-tsai_MEs[,c(9,11,12,j)]
  names(EG)[4]<-"EG"
  tmp.avg<-ddply(EG, .(genotype_s,tissue_s,treatment_s), summarise, tmp.avg=mean(EG), tmp.se=sd(EG)/sqrt(length(EG)))
  
  tmp.avg$genotype_s<-factor(tmp.avg$genotype_s, levels = c("wild type","SUT4-RNAi transgenic"))
  
  #relevel genotype
  tmp.avg$tissue_s<-factor(tmp.avg$tissue_s, levels = c("leaf LPI-15","root","bark","xylem"))
  
  
  allME<-rbind(allME,cbind(names(tsai_MEs)[j],tmp.avg))
}

#plot allME
names(allME)[1]<-"module"
allME$module<-factor(allME$module, levels = c("MEyellow","MEblue","MEturquoise","MEgrey","MEbrown"))

plot2<-ggplot(allME,aes(x=tissue_s,y=tmp.avg, group=treatment_s, colour=treatment_s)) + geom_line() + geom_errorbar(aes(ymax=tmp.avg+tmp.se, ymin=tmp.avg-tmp.se), width=0.25,size=0.5)+theme_bw() + ylab("Gene expression (rpkm)") +xlab("wood type") + theme(text = element_text(size=8)) + theme(plot.background = element_blank() ,panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border = element_blank()) + theme(axis.line = element_line(color = 'black',size=0.5))+ scale_fill_grey( start=0.8,end=0.3) + scale_y_continuous(limits=c(-0.45,0.85)) + facet_grid(module ~ genotype_s)

#print plot
plot2
ggsave(file="Data/results/Drought_allMEline.pdf",plot=plot2,width=4.5,height=6)

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
allME$module<-factor(allME$module, levels = c("MEyellow","MEblue","MEturquoise","MEgrey","MEbrown"))

plot3<-ggplot(allME,aes(x=tis,y=tmp.avg,group=1)) + geom_line() + geom_errorbar(aes(ymax=tmp.avg+tmp.se, ymin=tmp.avg-tmp.se), width=0.25,size=0.5)+theme_bw() + ylab("Gene expression (rpkm)") +xlab("P. trichocarpa tissue") + theme(text = element_text(size=8)) + theme(plot.background = element_blank() ,panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border = element_blank()) + theme(axis.line = element_line(color = 'black',size=0.5))+ scale_fill_grey( start=0.8,end=0.3) + scale_y_continuous(limits=c(-0.45,0.85)) + facet_grid(module ~ .)

#print plot
plot3
ggsave(file="Data/results/WoodyTissue_allMEline.pdf",plot=plot3,width=2.5,height=6)

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
allME$module<-factor(allME$module, levels = c("MEyellow","MEblue","MEturquoise","MEgrey","MEbrown"))

plot4<-ggplot(allME,aes(x=coordinate,y=EG,group=meassure,colour=meassure)) + geom_point() +theme_bw() + stat_smooth(method="lm", se=FALSE, colour="black")+ ylab("Gene expression (rpkm)") +xlab("P. trichocarpa Provenance") + theme(text = element_text(size=8)) + theme(plot.background = element_blank() ,panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border = element_blank()) + theme(axis.line = element_line(color = 'black',size=0.5))+ scale_fill_grey( start=0.8,end=0.3) + facet_grid(module ~ meassure,scale="free_x")+ scale_y_continuous(limits=c(-0.45,0.85))

#print plot
plot4
ggsave(file="Data/results/Provenance_allMEline.pdf",plot4,width=4,height=6)







###############################
#
# GO analysis of Populus genes using AT best blast hit
# Potri to AT is obtained from the phytozome annotation Ptrichocarpa_210
#
###############################

source("R/GO_populus.R")
GOtems<-read.table("Data/GO_terms.txt",sep="\t")
pt<-read.table("Data/Ptrichocarpa_210_annotation_primary.txt",sep="\t")
consColors = labels2colors(labels[mods$goodGenes])

out.anno<-data.frame(cbind(row.names(data),consColors))
out.anno2<-merge(out.anno,pt,by.x="V1",by.y="V2")
 
Iblue<-which(out.anno2$consColors=="blue")
GOblue<-atGOanalysis(out.anno2$V10[Iblue])
#View(summary(GOblue$BP))

Ibrown<-which(out.anno2$consColors=="brown")
GObrown<-atGOanalysis(out.anno2$V10[Ibrown])
#View(summary(GObrown$BP))
 
Iyellow<-which(out.anno2$consColors=="yellow")
GOyellow<-atGOanalysis(out.anno2$V10[Iyellow])
 
Iturquoise<-which(out.anno2$consColors=="turquoise")
GOturquoise<-atGOanalysis(out.anno2$V10[Iturquoise])
 
# Igrey<-which(out.anno2$consColors=="grey")
# GOgrey<-atGOanalysis(out.anno2$V10[Igrey])
 
save(GOblue,GObrown,GOyellow,GOturquoise,file="Data/results/Module_GO.rdata")

load("Data/results/Module_GO.rdata")



###############
#
#Plot Figure 5: Functional enrichment of Consensus module using GO analysis
#buid heatmap showing differences in GO enrichment for auxin, cell-wall, hormone and meristem
#
################


key<-"auxin|hormone|gibberelli|brassino|cell wall|meristem|localization|replication|methylation|histone|chromatin"

blueGO<-summary(GOblue$BP)[grep(key,summary(GOblue$BP)[,7],perl=TRUE),1:2]
names(blueGO)[2]<-"blue"

brownGO<-summary(GObrown$BP)[grep(key,summary(GObrown$BP)[,7],perl=TRUE),1:2]
names(brownGO)[2]<-"brown"

yellowGO<-summary(GOyellow$BP)[grep(key,summary(GOyellow$BP)[,7],perl=TRUE),1:2]
names(yellowGO)[2]<-"yellow"

turquoiseGO<-summary(GOturquoise$BP)[grep(key,summary(GOturquoise$BP)[,7],perl=TRUE),1:2]
names(turquoiseGO)[2]<-"turquoise"


GOtable<-merge(blueGO,brownGO,by.x="GOBPID",by.y="GOBPID",all=T)
GOtable<-merge(GOtable,yellowGO,by.x="GOBPID",by.y="GOBPID",all=T)
GOtable<-merge(GOtable,turquoiseGO,by.x="GOBPID",by.y="GOBPID",all=T)
GOtable<-merge(GOtable,GOtems,by.x="GOBPID",by.y="V1",all.x=T)

#generate order
o1<-grep("auxin",GOtable$V3)
o3<-grep("hormone|gibberelli|brassino",GOtable$V3)
o4<-grep("cell wall",GOtable$V3)
o5<-grep("meristem",GOtable$V3)
o2<-grep("localization",GOtable$V3)
o6<-grep("replication",GOtable$V3)
o7<-grep("methylation|histone|chromatin",GOtable$V3)

o<-c(o1,o3,o4,o5,o2,o6,o7)

#convert to -log10(pvalue)
GOtable[,2:5]<--log10(GOtable[,2:5])
GOtable[is.na(GOtable)] <- 0

results<-data.frame(t(GOtable[o,2:5]))
names(results)<-GOtable[o,1]

#reorder modules
results<-results[c("yellow","blue","turquoise","brown"),]

my_palette <- colorRampPalette(c("white", "red"))(n = 8)

pdf(file="Data/results/Module_Go_enrichment.pdf",width=8,height=4)
par(mar = c(6, 8.5, 3, 3));
labeledHeatmap(Matrix = results,
               xLabels = names(results),
               yLabels = paste("ME",row.names(results),sep=""),
               yColorLabels = TRUE,
               yColorWidth = 0.05,
               ySymbols = paste("ME",row.names(results),sep=""),
               colors = my_palette,
               #textMatrix = txt,
               setStdMargins = FALSE,
               cex.text = 0.5,
               cex.lab.y = 0.8,
               cex.lab.x = 0.5,
               xLabelsAngle = 90,
               zlim = c(0,15),
               main = paste("GO Enrichment of Co-expression Modules")
)

dev.off()



###############
#
#Plot Figure 6:
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
    
    tf_binding[count,1]<-i
    tf_binding[count,j]<-p
    count<-count+1
  }
}

row.names(tf_binding)<-paste("ME",tf_binding[,1],sep="")
tf_binding<-tf_binding[c("MEyellow","MEblue","MEturquoise","MEgrey","MEbrown"),]

pdf(file="Data/results/consensus_modules_TFbinding_heatmap.pdf",w=8,h=6)
par(mar=c(5,10,5,1))
labeledHeatmap(Matrix = -log10(tf_binding[,2:ncol(tf_binding)]),
               xLabels = names(tf_binding)[2:ncol(tf_binding)],
               yLabels = paste("ME",tf_binding[,1],sep=""),
               ySymbols = paste("ME",tf_binding[,1],sep=""),
               colorLabels = TRUE,
               colors = colorRampPalette(c("white", "red"))(n = 8),
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
fp_binding<-fp_binding[c("MEyellow","MEblue","MEturquoise","MEgrey","MEbrown"),c("Module","fpuS","fpOS","fpin","fpOE","fpdS","fpOA")]

pdf(file="Data/results/consensus_modules_DNase_enrich_heatmap.pdf",w=8,h=6)
par(mar=c(5,10,5,1))
labeledHeatmap(Matrix = -log10(fp_binding[,2:ncol(fp_binding)]),
               xLabels = names(fp_binding)[2:ncol(fp_binding)],
               yLabels = paste("ME",fp_binding[,1],sep=""),
               ySymbols = paste("ME",fp_binding[,1],sep=""),
               colorLabels = TRUE,
               colors = colorRampPalette(c("white", "red"))(n = 8),
               #textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(0,15),
               main = paste("Module-DNase relationships"))
dev.off()


###############
#
#Plot Figure 7:
#
################

row.names(out.anno)<-out.anno[,1]
para<-read.table("Data/version3_duplications.out",sep=",",header=T)

para[,9]<-out.anno[para[,1],2]
para[,10]<-out.anno[para[,2],2]
para[is.na(para)]<-"grey"

tmp<-chisq.test(para[,9],para[,10])

resid<-as.matrix(tmp$stdres)

pdf(file="Data/results/Paralog_module.pdf",w=8,h=6)
labeledHeatmap(Matrix = resid,
               xLabels =dimnames(resid)[[1]],
               yLabels =dimnames(resid)[[2]],
               ySymbols =dimnames(resid)[[2]],
               colorLabels = FALSE,
               colors= greenWhiteRed(50),
               setStdMargins = FALSE,
               cex.text= 0.5,
               zlim =c(-25,25))
dev.off()