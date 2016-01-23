#################
#Consensus CoExpression analysis of Populus RNA-seq 
library(edgeR)
library(WGCNA) 
library(RColorBrewer)
library(preprocessCore)
source("~/R/tensionwood/module_GO_enrichment.R")
source("/biodata/pipeline/Rscripts/GOstat//GO_populus.R")
source("../networkFunctions-extras-05.R");
allowWGCNAThreads(n=8);
options(stringsAsFactors = FALSE);

# #gravitropism experiment
# Data1<-read.table("Data/Gravitropism_normALLrpkm.txt",sep="\t",stringsAsFactors=F)
# 
# #xylem/phloem named PT_datExpr0
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
# save(data,TW_mods,PT_mods,man_mods,tsai_mods,file="Grand_analysis_expressions_mods.rdata")

#load expression and modules
load("Grand_analysis_expressions_mods.rdata")

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
save(powerTables,file="Data/scaleFreeAnalysis-powerTables.RData");
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

pdf(file = "Data/scaleFreeTopologyAnalysis.pdf", wi = 10, h=8);
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
    saveTOMs = TRUE,
    saveConsensusTOMs = TRUE,
    getTOMScalingSamples = TRUE,
    verbose = 3, indent = 2);
} ) );

collectGarbage();

#Save the results for future use
#save(mods,TOMinfo, file="Data/Grand_mods.RData");
load("Data/Grand_mods.RData")

mergeCut = 0.4
merge= mergeCloseModules(multiExpr, mods$unmergedColors, cutHeight = mergeCut, getNewUnassdME = TRUE, relabel = TRUE, consensusQuantile = 0.25);
labels=merge$colors;

#Save the resulting labels for future use
#save(merge, labels, file="Data/merge-labels.RData");
load("Data/merge-labels.RData")

pdf(file="Data/consensusDendro.pdf",wi=10,5);
plotDendroAndColors(mods$dendrograms[[1]], labels2colors(labels[mods$goodGenes]),
                    "Consensus",
                    main ="Consensus modules",
                    dendroLabels = FALSE,
                    addGuide = TRUE, hang = 0.01, abHeight = 0.99,
                    guideHang = 0.03, marAll =c(1,5,1,1));
#If plotting into a file, close it 
dev.off();

#plot eigengenes
MEs0 = multiSetMEs(multiExpr, universalColors = labels2colors(labels))

MEs<-orderMEs(rbind(MEs0[[1]]$data,MEs0[[2]]$data,MEs0[[3]]$data,MEs0[[4]]$data))

MEDiss = 1-cor(MEs);
# Cluster module eigengenes
METree = hclust(as.dist(MEDiss), method = "average");
# Plot the result
sizeGrWindow(7, 6)
pdf(file="Data/ME_Dendro.pdf",10,10)

plot(as.dendrogram(METree), horiz=TRUE, main = "Clustering of module eigengenes", xlab = "", sub = "",cex = 0.8)
dev.off()

#relate consensus modules to original modules

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
  
  
  #=====================================================================================
  #
  #  Code chunk 5
  #
  #=====================================================================================
  
  
  # Truncate p values smaller than 10^{-50} to 10^{-50} 
  pTable[is.infinite(pTable)] = 1.3*max(pTable[is.finite(pTable)]);
  pTable[pTable>100 ] = 100 ;
  # Marginal counts (really module sizes)
  TWModTotals = apply(CountTbl, 1, sum)
  consModTotals = apply(CountTbl, 2, sum)
  # Actual plotting
  sizeGrWindow(10,7 );
  par(mfrow=c(1,1));
  par(cex = 1.0);
  #par(mar=c(20, 15, 2.7, 1)+0.3);
  # Use function labeledHeatmap to produce the color-coded table with all the trimmings
  pdf(file=paste("Data/",names(multiColor)[l],"_connsensus_module_comp_heatmap.pdf",sep=""),w=7,h=10)
  par(mar=c(10,15,5,1))
  labeledHeatmap(Matrix = t(pTable),
                 yLabels = paste(" ", consModules),
                 xLabels = paste(" ", TWModules),
                 colorLabels = TRUE,
                 #ySymbols = paste("Cons ", consModules, ": ", consModTotals, sep=""),
                 #xSymbols = paste(setLabels[l]," ", TWModules, ": ", TWModTotals, sep=""),
                 #textMatrix = CountTbl,
                 ySymbols = paste(consModules, ": ", consModTotals, sep=""),
                 xSymbols = TWModules,
                 colors = greenWhiteRed(100)[50:100],
                 main = paste("Correspondence of ",setLabels[l]," set-specific and consensus modules",sep=""),
                 cex.text = 1.0, cex.lab = 1.0, setStdMargins = FALSE);
  dev.off()
  
}

#left off on line 281