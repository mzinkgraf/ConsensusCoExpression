#################
#Consensus CoExpression analysis of Populus RNA-seq 
library(edgeR)
library(WGCNA) 
library(RColorBrewer)
library(preprocessCore)
#source("~/R/tensionwood/module_GO_enrichment.R")
#source("/biodata/pipeline/Rscripts/GOstat//GO_populus.R")
source("R/networkFunctions-extras-05.R")
source("R/NetworkFunctions-Mouse.R")
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
load("Data/scaleFreeAnalysis-powerTables.RData")
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

###############
#
#Plot Figure 2A: Dendrogram and Consensus module identification
#
################

pdf(file="Data/consensusDendro.pdf",wi=10,5);
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

MEDiss = 1-cor(MEs);
# Cluster module eigengenes
METree = hclust(as.dist(MEDiss), method = "average");

###############
#
#Plot Figure 2B: Consensus and Individual Eigenegene dendrogram
#
################

pdf(file="Data/Experiments_ME_Dendro.pdf",w=10,h=2)
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
pdf(file="Data/Experiment_MEcor.pdf",w=10,h=2)
PlotExpPCsCor(oMEs,dataset_labels,IncludeSign=TRUE,setMargins=TRUE,plotConsensus=TRUE)
dev.off()


###############
#
#Plot Figure 2D: Identify Global and Experiment specific co-expression relationships
#
################

pdf(file="Data/Consensus_module_ExperMod_heatmap.pdf",w=10,h=3)
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
  
  #par(mar=c(10,15,5,1))
  labeledHeatmap(Matrix = t(pTable),
                 yLabels = paste(" ", consModules),
                 xLabels = paste(" ", TWModules),
                 colorLabels = TRUE,
                 #ySymbols = paste("Cons ", consModules, ": ", consModTotals, sep=""),
                 #xSymbols = paste(setLabels[l]," ", TWModules, ": ", TWModTotals, sep=""),
                 #textMatrix = CountTbl,
                 ySymbols = consModules,
                 xSymbols = TWModules,
                 colors = greenWhiteRed(100)[50:100],
                 main = dataset_labels[l],
                 cex.text = 0.8, cex.lab = 0.8, setStdMargins = FALSE);
  
  
}
dev.off()


#relate MEs to traits
row.names(MEs)<-lib_names

#TW 
samples<-read.table("Data/GA_RNAandLibraries.txt",sep="\t",header=T,stringsAsFactors=F)
samples<-samples[-28,]
exper1<-read.table("Data/sample_and_library_names.csv",sep=",", row.names=1,stringsAsFactors=F,header=T)
TW_traits<-data.frame(cbind(c(exper1$Genotype,samples$Genotype),c(exper1$Sample1,samples$Sample1),c(exper1$sample2,samples$Treatment),c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2)))
names(TW_traits)<-c("genotype","woodtype","GA","year")
row.names(TW_traits)<-c(row.names(exper1),samples$Name)

TW_MEs<-merge(TW_traits,MEs,by.x="row.names",by.y="row.names")
row.names(TW_MEs)<-TW_MEs[,1]
TW_MEs<-TW_MEs[,-1]

TW_ME_cor<-data.frame(matrix(vector(), 0, 4, dimnames=list(c(), c("genotype","woodtype","GA","year"))), stringsAsFactors=F)
for (e in 1:4)
{
  count=1
  for(k in 5:15)
  {
    TW_ME_cor[count,1:4]<-anova(lm(TW_MEs[,k]~genotype+woodtype+GA+year,data=TW_MEs))[1:4,5]
    row.names(TW_ME_cor)[count]<-names(TW_MEs)[k]
    count=count+1
  }
}

my_palette <- colorRampPalette(c("white", "red"))(n = 12)

pdf(file="Data/consensus_modules_TWtrait_heatmap.pdf",w=8,h=8)
par(mar=c(5,10,5,1))
labeledHeatmap(Matrix = -log10(TW_ME_cor),
               xLabels = names(TW_ME_cor),
               yLabels = row.names(TW_ME_cor),
               ySymbols = row.names(TW_ME_cor),
               colorLabels = TRUE,
               colors = my_palette,
               #textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(0,50),
               main = paste("Module-trait relationships"))
dev.off()

#########################
#PT association

tis<-c("phloem","phloem","xylem","phloem","xylem","phloem","xylem","phloem","xylem","phloem","xylem","phloem","xylem","phloem","xylem")
tree<-c("t0","t0","t0","t1","t1","t2","t2","t3","t3","t4","t4","t5","t5","t6","t6")
pt_lib<-names(data)[57:71]
PT_traits<-data.frame(cbind(pt_lib,tis,tree))

PT_MEs<-merge(PT_traits,MEs,by.x="pt_lib",by.y="row.names")
row.names(PT_MEs)<-PT_MEs[,1]
PT_MEs<-PT_MEs[,-1]

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

my_palette <- colorRampPalette(c("white", "red"))(n = 12)

pdf(file="Data/consensus_modules_PTtrait_heatmap.pdf",w=8,h=8)
par(mar=c(5,10,5,1))
labeledHeatmap(Matrix = -log10(PT_ME_cor),
               xLabels = names(PT_ME_cor),
               yLabels = row.names(PT_ME_cor),
               ySymbols = row.names(PT_ME_cor),
               colorLabels = TRUE,
               colors = my_palette,
               #textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(0,20),
               main = paste("Module-trait relationships"))
dev.off()


#########
#Provenance Module associations
sample_loc<-read.table("Data/provenance_sample_locations.txt",sep="\t",header=T)

man_MEs<-merge(sample_loc,MEs,by.x="row.names",by.y="row.names")
row.names(man_MEs)<-man_MEs[,1]
man_MEs<-man_MEs[,-c(1,2)]

man_ME_cor<-data.frame(matrix(vector(), 0, 2, dimnames=list(c(), c("longitude","latitude"))), stringsAsFactors=F)
for (e in 1:2)
{
  count=1
  for(k in 5:ncol(man_MEs))
  {
    man_ME_cor[count,1:2]<-anova(lm(man_MEs[,k]~Longitude+Latitude,data=man_MEs))[1:2,5]
    row.names(man_ME_cor)[count]<-names(man_MEs)[k]
    count=count+1
  }
}

my_palette <- colorRampPalette(c("white", "red"))(n = 12)

pdf(file="Data/consensus_modules_ProvenanceTrait_heatmap.pdf",w=8,h=8)
par(mar=c(5,10,5,1))
labeledHeatmap(Matrix = -log10(man_ME_cor),
               xLabels = names(man_ME_cor),
               yLabels = row.names(man_ME_cor),
               ySymbols = row.names(man_ME_cor),
               colorLabels = TRUE,
               colors = my_palette,
               #textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(0,20),
               main = paste("Module-trait relationships"))
dev.off()

#########
#mtsai Module associations
sample_tsai<-read.table("Data/Tsai_SUT4_RNAseq.txt",sep="\t",header=T)

tsai_MEs<-merge(sample_tsai,MEs,by.x="Run_s",by.y="row.names")
row.names(tsai_MEs)<-tsai_MEs[,1]

#relevel genotype
tsai_MEs$genotype_s<-factor(tsai_MEs$genotype_s, levels = c("wild type","SUT4-RNAi transgenic"))

#relevel genotype
tsai_MEs$tissue_s<-factor(tsai_MEs$tissue_s, levels = c("leaf LPI-15","root","bark","xylem"))


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

my_palette <- colorRampPalette(c("white", "red"))(n = 12)

pdf(file="Data/consensus_modules_tsaitrait_heatmap.pdf",w=8,h=8)
par(mar=c(5,10,5,1))
labeledHeatmap(Matrix = -log10(tsai_ME_cor),
               xLabels = names(tsai_ME_cor),
               yLabels = row.names(tsai_ME_cor),
               ySymbols = row.names(tsai_ME_cor),
               colorLabels = TRUE,
               colors = my_palette,
               #textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(0,20),
               main = paste("Module-trait relationships"))
dev.off()

########################################
#chip enrichment of modules
out.anno<-data.frame(cbind(row.names(data),consColors))

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

pdf(file="Data/consensus_modules_binding_heatmap.pdf",w=8,h=8)
par(mar=c(5,10,5,1))
labeledHeatmap(Matrix = -log10(tf_binding[,2:ncol(tf_binding)]),
               xLabels = names(tf_binding)[2:ncol(tf_binding)],
               yLabels = paste("ME",tf_binding[,1],sep=""),
               ySymbols = paste("ME",tf_binding[,1],sep=""),
               colorLabels = TRUE,
               colors = my_palette,
               #textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(0,20),
               main = paste("Module-binding relationships"))
dev.off()


###################
#combine all results

comb<-merge(TW_ME_cor,PT_ME_cor,by.x="row.names",by.y="row.names")
comb<-merge(comb,man_ME_cor,by.x="Row.names",by.y="row.names")
comb<-merge(comb,tsai_ME_cor,by.x="Row.names",by.y="row.names")
comb<-merge(comb,tf_binding[,-1],by.x="Row.names",by.y="row.names")
row.names(comb)<-comb[,1]
comb<-comb[rev(METree$labels[METree$order]),]

#remove grey
comb<-comb[-which(row.names(comb)=="MEgrey"),]

pdf(file="Data/consensus_modules_allTraits_heatmap.pdf",w=10,h=8)
par(mar=c(5,10,5,1))
labeledHeatmap(Matrix = -log10(comb[,-c(1)]),
               xLabels = names(comb)[-c(1)],
               yLabels = comb[,1],
               ySymbols = comb[,1],
               colorLabels = TRUE,
               colors = my_palette,
               #textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(0,20),
               xLabelsAngle = 90,
               main = paste("Module-binding relationships"))
dev.off()

