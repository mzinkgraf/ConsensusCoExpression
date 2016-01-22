#################
#Consensus CoExpression analysis of Populus RNA-seq 
library(edgeR)
library(WGCNA)          ###used for topological overlap calculation and clustering steps
library(RColorBrewer)   ###used to create nicer colour palettes
library(preprocessCore) ###used by the quantile normalization function
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
