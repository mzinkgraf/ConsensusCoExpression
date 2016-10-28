
GOanalysis<-function (genes,pv=0.01)
{
  
  # List of packages for function
  .packages = c("GOstats", "GSEABase")
  source("http://bioconductor.org/biocLite.R")
  # Install BioC packages if not already installed
  .inst <- .packages %in% installed.packages()
  if(length(.packages[!.inst]) > 0) biocLite(.packages[!.inst], quietly=TRUE)
  
  require("GOstats",quietly = TRUE)
  require("GSEABase",quietly = TRUE)
  
  data<-read.table("Data/GO_Universe.txt", sep="\t",header=TRUE)
  goFrame<-GOFrame(data,organism="Populus trichocarpa")
  goAllFrame<-GOAllFrame(goFrame)
  
  
  gsc <- GeneSetCollection(goAllFrame, setType = GOCollection())
  universe<-unique(data$frame.gene_id)
  
  g_id<-as.factor(intersect(data$frame.gene_id,genes))
  
  params_MF <- GSEAGOHyperGParams(name="Annotation Params MF",
                                  geneSetCollection=gsc,
                                  geneIds =g_id,
                                  universeGeneIds = universe,
                                  ontology = "MF",
                                  pvalueCutoff = pv,
                                  conditional = FALSE,
                                  testDirection = "over")
  
  Over_MF <- hyperGTest(params_MF)
  
  params_BP <- GSEAGOHyperGParams(name="Annotation Params BP",
                                  geneSetCollection=gsc,
                                  geneIds =g_id,
                                  universeGeneIds = universe,
                                  ontology = "BP",
                                  pvalueCutoff = pv,
                                  conditional = FALSE,
                                  testDirection = "over")
  
  Over_BP <- hyperGTest(params_BP)
  params_CC <- GSEAGOHyperGParams(name="Annotation Params CC",
                                  geneSetCollection=gsc,
                                  geneIds =g_id,
                                  universeGeneIds = universe,
                                  ontology = "CC",
                                  pvalueCutoff = pv,
                                  conditional = FALSE,
                                  testDirection = "over")
  
  Over_CC <- hyperGTest(params_CC)
  return(list(Over_MF,Over_BP,Over_CC))
}

atGOanalysis<-function (genes,pv=0.01)
{
  # List of packages for function
  .packages = c("GOstats", "GSEABase")

  # Install BioC packages if not already installed
  .inst <- .packages %in% installed.packages()
  if(length(.packages[!.inst]) > 0) 
    {
    source("http://bioconductor.org/biocLite.R")
    biocLite(.packages[!.inst], quietly=TRUE)
  }
  require("GOstats",quietly = TRUE)
  require("GSEABase",quietly = TRUE)
  
  #remove transcript numbers from AT names
  genes1<-sub("(AT\\d+G\\d+)\\.\\d+","\\1",genes,perl=TRUE)
  
  
  data<-read.table("Data/AT_GO_Universe.txt", sep="\t",header=TRUE)
  goFrame<-GOFrame(data,organism="Arabidopsis thaliana")
  goAllFrame<-GOAllFrame(goFrame)
  
  gsc <- GeneSetCollection(goAllFrame, setType = GOCollection())
  universe<-unique(data$frame.gene_id)
  
  g_id<-as.character(intersect(data$frame.gene_id,genes1))
  
  params_MF <- GSEAGOHyperGParams(name="Annotation Params MF",
                                  geneSetCollection=gsc,
                                  geneIds =g_id,
                                  universeGeneIds = universe,
                                  ontology = "MF",
                                  pvalueCutoff = pv,
                                  conditional = FALSE,
                                  testDirection = "over")
  
  Over_MF <- hyperGTest(params_MF)
  
  params_BP <- GSEAGOHyperGParams(name="Annotation Params BP",
                                  geneSetCollection=gsc,
                                  geneIds =g_id,
                                  universeGeneIds = universe,
                                  ontology = "BP",
                                  pvalueCutoff = pv,
                                  conditional = FALSE,
                                  testDirection = "over")
  
  Over_BP <- hyperGTest(params_BP)
  params_CC <- GSEAGOHyperGParams(name="Annotation Params CC",
                                  geneSetCollection=gsc,
                                  geneIds =g_id,
                                  universeGeneIds = universe,
                                  ontology = "CC",
                                  pvalueCutoff = pv,
                                  conditional = FALSE,
                                  testDirection = "over")
  
  Over_CC <- hyperGTest(params_CC)
  return(list(MF=Over_MF,BP=Over_BP,CC=Over_CC))
  
}

atGOanalysisBP<-function (genes,pv=0.01)
{
  # List of packages for function
  .packages = c("GOstats", "GSEABase")
  
  # Install BioC packages if not already installed
  .inst <- .packages %in% installed.packages()
  if(length(.packages[!.inst]) > 0) 
  {
    source("http://bioconductor.org/biocLite.R")
    biocLite(.packages[!.inst], quietly=TRUE)
  }
  require("GOstats",quietly = TRUE)
  require("GSEABase",quietly = TRUE)
  
  #remove transcript numbers from AT names
  genes1<-sub("(AT\\d+G\\d+)\\.\\d+","\\1",genes,perl=TRUE)
  
  
  data<-read.table("Data/AT_GO_Universe.txt", sep="\t",header=TRUE)
  goFrame<-GOFrame(data,organism="Arabidopsis thaliana")
  goAllFrame<-GOAllFrame(goFrame)
  
  gsc <- GeneSetCollection(goAllFrame, setType = GOCollection())
  universe<-unique(data$frame.gene_id)
  
  g_id<-as.character(intersect(data$frame.gene_id,genes1))
  
  params_BP <- GSEAGOHyperGParams(name="Annotation Params BP",
                                  geneSetCollection=gsc,
                                  geneIds =g_id,
                                  universeGeneIds = universe,
                                  ontology = "BP",
                                  pvalueCutoff = pv,
                                  conditional = FALSE,
                                  testDirection = "over")
  
  Over_BP <- hyperGTest(params_BP)
  return(list(BP=Over_BP))
  
}


atGO2Potri<-function(go_id) 
{
  mapping<-read.table("Data/MAPPINGatGO2PotriV3.txt",sep="\t",header=T)
  if(length(go_id)==0)
  {
    return(print("NO arabidopsis GO term was recieved"))
  } else if (length(go_id)>1) {
    return(print("More than one GO term was recieved"))
  } else if (length(grep("GO\\:\\d+",go_id,perl=TRUE))!=1) {
    return(print("GO term format is incorrect and shoulf be GO:0000000"))
  } else {
    return(unique(mapping$V1[which(mapping$frame.go_id==go_id)]))
  }
  
}

atGO2PotriAnno<-function(go_id) 
{
  mapping<-read.table("Data/MAPPINGatGO2PotriV3.txt",sep="\t",header=T)
  pt<-read.table("Data/Ptrichocarpa_210_annotation_primary.txt",sep="\t")
  if(length(go_id)==0)
  {
    return(print("NO arabidopsis GO term was recieved"))
  } else if (length(grep("GO\\:\\d+",go_id,perl=TRUE))!=length(go_id)) {
    return(print("GO term format is incorrect and shoulf be GO:0000000"))
  } else if (length(go_id)>1) {
    potri<-unique(mapping$V1[which(mapping$frame.go_id%in%go_id)])
    return(pt[which(pt$V2 %in% potri),])
  } else {
    potri<-unique(mapping$V1[which(mapping$frame.go_id==go_id)])
    return(pt[which(pt$V2 %in% potri),])
  }
  
}