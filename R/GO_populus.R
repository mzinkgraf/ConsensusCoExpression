
GOanalysis<-function (genes,pv=0.01)
{
  # List of packages for function
  .packages = c("GOstats", "GSEABase")
  
  # Install CRAN packages (if not already installed)
  .inst <- .packages %in% installed.packages()
  if(length(.packages[!.inst]) > 0) install.packages(.packages[!.inst], quietly=TRUE)
  
  require("GOstats",quietly = TRUE)
  require("GSEABase",quietly = TRUE)
  
    data<-read.table("~/R/GOstat/pt210_GO_Universe.txt", sep="\t",header=TRUE)
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
  
  # Install CRAN packages (if not already installed)
  .inst <- .packages %in% installed.packages()
  if(length(.packages[!.inst]) > 0) install.packages(.packages[!.inst], quietly=TRUE)
  
  require("GOstats",quietly = TRUE)
  require("GSEABase",quietly = TRUE)
  
    data<-read.table("~/R/GOstat/PTat_GO_Universe.txt", sep="\t",header=TRUE)
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