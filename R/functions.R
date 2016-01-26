#Consensus Co-expression paper functions

#get values not in vector
`%ni%`=Negate(`%in%`)

PeakTypeByDiff<-function(gene_strand,SS, ES, SE, EE)
{
  
  #omitting this NA check speeds things up
  #also you return -1 below if there the arguments are not satisfied
  #   if ( is.na(c(SS, ES, SE, EE)) ) {
  #     print("DATA MISSING:")
  #     print(peakS)
  #     print(peakE)
  #     print(geneS)
  #     print(geneE)
  #     return(-2)
  #   }
  gene_strand=as.integer(gene_strand)
  pSgS = SS > 0
  pSgE = SE > 0
  pEgS = ES > 0
  pEgE = EE > 0
  
  if ( gene_strand==1 & !pSgS & !pEgS & !pSgE & !pEgE) {
    
    return("upstream")
    
  } else if( gene_strand==1 & !pSgS & pEgS & !pSgE & !pEgE) {
    
    return("OverlapStart") 
    
  } else if ( gene_strand==1 & pSgS & pEgS & !pSgE & !pEgE) {
    
    return("inside") 
    
  } else if ( gene_strand==1 & pSgS & pEgS & !pSgE & pEgE) {
    
    return("OverlapEnd") 
    
  } else if ( gene_strand==1 & pSgS & pEgS & pSgE & pEgE) {
    
    return("downstream") 
    
  } else if ( gene_strand==1 & !pSgS & pEgS & !pSgE & pEgE) {
    
    return("OverlapAll") 
    
  } else if( gene_strand==-1 & pSgS & pEgS & pSgE & pEgE) {
    
    return("upstream") # negative strand
    
  } else if( gene_strand==-1 & pSgS & pEgS & !pSgE & pEgE) {
    
    return("OverlapStart") # negative strand
    
  } else if( gene_strand==-1 & !pSgS & pEgS & !pSgE & !pEgE) {
    
    return("OverlapEnd") # negative strand
    
  } else if( gene_strand==-1 & !pSgS & !pEgS & !pSgE & !pEgE) {
    
    return("downstream") # negative strand
    
  } else if( gene_strand==-1 & pSgS & pEgS & !pSgE & !pEgE) {
    
    return("inside") # negative strand
    
  } else if( gene_strand==-1 & !pSgS & pEgS & !pSgE & pEgE) {
    
    return("OverlapAll") # negative strand
    
  } else { 
    
    return("NA") 
    
  }
}

peaks2genes<-function(peaks,genes=NA,n=3,max_distance=NA) 
{
  if(is.na(genes)) genes<-read.table("Data/pt210_gene_features.txt",sep="\t",header=T,stringsAsFactors=F)
  
  out<-list()
  
  weight<-c(1,0.8,0.6,0.4,0.2)
  #   g2 = genes
  #   #flip the gene start and end based on strand
  #   for ( j in 1:nrow(g2)) {
  #     gS = g2[j,4]
  #     gE = g2[j,5]
  #     if ( g2[j,3] < 0) {
  #       g2[j,5] = gS
  #       g2[j,4] = gE
  #     }
  #   }
  
  for (i in 1:nrow(peaks)) 
  {
    #limit to only genes on the peak's chromosome
    #g<-g2[which(g2[,2] == peaks[i,1]),]
    g<-genes[which(genes[,2] == peaks[i,1]),]
    
    if ( nrow(g) > 0 ) {
      
      #calculate distance from the peak to each feature on chrom
      d<-as.data.frame(matrix(data = NA, nrow = nrow(g), ncol = 9))
      d[,1]<-g[,1]
      d[,2]<-peaks[i,2]-g[,4] #SS
      d[,3]<-peaks[i,3]-g[,4] #ES
      d[,4]<-peaks[i,2]-g[,5] #SE
      d[,5]<-peaks[i,3]-g[,5] #EE
      d[,6]<-apply(d[,2:5], 1, function(x) min(abs(x)) )
      d[,7]<-g[,3]
      d[,8]<-apply(d[,c(7,2,3,4,5)], 1, function(x) PeakTypeByDiff(x[1], x[2], x[3], x[4], x[5]))
      #set distance for Overlapping features to zero
      d[which(d[,8]=="OverlapStart" | d[,8]=="OverlapEnd" | d[,8]=="inside" | d[,8]=="OverlapAll"),6]<-0
      #calculate weight of binding
      d[which(d[,6]<=1000),9]<-weight[1]
      d[which(d[,6]>1000 & d[,6]<=2000),9]<-weight[2]
      d[which(d[,6]>2000 & d[,6]<=4000),9]<-weight[3]
      d[which(d[,6]>4000 & d[,6]<=6000),9]<-weight[4]
      d[which(d[,6]>6000),9]<-weight[5]
      d[which(is.na(d[,6])),9]<-NA
      d<-d[order(d[,6]),]
      
      #output either the number of genes or max_distance to gene
      if(!is.na(n) & is.na(max_distance)) {
        out[i]<-list(cbind(peaks[i,4],d[1:n,]))
      } else if (!is.na(max_distance) & is.na(n) & nrow(d[which(d[,6]<=max_distance),])>0) {
        out[i]<-list(cbind(peaks[i,4],d[which(d[,6]<=max_distance),]))
      } else if (!is.na(n) & !is.na(max_distance) ) {
        q<-d[which(d[,6]<=max_distance),]
        if(nrow(q>0)) {
          N=min(nrow(q),n)
          out[i]<-list(cbind(peaks[i,4],q[1:N,]))
        }
      } 
    }
    print(i) 
  }
  return(do.call(rbind.data.frame,out)) 
}


#--------------------------------------------------------------------------------------
#
# PlotExpPCsCor
#
#--------------------------------------------------------------------------------------
# Modified PlotCorPCs 

PlotExpPCsCor = function(PCs, Titles, ColorLabels = FALSE, colors = NULL, IncludeSign = FALSE,  IncludeGrey = TRUE, setMargins = FALSE, PlotDiagAdj = FALSE, plotConsensus=TRUE,...)
{
  
  #print(PlotDiagAdj);
  
  if (is.null(colors)) 
    if (IncludeSign)
    {
      colors = RedWhiteGreen(50);
    } else {
      colors = heat.colors(30);
    }
  NSets = length(PCs);
  #cex = par("cex");
  cex=0.5
  mar = par("mar");
  
  #par(cex = cex);
  if (!IncludeGrey)
  {
    for (set in 1:NSets)
      PCs[[set]]$data = PCs[[set]]$data[ , substring(names(PCs[[set]]$data),3)!="grey"]
  }
  
  if(plotConsensus)
  {
    par(mfrow = c(1, (NSets+1)));
    
    Cons<-list()
    Cons[[1]]<-PCs[[1]]
    #Cons[[1]]$data<-PCs[[1]]$data
    for (set in 2:NSets)
    {
      Cons[[1]]$data<-rbind(Cons[[1]]$data,PCs[[set]]$data)
    }
    for(sets in 1:NSets+1)
    {
      Cons[[sets]]<-PCs[[sets-1]]
    }
    Titles=c("Consensus",Titles)
    No.Sets<-NSets+1
    PCs<-Cons
  } else { 
    No.Sets<-NSets
    par(mfrow = c(1, No.Sets));
  }
  
  
  
  for (i.col in (1:No.Sets))
  {
    par(cex = cex);
    if (setMargins) {
      if (ColorLabels) {
        par(mar = c(1,2,3,4)+0.2);
      } else {
        par(mar = c(6,7,3,5)+0.2);
      }
    } else {
      par(mar = mar);
    }
    No.Modules = dim(PCs[[i.col]]$data)[2]
    
    corPC = cor(PCs[[i.col]]$data, use="p") 
    if (IncludeSign)
    {
      if (PlotDiagAdj) {
        HeatmapWithTextLabels((1+corPC)/2, names(PCs[[i.col]]$data), names(PCs[[i.col]]$data),
                              main=Titles[[i.col]], InvertColors=TRUE, zlim=c(0,1.0),
                              ColorLabels = ColorLabels, colors = colors, SetMargins = FALSE, ...);
      } else {
        HeatmapWithTextLabels(corPC, names(PCs[[i.col]]$data), names(PCs[[i.col]]$data),
                              main=Titles[[i.col]], InvertColors=TRUE, zlim=c(-1,1.0),
                              ColorLabels = ColorLabels, colors = colors, SetMargins = FALSE, ...);
      }
    } else {
      HeatmapWithTextLabels(abs(corPC), names(PCs[[i.col]]$data), names(PCs[[i.col]]$data),
                            main=Titles[[i.col]], InvertColors=TRUE, zlim=c(0,1.0),
                            ColorLabels = ColorLabels, colors = colors, SetMargins = FALSE, ...);
    }
    
  }
  
}
