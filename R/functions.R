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

#---------------------------------------------------------------------------------------------------------
# HeatmapWithTextLabels.R
#---------------------------------------------------------------------------------------------------------

#--------------------------------------------------------------------------
#
# ReverseRows = function(Matrix)
#
#--------------------------------------------------------------------------
#


ReverseRows = function(Matrix)
{
  ind = seq(from=dim(Matrix)[1], to=1, by=-1);
  Matrix[ind,];
  #Matrix
}

ReverseVector = function(Vector)
{
  ind = seq(from=length(Vector), to=1, by=-1);
  Vector[ind];
  #Vector
}

#--------------------------------------------------------------------------
#
# HeatmapWithTextLabels = function ( Matrix, xLabels, yLabels, ... ) { 
#
#--------------------------------------------------------------------------
# This function plots a heatmap of the specified matrix 
# and labels the x and y axes wit the given labels.
# It is assumed that the number of entries in xLabels and yLabels is consistent 
# with the dimensions in.
# If ColorLabels==TRUE, the labels are not printed and instead interpreted as colors --
#  -- a simple symbol with the appropriate color is printed instead of the label.
# The x,yLabels are expected to have the form "..color" as in "MEgrey" or "PCturquoise".
# xSymbol, ySymbols are additional markers that can be placed next to color labels

HeatmapWithTextLabels = function ( Matrix, xLabels, yLabels = NULL, xSymbols = NULL, ySymbols = NULL, 
                                   InvertColors=FALSE, ColorLabels = NULL, xColorLabels = FALSE,
                                   yColorLabels = FALSE,
                                   SetMargins = TRUE,
                                   colors = NULL, NumMatrix = NULL, cex.Num = NULL, cex.lab = NULL, 
                                   plotLegend = TRUE, ... ) 
{
  if (!is.null(ColorLabels)) {xColorLabels = ColorLabels; yColorLabels = ColorLabels; }
  
  if (is.null(yLabels) & (!is.null(xLabels)) & (dim(Matrix)[1]==dim(Matrix)[2])) 
    yLabels = xLabels; 
  if (SetMargins)
  {
    if (ColorLabels)
    {
      par(mar=c(2,2,3,5)+0.2);
    } else {
      par(mar = c(7,7,3,5)+0.2);
    }
  }
  if (is.null(colors)) 
    #if (IncludeSign)
    #{
    #colors = GreenRedWhite(50);
    #} else {
    colors = heat.colors(30);
  #}
  if (InvertColors) colors = ReverseVector(colors);
  if (plotLegend)
  {
    image.plot(t(ReverseRows(Matrix)), xaxt = "n", xlab="", yaxt="n", ylab="", col=colors, ...);
  } else {
    image(z = t(ReverseRows(Matrix)), xaxt = "n", xlab="", yaxt="n", ylab="", col=colors, ...);
  }
  axis(1, labels = FALSE, tick = FALSE)
  axis(2, labels = FALSE, tick = FALSE)
  axis(3, labels = FALSE, tick = FALSE)
  axis(4, labels = FALSE, tick = FALSE)
  nxlabels = length(xLabels)
  #   plot x axis labels using:
  #   par("usr")[3] - 0.25 as the vertical placement
  #   srt = 45 as text rotation angle
  #   adj = 1 to place right end of text at tick mark
  #   pd = TRUE to allow for text outside the plot region
  plotbox = par("usr");
  xmin = plotbox[1]; xmax = plotbox[2]; ymin = plotbox[3]; yrange = plotbox[4]-ymin;
  ymax = plotbox[4]; xrange = xmax - xmin;
  
  # print(paste("plotbox:", plotbox[1], plotbox[2], plotbox[3], plotbox[4]));
  
  nylabels = length(yLabels)
  axis(2, labels = FALSE, tick = FALSE)
  xspacing = 1/(nxlabels-1); yspacing = 1/(nylabels-1)
  # print(paste("nxlabels:", nxlabels));
  if (!xColorLabels)
  {
    if (is.null(cex.lab)) cex.lab = 1;
    text(((1:nxlabels)-1)*xspacing , ymin - 0.02, srt = 45, 
         adj = 1, labels = xLabels, xpd = TRUE, cex = cex.lab)
  } else {
    rect(((1:nxlabels)-1)*xspacing - xspacing/2, ymin-xspacing*1.2,
         ((1:nxlabels)-1)*xspacing + xspacing/2, ymin-xspacing*0.2,
         density = -1,  col = substring(xLabels, 3), border = substring(xLabels, 3), xpd = TRUE)
    if (!is.null(xSymbols))
      text ( ((1:nxlabels)-1)*xspacing, ymin-xspacing*1.3, xSymbols, adj = c(0.5, 0), xpd = TRUE);
  }
  if (!yColorLabels)
  {
    if (is.null(cex.lab)) cex.lab = 1;
    text(xmin - 0.01*xrange, ((1:nylabels)-1)/(nylabels-1), srt = 0, 
         adj = 1, labels = ReverseVector(yLabels), xpd = TRUE, cex = cex.lab )
  } else {
    rect(xmin-yspacing*1.2, ((nylabels:1)-1)*yspacing - yspacing/2,
         xmin-yspacing*0.2, ((nylabels:1)-1)*yspacing + yspacing/2, 
         density = -1,  col = substring(yLabels, 3), border = substring(yLabels, 3), xpd = TRUE)
    if (!is.null(ySymbols))
      text (xmin-yspacing*1.2, ((nylabels:1)-1)*yspacing, ySymbols, adj = c(1, 0.5), xpd = TRUE);
  }
  
  if (!is.null(NumMatrix))
  {
    if (is.null(cex.Num)) cex.Num = par("cex");
    #if (dim(NumMatrix)!=dim(Matrix))
    #  stop("HeatmapWithTextLabels: NumMatrix was given, but has dimensions incompatible with Matrix.");
    for (rw in 1:dim(Matrix)[1])
      for (cl in 1:dim(Matrix)[2])
      {
        text((cl-1)*xspacing, (dim(Matrix)[1]-rw)*yspacing, 
             as.character(NumMatrix[rw,cl]), xpd = TRUE, cex = cex.Num, adj = c(0.5, 0.5));
      }
  }
  axis(1, labels = FALSE, tick = FALSE)
  axis(2, labels = FALSE, tick = FALSE)
  axis(3, labels = FALSE, tick = FALSE)
  axis(4, labels = FALSE, tick = FALSE)
}

RedWhiteGreen = function(n)
{
  half = as.integer(n/2);
  green = c(seq(from=0, to=1, length.out = half), rep(1, times = half+1));
  red = c(rep(1, times = half+1), seq(from=1, to=0, length.out = half));
  blue = c(seq(from=0, to=1, length.out = half), 1, seq(from=1, to=0, length.out = half));
  col = rgb(red, green, blue, maxColorValue = 1);
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

#####
#predict values from R cookbook

predictvals <- function(model, xvar, yvar, xrange=NULL, samples=100, ...) {
  # If xrange isn't passed in, determine xrange from the models.
  # Different ways of extracting the x range, depending on model type 
  if (is.null(xrange)) {
  if (any(class(model) %in% c("lm", "glm"))) xrange <- range(model$model[[xvar]])
  else if (any(class(model) %in% "loess")) xrange <- range(model$x)
}
newdata <- data.frame(x = seq(xrange[1], xrange[2], length.out = samples))
names(newdata) <- xvar
newdata[[yvar]] <- predict(model, newdata = newdata, ...)
newdata
}