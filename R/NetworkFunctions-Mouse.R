# This document contains functions that we find useful for gene co-expression network analysis
# We recommend to read the tutorial first.
# Steve Horvath, Bin Zhang, Jun Dong, Andy Yip
# To cite this code or the statistical methods please use
# Zhang B, Horvath S (2005) A General Framework for Weighted Gene Co-Expression Network Analysis. 
# Statistical Applications in Genetics and Molecular Biology. In Press.
# Technical Report and software code at: www.genetics.ucla.edu/labs/horvath/CoexpressionNetwork.

# Modifications by Peter langfelder: in function ScaleFreePlot1: changed min(kk) and max(kk) to min(kk, rm.na=T) and
# max(kk, na.rm=T).  
# in function ModulePrinComps1 made the printed statement a bit more informative and use print.flush.
# Changed titling in ScaleFreePlot1

# Added lots of new functions for "higher level" analysis of several datasets at once

# CONTENTS  
# This document contains function for carrying out the following tasks
# A) Assessing scale free topology and choosing the parameters of the adjacency function
#    using the scale free topology criterion (Zhang and Horvath 05)
# B) Computing the topological overlap matrix 
# C) Defining gene modules using clustering procedures
# D) Summing up modules by their first principal component (first eigengene)
# E) Relating a measure of gene significance to the modules 
# F) Carrying out a within module analysis (computing intramodular connectivity) 
#    and relating intramodular connectivity to gene significance.
# G) Miscellaneous other functions, e.g. for computing the cluster coefficient.
# H) Network functions from Bin Zhang for dynamic tree cutting of a hierarchical clustering tree (dendrogram)
# I) General statistical functions, e.g. for scatterplots


#--------------------------------------------------------------------------------------------
# Set global parameters.

# This parameter controls whether my implementation of pmax should use a call to an external function
# (recommended if the requisite library is available). Otherwise an R-only implementation will be used,
# which is significantly slower (but available and stable on all R platforms).

UseCpmax = FALSE;

#--------------------------------------------------------------------------------------------

# Load the requisite libraries

WorkingDirectory = getwd();

#if (!library(MASS, logical.return=TRUE)) { # standard, no need to install
#For some reason, MASS does not seem to be installed on Titan, so we'll try to load it 
# from my own library. If that fails as well, stop.
#  if (!library(MASS, logical.return=TRUE, lib.loc="M:/Work/RLibrary")) stop()
#}   

library(MASS);
library(class)  # standard, no need to install
library(cluster)	
#library(sma)	# install it for the function plot.mat 
library(impute)# install it for imputing missing value
library(Hmisc)	# install it for the C-index calculations
library(survival)
library(fields);

#oldwd = getwd();

if (UseCpmax) source("../ComputerDefinition/ComputerDefinition.R");

if (exists("memory.limit"))
{
  # increase the available memory 
  memory.limit(size=4000)   
}

#setwd(oldwd);


#####################################################################################################
# A) Assessing scale free topology and choosing the parameters of the adjacency function
#    using the scale free topology criterion (Zhang and Horvath 05)
########################################################################################################


# ===================================================
#For hard thresholding, we use the signum (step) function
if(exists("signum1") ) rm(signum1); 
signum1=function(corhelp,tau1)  {
  adjmat1= as.matrix(abs(corhelp)>=tau1)
  dimnames(adjmat1) <- dimnames(corhelp)
  diag(adjmat1) <- 0
  adjmat1}

# ===================================================
# For soft thresholding, one can use the sigmoid function 
# But we have focused on the power adjacency function in the tutorial...
if (exists("sigmoid1") ) rm(sigmoid1); sigmoid1=function(ss,mu1=0.8,alpha1=20) {
  1/(1+exp(-alpha1*(ss-mu1)))}





#This function is useful for speeding up the connectivity calculation.
#The idea is to partition the adjacency matrix into consecutive baches of a #given size.
#In principle, the larger the batch size the faster is the calculation. But #smaller batchsizes require #less memory...
# Input: gene expression data set where *rows* correspond to microarray samples #and columns correspond to genes. 
# If fewer than MinimumNoSamples contain gene expression information for a given
# gene, then its connectivity is set to missing. 
if(exists("SoftConnectivity")) rm(SoftConnectivity);
SoftConnectivity=function(datE, power=6, batchsize=1500, MinimumNoSamples=10) {
  no.genes=dim(datE)[[2]]
  no.samples=dim(datE)[[1]]
  if (no.genes<no.samples | no.genes<10 | no.samples<5 ) {stop("Error: Something seems to be wrong. Make sure that the input data frame has genes as rows and array samples as columns. Alternatively, there could be fewer than 10 genes or fewer than 5 samples. ") } else {
    sum1=function(x) sum(x,na.rm=T)
    k=rep(NA,no.genes)
    no.batches=as.integer(no.genes/ batchsize)
    if (no.batches>0) {
      for (i in 1:no.batches) {
        print(paste("batch number = ", i))
        index1=c(1:batchsize)+(i-1)* batchsize
        ad1=abs(cor(datE[,index1], datE,use="p"))^power
        ad1[is.na(ad1)]=0
        k[index1]=apply(ad1,1,sum1)
        # If fewer than MinimumNoSamples contain gene expression information for a given
        # gene, then we set its connectivity to 0.
        NoSamplesAvailable=apply(!is.na(datE[,index1]),2,sum)
        k[index1][NoSamplesAvailable< MinimumNoSamples]=NA
      } # end of for (i in 1:no.batches
    } # end of if (no.batches>0)...
    if (no.genes-no.batches*batchsize>0 ) {
      restindex=c((no.batches*batchsize+1):no.genes)
      ad1=abs(cor(datE[,restindex], datE,use="p"))^power
      ad1[is.na(ad1)]=0
      k[restindex]=apply(ad1,1,sum1)
      NoSamplesAvailable=apply(!is.na(datE[,restindex]),2,sum)
      k[restindex][NoSamplesAvailable< MinimumNoSamples]=NA
    } # end of if
  } # end of else statement
  k
} # end of function




# ===================================================
# The function PickHardThreshold can help one to estimate the cut-off value 
# when using the signum (step) function.
# The first column lists the threshold ("cut"),
# the second column lists the corresponding p-value based on the Fisher Transform 
# of the correlation. 
# The third column reports the resulting scale free topology fitting index R^2.
# The fourth column reports the slope of the fitting line, it shoud be negative for 
# biologically meaningul networks.
# The fifth column reports the fitting index for the truncated exponential model. 
# Usually we ignore it.
# The remaining columns list the mean, median and maximum resulting connectivity.
# To pick a hard threshold (cut) with the scale free topology criterion:
# aim for high scale free R^2 (column 3), high connectivity (col 6) and negative slope 
# (around -1, col 4).
# The output is a list with 2 components. The first component lists a sugggested cut-off
# while the second component contains the whole table.
# The removeFirst option removes the first point (k=0, P(k=0)) from the regression fit.
# no.breaks specifies how many intervals used to estimate the frequency p(k) i.e. the no. of points in the 
# scale free topology plot.
if (exists("PickHardThreshold")) rm(PickHardThreshold);
PickHardThreshold=function(datExpr1,RsquaredCut=0.85, cutvector=seq(0.1,0.9,by=0.05) ,removeFirst=FALSE,no.breaks=10) {
  no.genes   <- dim(datExpr1)[[2]]
  no.genes <- dim(datExpr1)[[2]]
  no.samples= dim(datExpr1)[[1]]
  colname1=c("Cut","p-value", "scale law R^2", "slope="  ,"truncated R^2","mean(k)","median(k)","max(k)")
  datout=data.frame(matrix(666,nrow=length(cutvector),ncol=length(colname1) ))
  names(datout)=colname1
  datout[,1]=cutvector
  for (i in c(1:length(cutvector) ) ){
    cut1=cutvector[i]
    datout[i,2]=2*(1-pt(sqrt(no.samples-1)*cut1/sqrt(1-cut1^2),no.samples-1))}
  if(exists("fun1")) rm(fun1)
  fun1=function(x) {
    corx=abs(cor(x,datExpr1,use="p"))
    out1=rep(NA, length(cutvector) )
    for (j in c(1:length(cutvector))) {out1[j]=sum(corx>cutvector[j])}
    out1
  } # end of fun1
  datk=t(apply(datExpr1,2,fun1))
  for (i in c(1:length(cutvector) ) ){
    nolinkshelp <- datk[,i]-1
    cut2=cut(nolinkshelp,no.breaks)
    binned.k=tapply(nolinkshelp,cut2,mean)
    freq1=as.vector(tapply(nolinkshelp,cut2,length)/length(nolinkshelp))
    # The following code corrects for missing values etc
    breaks1=seq(from=min(nolinkshelp),to=max(nolinkshelp),length=no.breaks+1)
    hist1=hist(nolinkshelp,breaks=breaks1,equidist=F,plot=FALSE,right=TRUE)
    binned.k2=hist1$mids
    binned.k=ifelse(is.na(binned.k),binned.k2,binned.k)
    binned.k=ifelse(binned.k==0,binned.k2,binned.k)
    freq1=ifelse(is.na(freq1),0,freq1)
    xx= as.vector(log10(binned.k))
    if(removeFirst) {freq1=freq1[-1]; xx=xx[-1]}
    plot(xx,log10(freq1+.000000001),xlab="log10(k)",ylab="log10(p(k))" )
    lm1= lm(as.numeric(log10(freq1+.000000001))~ xx )
    lm2=lm(as.numeric(log10(freq1+.000000001))~ xx+I(10^xx) )
    datout[i,3]=summary(lm1)$adj.r.squared 
    datout[i,4]=summary(lm1)$coefficients[2,1]  
    datout[i,5]=summary(lm2)$adj.r.squared
    datout[i,6]=mean(nolinkshelp)
    datout[i,7]= median(nolinkshelp)
    datout[i,8]= max(nolinkshelp) 
  } 
  datout=signif(datout,3) 
  print(data.frame(datout));
  # the cut-off is chosen as smallest cut with R^2>RsquaredCut 
  ind1=datout[,3]>RsquaredCut
  indcut=NA
  indcut=ifelse(sum(ind1)>0,min(c(1:length(ind1))[ind1]),indcut)
  # this estimates the cut-off value that should be used. 
  # Don't trust it. You need to consider slope and mean connectivity as well!
  cut.estimate=cutvector[indcut][[1]]
  list(cut.estimate, data.frame(datout));
} # end of function











# ===========================================================
# The function PickSoftThreshold allows one to estimate the power parameter when using
# a soft thresholding approach with the use of the power function AF(s)=s^Power
# The function PickSoftThreshold allows one to estimate the power parameter when using
# a soft thresholding approach with the use of the power function AF(s)=s^Power
# The removeFirst option removes the first point (k=1, P(k=1)) from the regression fit.
if (exists("PickSoftThreshold")) rm(PickSoftThreshold);
PickSoftThreshold=function(datExpr1,RsquaredCut=0.85, powervector=c(seq(1,10,by=1),seq(12,20,by=2)),
                           removeFirst=FALSE,no.breaks=10) {
  no.genes <- dim(datExpr1)[[2]]
  no.samples= dim(datExpr1)[[1]]
  colname1=c("Power", "scale law R^2" ,"slope", "truncated R^2","mean(k)","median(k)","max(k)")
  datout=data.frame(matrix(666,nrow=length(powervector),ncol=length(colname1) ))
  names(datout)=colname1
  datout[,1]=powervector
  if(exists("fun1")) rm(fun1)
  fun1=function(x) {
    corx=abs(cor(x,datExpr1,use="p"))
    out1=rep(NA, length(powervector) )
    for (j in c(1:length(powervector))) {out1[j]=sum(corx^powervector[j])}
    out1
  } # end of fun1
  datk=t(apply(datExpr1,2,fun1))
  for (i in c(1:length(powervector) ) ){
    nolinkshelp <- datk[,i]-1
    cut2=cut(nolinkshelp,no.breaks)
    binned.k=tapply(nolinkshelp,cut2,mean)
    freq1=as.vector(tapply(nolinkshelp,cut2,length)/length(nolinkshelp))
    # The following code corrects for missing values etc
    breaks1=seq(from=min(nolinkshelp),to=max(nolinkshelp),length=no.breaks+1)
    hist1=hist(nolinkshelp,breaks=breaks1,equidist=F,plot=FALSE,right=TRUE)
    binned.k2=hist1$mids
    binned.k=ifelse(is.na(binned.k),binned.k2,binned.k)
    binned.k=ifelse(binned.k==0,binned.k2,binned.k)
    freq1=ifelse(is.na(freq1),0,freq1)
    
    xx= as.vector(log10(binned.k))
    if(removeFirst) {freq1=freq1[-1]; xx=xx[-1]}
    plot(xx,log10(freq1+.000000001),xlab="log10(k)",ylab="log10(p(k))" )
    lm1= lm(as.numeric(log10(freq1+.000000001))~ xx )
    lm2=lm(as.numeric(log10(freq1+.000000001))~ xx+I(10^xx) )
    datout[i,2]=summary(lm1)$adj.r.squared 
    datout[i,3]=summary(lm1)$coefficients[2,1]  
    datout[i,4]=summary(lm2)$adj.r.squared
    datout[i,5]=mean(nolinkshelp)
    datout[i,6]= median(nolinkshelp)
    datout[i,7]= max(nolinkshelp) 
  } 
  datout=signif(datout,3) 
  print(data.frame(datout));
  # the cut-off is chosen as smallest cut with R^2>RsquaredCut 
  ind1=datout[,2]>RsquaredCut
  indcut=NA
  indcut=ifelse(sum(ind1)>0,min(c(1:length(ind1))[ind1]),indcut)
  # this estimates the power value that should be used. 
  # Don't trust it. You need to consider slope and mean connectivity as well!
  power.estimate=powervector[indcut][[1]]
  list(power.estimate, data.frame(datout));
}






# ===================================================
# The function ScaleFreePlot1 creates a plot for checking scale free topology
# when truncated1=T is specificed, it provides the R^2 measures for the following
# degree distributions: a) scale free topology, b) log-log R^2 and c) truncated exponential R^2

# The function ScaleFreePlot1 creates a plot for checking scale free topology
if(exists("ScaleFreePlot1")) rm(ScaleFreePlot1) ; 

ScaleFreePlot1=function(kk,no.breaks=10,AF1="" ,truncated1=FALSE, removeFirst=FALSE,cex.lab1=1){
  
  #bin data into no.breaks bins: first create the factor cut1 and code the values
  cut1=cut(kk,no.breaks)
  #now calculate the mean of each bin
  binned.k=tapply(kk,cut1,mean)
  freq1=tapply(kk,cut1,length)/length(kk)
  # The following code corrects for missing values etc
  breaks1=seq(from=min(kk, na.rm=T),to=max(kk, na.rm=T),length=no.breaks+1)
  hist1=hist(kk,breaks=breaks1,equidist=F,plot=FALSE,right=TRUE)
  binned.k2=hist1$mids
  binned.k=ifelse(is.na(binned.k),binned.k2,binned.k)
  binned.k=ifelse(binned.k==0,binned.k2,binned.k)
  freq1=ifelse(is.na(freq1),0,freq1)
  plot(log10(binned.k),log10(freq1+.000000001),xlab=paste(AF1,"log10(k)"),ylab="log10(p(k))",cex.lab=cex.lab1 )
  xx= as.vector(log10(binned.k))
  if(removeFirst) {freq1=freq1[-1]; xx=xx[-1]}
  lm1=lm(as.numeric(log10(freq1+.000000001))~ xx )
  lines(xx,predict(lm1),col=1)
  OUTPUT=data.frame(ScaleFreeRsquared=round(summary(lm1)$adj.r.squared,2),Slope=round(lm1$coefficients[[2]],2))
  if (truncated1==TRUE) { 
    lm2=lm(as.numeric(log10(freq1+.000000001))~ xx+I(10^xx) );
    OUTPUT=data.frame(ScaleFreeRsquared=round(summary(lm1)$adj.r.squared,2),Slope=round(lm1$coefficients[[2]],2),
                      TruncatedRsquared=round(summary(lm2)$adj.r.squared,2))
    print("the red line corresponds to the truncated exponential fit")
    lines(xx,predict(lm2),col=2);
    title(paste( 
      "scale free R^2=",as.character(round(summary(lm1)$adj.r.squared,2)),
      ", slope=", round(lm1$coefficients[[2]],2),
      ", trunc.R^2=",as.character(round(summary(lm2)$adj.r.squared,2))))} else { 
        title(paste("R^2=",as.character(round(summary(lm1)$adj.r.squared,2)) , 
                    " sl=", round(lm1$coefficients[[2]],2)), cex=0.4)
      }
  OUTPUT
} # end of function 









#################################################################################################################
################################################################################################################################
# B) Computing the topological overlap matrix 
#################################################################################################################
#################################################################################################################



# ===================================================
#The function TOMdist1 computes a dissimilarity 
# based on the topological overlap matrix (Ravasz et al)
# Input: an Adjacency matrix with entries in [0,1]
if(exists("TOMdist1")) rm(TOMdist1);

TOMdist1=function(adjmat1, maxADJ=FALSE) {
  diag(adjmat1)=0;
  adjmat1[is.na(adjmat1)]=0;
  maxh1=max(as.dist(adjmat1) ); minh1=min(as.dist(adjmat1) ); 
  if (maxh1>1 | minh1 < 0 ) {print(paste("ERROR: the adjacency matrix contains entries that are larger than 1 or smaller than 0!!!, max=",maxh1,", min=",minh1)) } else { 
    if (  max(c(as.dist(abs(adjmat1-t(adjmat1)))))>0   ) {print("ERROR: non-symmetric adjacency matrix!!!") } else { 
      kk=apply(adjmat1,2,sum)
      maxADJconst=1
      if (maxADJ==TRUE) maxADJconst=max(c(as.dist(adjmat1 ))) 
      Dhelp1=matrix(kk,ncol=length(kk),nrow=length(kk))
      denomTOM= pmin(as.dist(Dhelp1),as.dist(t(Dhelp1)))   +as.dist(maxADJconst-adjmat1); 
      gc();gc();
      numTOM=as.dist(adjmat1 %*% adjmat1 +adjmat1);
      #TOMmatrix=numTOM/denomTOM
      # this turns the TOM matrix into a dissimilarity 
      out1=1-as.matrix(numTOM/denomTOM) 
      diag(out1)=1
      out1
    }}
}

#---------------------------------------------------------------------------
# This is a somewhat modified TOMdist1.

SignedTOMdist = function(adjmat1, maxADJ=FALSE)
{
  diag(adjmat1)=0;
  adjmat1[is.na(adjmat1)]=0;
  collect_garbage();
  kk=apply(adjmat1,2,sum)
  collect_garbage();
  maxADJconst=1
  if (maxADJ==TRUE) maxADJconst=max(c(as.dist(adjmat1 ))) 
  collect_garbage();
  Dhelp1 = matrix(kk,ncol=length(kk),nrow=length(kk))
  collect_garbage();
  denomTOM = pmin(as.dist(Dhelp1),as.dist(t(Dhelp1))) + as.dist(maxADJconst-adjmat1); 
  rm(Dhelp1);
  collect_garbage();
  gc(); gc();
  numTOM=as.dist(adjmat1 %*% adjmat1 +adjmat1);
  collect_garbage();
  #TOMmatrix=numTOM/denomTOM
  # this turns the TOM matrix into a dissimilarity 
  out1=1-as.matrix(numTOM/denomTOM) 
  rm(numTOM); rm(denomTOM);
  collect_garbage();
  diag(out1)=1
  out1
}



# ===================================================
# This function computes a TOMk dissimilarity
# which generalizes the topological overlap matrix (Ravasz et al)
# Input: an Adjacency matrix with entries in [0,1]
# WARNING:  ONLY FOR UNWEIGHTED NETWORKS, i.e. the adjacency matrix contains binary entries...
# This function is explained in Yip and Horvath (2005)
# http://www.genetics.ucla.edu/labs/horvath/TOM/
if(exists("TOMkdist1")) rm(TOMkdist1);
TOMkdist1 = function(adjmat1,k=1){
  maxh1=max(as.dist(adjmat1) ); minh1=min(as.dist(adjmat1) );
  if (k!=round(abs(k))) {
    stop("k must be a positive integer!!!", call.=TRUE);}
  if (maxh1>1 | minh1 < 0 ){
    print(paste("ERROR: entries of the adjacency matrix must be between inclusively 0 and 1!!!, max=",maxh1,", min=",minh1))}
  else {
    if (  max(c(as.dist(abs(adjmat1-t(adjmat1)))))>0   ) {print("ERROR: non-symmetric adjacency matrix!!!") } else { 
      
      B <- adjmat1;
      if (k>=2) {
        for (i in 2:k) {
          diag(B) <- diag(B) + 1;
          B = B %*% adjmat1;}}   # this gives the number of paths with length at most k connecting a pair
      B <- (B>0);   # this gives the k-step reachability from a node to another
      diag(B) <- 0;   # exclude each node being its own neighbor
      B <- B %*% B   # this gives the number of common k-step-neighbor that a pair of nodes share
      
      Nk <- diag(B);
      B <- B +adjmat1;   # numerator
      diag(B) <- 1;
      denomTOM=outer(Nk,Nk,FUN="pmin")+1-adjmat1;
      diag(denomTOM) <- 1;
      1 - B/denomTOM   # this turns the TOM matrix into a dissimilarity
    }}
}


# IGNORE THIS function...
# The function TOMdistROW computes the TOM distance of a gene (node)
# with that of all other genes in the network.
# WhichRow is an integer that specifies which row of the adjacency matrix
# corresponds to the gene of interest.
# Output=vector of TOM distances.
if (exists("TOMdistROW") ) rm(TOMdistROW) 
TOMdistROW=function(WhichRow=1, adjmat1, maxADJ=FALSE) {
  diag(adjmat1)=0;
  maxh1=max(as.dist(adjmat1) ); minh1=min(as.dist(adjmat1) ); 
  if (maxh1>1 | minh1 < 0 ) {print(paste("ERROR: the adjacency matrix contains entries that are larger than 1 or smaller than 0!!!, max=",maxh1,", min=",minh1)) } else { 
    kk=apply(adjmat1,2,sum)
    numTOM=adjmat1[WhichRow,] %*% adjmat1 +adjmat1[WhichRow,]; 
    numTOM[WhichRow]=1
    maxADJconst=1
    if (maxADJ==TRUE) maxADJconst=max(c(as.dist(adjmat1 ))) 
    denomTOM=pmin(kk[WhichRow],kk)+maxADJconst-adjmat1[WhichRow,]; denomTOM[WhichRow]=1
    #TOMmatrix=numTOM/denomTOM
    # this turns the TOM matrix into a dissimilarity 
    1-numTOM/denomTOM 
  }
}


#####################################################################################################
################################################################################################################################
# C) Defining gene modules using clustering procedures
#####################################################################################################
################################################################################################################################

# ===================================================
#The function modulecolor2 function assigns colors to the observations 
# in the branches of a dendrogram
# we use it to define modules....
if (exists("modulecolor2")) rm(modulecolor2);
modulecolor2=function(hier1, h1=0.9,minsize1=50) {
  # here we define modules by using a height cut-off for the branches
  labelpred= cutree(hier1,h=h1)
  sort1=-sort(-table(labelpred))
  modulename= as.numeric(names(sort1))
  modulebranch= sort1>minsize1
  no.modules=sum(modulebranch)
  # now we assume that there are fewer than a certain number of colors
  #colorcode=c("turquoise","blue","brown","yellow","green","red","black","purple","orange","pink",
  #"greenyellow","lightcyan","salmon","midnightblue","lightyellow")
  colorcode=c("turquoise","blue","brown","yellow","green","red","black","pink","magenta",
              "purple","greenyellow","tan","salmon", "midnightblue", "lightcyan","grey60",
              "lightgreen", "lightyellow", "royalblue", "darkred", "darkgreen", "darkturquoise",
              "darkgrey", "orange", "darkorange", "white" )
  
  # "grey" means not in any module;
  colorhelp=rep("grey",length(labelpred))
  if ( no.modules==0 | no.modules >length(colorcode)){ print(paste("The number of modules is problematic. Number of modules = ", as.character(no.modules)))} else { for (i in c(1:no.modules)) {colorhelp=ifelse(labelpred==modulename[i],colorcode[i],colorhelp)};
                                                                                                                                                                    colorhelp=factor(colorhelp,levels=c(colorcode[1:no.modules],"grey"))
  }
  factor(colorhelp, levels=unique(colorhelp[hier1$order] ))
}


#---------------------------------------------------------------------------------
#
# ModuleNumber
#
#---------------------------------------------------------------------------------
# Similar to modulecolor2 above, but returns numbers instead of colors, which is oftentimes more useful.
# 0 means unassigned.
# Return value is a simple vector, not a factor.
# Caution: the module numbers are neither sorted nor sequential, the only guarranteed fact is that grey
# probes are labeled by 0 and all probes belonging to the same module have the same number.
# If size-sorted sequential labels are required, "normalize" the result by calling NormalizeLabels
# (below).

ModuleNumber = function(HierTree, CutHeight = 0.9, MinSize = 50)
{
  Branches = cutree(HierTree, h = CutHeight);
  NOnBranches = table(Branches);
  #NOnBranches[i][[1]] somehow gives the number of elements on branch i.
  TrueBranch = NOnBranches >= MinSize;
  Branches[!TrueBranch[Branches]] = 0;
  
  #NewLabels = levels(factor(Branches));
  #for (lab in 1:length(NewLabels)) if (NewLabels[lab]!=0)
  #  Branches[Branches==NewLabels[lab]] = lab;
  
  Branches;
  
}


# The function hclustplot1 creates a barplot where the colors of the bars are sorted according to 
# a hierarchical clustering tree (hclust object)
#if (exists("hclustplot1")) rm(hclustplot1);
#hclustplot1=function(hier1,couleur,title1="Colors sorted by hierarchical clustering") 
#{
#if (length(hier1$order) != length(couleur) ) {print("ERROR: length of color vector not compatible with no. of objects in the hierarchical tree")};
#if (length(hier1$order) == length(couleur) ) {
#barplot(height=rep(1, length(couleur)), col= as.character(couleur[hier1$order]),border=F, main=title1,space=0, axes=F)}
#}

if (exists("hclustplot1")) rm(hclustplot1);
hclustplot1=function(hier1,Color1, Color2=NULL,title1="Colors sorted by hierarchical clustering") 
{
  options(stringsAsFactors=FALSE);
  if (length(hier1$order) != length(Color1) ) 
  { 
    stop("ERROR: length of color vector not compatible with no. of objects in the hierarchical tree");
  } else {
    if (is.null(Color2))
    {
      barplot(height=rep(1, length(Color1)), col= as.character(Color1[hier1$order]),
              border=F, main=title1,space=0, axes=F)
    } else if (length(Color1)==length(Color2)) {
      # height = matrix(0.5, nrow = 2, ncol = length(Color1));
      C1 = Color1[hier1$order]; C2 = Color2[hier1$order]
      step = 1/length(Color1);
      barplot(height=1, col = "white", border=F, main=title1,space=0, axes=F)
      for (i in 1:(length(Color1)))
      {
        lines(x=rep((i*step), times=2), y=c(0,0.5),  col = as.character(C1[i]));
        lines(x=rep((i*step), times=2), y=c(0.5,1),  col = as.character(C2[i]));
      } 
    }
  }
}

if (exists("hclustplotn")) rm(hclustplotn);
hclustplotn=function(hier1, Color, RowLabels=NULL, cex.RowLabels = 0.9, ...) 
{
  options(stringsAsFactors=FALSE);
  if (length(hier1$order) != dim(Color)[[1]] ) 
  { 
    stop("ERROR: length of color vector not compatible with no. of objects in the hierarchical tree");
  } else {
    No.Sets = dim(Color)[[2]];
    C = Color[hier1$order, ]; 
    step = 1/dim(Color)[[1]];
    ystep = 1/No.Sets;
    barplot(height=1, col = "white", border=F,space=0, axes=F, ...)
    for (j in 1:No.Sets)
    {
      ind = (1:(dim(C)[1]));
      xl = (ind-1) * step; xr = ind * step; 
      yb = rep(ystep*(j-1), dim(C)[1]); yt = rep(ystep*j, dim(C)[1]);
      rect(xl, yb, xr, yt, col = as.character(C[,j]), border = as.character(C[,j]));
      if (is.null(RowLabels))
      {
        text(as.character(j), pos=2, x=0, y=ystep*(j-0.5), cex=cex.RowLabels, xpd = TRUE);
      } else {
        text(RowLabels[j], pos=2, x=0, y=ystep*(j-0.5), cex=cex.RowLabels, xpd = TRUE);
      }
    }
    for (j in 1:No.Sets) lines(x=c(0,1), y=c(ystep*j,ystep*j));
  }
}


# ===================================================
# The function TOMplotn creates a TOM plot
# Inputs:  distance measure, hierarchical (hclust) object, color label=couleur
if (exists("TOMplotn")) rm(TOMplotn);
TOMplot1=function(disttom,hier1, couleur,terrainColors=FALSE) {
  no.nodes=length(couleur)
  if (no.nodes != dim(disttom)[[1]] ) {print("ERROR: number of color labels does not equal number of nodes in disttom")} else {
    labeltree=as.character(couleur)
    labelrow  = labeltree
    labelrow[hier1$order[length(labeltree):1]]=labelrow[hier1$order]
    options(expressions = 10000)
    if (terrainColors) heatmap(as.matrix(disttom),Rowv=as.dendrogram(hier1),Colv= as.dendrogram(hier1), scale="none",revC=T, ColSideColors=as.character(labeltree),RowSideColors=as.character(labelrow),
                               labRow=F, labCol=F, col = terrain.colors(1000)) else heatmap(as.matrix(disttom),Rowv=as.dendrogram(hier1),Colv= as.dendrogram(hier1), scale="none",revC=T, ColSideColors=as.character(labeltree),RowSideColors=as.character(labelrow),
                                                                                            labRow=F, labCol=F)
  }
} #end of function


# ===================================================
# The function TOMplot2 creates a TOM plot where the top and left color bars can be different
# Inputs:  distance measure, hierarchical (hclust) object, color label=couleurTop, couleurLeft
if (exists("TOMplot2")) rm(TOMplot2);
TOMplot2=function(disttom,hier1, couleurTop, couleurLeft) {
  no.nodes=length(couleurTop)
  if (no.nodes != length(couleurLeft)) {stop("ERROR: number of top color labels does not equal number of left color labels")}
  if (no.nodes != dim(disttom)[[1]] ) {stop("ERROR: number of color labels does not equal number of nodes in disttom")} else {
    labeltree = as.character(couleurTop)
    labelrow  = as.character(couleurLeft)
    labelrow[hier1$order[length(labeltree):1]]=labelrow[hier1$order]
    options(expressions = 10000)
    heatmap(as.matrix(disttom),Rowv=as.dendrogram(hier1),Colv= as.dendrogram(hier1), scale="none", revC=T, ColSideColors=as.character(labeltree),RowSideColors=as.character(labelrow),
            labRow=F, labCol=F)
  }
} #end of function



# IGNORE THIS FUNCTION...
# The function "BestHeightCut" allows one to find the best height cut-off
# for a hierarchical clustering tree when external gene information is available
# It computes a Kruskal Wallis-test p-value for each height cut-off
# based on determining whether gene significance differs across branch membership.
if(exists("BestHeightCut")) rm(BestHeightCut);
BestHeightCut=function(hier1, GeneSignif, hcut=seq(0.1,.95,by=0.01) ) {
  pvalues=rep(NA, length(hcut))
  for (i in c(1:length(hcut))) {
    colorhelp=modulecolor2(hier1,hcut[i])
    if (length(unique(colorhelp))>1 ) {pvalues[i]=kruskal.test(GeneSignif, colorhelp)$p.value}
    data.frame(hcut,pvalues)
  }}




#####################################################################################################
################################################################################################################################
# D) Summing up modules using their first principal components (first eigengene)
#####################################################################################################
################################################################################################################################

# ===================================================
#The function ModulePrinComps1 finds the first principal component (eigengene) in each 
# module defined by the colors of the input vector "couleur" (Pardon my French).
# It also reports the variances explained by the first 5 principal components.
# And it yields a measure of module conformity for each gene,
# which is highly correlated to the within module connectivity.
# The theoretical underpinnings are described in Horvath, Dong, Yip (2005)
# http://www.genetics.ucla.edu/labs/horvath/ModuleConformity/
# This requires the R library impute
# The output is a list with 3 components: 
# 1) a data frame of module eigengenes (MEs), 
# 2) a data frame that lists the percent variance explained by the first 5 MEs of a module
# 3) a data frame that lists the module conformity for each gene. 
# The be used as alternative connectivity measure....
if(exists("ModulePrinComps1")) rm(ModulePrinComps1);
ModulePrinComps1=function(datexpr, couleur, verbose = 1, print.level = 0, Impute = FALSE,
                          GetConformity = TRUE) {
  if (is.null(datexpr))
  {  
    print("ModulePrinComps1: Error: datexpr is NULL. ");
    stop();
  }
  if (is.null(couleur))
  {  
    print("ModulePrinComps1: Error: couleur is NULL. ");
    stop()
  }
  MaxVectors = 5;
  #print(paste("datexpr dimensions:", as.character(dim(datexpr))));
  spaces = PrintSpaces(print.level);
  modlevels=levels(factor(couleur))
  PrinComps=data.frame(matrix(666,nrow=dim(datexpr)[[1]],ncol= length(modlevels))) 
  varexplained= data.frame(matrix(666,nrow= 5,ncol= length(modlevels)))
  names(PrinComps)=paste("ME",modlevels,sep="")
  for(i in c(1:length(modlevels)) )
  {
    if (verbose>0) 
      print.flush(paste(spaces, "ModulePrinComps1 : Working on ME for module ", modlevels[i], sep = ""));
    modulename    = modlevels[i]
    restrict1= as.character(couleur)== modulename
    #print(paste("length of couleur:", length(couleur), "; length of restric1:", length(restrict1)));
    # in the following, rows are genes and columns are samples     
    datModule=t(datexpr[, restrict1])
    if (Impute)
    {
      saved.seed = .Random.seed;
      datModule=impute.knn(as.matrix(datModule))
      datModule=t(scale(t(datModule)));
      .Random.seed = saved.seed;
    }
    n = dim(datModule)[1]; p = dim(datModule)[2];
    svd1=svd(datModule, nu = min(n, p, MaxVectors), nv = min(n, p, MaxVectors));
    mtitle=paste("MEs of ", modulename," module", sep="");
    varexplained[,i]= (svd1$d[1:5])^2/sum(svd1$d^2)
    # this is the first principal component
    pc1=svd1$v[,1]
    # signh1=sign(sum(cor(pc1,  t(datModule))))
    # if (signh1 != 0)  pc1=signh1* pc1
    PrinComps[,i]= pc1
  }
  ModuleConformity= rep(666,length=dim(datexpr)[[2]])
  if (GetConformity)
  {
    for(i in 1:(dim(datexpr)[[2]])) 
      ModuleConformity[i] = abs(cor(datexpr[,i], PrinComps[,match(couleur[i], modlevels)], 
                                    use="pairwise.complete.obs"))
  } else
  {
    ModuleConformity = NULL;
  }
  
  list(PrinComps=PrinComps, varexplained=varexplained, ModuleConformity=ModuleConformity)
}



#####################################################################################################
################################################################################################################################
# E) Relating a measure of gene significance to the modules 
#####################################################################################################
################################################################################################################################

# ===================================================
# The function ModuleEnrichment1 creates a bar plot that shows whether modules are enriched with
# significant genes.
# More specifically, it reports the mean gene significance for each module.
# The gene significance can be a binary variable or a quantitative variable. 
# It also plots the 95% confidence interval of the mean (CI=mean +/- 1.96* standard error).
# It also reports a Kruskal Wallis P-value.
if( exists("ModuleEnrichment1") ) rm(ModuleEnrichment1);
ModuleEnrichment1=function(genesignif1,couleur,title1="gene significance across modules",labely="Gene Significance",boxplot=F) {
  if (length(genesignif1) != length(couleur) ) print("Error: vectors don\'t have the same lengths") else {
    if (boxplot != TRUE) {
      mean1=function(x) mean(x,na.rm=T) 
      means1=as.vector(tapply(genesignif1,couleur,mean1));
      se1= as.vector(tapply(genesignif1,couleur,stderr1))
      #par(mfrow=c(1,1))
      barplot(means1,
              names.arg=names(table(couleur) ),col= names(table(couleur) )
              ,ylab=labely)
      err.bp(as.vector(means1), as.vector(1.96*se1), two.side=T)} else {
        boxplot(split(genesignif1,couleur),notch=T,varwidth=T, col= names(table(couleur) ),ylab=labely)}
    
    title(paste(title1,", p-value=", signif(kruskal.test(genesignif1,factor(couleur))$p.value,2)))
  }
} # end of function


# IGNORE THIS...
# ===================================================
#The function fisherPvector allows one to compute Fisher exact p-values
# Thus it allows one to carry out an EASE analysis
# Output: a table of Fisher��s exact p-values
# Input: annotation1 is a vector of gene annotations
# Input: couleur (French for color) denotes the module color of each gene
# Only those gene functions (levels of annotation1) that occur a certain mininum number of times
# (parameter= minNumberAnnotation) in the data set will be considered.  
if (exists("fisherPvector" ) ) rm(fisherPvector);
fisherPvector=function(couleur,annotation1,minNumberAnnotation=50) {
  levelsannotation1=levels(annotation1)
  levelscouleur=levels(factor(couleur))
  no.couleur=length(levelscouleur)
  restannotation1=table(annotation1)>minNumberAnnotation
  no.annotation=sum( restannotation1)
  datoutP=data.frame(matrix(666,nrow=no.annotation,ncol=no.couleur) )
  #datoutProp=data.frame(matrix(666,nrow=no.annotation,ncol=2*no.couleur) )
  #names(datoutProp)=paste("Prop",paste( rep(levelscouleur ,rep(2, length(levelscouleur))) ) , c("Y","N")  ,sep=".")
  datoutProp=data.frame(matrix(666,nrow=no.annotation,ncol=no.couleur) )
  names(datoutProp)=paste("Perc",levelscouleur , sep=".")
  names(datoutP)=paste("p",levelscouleur,sep=".")
  restlevelsannotation1= levelsannotation1[restannotation1]
  row.names(datoutP)= restlevelsannotation1
  for (i in c(1:no.annotation) ) {
    for (j in c(1:no.couleur) ){
      tab1=table( annotation1 !=restlevelsannotation1[i], couleur !=levelscouleur[j])
      datoutP[i,j]=signif(fisher.test(tab1)$p.value,2) 
      #datoutProp[i,2*j-1]=signif(tab1[1,1]/sum(tab1[,1] ),2)
      #datoutProp[i,2*j]= signif(tab1[1,2]/sum(tab1[,2]) ,2)
    } 
    table2=table(annotation1 !=restlevelsannotation1[i], couleur)
    datoutProp[i,]= signif(table2[1,]/apply(table2,2,sum),2)
  }
  data.frame(datoutP,datoutProp)
} # end of function fisherPvector



#####################################################################################################
################################################################################################################################
# F) Carrying out a within module analysis (computing intramodular connectivity etc) 
#####################################################################################################
################################################################################################################################

# ===================================================
#The function DegreeInOut computes for each gene 
#a) the total number of connections, 
#b) the number of connections with genes within its module, 
#c) the number of connections with genes outside its module
# When scaledToOne=TRUE, the within module connectivities are scaled to 1, i.e. the max(K.Within)=1 for each module
if (exists("DegreeInOut")) rm(DegreeInOut); DegreeInOut =function(adj1, couleur,scaledToOne=FALSE) {
  no.nodes=length(couleur)
  couleurlevels=levels(factor(couleur))
  no.levels=length(couleurlevels)
  kWithin=rep(-666,no.nodes )
  diag(adj1)=0
  for (i in c(1:no.levels) ) {
    rest1=couleur==couleurlevels[i];
    if (sum(rest1) <3 ) { kWithin[rest1]=0 } else {
      kWithin[rest1]=apply(adj1[rest1,rest1],2,sum)
      if (scaledToOne) kWithin[rest1]=kWithin[rest1]/max(kWithin[rest1])}
  }
  kTotal= apply(adj1,2,sum) 
  kOut=kTotal-kWithin
  if (scaledToOne) kOut=NA
  kDiff=kWithin-kOut
  data.frame(kTotal,kWithin,kOut,kDiff)
}


# =======================================================================
# The function WithinModuleCindex1 relates the node measures (e.g. connectivities) 
# to "external" node significance information within each  module,
# i.e. it  carries out a by module analysis.
# Output: first column reports the spearman correlation p-value between the network variable and the 
# node significance. The next columns contain the Spearman correlations between the variables.
if (exists("WithinModuleAnalysis1")) rm(WithinModuleAnalysis1);
WithinModuleAnalysis1=function(datnetwork,nodesignif, couleur) 
{
  cortesthelp=function( x ) {
    len1=dim(x)[[2]]-1
    out1=rep(666, len1);
    for (i in c(1:len1) ) {out1[i]= signif( cor.test(x[,i+1], x[,1], method="s",use="p" )$p.value ,2) }
    data.frame( variable=names(x)[-1] , NS.CorPval=out1, NS.cor=t(signif(cor (x[,1], x[,-1],use="p",method="s"),2)), 
                signif(cor(x[,-1],use="p",method="s"),2) )
  } #end of function cortesthelp
  print("IGNORE  the warnings...");
  by( data.frame(nodesignif, datnetwork), couleur, cortesthelp);
} #end of function WithinModuleAnalysis


# =======================================================================
# The function WithinModuleCindex1 relates the node measures (e.g. connectivities) 
# to "external" node significance information within each  module, 
# i.e. it  carries out a by module analysis.
# BUT it focuses on the C-index also known as area under the ROC curve
# This measure is related to Kendall's Tau statistic and Somer's D, 
# see F. Harrel (Regression Modeling Strategies). Springer. 
# It requires the following library
library(Hmisc)
# Output: the first column reports the C-index and the second, p-value 
if (exists("WithinModuleCindex1")) rm(WithinModuleCindex1);
WithinModuleCindex1=function(datnetwork,nodesignif, couleur) {
  CindexFun=function( x ) {
    len1=dim(x)[[2]]-1
    outC=rep(666, len1);
    outP=rep(666, len1);
    for (i in c(1:len1) ) {rcor1=rcorr.cens(x[,i+1], x[,1])
                           outC[i]=rcor1[[1]] 
                           outP[i]=1- pnorm(abs(rcor1[[2]]/rcor1[[3]]))
    }
    data.frame( variable=names(x)[-1] , C.index=outC, p.value=outP)
  } #end of function CindexFun
  #print("IGNORE  the warnings...");
  by( data.frame(nodesignif, datnetwork),couleur,CindexFun);
} #end of function WithinModuleAnalysis


# The following function allows on to plot a gene (node) significance measure versus
# connectivity.
if(exists("plotConnectivityGeneSignif1") ) rm( plotConnectivityGeneSignif1);
plotConnectivityGeneSignif1=function(degree1,genesignif1,color1="black", 
                                     title1="Gene Significance vs Connectivity" , xlab1="Connectivity", ylab1="GeneSignificance") {
  lm1=lm(genesignif1~degree1 ,na.action="na.omit")
  plot(degree1, genesignif1, col=color1,ylab=ylab1,xlab=xlab1,main=paste(title1, ", cor=",  
                                                                         signif(cor( genesignif1,degree1, method="s",use="p" )   ,2) ))
  abline(lm1)
}




#####################################################################################################
################################################################################################################################
# G) Miscellaneous other functions, e.g. for computing the cluster coefficient.
#####################################################################################################
################################################################################################################################



# ===================================================
# The function ClusterCoef.fun computes the cluster coefficients.
# Input is an adjacency matrix 
if(exists("ClusterCoef.fun")) rm(ClusterCoef.fun) ; ClusterCoef.fun=function(adjmat1) {
  diag(adjmat1)=0
  no.nodes=dim(adjmat1)[[1]]
  computeLinksInNeighbors <- function(x, imatrix){x %*% imatrix %*% x}
  nolinksNeighbors <- c(rep(-666,no.nodes))
  total.edge <- c(rep(-666,no.nodes))
  maxh1=max(as.dist(adjmat1) ); minh1=min(as.dist(adjmat1) ); 
  if (maxh1>1 | minh1 < 0 ) {print(paste("ERROR: the adjacency matrix contains entries that are larger than 1 or smaller than 0!!!, max=",maxh1,", min=",minh1)) } else { 
    nolinksNeighbors <- apply(adjmat1, 1, computeLinksInNeighbors, imatrix=adjmat1)
    plainsum  <- apply(adjmat1, 1, sum)
    squaresum <- apply(adjmat1^2, 1, sum)
    total.edge = plainsum^2 - squaresum
    CChelp=rep(-666, no.nodes)
    CChelp=ifelse(total.edge==0,0, nolinksNeighbors/total.edge)
    CChelp}
} # end of function



# ===================================================
# The function err.bp  is used to create error bars in a barplot
# usage: err.bp(as.vector(means), as.vector(stderrs), two.side=F)

err.bp<-function(daten,error,two.side=F){
  if(!is.numeric(daten)) {
    stop("All arguments must be numeric")}
  if(is.vector(daten)){ 
    xval<-(cumsum(c(0.7,rep(1.2,length(daten)-1)))) 
  }else{
    if (is.matrix(daten)){
      xval<-cumsum(array(c(1,rep(0,dim(daten)[1]-1)),
                         dim=c(1,length(daten))))+0:(length(daten)-1)+.5
    }else{
      stop("First argument must either be a vector or a matrix") }
  }
  MW<-0.25*(max(xval)/length(xval)) 
  ERR1<-daten+error 
  ERR2<-daten-error
  for(i in 1:length(daten)){
    segments(xval[i],daten[i],xval[i],ERR1[i])
    segments(xval[i]-MW,ERR1[i],xval[i]+MW,ERR1[i])
    if(two.side){
      segments(xval[i],daten[i],xval[i],ERR2[i])
      segments(xval[i]-MW,ERR2[i],xval[i]+MW,ERR2[i])
    } 
  } 
} 

# ===================================================
# this function computes the standard error
if (exists("stderr1")) rm(stderr1)
stderr1 <- function(x){ sqrt( var(x,na.rm=T)/sum(!is.na(x))   ) }



# ===================================================
# The following two functions are for displaying the pair-wise correlation in a panel when using the command "pairs()"
# Typically, we use "pairs(DATA, upper.panel=panel.smooth, lower.panel=panel.cor, diag.panel=panel.hist)" to
# put the correlation coefficients on the lower panel.
panel.cor <- function(x, y, digits=2, prefix="", cex.cor){
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- abs(cor(x, y))
  txt <- format(c(r, 0.123456789), digits=digits)[1]
  txt <- paste(prefix, txt, sep="")
  if(missing(cex.cor)) cex <- 0.8/strwidth(txt)
  text(0.5, 0.5, txt, cex = cex * r)
}
panel.hist <- function(x, ...){
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(usr[1:2], 0, 1.5) )
  h <- hist(x, plot = FALSE)
  breaks <- h$breaks; nB <- length(breaks)
  y <- h$counts; y <- y/max(y)
  rect(breaks[-nB], 0, breaks[-1], y, col="cyan", ...)
}



# ===================================================
# this function computes the standard error
if (exists("stderr1")) rm(stderr1);
stderr1 <- function(x){ sqrt( var(x,na.rm=T)/sum(!is.na(x))   ) }




# ===================================================
# This function collects garbage
if (exists("collect_garbage")) rm(collect_garbage);
collect_garbage=function(){while (gc()[2,4] != gc()[2,4] | gc()[1,4] != gc()[1,4]){}}
collect_garbage()


# this function is used for computing the Rand index below...
# ===================================================
if (exists("choosenew") ) rm(choosenew)
choosenew <- function(n,k){
  n <- c(n)
  out1 <- rep(0,length(n))
  for (i in c(1:length(n)) ){
    if (n[i]<k) {out1[i] <- 0}
    else {out1[i] <- choose(n[i], k)}}
  out1	
}


# ===================================================
# the following function computes the Rand index between 2 clusterings
# assesses how similar two clusterings are
if (exists("Rand1") ) rm(Rand1)
Rand2 <- function(tab,adjust=T) {
  a <- 0; b <- 0; c <- 0; d <- 0; nn <- 0
  m <- nrow(tab);
  n <- ncol(tab);
  for (i in 1:m) {
    c<-0
    for(j in 1:n) {
      a <- a+choosenew(tab[i,j],2)
      nj <- sum(tab[,j])
      c <- c+choosenew(nj,2)
    }
    ni <- sum(tab[i,])
    b <- b+choosenew(ni,2)
    nn <- nn+ni
  }
  if(adjust==T) {
    d <- choosenew(nn,2)
    adrand <- (a-(b*c)/d)/(0.5*(b+c)-(b*c)/d)
    adrand
  } else {
    b <- b-a
    c <- c-a
    d <- choosenew(nn,2)-a-b-c
    rand <- (a+d)/(a+b+c+d)
    rand
  }
}

# ===================================================
# This function is used in "pairs()" function. The problem of the original  panel.cor is that 
# when the correlation coefficient is very small, the lower panel will have a large font 
# instead of a mini-font in a saved .ps file. This new function uses a format for corr=0.2 
# when corr<0.2, but it still reports the original value of corr, with a minimum format.

panel.cor1=function(x, y, digits=2, prefix="", cex.cor){
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- abs(cor(x, y))
  txt <- format(c(r, 0.123456789), digits=digits)[1]
  txt <- paste(prefix, txt, sep="")
  txt1=txt
  r1=r
  if (r<0.2) {
    r1=0.2
    txt1 <- format(c(r1, 0.123456789), digits=digits)[1]
    txt1 <- paste(prefix, txt1, sep="")
  }
  if(missing(cex.cor)) cex <- 0.8/strwidth(txt1)
  cex = cex * r1
  r <- round(r, digits)
  txt <- format(c(r, 0.123456789), digits=digits)[1]
  txt <- paste(prefix, txt, sep="")
  text(0.5, 0.5, txt, cex=cex)
}

merge2Clusters = function(mycolorcode, mainclusterColor, minorclusterColor){
  mycolorcode2 = ifelse(as.character(mycolorcode)==minorclusterColor, mainclusterColor, as.character(mycolorcode) )
  fcolorcode   =factor(mycolorcode2)
  fcolorcode
}



###############################################################################################################################
# I) GENERAL STATISTICAL FUNCTIONS

mean1=function(x) mean(x,na.rm=T)
var1=function(x) var(x,na.rm=T)



if (exists("scatterplot1") ) rm(scatterplot1);
scatterplot1=function(x,y, title = "", ... ){
  cor1=signif(cor(x,y,use="p",method="s"),2)
  corp=signif(cor.test(x,y,use="p",method="s")$p.value,2)
  if (corp<10^(-20) ) corp="<10^{-20}"
  plot(x,y, main=paste(title, "cor=", cor1,"p=",corp),...) 
}


# Functions for Network Screening



if(exists("HubGeneSignificance") ) rm(HubGeneSignificance)
HubGeneSignificance=function(k,GS,module,NumberHubs=10){
  colorlevels=levels(factor(module))
  HGS=rep(NA, length(colorlevels) )
  for (i in c(1:length(colorlevels) ) ) {
    restmodule= module==colorlevels[i]
    GSModule= GS[restmodule]
    kModule=k[restmodule] 
    HGS[i]=mean(GSModule[rank(-kModule,  ties.method="first")<=NumberHubs],na.rm=T)
  } # end of for loop
  barplot(HGS,col=colorlevels, ylab="Mean Hub Gene Significance", main=paste("Hub Gene Significance, top", NumberHubs, "hubs"))
  datout=data.frame(matrix(HGS,nrow=1,ncol=length(colorlevels)));
  names(datout)=colorlevels;
  datout
} # end of function




if(exists("StandardScreening1") ) rm(StandardScreening1); 
StandardScreening1=function( GS, LN=c(5,10)  ) {
  datout=data.frame(matrix(F, nrow=length(GS),ncol=length(LN) ))
  names(datout)=paste("List", LN,sep="")
  for ( j in c(1:length(LN))  ) {
    datout[,j]=rank(-GS,ties.method="first")<=LN[j]
  }
  datout
} # end of function StandardScreening1




if(exists("NetworkScreening1") ) rm(NetworkScreening1)
NetworkScreening1=function(datE, GS, MLN=1000, LN=10,beta=6,  minModuleSize=100,powerAllocation=3, NumberHubsHGS=50 , excludegrey=F,excludeturquoise=F, consider.sign.corKGS=F) {
  if (length(GS) != dim(datE)[[2]] ) print("Error: length(GS) not compatible with datE. Please check whether length(GS) = dim(datE)[[2]]?") 
  if ( max(LN)>MLN ) print("Error: requested list number bigger than maximum list number. Please increase MLN or decrease LN") 
  if ( MLN<minModuleSize) print("Warning: maximum list number smaller than minModuleSize-->every gene is grouped into the turquoise module") 
  if (length(GS) == dim(datE)[[2]]  &  max(LN) <= MLN  ){
    GS=abs(GS)
    restMLN=rank(-GS, ties.method="first")<=MLN
    GSrest=GS[restMLN]
    ADJ= abs(cor(datE[, restMLN], use="p"))^beta
    h1=hclust( as.dist(TOMdist1(ADJ)) , method="average")
    colorGS=rep("turquoise", dim(ADJ)[[2]] )
    if (MLN > minModuleSize ) {colorGS= factor(cutreeDynamic(h1, minModuleSize= minModuleSize)) }
    colorGSlevels=levels(factor(colorGS))
    if ( length(colorGSlevels)==1   ) { HGS=1; k=apply(ADJ,2,sum) } else {
      k= DegreeInOut(ADJ ,colorGS,scaledToOne=FALSE)$kWithin
      HGS=as.vector(as.matrix(HubGeneSignificance(k=k, GS=GSrest,module=colorGS,NumberHubs=NumberHubsHGS)[1,])) 
      # the following definition of HGS uses the correlation between k and GS...
      if ( consider.sign.corKGS==T ) {
        for (i in c(1:length(colorGSlevels))){
          HGS[i]=(cor( k[colorGS==colorGSlevels[i]], GSrest[colorGS==colorGSlevels[i]], use="p" ))}
      } # end of for loop
    } # end of if (consider.sign  )
    if (excludegrey) HGS[colorGSlevels=="grey"]=0 
    if (excludeturquoise==T) HGS[colorGSlevels=="turquoise"]=0 
    weightHGS=abs(HGS)^powerAllocation
    weightHGS=weightHGS/sum(weightHGS)
    datout=data.frame(matrix(F, nrow=length(GS),ncol=length(LN) ))
    names(datout)=paste("List", LN,sep="")
    for ( j in c(1:length(LN))  ) {
      LNcolor=round(LN[j]*weightHGS)
      maxNoColorGS=tapply(colorGS,colorGS,length)
      # if there are more hubs than the module size, then we pick genes from the next module
      # with highest hub gene significance
      for (ii in c(1:length(colorGSlevels))) { 
        ExcessNo=sum(c(LNcolor-maxNoColorGS)[LNcolor>maxNoColorGS])
        indexHighestWeightHGS=rank(-weightHGS, ties.method="first")==ii
        LNcolor[LNcolor>maxNoColorGS]=maxNoColorGS[LNcolor>maxNoColorGS]
        LNcolor[indexHighestWeightHGS  ]= LNcolor[indexHighestWeightHGS  ]+ExcessNo
        LNcolor[indexHighestWeightHGS  ]=LN[j]- sum(LNcolor[!indexHighestWeightHGS])
      }
      LNcolor[LNcolor>maxNoColorGS]=maxNoColorGS[LNcolor>maxNoColorGS]
      PickHubs=rep(F,length(GS))
      for (i in c(1:length(colorGSlevels))){
        # the following takes the sign between K and GS into account
        signKGS=1
        if ( consider.sign.corKGS==T ) {signKGS=sign(cor( k[colorGS==colorGSlevels[i]], GSrest[colorGS==colorGSlevels[i]], use="p" ))}
        PickHubs[restMLN][ colorGS==colorGSlevels[i]] =rank(-signKGS*k[colorGS==colorGSlevels[i]], ties.method="first")<=LNcolor[i] 
      } # end of for loop
      datout[,j]=PickHubs
      print(paste("For list", j, "with LN=",LN[j], "the algorithm picks the following number of hubs in each module")) 
      print( data.frame(colorGSlevels, modulesize=table(colorGS), HubGeneSignif=signif(HGS,2), weightHGS=signif(weightHGS,2),  LNcolor    ))
    } # end of for ( j in c(1:length(LN))  ) 
    datout
  } # end of if if (length(GS) == dim(datE)[[2]]  &  max(LN) < MLN
} # end of function NetworScreening1


#==========================================================================================================
#
# Peter Langfelder's additions
#
#==========================================================================================================

# PrintFlush.R

if (exists("print.flush")) { remove(print.flush); collect_garbage(); }
print.flush = function(...)
{
  x = print(...)
  if (exists("flush.console")) x=flush.console();
}

if (exists("PrintSpaces")) { remove(PrintSpaces); collect_garbage(); }
PrintSpaces = function(print.level)
{
  if (print.level>0) 
  {
    spaces = paste(" ",rep("  ", times=print.level-1), collapse="");
  } else
  {
    spaces = "";
  }
  spaces;
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

#--------------------------------------------------------------------------
#
# BarplotWithTextLabels = function ( Matrix, Labels, ... ) { 
#
#--------------------------------------------------------------------------
#
# Plots a barplot of the Matrix and writes the Labels underneath such that they are readable.

BarplotWithTextLabels = function ( Matrix, Labels, ColorLabels = FALSE, Colored = TRUE, 
                                   SetMargins = TRUE, Errors = NULL, ... ) 
{ 
  if (SetMargins) par(mar=c(3,3,2,2)+0.2)
  
  if (Colored)
  {
    colors = substring(Labels, 3);
  } else {
    colors = rep("grey", times = ifelse(length(dim(Matrix))<2, length(Matrix), dim(Matrix)[[2]]));
  }
  
  mp = barplot(Matrix, col = colors, xaxt = "n", xlab="", yaxt="n", ylab="", ...)
  
  if (length(dim(Matrix))==2) {
    means = apply(Matrix, 2, sum);
  } else {
    means = Matrix;  
  }
  
  err.bp(means, 1.96*Errors, two.side = T);
  
  # axis(1, labels = FALSE)
  nlabels = length(Labels)
  plotbox = par("usr");
  xmin = plotbox[1]; xmax = plotbox[2]; ymin = plotbox[3]; yrange = plotbox[4]-ymin;
  ymax = plotbox[4];
  # print(paste("yrange:", yrange));
  if (nlabels>1)
  {
    spacing = (mp[length(mp)] - mp[1])/(nlabels-1);
  } else {
    spacing = (xmax-xmin);
  }
  yoffset = yrange/30
  xshift = spacing/2;
  xrange = spacing * nlabels;
  if (ColorLabels)
  {
    #rect(xshift + ((1:nlabels)-1)*spacing - spacing/2.1, ymin - spacing/2.1 - spacing/8,
    #     xshift + ((1:nlabels)-1)*spacing + spacing/2.1, ymin - spacing/8,
    #     density = -1,  col = substring(Labels, 3), border = substring(Labels, 3), xpd = TRUE)
    rect(mp - spacing/2.1, ymin - 2*spacing/2.1 * yrange/xrange - yoffset,
         mp + spacing/2.1, ymin - yoffset,
         density = -1,  col = substring(Labels, 3), border = substring(Labels, 3), xpd = TRUE)
  } else {
    text(((1:nlabels)-1)*spacing +spacing/2 , ymin - 0.02*yrange, srt = 45, 
         adj = 1, labels = Labels, xpd = TRUE)
  }
  axis(2, labels = T)
}

#--------------------------------------------------------------------------
#
# SizeWindow
#
#--------------------------------------------------------------------------
# if the current device isn't of the required dimensions, close it and open a new one.

SizeWindow = function(width, height)
{
  din = par("din");
  if ( (din[1]!=width) | (din[2]!=height) )
  {
    dev.off();
    X11(width = width, height=height);
  }
}

#======================================================================================================
# GreenToRed.R
#======================================================================================================

GreenBlackRed = function(n)
{
  half = as.integer(n/2);
  red = c(rep(0, times = half), 0, seq(from=0, to=1, length.out = half));
  green = c(seq(from=1, to=0, length.out = half), rep(0, times = half+1));
  blue = rep(0, times = 2*half+1);
  col = rgb(red, green, blue, maxColorValue = 1);
}

GreenWhiteRed = function(n)
{
  half = as.integer(n/2);
  red = c(seq(from=0, to=1, length.out = half), rep(1, times = half+1));
  green = c(rep(1, times = half+1), seq(from=1, to=0, length.out = half));
  blue = c(seq(from=0, to=1, length.out = half), 1, seq(from=1, to=0, length.out = half));
  col = rgb(red, green, blue, maxColorValue = 1);
}

RedWhiteGreen = function(n)
{
  half = as.integer(n/2);
  green = c(seq(from=0, to=1, length.out = half), rep(1, times = half+1));
  red = c(rep(1, times = half+1), seq(from=1, to=0, length.out = half));
  blue = c(seq(from=0, to=1, length.out = half), 1, seq(from=1, to=0, length.out = half));
  col = rgb(red, green, blue, maxColorValue = 1);
}

#======================================================================================================
# AverageExpression.R
#======================================================================================================

# Getting a mean of expression data classified by a factor: 
# The colors are stored in a vector colorhdataOne. 
# Getting a mean of columns of a dataframe juts requires a strightforward application of mean.
# or apply(frame, 2, mean).
# Then can do by(frame, factor, function)
# but it returns an object of class "by" which seems to be a list with elements named by the factors.

AverageExprMatrix = function(NormExprData, colors) {
  no.genes = dim(NormExprData)[2]
  no.samples = dim(NormExprData)[1]
  colorsf = as.factor(colors)
  AverageExpr = matrix(ncol=nlevels(colorsf), nrow = no.samples)
  ExprDataMtrx = as.matrix(NormExprData)
  for (i in (1:no.samples)) AverageExpr[i,] = tapply(ExprDataMtrx[i,], colorsf, mean)
  AverageExpr
}

# A wrapper that turns the average expression matrix into a data frame 
# with appropriate column and row names 
AverageExprFrame = function(NormExprData, colors) 
{
  tmp = as.data.frame(AverageExprMatrix(NormExprData, colors), 
                      row.names = row.names(NormExprData) )
  names(tmp) = paste("AE", levels(as.factor(colors)), sep="")
  tmp
}

# Calculate correlation of Matrix1[,i] with Matrix2[,i] and return as a vector
ColumnCor = function ( Mtrx1, Mtrx2 ) 
{
  if ( var(dim(Mtrx1) - dim(Mtrx2))!=0 ) { 
    stop("ColumnCorr: Error: Given matrices have different dimensions.")
  }
  corrs = vector(mode="numeric", length=dim(Mtrx1)[2] )
  for (i in (1:dim(Mtrx1)[2]) ) corrs[i] = cor(Mtrx1[,i], Mtrx2[,i], use="p")
  names(corrs) = paste(names(Mtrx1), names(Mtrx2))
  corrs
}


#======================================================================================================
# NetworkAndModules
#======================================================================================================

MaxSets = 20;

#-------------------------------------------------------------------------------------------
#
# CheckSets
#
#-------------------------------------------------------------------------------------------
# Checks sets for consistency and returns some diagnostics.

CheckSets = function(ExprData)
{
  No.Sets = length(ExprData);
  if (No.Sets<=0) stop("No expression data given!");
  if (No.Sets>MaxSets) stop(paste("Number of datasets is excessive. To use more than", MaxSets, "sets,",
                                  "change the line MaxSets =", MaxSets, 
                                  "to your needed value and run the script again."));
  
  No.Samples = vector(length = No.Sets);
  
  No.Genes = dim(ExprData[[1]]$data)[2];
  for (set in 1:No.Sets) 
  {
    if (No.Genes!=dim(ExprData[[set]]$data)[2])
      stop(paste("Incompatible number of genes in set 1 and", set));
    No.Samples[set] = dim(ExprData[[set]]$data)[1];
  }
  
  list(ngenes = No.Genes, nsamples = No.Samples);
}

#-------------------------------------------------------------------------------------------
#
# KeepCommonProbes
#
#-------------------------------------------------------------------------------------------
# Filters out probes that are not common to all datasets, and puts probes into the same order in each
# set. Works by creating dataframes of probe names and their indices and merging them all.

KeepCommonProbes = function(ExprData, OrderBy = 1)
{
  No.Sets = length(ExprData);
  if (No.Sets<=0) stop("No expression data given!");
  if (No.Sets>MaxSets) stop(paste("Number of datasets is excessive. To use more than", MaxSets, "sets,",
                                  "change the line MaxSets =", MaxSets, 
                                  "to your needed value and run the script again."));
  
  Names = data.frame(Names = names(ExprData[[OrderBy]]$data));
  
  if (No.Sets>1) for (set in (1:No.Sets))
  {
    SetNames = data.frame(Names = names(ExprData[[set]]$data), index = c(1:dim(ExprData[[set]]$data)[2]));
    Names = merge(Names, SetNames, by.x = "Names", by.y = "Names", all = FALSE, sort = FALSE);
  }
  
  for (set in 1:No.Sets)
    ExprData[[set]]$data = ExprData[[set]]$data[, Names[, set+1]];
  
  ExprData;
}

#-------------------------------------------------------------------------------------------
#
# GetConnectivity
#
#-------------------------------------------------------------------------------------------
# This function takes expression data (rows=samples, colummns=genes), 
# and the SoftPower exponent used in weighting the
# correlations to get the network adjacency matrix, and returns an array of dimensions
# No.Genes * No.Sets containing the connectivities of each gene in each subset.

# Caution: KeepOverlapSign and negative powers are not functional.

GetConnectivity = function(ExprData, SoftPower=6, KeepOverlapSign = FALSE, 
                           verbose=1, print.level=0, BatchSize = 1500)
{
  No.Sets = length(ExprData);
  spaces = PrintSpaces(print.level);
  setsize = CheckSets(ExprData);
  No.Genes = setsize$ngenes;
  No.Samples = setsize$nsamples;
  
  if (verbose>1) print.flush(paste(spaces, "GetConnectivity: received", No.Sets, "datasets with No.Genes =", 
                                   as.character(No.Genes)));
  if (verbose>1) print.flush(paste(spaces, "  Received No.Sets =", as.character(No.Sets)));
  
  Connectivity = matrix(nrow = No.Genes, ncol = No.Sets);
  
  for (set in 1:No.Sets) 
  {
    if (verbose>1) print.flush(paste(spaces, "  Working on set", set));
    No.Batches = as.integer((No.Genes-1)/BatchSize);
    SetRestrConn = NULL;
    for (batch in 1:(No.Batches+1))
    {
      if (batch<=No.Batches)
      {
        if (verbose>2) print.flush(paste(spaces, "    Working on batch", batch));
        BatchIndex = c(1:BatchSize) + (batch-1)*BatchSize;
      } else {
        BatchIndex = c( (BatchSize*(batch-1)+1):No.Genes)
      }
      
      #      adj_mat = AdjacencyMatrix(ExprData[[set]]$data, ExprData[[set]]$data[, BatchIndex], SoftPower = SoftPower,
      #                                KeepSign = KeepOverlapSign, 
      #                                verbose = verbose-1, print.level = print.level+1);
      adj_mat = AdjacencyMatrixR(ExprData[[set]]$data, ExprData[[set]]$data[, BatchIndex], 
                                 SoftPower = SoftPower,
                                 verbose = verbose-1, print.level = print.level+1);
      adj_mat[is.na(adj_mat)] = 0;
      Connectivity[BatchIndex, set] = apply(adj_mat, 2, sum)-1;
    }
  }
  Connectivity;
}

#-------------------------------------------------------------------------------
#
# SelectGenesByConnectivity
#
#-------------------------------------------------------------------------------
# 
# Here network genes are selected based on 1. nonzero variance in all subsets, 2.
# ranking in a given subset. The retrun value is a boolean vector of length No.Genes
# signifying whether the corresponding gene is to be included in the network. 

SelectGenesByConnectivity = function(ExprData, Connectivity, Subs.Ind=1,
                                     DegreeCut=3600, verbose=1, print.level=0)
{
  spaces = PrintSpaces(print.level);
  
  No.Sets = length(ExprData);
  setsize = CheckSets(ExprData);
  No.Genes = setsize$ngenes;
  No.Samples = setsize$nsamples;
  
  variance = matrix(ncol = No.Sets, nrow = No.Genes)
  
  for (set in (1:No.Sets)) 
  {
    variance[ ,set] = as.vector(apply(ExprData[[set]]$data, 2, var, na.rm=T))
  }
  VarProduct = as.vector(apply(variance, 1, prod))
  if ( (DegreeCut>0) & (DegreeCut<No.Genes))
  {
    DegreeRank = rank(-Connectivity[, Subs.Ind], ties.method = "first");	
    SelectedGenes = DegreeRank <= DegreeCut & VarProduct > 0
  } else {
    SelectedGenes = VarProduct > 0;
  }
  if (verbose>0) 
    print.flush(paste(spaces, "SelectGenesByConnectivity:", as.character(sum(VarProduct==0)), 
                      "genes have zero variance in at least one subset, selected",
                      sum(SelectedGenes), "probes."));
  SelectedGenes;
}


#-------------------------------------------------------------------------------
#
# SelectGenesByMinConnectivity
#
#-------------------------------------------------------------------------------
# 
# Here network genes are selected based on 1. nonzero variance in all subsets, 2.
# ranking of the minimum (over all subsets) connectivity.

SelectGenesByMinConnectivity = function(ExprData, Connectivity, Subs.Ind=1,
                                        DegreeCut=3600, verbose=1, print.level=0)
{
  spaces = PrintSpaces(print.level);
  
  No.Sets = length(ExprData);
  setsize = CheckSets(ExprData);
  No.Genes = setsize$ngenes;
  No.Samples = setsize$nsamples;
  
  variance = matrix(ncol = No.Sets, nrow = No.Genes)
  
  for (set in 1:No.Sets) 
  {
    variance[ ,set] = as.vector(apply(ExprData[[set]]$data, 2, var, na.rm=T))
  }
  VarProduct = as.vector(apply(variance, 1, prod));
  if ( (DegreeCut>0) & (DegreeCut<No.Genes))
  {
    MinConn = as.vector(apply(Connectivity, 1, min));
    ConnRank = rank(-MinConn, ties.method = "first");
    SelectedGenes = ConnRank <= DegreeCut & VarProduct > 0;
  } else {
    SelectedGenes = VarProduct > 0;
  }
  
  if (verbose>0) 
    print.flush(paste(spaces, "SelectGenesByMinConnectivity:", as.character(sum(VarProduct==0)), 
                      "genes have zero variance in at least one subset, selected",
                      sum(SelectedGenes), "probes."));
  SelectedGenes;
}

#-------------------------------------------------------------------------------
#
# SelectGenesByRestrConnectivity
#
#-------------------------------------------------------------------------------
# 
# Here network genes are selected based on 1. nonzero variance in all subsets, 2.
# calculating the "N best" gene-gene connectivities, 3. ranking their sum
# Assume NAs have been imputed in the data.

SelectGenesByRestrConnectivity = function(ExprData, NBest = 20, SoftPower = 6,
                                          DegreeCut=3600, BatchSize = 1500, verbose=1, print.level=0, 
                                          ReturnDiags = FALSE)
{
  spaces = PrintSpaces(print.level);
  
  No.Sets = length(ExprData);
  setsize = CheckSets(ExprData);
  No.Genes = setsize$ngenes;
  No.Samples = setsize$nsamples;
  
  if (verbose>1) 
    print.flush(paste(spaces, "SelectGenesByRestrConnectivity: received a dataset with No.Genes =", 
                      as.character(No.Genes), " and No.Samples =", as.character(sum(No.Samples))));
  
  if (verbose>1) print.flush(paste(spaces, "  Received No.Sets =", as.character(No.Sets)));
  
  RestrConnectivity = matrix(nrow = No.Genes, ncol = No.Sets);
  
  for (set in 1:No.Sets)
  {
    if (verbose>1) print.flush(paste(spaces, "  Working on set", set));
    SetExprData = ExprData[[set]]$data;
    No.Batches = as.integer(No.Genes/BatchSize);
    SetRestrConn = NULL;
    for (batch in 1:(No.Batches+1))
    {
      if (verbose>2) print.flush(paste(spaces, "    Working on batch", batch));
      if (batch<=No.Batches)
      {
        BatchIndex = c(1:BatchSize) + (batch-1)*BatchSize;
      } else {
        BatchIndex = c( (BatchSize*(batch-1)+1):No.Genes)
      }
      if (verbose>2) print.flush(paste(spaces, "      Calculating adjacencies..."));
      BatchExprData = SetExprData[, BatchIndex];
      # In the following BatchConn[i,j] = cor(BatchExprData[,i], ExprData[,j])
      BatchConn = AdjacencyMatrixR(SetExprData, BatchExprData, SoftPower, 
                                   verbose = verbose - 1, print.level = print.level +1);
      BatchConn[is.na(BatchConn)] = 0;
      if (verbose>2) print.flush(paste(spaces, "      Sorting adjacencies..."));
      SortedBatchConn = apply(BatchConn, 2, sort, decreasing = TRUE);
      if (is.null(SetRestrConn)) 
      {
        SetRestrConn = SortedBatchConn[c(1:NBest), ];
      } else {
        SetRestrConn = cbind(SetRestrConn, SortedBatchConn[c(1:NBest),]);
      }
    }
    if (verbose>2) print.flush(paste(spaces, "    Summing the highest", NBest, "adjacencies..."));
    RestrConnectivity[, set] = apply(SetRestrConn, 2, sum);
    collect_garbage();
  }
  
  
  # The variance checking should strictly spekaing not be necessary, since zero variance would lead to a
  # correlation of NA and hence a zero minimum connectivity... but keep it in to remove such probes in a
  # case when the number of non-zero-variance probes is less than DegreeCut.
  
  if (verbose>2) print.flush(paste(spaces, " Checking non-zero variance..."));
  variance = matrix(ncol = No.Sets, nrow = No.Genes)
  
  for (set in 1:No.Sets) {
    variance[ ,set] = as.vector(apply(ExprData[[set]]$data, 2, var, na.rm=T))
  }
  
  VarProduct = as.vector(apply(variance, 1, prod));
  
  # This removes zero-variance probes from the ranking (again).
  
  RestrConnectivity[VarProduct==0, ] = 0;
  
  # Take the minimum of the restricted sums of connectivities...
  
  MinConn = as.vector(apply(RestrConnectivity, 1, min));
  ConnRank = rank(-MinConn, ties.method = "first");
  
  # Select genes based on the rank of the minima.
  
  if (verbose>2) print.flush(paste(spaces, " Selecting genes..."));
  SelectedGenes = ConnRank <= DegreeCut & VarProduct > 0
  if (verbose>0) 
    print.flush(paste(spaces, "  SelectGenesByRestrConnectivity:", as.character(sum(VarProduct==0)), 
                      "genes have zero variance in at least one subset, selected",
                      sum(SelectedGenes), "probes."));
  collect_garbage();
  if (ReturnDiags)
  {
    RetVal = list(SelectedGenes = SelectedGenes, RestrConnectivity = RestrConnectivity);
  } else {
    RetVal = SelectedGenes;
  }
  
  RetVal;
}

#---------------------------------------------------------------------
#
# TrafoCorMatrix
#
#---------------------------------------------------------------------
# Transforms correlation matrix to filter noise. Caution, this function discards sign of the correlations
# in the matrix.
# 

TrafoCorMatrix = function(Matrix, method, Power, No.Samples) 
{
  if (method=="power")
  {
    return (abs(Matrix)^Power); 
  } else if (method=="probability")
  { 
    abs_M = abs(Matrix);
    norm = pnorm(1, sd = 1/sqrt(No.Samples));
    raw_weight = pnorm(abs_M, sd = 1/sqrt(No.Samples));
    weight = (raw_weight - 0.5)/(norm - 0.5); 
    return (abs_M * weight^Power);
  } else stop("Unrecognized \'method\' given.");
}


#---------------------------------------------------------------------
#
# AdjacencyMatrix
#
#---------------------------------------------------------------------
# Computes the adjacency from the expression data: takes cor, transforms it as appropriate and possibly
# adds a sign if requested. No subselection on ExprData is performed.

AdjacencyMatrix = function(ExprData, ExprData2=NULL, SoftPower, KeepSign = FALSE,
                           verbose=1, print.level = 0)
{
  spaces = PrintSpaces(print.level);
  No.Samples = dim(ExprData)[1];
  
  if (SoftPower<=0)
  {
    method = "probability"; Power = -SoftPower;
  } else {
    method = "power"; Power = SoftPower;
  }
  if (verbose>2) print.flush(paste(spaces, "Transforming the correlation matrix using", method,
                                   "method with power", Power));
  cor_mat = cor(ExprData, ExprData2, use="p");
  # trafoed_cor = TrafoCorMatrix(cor_mat, method, SoftPower, No.Samples); collect_garbage();
  if (method=="power")
  {
    trafoed_cor = abs(cor_mat)^Power;
  } else 
  {
    abs_M = abs(cor_mat);
    norm = pnorm(1, sd = 1/sqrt(No.Samples));
    raw_weight = pnorm(abs_M, sd = 1/sqrt(No.Samples));
    weight = (raw_weight - 0.5)/(norm - 0.5); 
    trafoed_cor = abs_M * weight^Power;
  }
  if (KeepSign)
  {
    sign_cm = sign(cor_mat);
    adj_mat = sign_cm * trafoed_cor;
    #print("Keeping overlap sign");
    #print(paste("No. of negative entries: ", sum(gtom0<0)));
  } else {
    adj_mat = trafoed_cor;
  }
  rm(trafoed_cor, cor_mat);
  collect_garbage();
  adj_mat;
}

# A faster and less memory-intensive version, but it only works for method==power and it won't keep
# overlap sign

AdjacencyMatrixR = function(ExprData, ExprData2=NULL, SoftPower, 
                            verbose=1, print.level = 0)
{
  abs(cor(ExprData, ExprData2, use="p"))^SoftPower;
}

#---------------------------------------------------------------------
#
# GetModules
#
#---------------------------------------------------------------------

# Get the vector of module colors for each gene in the given dataset.
# The following code computes the topological overlap matrices for the selected genes
# in each subset: 

# CAUTION:  KeepOverlapSign and negative soft powers are not functional in this version.

GetModules = function(ExprData, SelectedGenes, SoftPower = 6, 
                      BranchHeightCutoff = 0.94, ModuleMinSize = 125, 
                      GetDissTOM = NULL, DissimilarityLevel = 1, KeepOverlapSign = FALSE,
                      ClusterOnlyOne = NULL, 
                      verbose = 1, print.level = 0 ) 
{
  spaces = PrintSpaces(print.level);
  
  if (verbose>0) print.flush(paste(spaces, "Entering GetModules:"));
  
  No.Sets = length(ExprData);
  setsize = CheckSets(ExprData);
  No.Genes = setsize$ngenes;
  No.Samples = setsize$nsamples;
  
  No.NetworkGenes = sum(SelectedGenes==TRUE)
  
  dissTOM = vector(mode="list", length = No.Sets)
  hierTOM = vector(mode="list", length = No.Sets);
  colordata = array(dim = c(No.NetworkGenes, No.Sets));
  
  if (!is.null(GetDissTOM))
  {
    if (!is.null(ClusterOnlyOne))
    {
      if (GetDissTOM!=ClusterOnlyOne)
      {
        stop(paste("GetModules: Error: Incompatible parameters given:",
                   "if both given, GetDissTOM and ClusterOnlyOne must equal. Their values are", 
                   GetDissTOM, ClusterOnlyOne));
      }
    } else 
    {
      ClusterOnlyOne = GetDissTOM;
    }
  }
  
  if (is.null(ClusterOnlyOne))
  {
    ClusterOnlyOne_n = 0;
  } else
  { 
    ClusterOnlyOne_n = ClusterOnlyOne;
  }
  
  if (is.null(GetDissTOM))
  {
    GetDissTOM_n = 0;
  } else
  { 
    GetDissTOM_n = ClusterOnlyOne;
  }
  
  if ((!is.null(GetDissTOM)) && (!is.element(GetDissTOM_n, c(1:No.Sets))))
  {
    if (verbose > 0) print.flush(paste(spaces, 
                                       "  Calculation of DissTOM is disabled by GetDissTOM out of valid range."));
  } else {
    if (verbose>0) print.flush(paste(spaces, 
                                     "  Calculating TOM and clustering trees..."))
  }
  for (set in 1:No.Sets)
  {
    if (is.null(GetDissTOM) | (GetDissTOM_n==set))
    {
      if (verbose>1) print.flush(paste(spaces, 
                                       "  Calculating TOM in subset", as.character(set)))
      #      gtom0 = AdjacencyMatrix(ExprData[[set]]$data[, SelectedGenes], SoftPower = SoftPower, 
      #                              KeepSign = KeepOverlapSign,
      #                              verbose = verbose, print.level = print.level+1);
      gtom0 = AdjacencyMatrixR(ExprData[[set]]$data[, SelectedGenes], SoftPower = SoftPower, 
                               verbose = verbose, print.level = print.level+1);
      if (DissimilarityLevel==0)
      {
        dissTOM[[set]] = list(data = 1-gtom0)
      } else 
      {
        dissTOM[[set]] = list(data = SignedTOMdist(gtom0))
      }
      rm(gtom0); collect_garbage();
    }
    collect_garbage()
  }
  if ((!is.null(ClusterOnlyOne)) && (!is.element(ClusterOnlyOne_n, c(1:No.Sets))))
  {
    if (verbose > 0) print.flush(paste(spaces, 
                                       "  Module detection is disabled by ClusterOnlyOne out of valid range."));
  } 
  for (set in 1:No.Sets) 
  {
    if (is.null(ClusterOnlyOne) | (ClusterOnlyOne_n==set))
    {
      if (verbose>1) print.flush(paste(spaces, 
                                       "  Calculating clustering tree in set", set))
      hierTOM[[set]] = 
        list(data = hclust(as.dist(dissTOM[[set]]$data),method="average"));
      collect_garbage()
      colordata[,set] = as.character(modulecolor2(hierTOM[[set]]$data, 
                                                  h1=BranchHeightCutoff, minsize1=ModuleMinSize))
      collect_garbage();
    }
  }
  if (!is.null(ClusterOnlyOne) & is.element(ClusterOnlyOne_n, c(1:No.Sets)))
  {
    for (set in 1:No.Sets) 
    {
      hierTOM[[set]] = hierTOM[[ClusterOnlyOne]];
      dissTOM[[set]] = dissTOM[[ClusterOnlyOne]];
      colordata[,set] = colordata[,ClusterOnlyOne];
    }
  }
  
  Modules = list( SelectedGenes = SelectedGenes, Dissimilarity = dissTOM, 
                  DissimilarityLevel = DissimilarityLevel,
                  ClusterTree = hierTOM, Colors = colordata, 
                  BranchHeightCutoff = BranchHeightCutoff, ModuleMinSize = ModuleMinSize, 
                  ClusterOnlyOne = ClusterOnlyOne);
  Modules;
}


#---------------------------------------------------------------------
#
#    GetNetwork
#
#---------------------------------------------------------------------

# ExprData should be a matrix or a data.frame; columns containing valid
# expression data can be specified in ExpressionColumns. By default all columns
# are taken. DegreeCut is the number of genes to be taken into the network
# construction (dissimilarity and clustering). If set to NULL or 0, all genes
# will be taken. 
# The return value is a list with the following components:

#network = list(SoftPower = SoftPower, Connectivity = Connectivity, IsInNetwork = SelectedGenes,
#               Dissimilarity = dissTOM, ClusterTree = hierTOM, Colors = colordata); 

# This function needs the functions in
# ../CommonFunctions/NetworkFunctions-PL.txt and ../CommonFunctions/PrintFlush.R to be loaded.

GetNetwork= function(ExprData, ProbeSelection = "set connectivity rank", 
                     DegreeCut = 3600, NBest = 20, SoftPower = 6, 
                     BranchHeightCutoff = 0.94, ModuleMinSize = 125,
                     GetDissTOM = NULL, DissimilarityLevel = 1, 
                     KeepOverlapSign = FALSE, ClusterOnlyOne = NULL, 
                     verbose = 1, print.level = 0 ) 
{
  spaces = PrintSpaces(print.level);
  No.Sets = length(ExprData);
  setsize = CheckSets(ExprData);
  No.Genes = setsize$ngenes;
  No.Samples = setsize$nsamples;
  
  if (verbose>1) print.flush(paste(spaces, "GetNetwork: received a dataset with No.Genes =", 
                                   as.character(No.Genes), " and No.Samples =", as.character(sum(No.Samples))));
  if (verbose>1) print.flush(paste(spaces, "  Received No.Sets =", as.character(No.Sets)));
  
  if (ProbeSelection=="set connectivity rank")
  {
    if ((DegreeCut>0) & (DegreeCut < No.Genes))
    {
      Connectivity = GetConnectivity(ExprData, SoftPower=SoftPower,
                                     verbose=verbose-1, print.level = print.level+1);
    } else {
      Connectivity = NULL;
    } 
    SelectedGenes = SelectGenesByConnectivity(ExprData, Connectivity, Subs.Ind=1,
                                              DegreeCut, verbose, print.level+1);
  } else if (ProbeSelection=="min connectivity rank")
  {
    if ((DegreeCut>0) & (DegreeCut < No.Genes))
    {
      Connectivity = GetConnectivity(ExprData, SoftPower=SoftPower,
                                     verbose=verbose-1, print.level = print.level+1);
    } else {
      Connectivity = NULL;
    } 
    SelectedGenes = SelectGenesByMinConnectivity(ExprData, Connectivity, Subs.Ind=1,
                                                 DegreeCut, verbose, print.level+1);
  } else if (ProbeSelection=="min restricted connectivity rank")
  {
    SelectedGenes =  SelectGenesByRestrConnectivity(ExprData, NBest = NBest, 
                                                    SoftPower, DegreeCut, verbose = verbose-1, print.level = print.level+1);
    Connectivity = NULL;
  } else if (ProbeSelection=="all")
  {
    SelectedGenes = rep(TRUE, times = No.Genes);
    Connectivity = NULL;
  } else
    stop(paste("GetNetwork: unrecognized ProbeSelection:", ProbeSelection));
  
  collect_garbage();collect_garbage();collect_garbage();
  
  No.NetworkGenes = sum(SelectedGenes==TRUE)
  
  # The following code computes the topological overlap matrices for the selected genes
  # in each subset: 
  
  Modules = GetModules(ExprData, SelectedGenes = SelectedGenes, SoftPower = SoftPower, 
                       BranchHeightCutoff, ModuleMinSize, GetDissTOM, DissimilarityLevel, KeepOverlapSign, 
                       ClusterOnlyOne, verbose, print.level);
  
  network = list(SoftPower = SoftPower, Connectivity = Connectivity, 
                 SelectedGenes = SelectedGenes, Dissimilarity = Modules$Dissimilarity,
                 DissimilarityLevel = Modules$DissimilarityLevel, 
                 ClusterTree = Modules$ClusterTree, Colors = Modules$Colors,
                 BranchHeightCutoff = Modules$BranchHeightCutoff, 
                 ModuleMinSize = Modules$ModuleMinSize, 
                 SignedDiss = KeepOverlapSign); 
  network
} 

kWithinModule = function(TOM, Colors, gene)
{
  if ((dim(TOM)[[1]]!=dim(TOM)[[2]]) || (length(Colors)!=dim(TOM)[[2]]))
  {
    stop("kWithinModule: Error: TOM is not a square matrix or the dimensions of Colors does not match.");
  }
  color = Colors[gene];
  GeneModuleTOM = TOM[Colors==color, gene];
  sum(GeneModuleTOM);
}

#---------------------------------------------------------------------
#
# PowerConnectivities
#
#---------------------------------------------------------------------
# Computes the connectivities for a given vector of SoftPowers. It is assumed that all entries
# of SoftPowers have the same sign for simplicity.
# No subselection is performed on ExprData

PowerConnectivities = function(SetExprData, SoftPowers, KeepSign = FALSE,
                               verbose=1, print.level = 0, BatchSize = 1536)
{
  spaces = PrintSpaces(print.level);
  
  No.Genes = dim(SetExprData)[2];
  No.Samples = dim(SetExprData)[1];
  
  No.Powers = length(SoftPowers);
  
  Connectivities = matrix(0, nrow = No.Genes, ncol = No.Powers);
  
  if (SoftPowers[1]<=0)
  {
    method = "probability"; Powers = -SoftPowers;
  } else {
    method = "power"; Powers = SoftPowers;
  }
  if (verbose>2) print.flush(paste(spaces, "Transforming the correlation matrix using", method,
                                   "method with powers", paste(Powers, collapse = ", ")));
  No.Batches = as.integer(No.Genes/BatchSize);
  for (batch in 1:(No.Batches+1))
  {
    if (verbose>2) print.flush(paste(spaces, "    Working on batch", batch));
    if (batch<=No.Batches)
    {
      BatchIndex = c(1:BatchSize) + (batch-1)*BatchSize;
    } else {
      BatchIndex = c( (BatchSize*(batch-1)+1):No.Genes)
    }
    if (verbose>2) print.flush(paste(spaces, "      Calculating adjacencies..."));
    BatchExprData = SetExprData[, BatchIndex];
    # In the following BatchConn[i,j] = cor(BatchExprData[,i], ExprData[,j])
    BatchCorr = cor(SetExprData, BatchExprData, use = "p")
    BatchCorr[is.na(BatchCorr)] = 0;
    
    for (power in 1:No.Powers) 
    {
      trafoed_cor = TrafoCorMatrix(BatchCorr, method, Powers[power], No.Samples); collect_garbage();
      if (KeepSign)
      {
        sign_cm = sign(cor_mat);
        adj_mat = sign_cm * trafoed_cor;
      } else {
        adj_mat = trafoed_cor;
      }
      # adj_mat[is.na(adj_mat)] = 0;
      Connectivities[BatchIndex, power] = apply(adj_mat, 2, sum) - 1;
    }
  }
  Connectivities;
}

#---------------------------------------------------------------------
#
# ScaleFreeAnalysis
#
#---------------------------------------------------------------------
# Analyzes the scale-free criterion for a given SetExprData.
# Here SetExprData is assumed to come from only one set and be in a plain matrix format.

ScaleFreeAnalysis = function(SetExprData, Powers, KeepSign = FALSE, 
                             No.HistBins = 10,
                             verbose=1, print.level = 0, BatchSize = 1536)
{
  spaces = PrintSpaces(print.level);
  if (verbose>0) print.flush(paste(spaces, "Calculating scale-free diagnostics"));
  
  # Get connectivities
  
  if (length(Powers[Powers>0]))
  {
    PositiveConnects = PowerConnectivities(SetExprData, Powers[Powers>0],  KeepSign, 
                                           verbose = verbose - 1, print.level = print.level + 1);
  } else
    PositiveConnects = NULL;
  
  if (length(Powers[Powers<=0]))
  {
    NegativeConnects = PowerConnectivities(SetExprData, Powers[Powers<=0],  KeepSign, 
                                           verbose = verbose - 1, print.level = print.level + 1);
  } else
    NegativeConnects = NULL;
  
  
  # Reorganize powers if necessary
  
  Powers = c(Powers[Powers>0], Powers[Powers<=0]);
  No.Powers = length(Powers);
  
  Connectivities = cbind(PositiveConnects, NegativeConnects);
  
  DiagNames = c("Power", "scale law R^2" ,"slope", "truncated R^2","mean(k)","median(k)","max(k)",
                "intercept", "a", "b", "c" );
  Diags = matrix(0, nrow = No.Powers, ncol = length(DiagNames));
  
  Diags[, 1] = Powers;
  
  Histos = vector(mode = "list", length = No.Powers); 
  
  # Calculate the fits and diagnostics for each power separately
  
  for (power in 1:No.Powers)
  {
    minC = min(Connectivities[, power]);
    maxC = max(Connectivities[, power]);
    #breaks = seq(from = 0, to = maxC * (1+1/(2*No.HistBins)) , length.out = No.HistBins +1);
    breaks = seq(from = log10(minC), to = log10(maxC), length.out = No.HistBins +1);
    h = hist(log10(Connectivities[, power]), breaks = breaks, plot = FALSE, 
             include.lowest = TRUE, right = TRUE);
    x = h$mids;
    y = log10(h$counts+0.1);
    fit = lm(y ~ x);
    xx = 10^x;
    fit2 = lm(y ~ x + xx);
    Diags[power, 2] = summary(fit)$adj.r.squared
    Diags[power, 3] = summary(fit)$coefficients[2,1]
    Diags[power, 4] = summary(fit2)$adj.r.squared
    Diags[power, 5] = mean(Connectivities[, power]);
    Diags[power, 6] = median(Connectivities[, power]);
    Diags[power, 7] = max(Connectivities[, power]);
    Diags[power, 8] = summary(fit)$coefficients[1,1];
    Diags[power, 9:11] = summary(fit2)$coefficients[1:3,1];
    Histos[[power]] = list(h = h);
  }
  
  # adjust the format...
  
  Diags = data.frame(Diags);
  names(Diags) = DiagNames;
  
  list(Diags = Diags, Histos = Histos);
} 


#======================================================================================================
# ModulePrincipalComponents-03.R
#======================================================================================================

#-------------------------------------------------------------------------------------
#
#  ModulePrincipalComponents
#
#-------------------------------------------------------------------------------------

# Calculates the principal components of modules of a given network.
#   - - should multiply them by -1 wherever appropriate.
# Input: Data: expression data, module colors. AlignPCs can take the values "", "along average".
# output : a dataframe of principal components.

ModulePrincipalComponents = function(Data, ModuleColors, AlignPCs = "along average", Impute = FALSE, 
                                     verbose = 1, print.level=0) 
{
  spaces = PrintSpaces(print.level);
  
  AlignPCsRecognizedValues =  c("", "along average");
  if (!is.element(AlignPCs, AlignPCsRecognizedValues)) {
    print.flush(paste("ModulePrincipalComponents: Error:",
                      "parameter AlignPCs has an unrecognised value:", 
                      AlignPCs, "; Recognized values are ", AlignPCsRecognizedValues));
    stop()
  }
  
  if (verbose>0) print.flush(paste(spaces, "ModulePrincipalComponents: Calculating PCs"));
  FullPCs = ModulePrinComps1(Data, ModuleColors, verbose = verbose-1, print.level = print.level+1,
                             GetConformity = FALSE, Impute = Impute);
  PCs = FullPCs$PrinComps;
  
  if (AlignPCs == "") AlignedPCs = PCs;
  if (AlignPCs == "along average") 
  {
    if (verbose>0) print.flush(paste(spaces,"ModulePrincipalComponents:", 
                                     "Aligning PCs with average expression for each module."))
    if (verbose>1) print.flush(paste(spaces,"  ++ Calculating averages..."));
    NormData = scale(Data);
    AverageModuleExpr = data.frame(AverageExprMatrix(NormData, ModuleColors));
    if (verbose>1) print.flush(paste(spaces,"  ++ Aligning principal components..."));
    AverageAndPCCor = diag(cor(PCs, AverageModuleExpr, use = "pairwise.complete.obs"));
    sign.matrix = matrix(0, nrow = dim(PCs)[2], ncol = dim(PCs)[2]);
    diag(sign.matrix) = sign(AverageAndPCCor);
    AlignedPCs = as.data.frame(as.matrix(PCs) %*% sign.matrix);
    names(AlignedPCs) = names(PCs);
    rownames(AlignedPCs) = rownames(PCs);
    names(AverageModuleExpr) = names(PCs);
    rownames(AverageModuleExpr) = rownames(PCs);
    if (verbose>1) print.flush(paste(spaces,"  ++ done."));
  }
  RetPCs = list(data = AlignedPCs, VarExplained = FullPCs$varexplained, 
                ModuleConformity = FullPCs$ModuleConformity, AverageExpr = AverageModuleExpr);
  RetPCs;
}

#--------------------------------------------------------------------------------------
#
# ModulePCs
# NetworkModulePCs
#
#--------------------------------------------------------------------------------------

ModulePCs = function(ExprData, SelectedGenes, ModuleColors, UniversalModuleColors = NULL,
                     OnlySet = NULL,
                     AlignPCs="along average", Impute = FALSE,  verbose=1, print.level=0)
{
  spaces = PrintSpaces(print.level)
  No.Sets = length(ExprData);
  setsize = CheckSets(ExprData);
  No.Genes = setsize$ngenes;
  No.Samples = setsize$nsamples;
  
  if (verbose>0) print.flush(paste(spaces,"ModulePCs: Looking for module PCs."));
  PCs = vector(mode="list", length=No.Sets);
  if (is.null(OnlySet))
  {
    CalculatedSubsets = c(1:No.Sets);
  } else {
    CalculatedSubsets = c(OnlySet);
  }
  for (set in CalculatedSubsets) {
    if (verbose>0) print.flush(paste(spaces,"  Working on subset", as.character(set), "...")); 
    if (is.null(UniversalModuleColors))
    {
      SubsetColors = ModuleColors[,set]; 
    } else { 
      SubsetColors = UniversalModuleColors; }
    SubsetPCs = ModulePrincipalComponents(Data = ExprData[[set]]$data[,SelectedGenes==1],
                                          ModuleColors = SubsetColors, AlignPCs = AlignPCs, Impute = Impute, 
                                          verbose = verbose-1, print.level = print.level+1);
    PCs[[set]] = list(data = SubsetPCs$data, AverageExpr = SubsetPCs$AverageExpr, 
                      ModuleConformity = SubsetPCs$ModuleConformity, 
                      VarExplained = SubsetPCs$VarExplained);
    rm(SubsetColors); rm(SubsetPCs); collect_garbage();
  }
  PCs;
}

# The Network is the same list that is returned by GetModules. It does not contain expression
# data, so those must be given separately. 

NetworkModulePCs = function(ExprData, Network, UniversalModuleColors = NULL, 
                            OnlySet = NULL,
                            AlignPCs="along average", Impute = FALSE, verbose=1, print.level=0)
{
  PCs = ModulePCs(ExprData, SelectedGenes = Network$SelectedGenes, 
                  ModuleColors = Network$Colors, UniversalModuleColors, OnlySet = OnlySet,
                  AlignPCs, Impute = Impute, verbose, print.level);
  
  PCs;
}

#--------------------------------------------------------------------------------------
#
# AddTraitToPCs
#
#--------------------------------------------------------------------------------------

# Adds a trait vector to a set of eigenvectors.
# Caution: Traits is assumed to be a vector of lists with each list having an entry data which is 
# a No.Samples x No.Traits data frame with an appropriate column name, not a vector.

AddTraitToPCs = function(PCs, Traits, verbose=0, print.level=0)
{
  spaces = PrintSpaces(print.level);
  
  No.Sets = length(Traits);
  setsize = CheckSets(Traits);
  No.Traits = setsize$ngenes;
  No.Samples = setsize$nsamples;
  
  if (length(PCs)!=No.Sets)
    stop("Numbers of sets encoded in the length of PCs and Traits parameters differ - must be the same.");
  
  if (verbose>0) print.flush(paste(spaces, "AddTraitToPCs: Adding traits to principal components."));
  
  PCTs = vector(mode="list", length=No.Sets);
  for (set in 1:No.Sets)
  {
    trait.subs = Traits[[set]]$data;
    PCT = cbind(PCs[[set]]$data, trait.subs);
    AET = cbind(PCs[[set]]$AverageExpr, trait.subs);
    names(PCT) = c(names(PCs[[set]]$data), names(trait.subs));
    names(AET) = c(names(PCs[[set]]$AverageExpr), names(trait.subs));
    PCTs[[set]] = list(data=PCT, AverageExpr = AET,
                       ModuleConformity = PCs[[set]]$ModuleConformity, 
                       VarExplained = PCs[[set]]$VarExplained);
  }
  PCTs;
}


#--------------------------------------------------------------------------------------
#
# OrderPCs
#
#--------------------------------------------------------------------------------------
#
# performs hierarchical clustering on PCs and returns the order suitable for plotting.

OrderNames = function(Names, Order)
{
  if (length(Names)!=length(Order))
  {
    cat("OrderNames: Error: Length of names is different from the length of the order vector."); 
    stop();
  }
  OrderedNames = Names;
  for (i in 1:length(Order))
  {
    OrderedNames[i] = Names[Order[i]];
  }
  OrderedNames;
}

OrderPCs = function(PCs, GreyLast = TRUE, GreyName = "PCgrey", OrderBy = 1, Order = NULL, OnlySet = NULL)
{
  if (!is.null(OnlySet)) OrderBy = OnlySet;
  
  if (is.null(Order))
  {
    print.flush(paste("OrderPCs: order not given, clustering given data in subset", OrderBy));
    corPC = cor(PCs[[OrderBy]]$data, use="p")
    disPC = as.dist(1-corPC);
    clust = hclust(disPC, method = "average");
    Order = clust$order;
  } 
  
  if (length(Order)!=dim(PCs[[OrderBy]]$data)[2])
    stop("OrderPCs: given PCs and Order have incompatible dimensions.");
  
  if (GreyLast)
  {
    print.flush("OrderPCs:: Putting grey module last");
    ind.grey = 0;
    PCNames = names(PCs[[OrderBy]]$data);  
    for (i in 1:length(PCNames))
    {
      if (PCNames[Order[i]]==GreyName) {order.grey = i; ind.grey = Order[i]; }
    }
    if (ind.grey==0)
    {
      print(paste("OrderPCs:: Error: The grey ME name", GreyName, 
                  "was not found among the names of the given MEs:", PCNames));
      stop();
    }
    
    if (order.grey<length(Order))
    { 
      for (i in order.grey:(length(Order)-1))
      {
        Order[i] = Order[i+1];
      }
    }
    Order[length(Order)] = ind.grey;
  }  
  No.Sets = length(PCs);
  OrderedPCs = vector(mode="list", length = No.Sets);
  if (is.null(OnlySet))
  {
    CalculatedSubsets = c(1:No.Sets);
  } else {
    CalculatedSubsets = c(OnlySet);
  }
  for (subset in CalculatedSubsets) 
  {
    OrderedData = PCs[[subset]]$data;
    OrderedAE = PCs[[subset]]$AverageExpr;
    for (col in (1:dim(OrderedData)[2]))
    {
      OrderedData[col] = PCs[[subset]]$data[Order[col]];
      OrderedAE[col] = PCs[[subset]]$AverageExpr[Order[col]];
    }
    names(OrderedData) = OrderNames(names(PCs[[subset]]$data), Order);
    names(OrderedAE) = names(OrderedData);
    OrderedPCs[[subset]] = list(data = OrderedData, VarExplained = PCs[[subset]]$VarExplained,
                                ModuleConformity = PCs[[subset]]$ModuleConformity, AverageExpr = OrderedAE, 
                                order = Order);
  }
  OrderedPCs;
}

#---------------------------------------------------------------------------------------------
#
# ConsensusOrderPCs
#
#---------------------------------------------------------------------------------------------
# Orders PCs by the dendrogram of their consensus dissimilarity.

ConsensusOrderPCs = function(PCs, UseAbs = FALSE, OnlySet = NULL, GreyLast = TRUE, GreyName = "MEgrey")
{
  Diss = ConsensusPCDissimilarity(PCs, UseAbs = UseAbs, OnlySet = OnlySet);
  h = hclust(as.dist(Diss));
  
  OrderPCs(PCs, GreyLast = GreyLast, GreyName = GreyName, Order = h$order);
} 

#--------------------------------------------------------------------------------------
#
# CorrelationPreservation
#
#--------------------------------------------------------------------------------------
#
# Given a set of PCs (or OrderedPCs), calculate the preservation values for each module in each pair
# of datasets and return them as a matrix

CorrelationPreservation = function(PCs, Set.Labels, LeaveGreyOut = TRUE)
{
  No.Sets = length(PCs);
  if (No.Sets!=length(Set.Labels)) stop("The lengths of PCs and Set.Labels must equal.");
  if (No.Sets<=1) stop("Something is wrong with argument PCs: its length is 0 or 1");
  Names = names(PCs[[1]]$data);
  if (LeaveGreyOut)
  {
    Use = substring(Names, 3)!="grey";
  } else {
    Use = rep(TRUE, times = length(Names));
  }
  No.Mods = ncol(PCs[[1]]$data[, Use]); 
  CP = matrix(0, nrow = No.Mods, ncol = No.Sets*(No.Sets-1)/2);
  diag(CP) = 1;
  CPInd = 1;
  CPNames = NULL;
  for (i in 1:(No.Sets-1))
    for (j in (i+1):No.Sets)
    {
      corPC1 = cor(PCs[[i]]$data[, Use], use="p");
      corPC2 = cor(PCs[[j]]$data[, Use], use="p");
      d = 1-abs(tanh((corPC1 - corPC2) / (abs(corPC1) + abs(corPC2))^2));
      CP[ ,CPInd] = apply(d, 1, sum)-1;
      CPNames = c(CPNames, paste(Set.Labels[i], "::", Set.Labels[j], collapse = ""));
      CPInd = CPInd + 1;
    }
  CPx = as.data.frame(CP);
  names(CPx) = CPNames;
  rownames(CPx) = Names[Use];
  CPx;
}


#--------------------------------------------------------------------------------------
#
# CorrelationPreservation
#
#--------------------------------------------------------------------------------------
#
# Given a set of PCs (or OrderedPCs), calculate the preservation values for each each pair
# of datasets and return them as a matrix.

SetCorrelationPreservation = function(PCs, Set.Labels, LeaveGreyOut = TRUE, method = "absolute")
{
  m = charmatch(method, c("absolute", "hyperbolic"));
  if (is.na(m))
  {
    stop("Unrecognized method given. Recognized methods are absolute, hyperbolic. ");
  }
  No.Sets = length(PCs);
  if (No.Sets!=length(Set.Labels)) stop("The lengths of PCs and Set.Labels must equal.");
  if (No.Sets<=1) stop("Something is wrong with argument PCs: its length is 0 or 1");
  Names = names(PCs[[1]]$data);
  if (LeaveGreyOut)
  {
    Use = substring(Names, 3)!="grey";
  } else {
    Use = rep(TRUE, times = length(Names));
  }
  No.Mods = ncol(PCs[[1]]$data[, Use]);
  SCP = matrix(0, nrow = No.Sets, ncol = No.Sets);
  diag(SCP) = 0;
  for (i in 1:(No.Sets-1))
    for (j in (i+1):No.Sets)
    {
      corPC1 = cor(PCs[[i]]$data[, Use], use="p");
      corPC2 = cor(PCs[[j]]$data[, Use], use="p");
      if (m==1) {
        d = 1 - abs(corPC1 - corPC2)/2;
      } else {
        d = 1-abs(tanh((corPC1 - corPC2) / (abs(corPC1) + abs(corPC2))^2));
      }
      SCP[i,j] = sum(d[upper.tri(d)])/sum(upper.tri(d));
      SCP[j,i] = SCP[i,j];
    }
  SCPx = as.data.frame(SCP);
  names(SCPx) = Set.Labels;
  rownames(SCPx) = Set.Labels;
  SCPx;
}

#--------------------------------------------------------------------------------------
#
# PlotCorPCs
#
#--------------------------------------------------------------------------------------
# Plots a matrix plot of the PC(T)s. On the diagonal the heatmaps show correlation of PCs in the
# particular subset; off-diagonal are differences in the correlation matrix. 
# Titles is a vector of titles for the diagonal diagrams; the off-diagonal will have no title
# for now.

# Now using the d = tanh( (C1 - C2)/(|C1| + |C2|)^2 ) measure.

#PlotCorPCs = function(PCs, Titles, Powers = c(1,2))
PlotCorPCs = function(PCs, Titles, ColorLabels = FALSE, colors = NULL, IncludeSign = FALSE, 
                      ColoredBarPlot = TRUE, LetterSubPlots = TRUE, Letters = NULL, IncludeGrey = TRUE, 
                      setMargins = TRUE, plotCPMeasure = TRUE, plotMeans = FALSE, CPzlim = c(0,1),
                      printCPVals = FALSE, CPcex = 0.9, PlotDiagAdj = FALSE,
                      ...)
{
  #letters = "abcdefghijklmnopqrstuvwxyz";
  if (is.null(Letters)) Letters = "ABCDEFGHIJKLMNOPQRSTUVWXYZ";
  
  print(PlotDiagAdj);
  
  if (is.null(colors)) 
    if (IncludeSign)
    {
      colors = RedWhiteGreen(50);
    } else {
      colors = heat.colors(30);
    }
  No.Sets = length(PCs);
  cex = par("cex");
  mar = par("mar");
  par(mfrow = c(No.Sets, No.Sets));
  par(cex = cex);
  if (!IncludeGrey)
  {
    for (set in 1:No.Sets)
      PCs[[set]]$data = PCs[[set]]$data[ , substring(names(PCs[[set]]$data),3)!="grey"]
  }
  for (i.row in (1:No.Sets))
  {
    for (i.col in (1:No.Sets))
    {
      letter.ind = (i.row-1) * No.Sets + i.col;
      if (LetterSubPlots) 
      {
        #letter = paste("(", substring(letters, first = letter.ind, last = letter.ind), ")", sep = "");
        letter = paste( substring(Letters, first = letter.ind, last = letter.ind), ".  ", sep = "");
      } else {
        letter = NULL;
      }
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
      if (i.row==i.col)
      {
        corPC = cor(PCs[[i.col]]$data, use="p") 
        if (IncludeSign)
        {
          if (PlotDiagAdj) {
            HeatmapWithTextLabels((1+corPC)/2, names(PCs[[i.col]]$data), names(PCs[[i.col]]$data),
                                  main=paste(letter, Titles[[i.col]]), InvertColors=TRUE, zlim=c(0,1.0),
                                  ColorLabels = ColorLabels, colors = colors, SetMargins = FALSE, ...);
          } else {
            HeatmapWithTextLabels(corPC, names(PCs[[i.col]]$data), names(PCs[[i.col]]$data),
                                  main=paste(letter, Titles[[i.col]]), InvertColors=TRUE, zlim=c(-1,1.0),
                                  ColorLabels = ColorLabels, colors = colors, SetMargins = FALSE, ...);
          }
        } else {
          HeatmapWithTextLabels(abs(corPC), names(PCs[[i.col]]$data), names(PCs[[i.col]]$data),
                                main=paste(letter, Titles[[i.col]]), InvertColors=TRUE, zlim=c(0,1.0),
                                ColorLabels = ColorLabels, colors = colors, SetMargins = FALSE, ...);
        }
      } else
      {
        corPC1 = cor(PCs[[i.col]]$data, use="p");
        corPC2 = cor(PCs[[i.row]]$data, use="p");
        cor.dif = (corPC1 - corPC2)/2;
        d = tanh((corPC1 - corPC2) / (abs(corPC1) + abs(corPC2))^2);
        # d = abs(corPC1 - corPC2) / (abs(corPC1) + abs(corPC2));
        dispd = cor.dif;
        if (plotCPMeasure) dispd[upper.tri(d)] = d[upper.tri(d)];
        if (i.row>i.col)
        {
          if (IncludeSign)
          {
            half = as.integer(length(colors)/2);
            halfColors = colors[1:half+1];
          } else {
            halfColors = colors;
          }
          if (printCPVals) {
            printMtx = matrix(paste(".", as.integer((1-abs(dispd))*100), sep = ""), 
                              nrow = nrow(dispd), ncol = ncol(dispd));
            printMtx[printMtx==".100"] = "1";
          } else { 
            printMtx = NULL; 
          }
          if (sum( (1-abs(dispd)<CPzlim[1]) | (1-abs(dispd)>CPzlim[2]) )>0)
            warning("PlotCorPCs: Correlation preservation data out of zlim range!");
          if (plotCPMeasure) {
            HeatmapWithTextLabels(1-abs(dispd), names(PCs[[i.col]]$data), names(PCs[[i.col]]$data), 
                                  main=paste(letter, "UT: Cor.Pres\nLT: 1-Cor.Diff"), InvertColors=TRUE, 
                                  ColorLabels = ColorLabels, zlim = CPzlim, colors = halfColors,
                                  SetMargins = FALSE, 
                                  NumMatrix = printMtx, cex.Num = CPcex, ...);
          } else {
            HeatmapWithTextLabels(1-abs(dispd), names(PCs[[i.col]]$data), names(PCs[[i.col]]$data), 
                                  main=paste(letter, "Preservation"), InvertColors=TRUE, 
                                  ColorLabels = ColorLabels, zlim = CPzlim, colors = halfColors,
                                  SetMargins = FALSE,  NumMatrix= printMtx, cex.Num = CPcex, ...);
          }
        } else {
          if (plotCPMeasure) {
            dp = 1-abs(d);
            method = "";
          } else {
            dp = 1-abs(cor.dif); 
            method = "Preservation:";
          }
          diag(dp) = 0;
          if (plotMeans) {
            sum_dp = mean(dp[upper.tri(dp)]);
            means = apply(dp, 2, sum)/(ncol(dp)-1);
            BarplotWithTextLabels(means, names(PCs[[i.col]]$data), 
                                  main=paste(letter, "D=", signif(sum_dp,2)), 
                                  ylim=c(0,1),
                                  ColorLabels = ColorLabels, Colored = ColoredBarPlot,
                                  SetMargins = FALSE, ... )
          } else {
            sum_dp = sum(dp[upper.tri(dp)]);
            BarplotWithTextLabels(dp, names(PCs[[i.col]]$data),
                                  main=paste(letter, method, "sum = ", signif(sum_dp,3)), 
                                  ylim=c(0,dim(dp)[[1]]),
                                  ColorLabels = ColorLabels, Colored = ColoredBarPlot, 
                                  SetMargins = FALSE, ... )
          }
        }
      }
    }
  }
}


#--------------------------------------------------------------------------------------
#
# PlotCorPCsAndDendros
#
#--------------------------------------------------------------------------------------
# Plots a matrix plot of the PC(T)s. On the diagonal the heatmaps show correlation of PCs in the
# particular subset; off-diagonal are differences in the correlation matrix. 
# Titles is a vector of titles for the diagonal diagrams; the off-diagonal will have no title
# for now.

# Now using the d = tanh( (C1 - C2)/(|C1| + |C2|)^2 ) measure.

#PlotCorPCs = function(PCs, Titles, Powers = c(1,2))
PlotCorPCsAndDendros = function(PCs, Titles, ColorLabels = FALSE, colors = NULL, IncludeSign = FALSE, 
                                ColoredBarPlot = TRUE, LetterSubPlots = TRUE, Letters = NULL, IncludeGrey = TRUE, 
                                setMargins = TRUE, plotCPMeasure = TRUE, plotMeans = FALSE, CPzlim = c(0,1),
                                printCPVals = FALSE, CPcex = 0.9, plotErrors = FALSE, marDendro = NULL,
                                marHeatmap = NULL, PlotDiagAdj = FALSE, 
                                ...)
{
  #Letters = "abcdefghijklmnopqrstuvwxyz";
  if (is.null(Letters)) Letters = "ABCDEFGHIJKLMNOPQRSTUVWXYZ";
  
  if (is.null(colors)) 
    if (IncludeSign)
    {
      colors = RedWhiteGreen(50);
    } else {
      colors = heat.colors(30);
    }
  No.Sets = length(PCs);
  cex = par("cex");
  mar = par("mar");
  par(mfrow = c(No.Sets+1, No.Sets));
  par(cex = cex);
  if (!IncludeGrey)
  {
    for (set in 1:No.Sets)
      PCs[[set]]$data = PCs[[set]]$data[ , substring(names(PCs[[set]]$data),3)!="grey"]
  }
  letter.ind = 1;
  for (set in 1:No.Sets)
  {
    #par(cex = StandardCex/1.4);
    par(mar = marDendro);
    labels = names(PCs[[set]]$data);
    uselabels = labels[substring(labels,3)!="grey"];
    corPC = cor(PCs[[set]]$data[substring(labels,3)!="grey",
                                substring(labels,3)!="grey"], use="p");
    disPC = as.dist(1-corPC);
    clust = hclust(disPC, method = "average");
    if (LetterSubPlots) {
      main = paste(substring(Letters, letter.ind, letter.ind), ". ", Set.Labels[set], sep="");
    } else {
      main = Set.Labels[set];
    }
    plot(clust, main = main, sub="", xlab="", 
         labels = substring(uselabels, 3), ylab="", ylim=c(0,1));
    letter.ind = letter.ind + 1;
  }
  
  for (i.row in (1:No.Sets))
  {
    for (i.col in (1:No.Sets))
    {
      letter.ind = i.row * No.Sets + i.col;
      if (LetterSubPlots) 
      {
        #letter = paste("(", substring(Letters, first = letter.ind, last = letter.ind), ")", sep = "");
        letter = paste( substring(Letters, first = letter.ind, last = letter.ind), ".  ", sep = "");
      } else {
        letter = NULL;
      }
      par(cex = cex);
      if (setMargins) {
        if (ColorLabels) {
          par(mar = c(1,2,3,4)+0.2);
        } else {
          par(mar = c(6,7,3,5)+0.2);
        }
      } else {
        par(mar = marHeatmap);
      }
      No.Modules = dim(PCs[[i.col]]$data)[2]
      if (i.row==i.col)
      {
        corPC = cor(PCs[[i.col]]$data, use="p") 
        if (IncludeSign)
        {
          if (PlotDiagAdj) {
            HeatmapWithTextLabels((1+corPC)/2, names(PCs[[i.col]]$data), names(PCs[[i.col]]$data),
                                  main=paste(letter, Titles[[i.col]]), InvertColors=TRUE, zlim=c(0,1.0),
                                  ColorLabels = ColorLabels, colors = colors, SetMargins = FALSE, ...);
          } else {
            HeatmapWithTextLabels(corPC, names(PCs[[i.col]]$data), names(PCs[[i.col]]$data),
                                  main=paste(letter, Titles[[i.col]]), InvertColors=TRUE, zlim=c(-1,1.0),
                                  ColorLabels = ColorLabels, colors = colors, SetMargins = FALSE, ...);
          }
        } else {
          HeatmapWithTextLabels(abs(corPC), names(PCs[[i.col]]$data), names(PCs[[i.col]]$data),
                                main=paste(letter, Titles[[i.col]]), InvertColors=TRUE, zlim=c(0,1.0),
                                ColorLabels = ColorLabels, colors = colors, SetMargins = FALSE, ...);
        }
      } else
      {
        corPC1 = cor(PCs[[i.col]]$data, use="p");
        corPC2 = cor(PCs[[i.row]]$data, use="p");
        cor.dif = (corPC1 - corPC2)/2;
        d = tanh((corPC1 - corPC2) / (abs(corPC1) + abs(corPC2))^2);
        # d = abs(corPC1 - corPC2) / (abs(corPC1) + abs(corPC2));
        dispd = cor.dif;
        if (plotCPMeasure) dispd[upper.tri(d)] = d[upper.tri(d)];
        if (i.row>i.col)
        {
          if (IncludeSign)
          {
            half = as.integer(length(colors)/2);
            halfColors = colors[1:half+1];
          } else {
            halfColors = colors;
          }
          if (printCPVals) {
            printMtx = matrix(paste(".", as.integer((1-abs(dispd))*100), sep = ""), 
                              nrow = nrow(dispd), ncol = ncol(dispd));
            printMtx[printMtx==".100"] = "1";
          } else { 
            printMtx = NULL; 
          }
          if (sum( (1-abs(dispd)<CPzlim[1]) | (1-abs(dispd)>CPzlim[2]) )>0)
            warning("PlotCorPCs: Correlation preservation data out of zlim range!");
          if (plotCPMeasure) {
            HeatmapWithTextLabels(1-abs(dispd), names(PCs[[i.col]]$data), names(PCs[[i.col]]$data), 
                                  main=paste(letter, "UT: Cor.Pres\nLT: 1-Cor.Diff"), InvertColors=TRUE, 
                                  ColorLabels = ColorLabels, zlim = CPzlim, colors = halfColors,
                                  SetMargins = FALSE, 
                                  NumMatrix = printMtx, cex.Num = CPcex, ...);
          } else {
            HeatmapWithTextLabels(1-abs(dispd), names(PCs[[i.col]]$data), names(PCs[[i.col]]$data), 
                                  main=paste(letter, "Preservation"), InvertColors=TRUE, 
                                  ColorLabels = ColorLabels, zlim = CPzlim, colors = halfColors,
                                  SetMargins = FALSE,  NumMatrix= printMtx, cex.Num = CPcex, ...);
          }
        } else {
          if (plotCPMeasure) {
            dp = 1-abs(d);
            method = "Cor.Pres.:";
          } else {
            dp = 1-abs(cor.dif); 
            method = "Preservation:";
          }
          diag(dp) = 0;
          if (plotMeans) {
            sum_dp = mean(dp[upper.tri(dp)]);
            means = apply(dp, 2, sum)/(ncol(dp)-1);
            if (plotErrors) {
              Errors = sqrt( (apply(dp^2, 2, sum)/(ncol(dp)-1) - means^2)/(ncol(dp)-2));
            } else {
              Errors = NULL; 
            }
            BarplotWithTextLabels(means, names(PCs[[i.col]]$data), 
                                  main=paste(letter, "D=", signif(sum_dp,2)), 
                                  ylim=c(0,1),
                                  ColorLabels = ColorLabels, Colored = ColoredBarPlot,
                                  SetMargins = FALSE, Errors = Errors, ... )
          } else {
            sum_dp = sum(dp[upper.tri(dp)]);
            BarplotWithTextLabels(dp, names(PCs[[i.col]]$data),
                                  main=paste(letter, method, "sum = ", signif(sum_dp,3)), 
                                  ylim=c(0,dim(dp)[[1]]),
                                  ColorLabels = ColorLabels, Colored = ColoredBarPlot, 
                                  SetMargins = FALSE, ... )
          }
        }
      }
    }
  }
}

#---------------------------------------------------------------------------------------------
#
# PlotCorPCs_1
#
#---------------------------------------------------------------------------------------------
# This function plots the correlations and differences/preservation only, no barplot

PlotCorPCs_1 = function(PCs, Titles, ColorLabels = FALSE, colors = NULL, IncludeSign = FALSE, 
                        ColoredBarPlot = TRUE, ...)
{
  if (is.null(colors)) 
    if (IncludeSign)
    {
      colors = GreenWhiteRed(50);
    } else {
      colors = heat.colors(30);
    }
  No.Sets = length(PCs);
  par(mfrow = c(No.Sets, No.Sets));
  for (i.row in (1:No.Sets))
  {
    for (i.col in (1:No.Sets))
    {
      No.Modules = dim(PCs[[i.col]]$data)[2]
      if (i.row==i.col)
      {
        corPC = cor(PCs[[i.col]]$data, use="p") 
        if (IncludeSign)
        {
          HeatmapWithTextLabels(corPC, names(PCs[[i.col]]$data), names(PCs[[i.col]]$data),
                                main=Titles[[i.col]], InvertColors=TRUE, zlim=c(-1,1.0),
                                ColorLabels = ColorLabels, colors = colors, ...);
        } else {
          HeatmapWithTextLabels(abs(corPC), names(PCs[[i.col]]$data), names(PCs[[i.col]]$data),
                                main=Titles[[i.col]], InvertColors=TRUE, zlim=c(0,1.0),
                                ColorLabels = ColorLabels, colors = colors, ...);
        }
      } else
      {
        corPC1 = cor(PCs[[i.col]]$data, use="p");
        corPC2 = cor(PCs[[i.row]]$data, use="p");
        cor.dif = corPC1 - corPC2;
        d = tanh((corPC1 - corPC2) / (abs(corPC1) + abs(corPC2))^2);
        # d = abs(corPC1 - corPC2) / (abs(corPC1) + abs(corPC2));
        if (IncludeSign)
        {
          half = as.integer(length(colors)/2);
          halfColors = colors[1:half+1];
        } else {
          halfColors = colors;
        }
        if (i.row>i.col)
        {
          HeatmapWithTextLabels(d, names(PCs[[i.col]]$data), names(PCs[[i.col]]$data), 
                                main="Correlation NONPreservation", InvertColors=TRUE, 
                                ColorLabels = ColorLabels, zlim=c(-1,1.0), colors = colors, ...);
        } else {
          HeatmapWithTextLabels(cor.dif, names(PCs[[i.col]]$data), names(PCs[[i.col]]$data), 
                                main="Correlation Difference", InvertColors=TRUE, 
                                ColorLabels = ColorLabels, zlim=c(-1,1.0), colors = colors, ...);
        }
      }
    }
  }
}

#---------------------------------------------------------------------------------------------
#
# kME
#
#---------------------------------------------------------------------------------------------
# This function calculates the ME-based connectivity.

kME = function(SetExprData, PCs, ColorF, gene)
{
  color.ind = as.integer(ColorF[gene]);
  #print(paste(gene, color.ind));
  #if (is.na(color.ind)) stop(paste("color.ind is NA for gene", gene));
  #if (is.null(PCs) | is.null(dim(PCs))) stop(paste("PCs is NULL for gene", gene));
  #
  #if (color.ind>ncol(PCs)) stop(paste("For gene", gene, "have color factor", ColorF[gene], 
  #"whose index", color.ind, "exceeds the dimension", ncol(PCs), 
  #"of PCs."));
  kme = cor(SetExprData[, gene], PCs[, color.ind], use="p");
  kme;
}


#---------------------------------------------------------------------------------------------
#
# ConsensusPCDissimilarity
#
#---------------------------------------------------------------------------------------------
# This function calcualtes a consensus dissimilarity (i.e., correlation) among sets of PCs (more generally,
# any sets of vectors). 
# CAUTION: when not using absolute value, the minimum similarity will favor the large negative values!

ConsensusPCDissimilarity = function(PCs, UseAbs = FALSE, OnlySet = NULL)
{
  MEDiss = vector(mode="list", length = No.Sets);
  if (is.null(OnlySet))
  {
    CalculatedSubsets = c(1:No.Sets);
  } else {
    CalculatedSubsets = c(OnlySet);
  }
  for (set in CalculatedSubsets)
  {
    if (UseAbs)
    {
      diss = 1-abs(cor(PCs[[set]]$data, use="p"));
    } else
    {
      diss = 1-cor(PCs[[set]]$data, use="p");
    }
    MEDiss[[set]] = list(Diss = diss);
  }
  
  if (is.null(OnlySet))
  {
    ConsDiss = (MEDiss[[1]]$Diss)
    if (No.Sets>1) for (set in 2:No.Sets)
      ConsDiss = pmax(ConsDiss, MEDiss[[set]]$Diss);
  } else {
    ConsDiss = MEDiss[[OnlySet]]$Diss;
  }
  
  ConsDiss = as.data.frame(ConsDiss);
  names(ConsDiss) = names(PCs[[1]]$data);
  rownames(ConsDiss) = names(PCs[[1]]$data);
  
  ConsDiss;
}


#---------------------------------------------------------------------------------------------
#
# MergeCloseModules
#
#---------------------------------------------------------------------------------------------
# This function merges modules whose PCs fall on one branch of a hierarchical clustering tree

MergeCloseModules = function(ExprData, Network, SmallModuleColors, CutHeight, OnlySet = NULL, 
                             StandardColors = NULL, 
                             OrderedPCs = NULL, UseAbs = FALSE, IncludeGrey = FALSE, 
                             Relabel = TRUE, Impute = FALSE, 
                             verbose = 1, print.level=0)
{
  spaces = PrintSpaces(print.level);
  if (verbose>0) print.flush(paste(spaces, 
                                   "MergeCloseModules: Merging modules whose distance is less than", CutHeight));
  
  # If ordered PCs were not given, calculate them
  
  if (is.null(OrderedPCs)) 
  {
    PCs = NetworkModulePCs(ExprData, Network, UniversalModuleColors = SmallModuleColors,
                           OnlySet = OnlySet, Impute = Impute,
                           verbose = verbose-1, print.level = print.level+1);
    OrderedPCs = OrderPCs(PCs, GreyLast=TRUE, GreyName = "MEgrey", OnlySet = OnlySet );     
  } else if (nlevels(as.factor(SmallModuleColors))!=dim(OrderedPCs[[1]]$data)[2])
  {
    if (verbose>0) print.flush(paste(spaces, "MergeCloseModules: Number of goven module colors", 
                                     "does not match number of given MEs => recalculating the MEs."))
    PCs = NetworkModulePCs(ExprData, Network, UniversalModuleColors = SmallModuleColors,
                           OnlySet = OnlySet, Impute = Impute,
                           verbose = verbose-1, print.level = print.level+1);
    OrderedPCs = OrderPCs(PCs, GreyLast=TRUE, GreyName = "MEgrey", OnlySet = OnlySet);     
  }
  
  # Cluster the found module eigengenes and merge ones that are too close to one another _in both sets_.
  
  No.Sets = length(ExprData);
  setsize = CheckSets(ExprData);
  No.Genes = setsize$ngenes;
  No.Samples = setsize$nsamples;
  
  MEDiss = vector(mode="list", length = No.Sets);
  if (is.null(OnlySet))
  {
    CalculatedSubsets = c(1:No.Sets);
  } else {
    CalculatedSubsets = c(OnlySet);
  }
  for (set in CalculatedSubsets)
  {
    if (IncludeGrey)
    {
      IndexRange = c(1:(nlevels(as.factor(SmallModuleColors))));
    } else {
      IndexRange = c(1:(nlevels(as.factor(SmallModuleColors))-1));
    }
    if (UseAbs)
    {
      diss = 1-abs(cor(OrderedPCs[[set]]$data[, IndexRange], use="p"));
    } else
    {
      diss = 1-cor(OrderedPCs[[set]]$data[, IndexRange], use="p");
    }
    MEDiss[[set]] = list(Diss = diss);
  }
  
  if (is.null(OnlySet))
  {
    ConsDiss = (MEDiss[[1]]$Diss)
    if (No.Sets>1) for (set in 2:No.Sets)
      ConsDiss = pmax(ConsDiss, MEDiss[[set]]$Diss);
  } else {
    ConsDiss = MEDiss[[OnlySet]]$Diss;
  }
  
  METree = hclust(as.dist(ConsDiss), method = "average");
  METreeBranches = as.factor(ModuleNumber(HierTree = METree, CutHeight = CutHeight, MinSize = 1));
  
  # Analyze the branches: look for the ones that contain more than one original module
  
  MEUniqueBranches = levels(METreeBranches);
  MENo.Branches = nlevels(METreeBranches)
  MENumberOnBranch = rep(0, times = MENo.Branches);
  for (branch in 1:MENo.Branches)
  {
    MENumberOnBranch[branch] = sum(METreeBranches==MEUniqueBranches[branch]);
  }
  
  MergedColors = SmallModuleColors;
  
  # Merge modules on the same branch
  
  for (branch in 1:MENo.Branches) if (MENumberOnBranch[branch]>1)
  {
    if (verbose>3) print.flush(paste(spaces, "   Working on branch", branch, "having", 
                                     MENumberOnBranch[branch], "original modules"));
    ModulesOnThisBranch = names(METreeBranches)[METreeBranches==MEUniqueBranches[branch]];
    ColorsOnThisBranch = substring(ModulesOnThisBranch, 3);
    if (verbose>3) print.flush(paste("    Original colors on this branch:", paste(ColorsOnThisBranch, 
                                                                                  collapse=", ")));
    for (color in 2:length(ColorsOnThisBranch))
      MergedColors[MergedColors==ColorsOnThisBranch[color]] = ColorsOnThisBranch[1];
  }
  
  No.Mods = nlevels(as.factor(MergedColors));
  RawModuleColors = levels(as.factor(MergedColors));
  
  # print(paste("No. of new modules: ", No.Mods));
  # print(paste("Merged module colors:"));
  # print(table(as.factor(MergedColors)));
  
  MergedNewColors = MergedColors;
  if (Relabel)
  { 
    # Relabel the merged colors to the usual order based on the number of genes in each module
    if (is.null(StandardColors))
    {
      StandardColors = c("turquoise","blue","brown","yellow","green","red","black","pink","magenta",
                         "purple","greenyellow","tan","salmon","cyan", "midnightblue", "lightcyan","grey60", 
                         "lightgreen", "lightyellow", "royalblue", "darkred", "darkgreen", "darkturquoise", 
                         "darkgrey", "orange", "darkorange", "white" );
      
    }
    No.GenesInModule = rep(0, No.Mods);
    for (mod in 1:No.Mods) No.GenesInModule[mod] = sum(MergedColors==RawModuleColors[mod]);
    
    SortedRawModuleColors = RawModuleColors[order(-No.GenesInModule)]
    
    # Change the color names to the standard sequence, but leave grey grey (that's why rank in general does
    # not equal color)
    
    if (verbose>3) print(paste(spaces, "   Changing original colors:"));
    rank = 0;
    for (color in 1:length(SortedRawModuleColors)) if (SortedRawModuleColors[color]!="grey")
    {
      rank = rank + 1;
      if (verbose>3) print(paste(spaces, "      ", SortedRawModuleColors[color], 
                                 "to ", StandardColors[rank]));
      MergedNewColors[MergedColors==SortedRawModuleColors[color]] = StandardColors[rank];
    }
  }
  
  list(Colors = MergedNewColors, ClustTree = METree, CutHeight = CutHeight);
}


#======================================================================================================
# ConsensusModules.R
#======================================================================================================

#

#-------------------------------------------------------------------------------------
#
# ConsensusModules
#
#-------------------------------------------------------------------------------------

# This function uses a call to a precompiled function written in C.

SetConsensusModules2 = function(Type = "consensus", 
                                ExprData, Network, PCs = NULL, ConsBranchHeightCut = 0.1, ConsModMinSize = 20, 
                                Impute = FALSE, 
                                verbose = 2, print.level = 0)
{
  if (Type != "consensus" & Type != "majority")
    stop("The Type parameter must be either \'consensus\' or \'majority\'.");
  
  No.Sets = length(ExprData);
  setsize = CheckSets(ExprData);
  No.Samples = setsize$nsamples;
  
  spaces = PrintSpaces(print.level);
  if (is.null(PCs))
  {
    print.flush(paste(spaces, "ConsensusModules: PCs not given, calculating them."));
    PCs = NetworkModulePCs(ExprData, Network, UniversalModuleColors = NULL, Impute = Impute,
                           verbose = verbose-1, print.level = print.level+1);
  }
  
  ModuleCor = vector(mode="list", length = No.Sets);
  No.Modules = vector(length = No.Sets);
  GreyIndex = vector(length = No.Sets);
  
  ModuleIndex = vector(mode="list", length = No.Sets);
  
  for (i in 1:No.Sets) 
  {
    No.Modules[i] = dim(PCs[[i]]$data)[[2]];
    ModuleCor[[i]] = list(cor = data.frame(cor(PCs[[i]]$data, use="p")), GreyIndex = NA, 
                          No.Modules = No.Modules[i]);
    names(ModuleCor[[i]]$cor) = names(PCs[[i]]$data);
    module_names = names(PCs[[i]]$data);
    
    # Find and penalize grey module
    for (j in (1:No.Modules[i])) if (substring(module_names[j], 3)=="grey") 
    {
      GreyIndex[i] = j;
      ModuleCor[[i]]$GreyIndex = j;
      ModuleCor[[i]]$cor[j,] = -1;
      ModuleCor[[i]]$cor[,j] = -1;
    }
    
    ColorF = factor(Network$Colors[,i])
    ModuleIndex[[i]] = list(Ind = as.integer(ColorF));
    rm(ColorF);
  }
  
  No.Genes = sum(Network$SelectedGenes);
  
  GeneModuleDiss = array(0, dim = c(No.Genes, No.Genes));
  GeneModuleDiss1 = array(0, dim = c(No.Genes, No.Genes));
  
  #void GeneModuleDiss(double * GeneDiss, int * NGenes, double * ModuleCor, int * NModules,
  #                    int * ModuleMembership)
  
  wd = getwd();
  setwd("../ConsensusModules");
  if (Computer=="genetics-Windows")
  {
    gmd.lib = "GeneModuleDiss.dll";
  } else {
    gmd.lib = "GeneModuleDiss.so";
  }
  
  if (is.loaded(gmd.lib)) dyn.unload(gmd.lib);
  
  dyn.load(gmd.lib);
  setwd(wd);
  
  if (verbose)
  {
    print.flush(paste(spaces, "Calculating module-based gene dissimilarity"));
  } 
  
  for (set in 1:No.Sets)
  {
    if (verbose>1) print.flush(paste(spaces, "++ Working on set ", set));
    ModuleMembership = as.integer(ModuleIndex[[set]]$Ind);
    ModuleCors = as.matrix(ModuleCor[[set]]$cor);
    GeneModuleDiss1 = .C("GeneModuleDiss", GeneModuleDiss1, as.integer(No.Genes), ModuleCors,
                         as.integer(No.Modules[set]), ModuleMembership)[[1]];
    dim(GeneModuleDiss1) = c(No.Genes, No.Genes);
    if (Type=="majority")
    {
      GeneModuleDiss = GeneModuleDiss + GeneModuleDiss1;
    } else {
      GeneModuleDiss = My_pmax(GeneModuleDiss, GeneModuleDiss1);
    }
    rm(ModuleMembership); rm(ModuleCors);
    collect_garbage();
  }
  
  rm(GeneModuleDiss1);
  rm(ModuleCor); rm(ModuleIndex); rm(GreyIndex); 
  collect_garbage();
  
  GeneModuleDiss = GeneModuleDiss + t(GeneModuleDiss);
  
  # Cluster genes based on the GeneModuleDiss dissimilarity measure
  
  if (verbose>0) 
  {
    print.flush(paste(spaces, "Clustering the gene module-based dissimilarity..."));
  }
  
  Cluster = hclust(as.dist(GeneModuleDiss), method = "average");
  
  if (!is.null(ConsBranchHeightCut))
  {
    ConsensusColor = as.character(modulecolor2(Cluster, h1 = ConsBranchHeightCut, minsize1 = ConsModMinSize));
  } else {
    ConsensusColor = NULL;
  }
  
  list(ClustTree = Cluster, Colors = ConsensusColor, ConsBranchHeightCut = ConsBranchHeightCut,
       ConsModMinSize = ConsModMinSize, Dissimilarity = GeneModuleDiss); 
}

# --------------------------------------------------------------------------------------------------------
#
# Same thing one more time, without the C call of GeneModuleDiss.
#
# --------------------------------------------------------------------------------------------------------

SetConsensusModules = function(Type = "consensus", 
                               ExprData, Network, PCs = NULL, ConsBranchHeightCut = 0.1, ConsModMinSize = 20, 
                               Impute = FALSE, verbose = 2, print.level = 0)
{
  if (Type != "consensus" & Type != "majority")
    stop("The Type parameter must be either \'consensus\' or \'majority\'.");
  
  No.Sets = length(ExprData);
  setsize = CheckSets(ExprData);
  No.Samples = setsize$nsamples;
  
  spaces = PrintSpaces(print.level);
  if (is.null(PCs))
  {
    print.flush(paste(spaces, "ConsensusModules: PCs not given, calculating them."));
    PCs = NetworkModulePCs(ExprData, Network, UniversalModuleColors = NULL, Impute = Impute,
                           verbose = verbose-1, print.level = print.level+1);
  }
  
  ModuleCor = vector(mode="list", length = No.Sets);
  No.Modules = vector(length = No.Sets);
  GreyIndex = vector(length = No.Sets);
  
  ModuleIndex = vector(mode="list", length = No.Sets);
  
  for (i in 1:No.Sets) 
  {
    No.Modules[i] = dim(PCs[[i]]$data)[[2]];
    ModuleCor[[i]] = list(cor = data.frame(cor(PCs[[i]]$data, use="p")), GreyIndex = NA, 
                          No.Modules = No.Modules[i]);
    names(ModuleCor[[i]]$cor) = names(PCs[[i]]$data);
    module_names = names(PCs[[i]]$data);
    
    # Find and penalize grey module
    for (j in (1:No.Modules[i])) if (substring(module_names[j], 3)=="grey") 
    {
      GreyIndex[i] = j;
      ModuleCor[[i]]$GreyIndex = j;
      ModuleCor[[i]]$cor[j,] = -1;
      ModuleCor[[i]]$cor[,j] = -1;
    }
    
    ColorF = factor(Network$Colors[,i])
    ModuleIndex[[i]] = list(Ind = as.integer(ColorF));
    rm(ColorF);
  }
  
  No.Genes = sum(Network$SelectedGenes);
  
  GeneModuleDiss = array(0, dim = c(No.Genes, No.Genes));
  GeneModuleDiss1 = array(0, dim = c(No.Genes, No.Genes));
  
  if (verbose)
  {
    print.flush(paste(spaces, "Calculating module-based gene dissimilarity"));
  } 
  
  Index = c(1:No.Genes)
  for (set in 1:No.Sets)
  {
    if (verbose>1) print.flush(paste(spaces, "++ Working on set ", set));
    ModuleMembership = as.integer(ModuleIndex[[set]]$Ind);
    ModuleCors = as.matrix(ModuleCor[[set]]$cor);
    GeneModuleDiss1 = ModuleCors[ Modulemebership[Index], ModuleMembership[Index] ];
    dim(GeneModuleDiss1) = c(No.Genes, No.Genes);
    if (Type=="majority")
    {
      GeneModuleDiss = GeneModuleDiss + GeneModuleDiss1/No.Sets;
    } else {
      GeneModuleDiss = My_pmax(GeneModuleDiss, GeneModuleDiss1);
    }
    rm(ModuleMembership); rm(ModuleCors);
    collect_garbage();
  }
  
  rm(GeneModuleDiss1);
  rm(ModuleCor); rm(ModuleIndex); rm(GreyIndex); 
  collect_garbage();
  
  # GeneModuleDiss = GeneModuleDiss + t(GeneModuleDiss);
  
  # Cluster genes based on the GeneModuleDiss dissimilarity measure
  
  if (verbose>0) 
  {
    print.flush(paste(spaces, "Clustering the gene module-based dissimilarity..."));
  }
  
  Cluster = hclust(as.dist(GeneModuleDiss), method = "average");
  
  if (!is.null(ConsBranchHeightCut))
  {
    ConsensusColor = as.character(modulecolor2(Cluster, h1 = ConsBranchHeightCut, minsize1 = ConsModMinSize));
  } else {
    ConsensusColor = NULL;
  }
  
  list(ClustTree = Cluster, Colors = ConsensusColor, ConsBranchHeightCut = ConsBranchHeightCut,
       ConsModMinSize = ConsModMinSize, Dissimilarity = GeneModuleDiss); 
}

ConsensusModules = function(ExprData, Network, PCs = NULL, ConsBranchHeightCut = 0.1, ConsModMinSize = 20, 
                            verbose = 2, print.level = 0)
{
  SetConsensusModules("majority", ExprData, Network, PCs, ConsBranchHeightCut, ConsModMinSize, 
                      verbose, print.level);
}
#-------------------------------------------------------------------------------------
#
# AverageModules
#
#-------------------------------------------------------------------------------------

AverageModules = function(Network, ConsBranchHeightCut = NULL, ConsModMinSize = 40, 
                          verbose = 2, print.level = 0)
{
  No.Sets = length(Network$Dissimilarity);
  No.Genes = ncol(Dissimilarity[[1]]$data);
  
  spaces = PrintSpaces(print.level);
  
  AverageDissTOM = matrix(0, nrow = No.Genes, ncol = No.Genes);
  
  for (set in (1:No.Sets))
    AverageDissTOM = AverageDissTOM + Network$Dissimilarity[[set]]$data/No.Sets;
  
  
  # Cluster genes based on the AverageDissTOM dissimilarity measure
  
  if (verbose>0) 
  {
    print.flush(paste(spaces, "Clustering the average gene dissimilarity..."));
  }
  
  Cluster = hclust(as.dist(AverageDissTOM), method = "average");
  
  if (!is.null(ConsBranchHeightCut))
  {
    ConsensusColor = as.character(modulecolor2(Cluster, h1 = ConsBranchHeightCut, minsize1 = ConsModMinSize));
  } else {
    ConsensusColor = NULL;
  }
  
  list(ClustTree = Cluster, Colors = ConsensusColor, ConsBranchHeightCut = ConsBranchHeightCut,
       ConsModMinSize = ConsModMinSize, Dissimilarity = AverageDissTOM); 
}


#-------------------------------------------------------------------------------------
#
# My_pmax
#
#-------------------------------------------------------------------------------------
# Re-implementation of R's pmax, because the latter eats up way too much memory.
# Two versions, one that uses an external function written in C; the other one uses only internal R calls,
# but is slower. Which one is chosen depends on the option UseCpmax ___at the time the function is
# loaded___.

if (UseCpmax)
{
  My_pmax = function(a,b)
  {
    if (Computer=="genetics-Windows")
    {
      pmax.lib = "pmax.dll";
    } else {
      pmax.lib = "pmax.so";
    }
    
    if (!is.loaded(pmax.lib)) 
    {
      path = getwd();
      setwd("../CommonFunctions");
      dyn.load(pmax.lib)
      setwd(path);
    }
    if (length(a)==length(b))
    {
      result = .C("pmax", as.double(a), as.double(b), as.integer(length(a)), DUP=FALSE);
      dim(result[[1]]) = dim(a);
    } else {
      stop(paste("My_pmax: the length of parameters a and b differ:", length(a), length(b)));
    } 
    result[[1]];
  }
} else
{
  My_pmax = function(a,b, batch=200000, verbose = 0)
  {
    dim = dim(a);
    a = as.vector(a); b = as.vector(b);
    if (length(a)==length(b))
    {
      start = 1; stop = min(start+batch-1, length(a));
      while (start <= length(a))
      {
        if (verbose>1) print.flush(paste("  My_pmax:", start, "through", stop, "of", length(a)));
        c = cbind(a[start:stop], b[start:stop])
        a[start:stop] = apply(c, 1, max);
        start = stop+1; stop = min(start+batch-1, length(a));
      }
      dim(a) = dim;
    } else {
      stop(paste("My_pmax: the length of parameters a and b differ:", length(a), length(b)));
    }
    a;
  }
}

#-------------------------------------------------------------------------------------
#
# IntersectModules
#
#-------------------------------------------------------------------------------------
IntersectModules = function(Network, ConsBranchHeightCut = NULL, ConsModMinSize = 40, 
                            verbose = 2, print.level = 0)
{
  No.Sets = length(Network$Dissimilarity);
  No.Genes = sum(Network$SelectedGenes);
  
  spaces = PrintSpaces(print.level);
  if (verbose>0) print.flush(paste(spaces, "IntersectModules: Calculating minimum", 
                                   "gene TOM overlap of given sets"));
  
  spaces = PrintSpaces(print.level);
  
  IntersectDissTOM = matrix(0, nrow = No.Genes, ncol = No.Genes);
  
  #print(paste(memory.size(), memory.size(TRUE)));
  
  #for (set in 1:No.Sets)
  #{
  #print.flush(paste("Checking on set", set, 
  #paste(dim(Network$Dissimilarity[[set]]$data), collapse=",")));
  #}
  
  #print(paste(memory.size(), memory.size(TRUE)));
  for (set in 1:No.Sets)
  {
    #print.flush(paste("Working on set", set, paste(dim(IntersectDissTOM), collapse=",") ,
    #paste(dim(Network$Dissimilarity[[set]]$data), collapse="")));
    IntersectDissTOM = My_pmax(IntersectDissTOM, Network$Dissimilarity[[set]]$data);
    collect_garbage();
  }
  
  #print(paste(memory.size(), memory.size(TRUE)));
  collect_garbage()
  #print(paste(memory.size(), memory.size(TRUE)));
  # Cluster genes based on the IntersectDissTOM dissimilarity measure
  
  if (verbose>1) 
  {
    print.flush(paste(spaces, "IntersectModules: Clustering the maximum gene dissimilarity..."));
  }
  
  Cluster = hclust(as.dist(IntersectDissTOM), method = "average");
  
  if (!is.null(ConsBranchHeightCut))
  {
    ConsensusColor = as.character(modulecolor2(Cluster, h1 = ConsBranchHeightCut, minsize1 = ConsModMinSize));
  } else {
    ConsensusColor = NULL;
  }
  collect_garbage()
  #print(paste(memory.size(), memory.size(TRUE)));
  
  list(ClustTree = Cluster, Colors = ConsensusColor, ConsBranchHeightCut = ConsBranchHeightCut,
       ConsModMinSize = ConsModMinSize, Dissimilarity = IntersectDissTOM); 
}

#-----------------------------------------------------------------------------------------
#
# ConsensusDiagnostics
#
#------------------------------------------------------------------------------------------

# a big WARNING: the order of PCs must be the same as the order of the colors when they are converted to
# a factor. This means PCs must NOT be the OrderedPCs!

ConsensusDiagnostics = function(ExprData, Network, ConsColors, PCs = NULL, Impute = FALSE, 
                                verbose=2, print.level=0)
{
  No.Sets = length(ExprData);
  setsize = CheckSets(ExprData);
  No.Genes = sum(Network$SelectedGenes)
  No.Samples = setsize$nsamples;
  
  
  spaces = PrintSpaces(print.level);
  if (verbose) print.flush(paste(spaces, "Calculating Consensus Module Diagnostics."));
  if (is.null(PCs))
  {
    PCs = NetworkModulePCs(ExprData, Network, UniversalModuleColors = ConsColors, Impute = Impute, 
                           verbose = verbose-1, print.level = print.level+1);
  }
  Modules = factor(ConsColors);
  No.Mods = nlevels(Modules);
  No.AssdGenes = sum(ConsColors!="grey");
  if (verbose>1) print.flush(paste(spaces, "  Calculating IM and ME-based connectivities")); 
  # k_me contains ME-based connectivities for every gene in each set
  k_me = matrix(0, ncol = No.Sets, nrow = No.Genes);
  # k_im contains intra-modular connectivities for every gene in each set
  k_im = matrix(0, ncol = No.Sets, nrow = No.Genes);
  for (set in 1:No.Sets)
  {
    TOM = 1-Network$Dissimilarity[[set]]$data;
    for (gene in 1:No.Genes)
    {
      k_me[gene, set] = kME(ExprData[[set]]$data[, Network$SelectedGenes], PCs[[set]]$data, 
                            Modules, gene);
      k_im[gene, set] = kWithinModule(TOM, ConsColors, gene);
    }
  }
  if (verbose>1) print.flush(paste(spaces, "  Calculating correlations of IM and ME-based connectivities")); 
  k_me_cor = vector(mode="list", length = No.Mods);
  k_im_cor = vector(mode="list", length = No.Mods);
  for (mod in 1:No.Mods)
  {
    color = levels(Modules)[mod];
    kme_c = matrix(0, nrow = No.Sets, ncol = No.Sets);
    kim_c = matrix(0, nrow = No.Sets, ncol = No.Sets);
    kme_p = matrix(0, nrow = No.Sets, ncol = No.Sets);
    kim_p = matrix(0, nrow = No.Sets, ncol = No.Sets);
    for (s1 in 1:No.Sets)
      for (s2 in 1:No.Sets)
      {
        kme1 = k_me[ConsColors==color, s1];
        kme2 = k_me[ConsColors==color, s2];
        kme_c[s1,s2] = cor(kme1, kme2, use="pairwise.complete.obs");
        kme_p[s1,s2] = cor.test(kme1, kme2, use="p", method="s")$p.value;
        kim1 = k_im[ConsColors==color, s1];
        kim2 = k_im[ConsColors==color, s2];
        kim_c[s1,s2] = cor(kim1, kim2, use="p");
        kim_p[s1,s2] = cor.test(kim1, kim2, use="p", method="s")$p.value;
      }
    k_me_cor[[mod]] = list(cor = kme_c, p.value = kme_p);
    k_im_cor[[mod]] = list(cor = kim_c, p.value = kim_p);
  }
  
  list(No.AssdGenes = No.AssdGenes, k_im = k_im, k_me = k_me, k_im_cor = k_im_cor,
       k_me_cor = k_me_cor, No.Mods = No.Mods);
  
}

#-----------------------------------------------------------------------------------------
#
# ProbePermutationTest
#
#------------------------------------------------------------------------------------------

ProbePermutationTest = function(#ExprData, 
  Network, No.Perms, method, ConsTreeCut, ConsModMinSize = 40, 
  No.Plots = 0,
  verbose = 2, print.level = 0)
{
  spaces = PrintSpaces(print.level);
  
  # No.Sets = length(ExprData);
  # setsize = CheckSets(ExprData);
  # No.Genes = setsize$ngenes;
  # No.Samples = setsize$nsamples;
  
  No.Sets = length(Network$Dissimilarity);
  
  # Initialize some useful numbers
  if (verbose)
  {
    print.flush(paste(spaces, "PermutationTest: Will generate", No.Perms, 
                      "permutations and consensus modules via", method, "method."));
    print.flush(paste(spaces, "  Using cuts:", paste(ConsTreeCut, collapse = ",")));
  }
  
  No.Probes = sum(Network$SelectedGenes);
  No.Cuts = length(ConsTreeCut);
  
  plot = 1;
  
  # Prepare some statistics
  
  No.Modules = matrix(NA, nrow = No.Cuts, ncol = No.Perms+1);
  ModuleSizes = vector(mode = "list", length = No.Perms+1);
  
  # Here we go...
  for (iperm in (1:(No.Perms+1)))
  {
    if ((verbose==2 && iperm/100 == as.integer(iperm/100)) || (verbose>2))
      print.flush(paste(spaces, "Working on permutation", iperm));
    # PermExprData = vector(mode="list", length = No.Sets);
    # PermExprData[[1]] = ExprData[[1]];
    
    PermNetwork = Network;
    # Permute expression data and the corresponding dissimilarity in Network
    
    for (set in 1:No.Sets)
    {
      if (iperm!=1)
      {
        perm = sample(No.Probes);
      } else {
        perm = c(1:No.Probes);
      }
      # copy original data
      # PermExprData[[set]] = ExprData[[set]];
      # permute selected samples
      # PermExprData[[set]]$data[, Network$SelectedGenes] 
      # = PermExprData[[set]]$data[,Network$SelectedGenes][, perm];
      # Permute dissimilarities in the Network
      PermNetwork$Dissimilarity[[set]]$data = PermNetwork$Dissimilarity[[set]]$data[perm, perm];
    }
    collect_garbage();
    
    # Get the consensus modules
    
    if (method=="Average")
    {
      Consensus = AverageModules(Network = PermNetwork, 
                                 verbose = verbose-1, print.level = print.level+1);
    } else if (method=="Intersection")
    {
      Consensus = IntersectModules(Network = PermNetwork, 
                                   verbose = verbose-1, print.level = print.level+1);
    } else {
      stop(paste("PermutationTest: The method parameter has an unrecognized value:", method));
    }
    if (plot<=No.Plots) 
    {
      if (plot==1) plot(Consensus$ClustTree, labels=FALSE, main = "Original consensus", sub = "", 
                        xlab = "", ylab = "");
      if (plot>1) plot(Consensus$ClustTree, labels=FALSE, main = paste("Permutated data", plot-1), sub = "", 
                       xlab = "", ylab = "");
      plot = plot + 1;
    }
    
    collect_garbage();
    
    # Count modules and their sizes
    
    PermModuleSizes = vector(mode = "list", length = No.Cuts);
    for (cut in 1:No.Cuts)
    {
      ConsensusColor = ModuleNumber(Consensus$ClustTree, CutHeight = ConsTreeCut[cut], 
                                    MinSize = ConsModMinSize);
      PermModuleSizes[[cut]] = list(Sizes = table(ConsensusColor));
      No.Modules[cut, iperm] = length(PermModuleSizes[[cut]]$Sizes);
      #print.flush(paste("Table of ConsensusColor for cut", ConsTreeCut[cut]));
      #print.flush(table(ConsensusColor));
    }
    rm(PermNetwork);
    collect_garbage();
    collect_garbage();
    ModuleSizes[[iperm]] = list(ForCut = PermModuleSizes);
    if (verbose>2) print.flush(paste(spaces, "Permutation", iperm, ": No.Modules for all cuts =", 
                                     paste(No.Modules[, iperm], collapse = ",")));
  }
  
  list(No.Modules = No.Modules, ModuleSizes = ModuleSizes);
  
}


#-----------------------------------------------------------------------------------------------
# Trait selection based on independence and significance
# Assumes that, just like with the PCs, the Traits have the same columns in each dataset (though the
# sample sets need not be the same).

SelectTraits = function(Traits, BranchCut = 0.25, SelectOnSignificance = FALSE, PCs = NULL, 
                        SignifThres = 0.03, Impute = FALSE, verbose = 1, print.level = 0)
{
  spaces = PrintSpaces(print.level);
  if (verbose>0) print.flush(paste(spaces, "SelectTraits: Selecting from ", dim(Traits[[1]]$data)[2],
                                   "traits."));
  No.Sets = length(Traits);
  TDiss = 1-cor(Traits[[1]]$data, use = "pairwise.complete.obs");
  if (No.Sets>1) for (set in 2:No.Sets)
  {
    TDiss = pmax(TDiss, 1-cor(Traits[[set]]$data, use = "pairwise.complete.obs"));
  }
  h = hclust(as.dist(TDiss), method = "average");
  TMods = ModuleNumber(h, CutHeight = BranchCut, MinSize = 1);
  No.TMods = nlevels(as.factor(TMods));
  SelTraits = vector(mode="list", length = No.Sets);
  for (set in 1:No.Sets)
  {
    TData = Traits[[set]]$data;
    TData[is.na(TData)] = 0;
    TPCs = ModulePrinComps1(TData, as.factor(TMods), Impute = Impute, verbose = 0, 
                            GetConformity = TRUE);
    SelTraits[[set]] = list(data = TPCs$PrinComps);
    for (tmod in 1:No.TMods)
    {
      if (sum(TMods==tmod)>1)
      {
        rnk = order(-TPCs$ModuleConformity[TMods==tmod]);
        SelTraits[[set]]$data[, tmod] = (TData[, TMods==tmod])[, rnk[1]];
        names(SelTraits[[set]]$data)[tmod] = (names(TData)[TMods==tmod])[rnk[1]];
      } else {
        SelTraits[[set]]$data[, tmod] = TData[, TMods==tmod];
        names(SelTraits[[set]]$data)[tmod] = names(TData)[TMods==tmod];
      }
    }
  }
  
  if (verbose>0) print.flush(paste(spaces, "SelectTraits: Clustering led to ", dim(SelTraits[[1]]$data)[2],
                                   "traits."));
  
  if (SelectOnSignificance)
  {
    # Reduce further: calculate cor.tests for each ME with each trait in each set; 
    # keep only traits that have at least one cor.test$p.value below a threshold
    
    if (is.null(PCs)) stop("PCs must be given when SelectOnSignificance is requested.");
    
    No.Mods = dim(PCs[[1]]$data)[2];
    if (is.null(No.Mods)) 
      stop("Given PCs do not appear to have the correct structure (vector of list",
           "with \'data\' component being a matrix whose columns are PC vectors");
    
    No.Traits = dim(SelTraits[[1]]$data)[2];
    
    SelectTrait = rep(FALSE, times = No.Traits);
    
    for (trait in 1:No.Traits)
      for (mod in 1:No.Mods)
      {
        Significant = TRUE;
        for (set in (1:No.Sets))
        {
          ct = cor.test(PCs[[set]]$data[, mod], SelTraits[[set]]$data[, trait]); 
          if (ct$p.value>SignifThres) Significant = FALSE;
        }
        if (Significant) SelectTrait[trait] = TRUE;
      }
    
    for (set in 1:No.Sets)
    {
      SelTraits[[set]]$data = SelTraits[[set]]$data[, SelectTrait];
    }
    
    # print(paste("No. of selected traits expected by chance:", No.Mods * No.Traits * TraitThres));
  }
  
  # Re-cluster the significant traits for diagnostic purposes
  if (sum(SelectTrait)>1)
  {
    TDiss = 1-cor(SelTraits[[1]]$data, use = "pairwise.complete.obs");
    if (No.Sets>1) for (set in 2:No.Sets)
    {
      TDiss = pmax(TDiss, 1-cor(SelTraits[[set]]$data, use = "pairwise.complete.obs"));
    }
    newh = hclust(as.dist(TDiss), method = "average");
  } else {
    newh = NULL;
  }
  
  if (verbose>0) 
  {
    print.flush(paste(spaces, "SelectTraits: Selected", sum(SelectTrait), "traits: "));
    print.flush(paste(spaces, paste(names(SelTraits[[1]]$data), collapse = ", ")));
  }
  
  
  list(No.SelectedTraits = sum(SelectTrait), Traits = SelTraits, ClusterTree = h, NewClusterTree = newh);
}

#======================================================================================================
# ColorHandler.R
#======================================================================================================

# A set of global variables and functions that should help handling color names for some 400+ modules.
# A vector called GlobalStandardColors is defined that holds color names with first few entries 
# being the well-known and -loved colors. The rest is randomly chosen from the color names of R,
# excluding grey colors.

#---------------------------------------------------------------------------------------------------------
#
# GlobalStandardColors 
#
#---------------------------------------------------------------------------------------------------------
# This code forms a vector of color names in which the first entries are given by BaseColors and the rest
# is "randomly" chosen from the rest of R color names that do not contain "grey" nor "gray".

BaseColors = c("turquoise","blue","brown","yellow","green","red","black","pink","magenta",
               "purple","greenyellow","tan","salmon","cyan", "midnightblue", "lightcyan",
               "grey60", "lightgreen", "lightyellow", "royalblue", "darkred", "darkgreen",
               "darkturquoise", "darkgrey",
               "orange", "darkorange", "white", "skyblue", "saddlebrown", "steelblue", 
               "paleturquoise", "violet", "darkolivegreen", "darkmagenta", "white" );

RColors = colors()[-grep("grey", colors())];
RColors = RColors[-grep("gray", RColors)];
InBase = match(BaseColors, RColors);
ExtraColors = RColors[-c(InBase[!is.na(InBase)])];
No.Extras = length(ExtraColors);

# Here is the vector of colors that should be used by all functions:

GlobalStandardColors = c(BaseColors, ExtraColors[rank(sin(13*c(1:No.Extras) +sin(13*c(1:No.Extras))) )] );

rm(BaseColors, RColors, ExtraColors, No.Extras);

#---------------------------------------------------------------------------------------------------------
#
# NormalizeLabels
#
#---------------------------------------------------------------------------------------------------------
# "Normalizes" numerical labels such that the largest group is labeled 1, the next largest 2 etc.
# If KeepZero == TRUE, label zero is preserved.

NormalizeLabels = function(Labels, KeepZero = TRUE)
{
  if (KeepZero)
  {
    NonZero = (Labels!=0);
  }
  else
  {
    NonZero = rep(TRUE, length(Labels));
  }
  f = as.numeric(factor(Labels[NonZero]));
  t = table(Labels[NonZero]);
  # print(t)
  r = rank(-as.vector(t));
  norm_labs = rep(0, times = length(Labels));
  norm_labs[NonZero] = r[f];
  norm_labs;
}

#---------------------------------------------------------------------------------------------------------
#
# ColorsFromLabels
#
#---------------------------------------------------------------------------------------------------------
# This function converts integer numerical labels Labels into color names in the order either given by
# StandardColors,
# or (if StandardColors==NULL) by GlobalStandardColors. If GreyIsZero == TRUE, labels 0 will be assigned
# the color grey; otherwise presence of labels below 1 will trigger an error.

ColorsFromLabels = function(Labels, ZeroIsGrey = TRUE, StandardColors = NULL)
{
  if (is.null(StandardColors)) StandardColors = GlobalStandardColors;
  if (ZeroIsGrey) MinLabel = 0 else MinLabel = 1
  if (sum( (Labels>=MinLabel) & (Labels <= length(StandardColors)) )!= length(Labels))
    stop(paste("Input error: something's wrong with Labels. Either they are not a numeric vector,",
               "or some values are below", MinLabel, 
               "or some values are above the maximum number of colors", length(StandardColors)));
  Colors = rep("grey", length(Labels));
  Colors[Labels!=0] = StandardColors[Labels[Labels!=0]];
  Colors;
}

#---------------------------------------------------------------------------------------------------------
#
# DisplayColors
#
#---------------------------------------------------------------------------------------------------------
# This function plots a barplot with all colors given. If Colors are not given, GlobalStandardColors are
# used, i.e. if you want to see the GlobalStandardColors, just call this function without parameters.

DisplayColors = function(Colors = NULL)
{
  if (is.null(Colors)) Colors = GlobalStandardColors;
  barplot(rep(1, length(Colors)), col = Colors, border = Colors);
}



#---------------------------------------------------------------------------------------------------------
#
# ImputeExprData
#
#---------------------------------------------------------------------------------------------------------

# Impute values for NAs. This handles all the t() and as.matrix() operations.

ImputeExprData = function(ExprData)
{
  imputed = data.frame(scale(t(impute.knn(t(as.matrix(ExprData))))))
  names(imputed) = names(ExprData);
  rownames(imputed) = rownames(ExprData);
  imputed;
}

#---------------------------------------------------------------------------------------------------------
#
# dynamicTreeCut and moduleColor functions
#
#---------------------------------------------------------------------------------------------------------

#------------------------------------------------------------------------------------------
#
# cutreeDynamic
#
#------------------------------------------------------------------------------------------
# Modification(s) by Peter Langfelder: returns numerical labels (not colors) in a vector (not a factor).
# Several rendundant blocks of code removed; duplicate function definitions removed; unused
# functions removed.

.minAttachModuleSize = 100;

if( exists("cutreeDynamicTree") ) rm(cutreeDynamicTree)

cutreeDynamicTree = function(dendro, maxTreeHeight=1, deepSplit=TRUE, minModuleSize=50)
{
  staticCutCluster = .cutTreeStatic(dendro=dendro, heightcutoff=maxTreeHeight, minsize1=minModuleSize)
  
  #get tree height for every singleton
  #node_index   tree_height
  demdroHeiAll= rbind( cbind(dendro$merge[,1], dendro$height), cbind(dendro$merge[,2], dendro$height) )
  
  #singletons will stand at the front of the list
  myorder = order(demdroHeiAll[,1])
  
  #get # of singletons
  no.singletons = length(dendro$order)
  demdroHeiAll.sort = demdroHeiAll[myorder, ]
  demdroHei.sort = demdroHeiAll.sort[c(1:no.singletons), ]
  demdroHei = demdroHei.sort[seq(no.singletons, 1, by=-1), ]
  demdroHei[,1]     = -demdroHei[,1]
  
  # combine with prelimilary cluster-cutoff results
  demdroHei         = cbind(demdroHei, as.integer(staticCutCluster))
  
  # re-order the order based on the dendrogram order dendro$order
  demdroHei.order = demdroHei[dendro$order, ]
  
  static.clupos = .locateCluster(demdroHei.order[, 3])
  
  if (is.null(static.clupos) ){
    module.assign   = rep(0, no.singletons)
    return ( module.assign )
  }
  
  static.no     = dim(static.clupos)[1]
  static.clupos2 =     static.clupos
  static.no2     =     static.no
  
  #split individual cluster if there are sub clusters embedded
  mcycle=1
  while(1==1){
    clupos = NULL
    for (i in c(1:static.no)){
      mydemdroHei.order = demdroHei.order[ c(static.clupos[i,1]:static.clupos[i,2]), ] #index to [1, clusterSize]
      mydemdroHei.order[, 1] = mydemdroHei.order[, 1] - static.clupos[i, 1] + 1
      
      #cat("Cycle ", as.character(mcycle), "cluster (", static.clupos[i,1], static.clupos[i,2], ")\n")
      #cat("i=", as.character(i), "\n")
      
      iclupos = .processIndividualCluster(mydemdroHei.order, 
                                          cminModuleSize       = minModuleSize, 
                                          cminAttachModuleSize = .minAttachModuleSize)
      
      iclupos[,1] = iclupos[,1] + static.clupos[i, 1] -1 #recover the original index
      iclupos[,2] = iclupos[,2] + static.clupos[i, 1] -1
      
      clupos  = rbind(clupos, iclupos) #put in the final output buffer
    }
    
    if(deepSplit==FALSE){
      break
    }
    
    if(dim(clupos)[1] != static.no) {
      static.clupos = clupos
      static.no     = dim(static.clupos)[1]
    }else{
      break
    }
    mcycle = mcycle + 1
    #static.clupos  
  }
  final.cnt = dim(clupos)[1]
  #assign colors for modules
  module.assign = rep(0, no.singletons)
  module.cnt=1
  for (i in c(1:final.cnt ))
  {
    sdx = clupos[i, 1] #module start point
    edx = clupos[i, 2] #module end point
    module.size = edx - sdx +1 
    if(module.size <minModuleSize){
      next
    }
    #assign module lable
    module.assign[sdx:edx] = rep(module.cnt, module.size)
    #update module label for the next module
    module.cnt = module.cnt + 1
  }
  colcode.reduced.order = .assignModuleNumber(module.assign, minsize1=minModuleSize)
  recov.order = order( demdroHei.order[,1])
  colcode.reduced = colcode.reduced.order[recov.order]
  colcode.reduced
}

#leftOrright >0 : running length (with same sign) to right, otherwise to the left
#mysign = -1: negative value, mysign = -1: positive value
.runlengthSign = function(mysequence, leftOrright=-1, mysign=-1){
  seqlen = length(mysequence)   
  if(leftOrright<0){
    pseq = rev(mysequence)
  }else{
    pseq = mysequence
  }
  
  if(mysign<0){ #see where the first POSITIVE number occurs
    nonezero.bool = (pseq > 0)
  }else{ #see where the first NEGATIVE number occur
    nonezero.bool = (pseq < 0)
  }
  if( sum(nonezero.bool) > 0){
    runlength = min( c(1:seqlen)[nonezero.bool] ) - 1
  }else{
    runlength = 0
  }
}

#"0" is for grey module
.assignModuleColor = function(labelpred, minsize1=50, anameallmodules=FALSE, auseblackwhite=FALSE) {
  # here we define modules by using a height cut-off for the branches
  #labelpred= cutree(dendro,h=heightcutoff)
  #cat(labelpred)
  
  #"0", grey module doesn't participate color assignment, directly assigned as "grey"
  labelpredNoZero = labelpred[ labelpred >0 ]
  sort1=-sort(-table(labelpredNoZero))
  # sort1
  modulename= as.numeric(names(sort1))
  modulebranch= sort1 >= minsize1
  no.modules = sum(modulebranch)
  
  colorcode=GlobalStandardColors;
  
  #"grey" means not in any module;
  colorhelp=rep("grey",length(labelpred))
  if ( no.modules==0){
    print("No module detected.")
  }
  else{
    if ( no.modules > length(colorcode)  ){
      print( paste("Too many modules:", as.character(no.modules)) )
    }
    
    if ( (anameallmodules==FALSE) || (no.modules <=length(colorcode)) ){
      labeledModules = min(no.modules, length(colorcode) )
      for (i in c(1:labeledModules)) {
        colorhelp=ifelse(labelpred==modulename[i],colorcode[i],colorhelp)
      }
      colorhelp=factor(colorhelp,levels=c(colorcode[1:labeledModules],"grey"))
    }else{#nameallmodules==TRUE and no.modules >length(colorcode)
      maxcolors=length(colorcode)
      labeledModules = no.modules
      extracolors=NULL
      blackwhite=c("red", "black")
      for(i in c((maxcolors+1):no.modules)){
        if(auseblackwhite==FALSE){
          icolor=paste("module", as.character(i), sep="")
        }else{#use balck white alternatively represent extra colors, for display only
          #here we use the ordered label to avoid put the same color for two neighboring clusters
          icolor=blackwhite[1+(as.integer(modulename[i])%%2) ]
        }
        extracolors=c(extracolors, icolor)
      }
      
      #combine the true-color code and the extra colorcode into a uniform colorcode for 
      #color assignment
      allcolorcode=c(colorcode, extracolors)
      
      for (i in c(1:labeledModules)) {
        colorhelp=ifelse(labelpred==modulename[i],allcolorcode[i],colorhelp)
      }
      colorhelp=factor(colorhelp,levels=c(allcolorcode[1:labeledModules],"grey"))
    }
  }
  
  
  colorhelp
}

# This function written by Peter Langfelder, based on .assignModuleColor above but simplified.
# Assigns module numbers, not colors. All modules are labeled.
#"0" is for grey module
.assignModuleNumber = function(labelpred, minsize1=50) {
  #"0", grey module doesn't participate color assignment, directly assigned as "grey"
  labelpredNoZero = labelpred[ labelpred >0 ]
  sort1=-sort(-table(labelpredNoZero))
  # sort1
  modulename= as.numeric(names(sort1))
  modulebranch= sort1 >= minsize1
  no.modules = sum(modulebranch)
  
  #"grey" means not in any module;
  colorhelp=rep(0,length(labelpred))
  
  for (i in c(1:no.modules)) {
    colorhelp=ifelse(labelpred==modulename[i],i,colorhelp)
  }
  
  colorhelp
}


#locate the start/end positions of each cluster in the ordered cluster label sequence 
#where "-1" indicating no cluster
#3-1 -1 1 1 1 1 2 2 2
#3 3 -1-1 1 1 1 1 2 2 2  (shift)
#---------------------------------
#0-4  0 2 0 0 0 1 0 0 0   (difference)
#       *     * @
.locateCluster = function(clusterlabels)
{
  no.nodes = length(clusterlabels)
  clusterlabels.shift = c(clusterlabels[1], c(clusterlabels[1:(no.nodes-1)]) )
  
  #a non-zero point is the start point of a cluster and it previous point is the end point of the previous
  #cluster 
  label.diff = abs(clusterlabels - clusterlabels.shift)
  
  #process the first and last positions as start/end points if they belong to a cluster instead of no
  # cluster "-1" 
  if(clusterlabels[1]       >0) {label.diff[1]=1} 
  if(clusterlabels[no.nodes]>0) {label.diff[no.nodes]=1} 
  
  flagpoints.bool = label.diff > 0
  if( sum(flagpoints.bool) ==0){
    return(NULL)
  }
  
  flagpoints = c(1:no.nodes)[flagpoints.bool]
  no.points  = length(flagpoints)
  
  myclupos=NULL
  for(i in c(1:(no.points-1)) ){
    idx = flagpoints[i]
    if(clusterlabels[idx]>0){
      myclupos = rbind(myclupos, c(idx, flagpoints[i+1]-1) )
    }
  }
  myclupos
}


#input is the cluster demdrogram of an individual cluster, we want to find its embbeded subclusters
#execution order: mean-height ==> (mean+max)/2 ==> (mean+min)/2
#useMean: =0 ~ use mean-height   as calibation line
#         =1 ~ use (mean+max)/2  as calibation line to detect relatively a small cluster sitting on the head of a bigger one,
#                      so mean-height is too low to detect the two modules.
#         =-1~ use (mean+min)/2  as calibation line to detect relatively a small cluster sitting on the tail of a bigger one,
#                      so mean-height & (mean+max)/2 are too high to detect the two modules

.processIndividualCluster = function(clusterDemdroHei, cminModuleSize=50, cminAttachModuleSize=100, minTailRunlength=12, useMean=0){
  #for debug: use all genes
  #clusterDemdroHei =demdroHei.order
  
  no.cnodes = dim(clusterDemdroHei)[1]
  
  cmaxhei   = max(clusterDemdroHei[, 2])
  cminhei   = min(clusterDemdroHei[, 2])
  
  cmeanhei  = mean(clusterDemdroHei[, 2])
  cmidhei = (cmeanhei + cmaxhei)/2.0
  cdwnhei = (cmeanhei + cminhei)/2.0
  
  if (useMean==1){
    comphei = cmidhei
  }else if (useMean==-1){
    comphei = cdwnhei
  }else{ #normal case
    comphei = cmeanhei
  }
  
  # compute height diffrence with mean height
  heidiff       = clusterDemdroHei[,2] - comphei
  heidiff.shift = .shiftSequence(heidiff, -1)
  
  # get cut positions
  # detect the end point of a cluster, whose height should be less than meanhei 
  #  and the node behind it is the start point of the next cluster which has a height above meanhei
  cuts.bool = (heidiff<0) & (heidiff.shift > 0)
  cuts.bool[1]         = TRUE
  cuts.bool[no.cnodes] = TRUE
  
  if(sum(cuts.bool)==2){
    if (useMean==0){
      new.clupos=.processIndividualCluster(clusterDemdroHei=clusterDemdroHei, cminModuleSize=cminModuleSize, 
                                           cminAttachModuleSize=cminAttachModuleSize,
                                           useMean=1)
    }else if(useMean==1){
      new.clupos=.processIndividualCluster(clusterDemdroHei=clusterDemdroHei, cminModuleSize=cminModuleSize, 
                                           cminAttachModuleSize=cminAttachModuleSize,
                                           useMean=-1)          
    }else{
      new.clupos = rbind(c(1, no.cnodes))
    }
    return (new.clupos)
  }
  
  #a good candidate cluster-end point should have significant # of ahead nodes with head < meanHei
  cutindex =c(1:no.cnodes)[cuts.bool]
  no.cutps = length(cutindex)
  runlens  = rep(999, no.cutps)
  cuts.bool2 = cuts.bool
  for(i in c(2:(no.cutps-1)) ){
    seq = c( (cutindex[i-1]+1):cutindex[i] )
    runlens[i] = .runlengthSign(heidiff[seq], leftOrright=-1, mysign=-1)
    
    if(runlens[i] < minTailRunlength){
      #cat("run length=", runlens[i], "\n")
      cuts.bool2[ cutindex[i] ] = FALSE
    }
  }
  
  #attach SMALL cluster to the left-side BIG cluster if the small one has smaller mean height
  cuts.bool3=cuts.bool2
  if(sum(cuts.bool2) > 3) {
    curj = 2
    while (1==1){
      cutindex2 =c(1:no.cnodes)[cuts.bool2]
      no.clus = length(cutindex2) -1
      if (curj>no.clus){
        break
      }
      pre.sdx = cutindex2[ curj-1 ]+1 #previous module start point
      pre.edx = cutindex2[ curj ] #previous module end   point
      pre.module.size = pre.edx - pre.sdx +1 
      pre.module.hei  = mean(clusterDemdroHei[c(pre.sdx:pre.edx) , 2])
      
      cur.sdx = cutindex2[ curj ]+1 #previous module start point
      cur.edx = cutindex2[ curj+1 ] #previous module end   point
      cur.module.size = cur.edx - cur.sdx +1 
      cur.module.hei  = mean(clusterDemdroHei[c(cur.sdx:cur.edx) , 2])
      
      #merge to the leftside major module, don't change the current index "curj"
      #if( (pre.module.size >minAttachModuleSize)&(cur.module.hei<pre.module.hei)&(cur.module.size<minAttachModuleSize) ){
      if( (cur.module.hei<pre.module.hei)&(cur.module.size<cminAttachModuleSize) ){
        cuts.bool2[ cutindex2[curj] ] = FALSE
      }else{ #consider next cluster
        curj = curj + 1
      }
    }#while
  }#if 
  
  cutindex2 =c(1:no.cnodes)[cuts.bool2]
  no.cutps = length(cutindex2)
  
  #we don't want to lose the small cluster at the tail, attch it to the previous big cluster
  #cat("Lclu= ", cutindex2[no.cutps]-cutindex2[no.cutps-1]+1, "\n")
  if(no.cutps > 2){
    if( (cutindex2[no.cutps] - cutindex2[no.cutps-1]+1) < cminModuleSize ){
      cuts.bool2[ cutindex2[no.cutps-1] ] =FALSE  
    }
  }
  
  cutindex2  = c(1:no.cnodes)[cuts.bool2]
  cutindex2[1]=cutindex2[1]-1 #the first 
  no.cutps2  = length(cutindex2)
  
  if(no.cutps2 > 2){
    new.clupos = cbind( cutindex2[c(1:(no.cutps2-1))]+1, cutindex2[c(2:no.cutps2)] )
  }else{
    new.clupos = cbind( 1, no.cnodes)
  }
  
  if ( dim(new.clupos)[1] == 1 ){   
    if (useMean==0){
      new.clupos=.processIndividualCluster(clusterDemdroHei=clusterDemdroHei, cminModuleSize=cminModuleSize, 
                                           cminAttachModuleSize=cminAttachModuleSize,
                                           useMean=1)
    }else if(useMean==1){
      new.clupos=.processIndividualCluster(clusterDemdroHei=clusterDemdroHei, cminModuleSize=cminModuleSize, 
                                           cminAttachModuleSize=cminAttachModuleSize,
                                           useMean=-1)          
    }   
  }
  new.clupos
}


#delta >0 : shift to right, otherwise to the left
.shiftSequence = function(mysequence, delta){
  seqlen = length(mysequence)
  if(delta>0){
    finalseq=c(mysequence[1:delta], mysequence[1:(seqlen-delta)])
  }else{
    posdelta = -delta 
    finalseq=c(mysequence[(posdelta+1):seqlen], mysequence[(seqlen-posdelta+1):seqlen])
  }
  finalseq
}

#use height cutoff to remove
.cutTreeStatic = function(dendro,heightcutoff=0.99, minsize1=50) {
  
  # here we define modules by using a height cut-off for the branches
  labelpred= cutree(dendro,h=heightcutoff)
  sort1=-sort(-table(labelpred))
  sort1
  modulename= as.numeric(names(sort1))
  modulebranch= sort1 >= minsize1
  no.modules=sum(modulebranch)
  
  colorhelp = rep(-1, length(labelpred) )
  if ( no.modules==0){
    print("No module detected")
  }
  else{
    for (i in c(1:no.modules)) {
      colorhelp=ifelse(labelpred==modulename[i],i ,colorhelp)
    }
  }
  colorhelp
}

if (exists("printFlush")) { remove(printFlush); collect_garbage(); }
printFlush = function(...)
{
  x = print(...)
  if (exists("flush.console")) x=flush.console();
}

if (exists("indentSpaces")) { remove(indentSpaces); collect_garbage(); }
indentSpaces = function(indent = 0)
{
  if (indent>0) 
  {
    spaces = paste(rep("  ", times=indent), collapse="");
  } else
  {
    spaces = "";
  }
  spaces;
}



#------------------------------------------------------------------------------------------------------
#
# Hybrid tree cut
#
#------------------------------------------------------------------------------------------------------

.NewBranch = function(IsBasic = TRUE, IsClosed = FALSE,
                      Content = NULL, MergingHeights = NULL,
                      RootHeight = 1, LastMerge = 0, Size = 0, MergedInto = 0, Singletons = NULL,
                      SingletonHeights = NULL, 
                      AttachHeight = NULL, FailSize = FALSE )
{
  list(IsBasic = IsBasic, IsClosed = IsClosed,
       Content = Content, MergingHeights = MergingHeights,
       RootHeight = RootHeight, LastMerge = LastMerge, Size = Size, MergedInto = MergedInto,
       Singletons = Singletons, SingletonHeights = SingletonHeights,
       AttachHeight = AttachHeight, FailSize = FailSize);
}

# The following are supporting function for GetClusters. 

.CoreSize = function(BranchSize, minClusterSize)
{
  BaseCoreSize = minClusterSize/2 + 1;
  if (BaseCoreSize < BranchSize)
  {
    CoreSize = as.integer(BaseCoreSize + sqrt(BranchSize - BaseCoreSize));
  } else CoreSize = BranchSize;
  CoreSize;
}

# This assumes the diagonal of the distance matrix
# is zero, BranchDist is a square matrix whose dimension is at least 2.

.CoreScatter = function(BranchDist, minClusterSize)
{
  nPoints = dim(BranchDist)[1];
  PointAverageDistances = apply(BranchDist, 2, sum) / (nPoints-1);
  CoreSize = minClusterSize/2 + 1;
  if (CoreSize < nPoints)
  {
    EffCoreSize = as.integer(CoreSize + sqrt(nPoints - CoreSize));
    ord = order(PointAverageDistances);
    Core = ord[c(1:EffCoreSize)];
  } else {
    Core = c(1:nPoints);
    EffCoreSize = nPoints;
  }
  CoreAverageDistances = apply(BranchDist[Core, Core], 2, sum) / (EffCoreSize-1);
  mean(CoreAverageDistances);
}

#-------------------------------------------------------------------------------------------
#
# cutreeHybrid
#
#-------------------------------------------------------------------------------------------
# Traverses a given clustering tree and detects branches whose size is at least minClusterSize, average
# singleton joining height is at most maxCoreScatter and split (attaching height minus average
# height) is at least minGap. If cutHeight is set, all clusters are cut at that height.

# clusterTrim is the fraction of the cluster gap that will be trimmed away.
# Objects whose joining height is above that will be (re-)assigned based on distance to medoids. If
# clusterTrim<=0, all assigments of stage 1 will be respected. 

cutreeHybrid = function(dendro, cutHeight = NULL, minClusterSize = 20, deepSplit = 1,
                        maxCoreScatter = NULL, minGap = NULL, 
                        maxAbsCoreScatter = NULL, minAbsGap = NULL, clusterTrim = 0,
                        labelUnlabeled = TRUE, distM = NULL,
                        useMedoids = FALSE, maxDistToLabel = cutHeight, 
                        respectSmallClusters = TRUE, 
                        verbose = 2, indent = 0)
{
  
  spaces = indentSpaces(indent);
  
  MxBranches = length(dendro$height)
  Branches = vector(mode="list", length = MxBranches);
  
  nBranches = 0;
  
  if (verbose>0) printFlush(paste(spaces, "Detecting clusters..."));
  
  if (is.null(distM))
  {
    if (verbose>0) 
      printFlush(paste(spaces, "..distM not given, will ignore distance information."));
    distM = matrix(0, nrow = length(dendro$height)+1, ncol = length(dendro$height)+1);
  } else {
    diag(distM) = 0;
  }
  
  # No. of merges in the tree
  nMerge = length(dendro$height);
  if (nMerge < 1)
    stop("The given dendrogram is suspicious: number of merges is zero.");
  
  if (is.null(cutHeight))
  {
    if (verbose>0) 
      printFlush(paste(spaces, 
                       "..cutHeight not given, setting it to the maximum height in given dendrogram."));
    cutHeight = max(dendro$height);
    nMergeBelowCut = nMerge;
  } else {
    nMergeBelowCut = sum(dendro$height<=cutHeight);
    if (nMergeBelowCut < 1) stop("cutHeight set too low: no merges below the cut.");
  }
  
  # Default values for maxCoreScatter and minGap:
  
  defMCS = c(0.64, 0.73, 0.82, 0.91);
  defMG = (1-defMCS)*3/4;
  nSplitDefaults = length(defMCS);
  
  # Convert deep split to range 1..4
  if (is.logical(deepSplit)) deepSplit = as.integer(as.integer(deepSplit)*nSplitDefaults + 1);
  if ((deepSplit<0) | (deepSplit>nSplitDefaults-1))
    stop(paste("Parameter deepSplit out of range: allowable range is 0 through", nSplitDefaults-1));
  
  deepSplit = deepSplit + 1;
  
  # If not set, set the cluster gap and core scatter according to deepSplit.
  if (is.null(maxCoreScatter)) maxCoreScatter = defMCS[deepSplit];
  if (is.null(minGap)) minGap = defMG[deepSplit];
  
  # If maxDistToLabel is not set, set it to cutHeight
  
  if (is.null(maxDistToLabel)) maxDistToLabel = cutHeight;
  
  # Convert (relative) minGap and maxCoreScatter to corresponding absolute quantities if the latter were
  # not given.
  
  refQuantile = 0.05;
  refMerge = as.integer(nMerge * refQuantile + 0.5);
  if (refMerge < 1) refMerge = 1;
  refHeight = dendro$height[refMerge]; 
  
  if (is.null(maxAbsCoreScatter))
    maxAbsCoreScatter = refHeight + maxCoreScatter * (cutHeight - refHeight);
  if (is.null(minAbsGap))
    minAbsGap = minGap * (cutHeight - refHeight);
  
  # For each merge, record the cluster that it belongs to
  IndMergeToBranch = rep(0, times = nMerge)
  
  # The root
  RootBranch = 0;
  if (verbose>2) printFlush(paste(spaces, "..Going through the merge tree"));
  
  for (merge in 1:nMerge) if (dendro$height[merge]<=cutHeight)
  {
    # are both merged objects sigletons?
    if (dendro$merge[merge,1]<0 & dendro$merge[merge,2]<0)
    {
      # Yes; start a new cluster.
      nBranches = nBranches + 1;
      Branches[[nBranches]] = 
        .NewBranch(Content = dendro$merge[merge,], MergingHeights = rep(dendro$height[merge], 2), 
                   LastMerge = merge, Size = 2, Singletons = -dendro$merge[merge,],
                   SingletonHeights = rep(dendro$height[merge], 2));
      IndMergeToBranch[merge] = nBranches;
      RootBranch = nBranches;
    } else if (dendro$merge[merge,1] * dendro$merge[merge,2] <0)
    {
      # merge the sigleton into the cluster
      clust = IndMergeToBranch[max(dendro$merge[merge,])];
      if (clust==0) stop("Internal error: a previous merge has no associated cluster. Sorry!");
      gene = min(dendro$merge[merge,]);
      Branches[[clust]]$Content = c(Branches[[clust]]$Content, gene);
      if (Branches[[clust]]$IsBasic) Branches[[clust]]$Singletons = c(Branches[[clust]]$Singletons, -gene);
      Branches[[clust]]$MergingHeights = c(Branches[[clust]]$MergingHeights, dendro$height[merge]);
      Branches[[clust]]$SingletonHeights = c(Branches[[clust]]$SingletonHeights, dendro$height[merge]);
      Branches[[clust]]$LastMerge = merge;
      Branches[[clust]]$Size = Branches[[clust]]$Size + 1;
      IndMergeToBranch[merge] = clust;
      RootBranch = clust;
    } else
    {
      # attempt to merge two clusters:
      clusts = IndMergeToBranch[dendro$merge[merge,]];
      sizes = c(Branches[[clusts[1]]]$Size, Branches[[clusts[2]]]$Size);
      rnk = rank(sizes, ties.method = "first");
      smaller = clusts[rnk[1]]; larger = clusts[rnk[2]];
      if (Branches[[smaller]]$IsBasic)
      {
        coresize = .CoreSize(length(Branches[[smaller]]$Singletons), minClusterSize);
        Core = Branches[[smaller]]$Singletons[c(1:coresize)];
        SmAveDist = mean(apply(distM[Core, Core], 2, sum)/(coresize-1)); 
      } else { SmAveDist = 0; }
      
      if (Branches[[larger]]$IsBasic)
      {
        coresize = .CoreSize(length(Branches[[larger]]$Singletons), minClusterSize);
        Core = Branches[[larger]]$Singletons[c(1:coresize)];
        LgAveDist = mean(apply(distM[Core, Core], 2, sum)/(coresize-1)); 
      } else { LgAveDist = 0; }
      
      # Is the smaller cluster small or shallow enough to be merged?
      SmallerScores = c(Branches[[smaller]]$IsBasic, 
                        (Branches[[smaller]]$Size < minClusterSize),
                        (SmAveDist > maxAbsCoreScatter), 
                        (dendro$height[merge] - SmAveDist < minAbsGap));
      
      if ( SmallerScores[1] * sum(SmallerScores[c(2:4)]) > 0 )
      {
        DoMerge = TRUE;
        SmallerFailSize = !(SmallerScores[3] | SmallerScores[4]);  # Smaller fails only due to size
      } else 
      {
        LargerScores = c(Branches[[larger]]$IsBasic, 
                         (Branches[[larger]]$Size < minClusterSize),
                         (LgAveDist > maxAbsCoreScatter), 
                         (dendro$height[merge] - LgAveDist < minAbsGap));
        if ( LargerScores[1] * sum(LargerScores[c(2:4)]) > 0 )
        { # Actually: the larger one is the one to be merged
          DoMerge = TRUE;
          SmallerFailSize = !(LargerScores[3] | LargerScores[4]);  # cluster fails only due to size
          x = smaller; smaller = larger; larger = x;
        } else {
          DoMerge = FALSE; # None of the two satisfies merging criteria
        }
      }
      if (DoMerge)
      {
        # merge the smaller into the larger cluster and close it.
        Branches[[smaller]]$IsClosed = TRUE;
        Branches[[smaller]]$FailSize = SmallerFailSize;
        Branches[[smaller]]$MergedInto = larger; 
        Branches[[smaller]]$AttachHeight = dendro$height[merge];
        Branches[[larger]]$Content = c(Branches[[larger]]$Content, smaller);
        if (Branches[[larger]]$IsBasic) 
        {
          Branches[[larger]]$Singletons = 
            c(Branches[[larger]]$Singletons, Branches[[smaller]]$Singletons);
          Branches[[larger]]$SingletonHeights = 
            c(Branches[[larger]]$SingletonHeights, Branches[[smaller]]$SingletonHeights);
        }
        Branches[[larger]]$MergingHeights = c(Branches[[larger]]$MergingHeights, dendro$height[merge]);
        Branches[[larger]]$Size = Branches[[larger]]$Size + Branches[[smaller]]$Size;
        Branches[[larger]]$LastMerge = merge;
        IndMergeToBranch[merge] = larger;
        RootBranch = larger;
      } else
      {
        # Close both clusters and start a composite cluster.
        nBranches = nBranches + 1;
        Branches[[smaller]]$IsClosed = TRUE;
        Branches[[larger]]$IsClosed = TRUE;
        Branches[[smaller]]$AttachHeight = dendro$height[merge];
        Branches[[larger]]$AttachHeight = dendro$height[merge];
        # print(paste("  Starting a composite cluster with number", nBranches));
        Branches[[nBranches]] = .NewBranch(IsBasic = FALSE, Content = clusts, 
                                           MergingHeights = rep(dendro$height[merge], 2), Size = sum(sizes),
                                           LastMerge = merge, 
                                           #Singletons = c(Branches[[larger]]$Singletons, Branches[[smaller]]$Singletons)
                                           Singletons = NULL
        );
        IndMergeToBranch[merge] = nBranches;
        RootBranch = nBranches;
      }
    }
  }
  
  if (verbose>2) printFlush(paste(spaces, "..Going through detected branches and marking clusters.."));
  
  nPoints = nMerge+1;
  
  IsBasic = rep(TRUE, times = nBranches);
  IsBranch = rep(FALSE, times = nBranches);
  SmallLabels = rep(0, times = nPoints);
  Trimmed = rep(0, times = nPoints);
  for (clust in 1:nBranches)
  {
    if (is.null(Branches[[clust]]$AttachHeight)) Branches[[clust]]$AttachHeight = cutHeight;
    IsBasic[clust] = Branches[[clust]]$IsBasic;
    if (Branches[[clust]]$IsBasic)
    {
      coresize = .CoreSize(length(Branches[[clust]]$Singletons), minClusterSize);
      Core = Branches[[clust]]$Singletons[c(1:coresize)];
      Branches[[clust]]$Core = Core;
      CoreScatter = mean(apply(distM[Core, Core], 2, sum)/(coresize-1)); 
    } else { CoreScatter = 0; }
    IsBranch[clust] = Branches[[clust]]$IsBasic & (Branches[[clust]]$Size >= minClusterSize) & 
      (CoreScatter < maxAbsCoreScatter) &
      (Branches[[clust]]$AttachHeight - CoreScatter > minAbsGap);
    if (Branches[[clust]]$FailSize) SmallLabels[Branches[[clust]]$Singletons] = clust;
  }
  if (!respectSmallClusters) SmallLabels = rep(0, times = nPoints);
  
  # Trim objects from clusters that are too close to the boundary.
  
  if (clusterTrim>0) for (clust in 1:nBranches) 
    if (Branches[[clust]]$IsBasic & (Branches[[clust]]$Size>2))
    {
      if (is.null(Branches[[clust]]$AttachHeight)) Branches[[clust]]$AttachHeight = cutHeight;
      #if (length(Branches[[clust]]$SingletonHeights)!=length(Branches[[clust]]$Singletons))
      #  stop("Internal error: length of SingletonHeights differs from length of Singletons. Sorry!");
      bottom = min(Branches[[clust]]$MergingHeights);
      top = Branches[[clust]]$AttachHeight;
      # First object in the cluster to be trimmed:
      FirstTrim = match( TRUE, (Branches[[clust]]$MergingHeights - bottom)/(top-bottom) > 1-clusterTrim);
      if (!is.na(FirstTrim))
      {
        FirstCont = Branches[[clust]]$Content[FirstTrim[1]];
        NSingls = Branches[[clust]]$Size;
        if (FirstCont<0)
        {
          FirstSingl = c(1:NSingls)[Branches[[clust]]$Singletons == -FirstCont];
          if (length(FirstSingl)==0) 
          { 
            print(paste("FirstCont:", FirstCont))
            print("Content:");
            print( Branches[[clust]]$Content);
            print(paste("FirstSingl:", FirstSingl))
            print("Singletons:");
            print( Branches[[clust]]$Singletons);
            stop(paste("Internal error: Trimming: First trimmed content",
                       "points to an invalid singleton. Sorry!"));
          }
          
        } else
        {
          FirstSingl = 
            c(1:NSingls)[Branches[[clust]]$Singletons == Branches[[FirstCont]]$Singletons[1]];
          if (length(FirstSingl)==0) stop(paste("Internal error: Trimming: First trimmed content",
                                                "points to an invalid cluster singleton. Sorry!"));
        }
        if ((FirstSingl < NSingls/2) | (FirstSingl < 3))
        {
          #printFlush(paste("Trimming cluster ", clust));
          #printFlush("Merging heights:")
          #print(Branches[[clust]]$MergingHeights);
          #printFlush(paste("FirstTrim:", FirstTrim));
          #printFlush("Trim Condition:");
          #printFlush((Branches[[clust]]$MergingHeights - bottom)/(top-bottom) <= clusterTrim);
          warning(paste("GetClusters: keeping a low proportion of a cluster:", FirstSingl,
                        "out of", Branches[[clust]]$Size, "objects.\n", 
                        "Increasing the proportion of kept objects to preserve a branch."));
          if (verbose>3)
          {
            printFlush(paste("GetClusters: keeping a low proportion of a cluster:", FirstSingl,
                             "out of", Branches[[clust]]$Size, "objects."));
            printFlush(paste(spaces, "  "));
            printFlush(paste(spaces, "  ..SingletonHeights:", 
                             paste(signif(Branches[[clust]]$SingletonHeights,2), collapse=", ")));
            printFlush(paste(spaces, "  ..bottom =", signif(bottom,2), ", top =", signif(top,2)));
          }
          # Make sure we keep a certain minimum of singletons
          FirstSingl = max(FirstSingl, as.integer(NSingls/3)+1, 3);
        }
        if (FirstSingl<=NSingls) 
        {
          Trimmed[Branches[[clust]]$Singletons[c(FirstSingl:NSingls)]] = clust;
          Branches[[clust]]$Singletons = Branches[[clust]]$Singletons[-c(FirstSingl:NSingls)];
          Branches[[clust]]$SingletonHeights = Branches[[clust]]$SingletonHeights[-c(FirstSingl:NSingls)];
          Branches[[clust]]$Size = length(Branches[[clust]]$Singletons);
        }
      }
    }
  
  # Here's where the original AssignLabel starts
  
  if (verbose>2) printFlush(paste(spaces, "..Assigning stage 1 labels.."));
  
  Colors = rep(0, times = nPoints);
  IsCore = rep(0, times = nPoints);
  BranchBranches = c(1:nBranches)[IsBranch];
  color = 0;
  for (clust in BranchBranches)
  {
    color = color+1;
    Colors[Branches[[clust]]$Singletons] = color;
    SmallLabels[Branches[[clust]]$Singletons] = 0;
    coresize = .CoreSize(length(Branches[[clust]]$Singletons), minClusterSize);
    Core = Branches[[clust]]$Singletons[c(1:coresize)];
    IsCore[Core] = color;
  } 
  
  Labeled = c(1:nPoints)[Colors!=0];
  Unlabeled = c(1:nPoints)[Colors==0];
  nUnlabeled = length(Unlabeled);
  UnlabeledExist = (nUnlabeled>0);
  if (length(Labeled)>0)
  {
    LabelFac = factor(Colors[Labeled]);
    nProperLabels = nlevels(LabelFac);
  } else
  {
    nProperLabels = 0;
  }
  if (labelUnlabeled & UnlabeledExist & nProperLabels>0)
  {
    if (verbose>1) printFlush(paste(spaces, "..Assigning stage 2 labels.."));
    # Assign some of the grey genes to the nearest module. Define nearest as the distance to the medoid,
    # that is the point in the cluster that has the lowest average distance to all other points in the
    # cluster. First get the medoids.
    if (is.null(distM)) 
      stop("When requesting assigning points based on distance, distM must be given.");
    diag(distM) = 0;
    if (useMedoids)
    {
      Medoids = rep(0, times = nProperLabels);
      ClusterRadii = rep(0, times = nProperLabels);	
      for (cluster in 1:nProperLabels)
      {
        InCluster = c(1:nPoints)[Colors==cluster];
        DistInCluster = distM[InCluster, InCluster];
        DistSums = apply(DistInCluster, 2, sum);
        Medoids[cluster] = InCluster[which.min(DistSums)];
        ClusterRadii[cluster] = max(DistInCluster[, which.min(DistSums)])
      }
      # If small clusters are to be respected, assign those first based on medoid-medoid distances.
      if (respectSmallClusters)
      {
        FSmallLabels = factor(SmallLabels);
        SmallLabLevs = as.numeric(levels(FSmallLabels));
        nSmallClusters = nlevels(FSmallLabels) - (SmallLabLevs[1]==0);
        if (nSmallClusters>0) for (sclust in SmallLabLevs[SmallLabLevs!=0])
        {
          InCluster = c(1:nPoints)[SmallLabels==sclust];
          # printFlush(paste("SmallCluster", sclust, "has", length(InCluster), "elements."));
          DistInCluster = distM[InCluster, InCluster];
          if (length(InCluster)>1)
          {
            DistSums = apply(DistInCluster, 2, sum);
            smed = InCluster[which.min(DistSums)];
            DistToMeds = distM[Medoids, smed];
            Nearest = which.min(DistToMeds);
            DistToNearest = DistToMeds[Nearest];
            if ( (DistToNearest < ClusterRadii[Nearest]) | (DistToNearest <  maxDistToLabel) )
            {
              Colors[InCluster] = Nearest;
            } else Colors[InCluster] = -1;  # This prevents individual points from being assigned later 
          }
        }
      }
      # Assign leftover unlabeled objects to clusters with nearest medoids
      Unlabeled = c(1:nPoints)[Colors==0];
      UnassdToMedoidDist = distM[Medoids, Unlabeled];
      if (nProperLabels>1)
      {
        NearestMedoids = apply(UnassdToMedoidDist, 2, which.min);
        NearestCenterDist = apply(UnassdToMedoidDist, 2, min);
      } else {
        NearestMedoids = rep(1, times = nUnlabeled);
        NearestCenterDist = UnassdToMedoidDist;
      }
      Colors[Unlabeled] = ifelse((NearestCenterDist < ClusterRadii[NearestMedoids]) | 
                                   (NearestCenterDist < maxDistToLabel) ,
                                 NearestMedoids, 0);
      UnlabeledExist = (sum(Colors==0)>0);
      if (verbose>1)
        printFlush(paste(spaces, "   ...assigned", 
                         sum((NearestCenterDist < ClusterRadii[NearestMedoids]) | 
                               (NearestCenterDist < maxDistToLabel)), 
                         "of", nUnlabeled, "previously unassigned points.")); 
    } else # Instead of medoids, use average distances
    {
      ClusterDiam = rep(0, times = nProperLabels);	
      for (cluster in 1:nProperLabels)
      {
        InCluster = c(1:nPoints)[Colors==cluster];
        nInCluster = length(InCluster)
        DistInCluster = distM[InCluster, InCluster];
        if (nInCluster>1) {
          AveDistInClust = apply(DistInCluster, 2, sum)/(nInCluster-1);
          ClusterDiam[cluster] = max(AveDistInClust);
        } else {
          ClusterDiam[cluster] = 0;
        }
      }
      # If small clusters are respected, assign them first based on average cluster-cluster distances.
      if (respectSmallClusters)
      {
        FSmallLabels = factor(SmallLabels);
        SmallLabLevs = as.numeric(levels(FSmallLabels));
        nSmallClusters = nlevels(FSmallLabels) - (SmallLabLevs[1]==0);
        if (nSmallClusters>0) for (sclust in SmallLabLevs[SmallLabLevs!=0])
        {
          InCluster = c(1:nPoints)[SmallLabels==sclust];
          # printFlush(paste("SmallCluster", sclust, "has", length(InCluster), "elements."));
          if (length(InCluster)>1)
          {
            DistSClustClust = distM[InCluster, Labeled];
            MeanDist = apply(DistSClustClust, 2, mean);
            MeanMeanDist = tapply(MeanDist, LabelFac, mean);
            Nearest = which.min(MeanMeanDist);
            NearestDist = MeanMeanDist[Nearest];
            if ( ((NearestDist < ClusterDiam[Nearest]) | (NearestDist <  maxDistToLabel)) )
            {
              Colors[InCluster] = Nearest;
            } else Colors[InCluster] = -1;  # This prevents individual points from being assigned later
          } 
        }
      }
      # Assign leftover unlabeled objects to clusters with nearest medoids
      Unlabeled = c(1:nPoints)[Colors==0];
      if (length(Unlabeled)>0)
      {
        UnassdToClustDist = apply(distM[Labeled, Unlabeled], 2, tapply, LabelFac, mean);
        if (nProperLabels>1)
        {
          NearestClusters = apply(UnassdToClustDist, 2, which.min);
          NearestClusterDist = apply(UnassdToClustDist, 2, min);
        } else
        {
          NearestClusters = rep(1, length(Unlabeled));
          NearestClusterDist = UnassdToClustDist;
        }
        Colors[Unlabeled] = ifelse((NearestClusterDist < ClusterDiam[NearestClusters]) | 
                                     (NearestClusterDist < maxDistToLabel),
                                   NearestClusters, 0);
        if (verbose>1)
          printFlush(paste(spaces, "   ...assigned", 
                           sum((NearestClusterDist < ClusterDiam[NearestClusters]) | 
                                 (NearestClusterDist < maxDistToLabel)), 
                           "of", nUnlabeled, "previously unassigned points.")); 
      }
    }
  }
  
  # Relabel labels such that 1 corresponds to the largest cluster etc.
  Colors[Colors<0] = 0;
  UnlabeledExist = (sum(Colors==0)>0);
  NumLabs = as.numeric(as.factor(Colors));
  Sizes = table(NumLabs);
  if (UnlabeledExist)
  {
    if (length(Sizes)>1)
    {
      SizeRank = c(1, rank(-Sizes[2:length(Sizes)], ties.method="first")+1);
    } else {
      SizeRank = 1;
    }
    OrdNumLabs = SizeRank[NumLabs];
  } else {
    SizeRank = rank(-Sizes[1:length(Sizes)], ties.method="first");
    OrdNumLabs = SizeRank[NumLabs];
  }
  OrdIsCore = OrdNumLabs-UnlabeledExist;
  OrdIsCore[IsCore==0] = 0; 
  
  list(labels = OrdNumLabs-UnlabeledExist,
       cores = OrdIsCore,
       smallLabels = SmallLabels,
       trimmed = as.numeric(factor(Trimmed))-1,
       branches  = list(nBranches = nBranches, Branches = Branches, 
                        IndMergeToBranch = IndMergeToBranch,
                        RootBranch = RootBranch, IsBasic = IsBasic, IsBranch = IsBranch, 
                        nPoints = nMerge+1));
} 


#----------------------------------------------------------------------------------------------
#
# cutreeDynamic
#
#----------------------------------------------------------------------------------------------
# A wrapper function for cutreeHybrid and cutreeDynamicTree.

cutreeDynamic = function(dendro, cutHeight = NULL, minClusterSize = 20, 
                         method = "hybrid", deepSplit = (ifelse(method=="dynamic", 1, FALSE)), 
                         maxCoreScatter = NULL, minGap = NULL,
                         maxAbsCoreScatter = NULL, minAbsGap = NULL, clusterTrim = 0,  
                         labelUnlabeled = TRUE, distM = NULL, 
                         useMedoids = FALSE, maxDistToLabel = cutHeight,
                         respectSmallClusters = TRUE, 
                         verbose = 2, indent = 0)
{
  
  if (class(dendro)!="hclust") stop("Argument dendro must have class hclust.");
  methods = c("hybrid", "tree");
  met = charmatch(method, methods);
  if (is.na(met))
  {
    stop(paste("Invalid method argument. Accepted values are (unique abbreviations of)", 
               paste(methods, collapse = ", ")));
  } else if (met==1)
  {
    return(cutreeHybrid(dendro = dendro, distM = distM, minClusterSize = minClusterSize, 
                        cutHeight = cutHeight, deepSplit = deepSplit,
                        maxCoreScatter = maxCoreScatter, minGap = minGap,
                        maxAbsCoreScatter = maxAbsCoreScatter, minAbsGap = minAbsGap,
                        labelUnlabeled = labelUnlabeled, useMedoids = useMedoids, 
                        maxDistToLabel = maxDistToLabel, clusterTrim = clusterTrim, 
                        respectSmallClusters = respectSmallClusters, 
                        verbose = verbose, indent = indent)$labels);
  } else
    return(cutreeDynamicTree(dendro = dendro, maxTreeHeight = cutHeight, deepSplit = deepSplit,
                             minModuleSize = minClusterSize));
}


# These functions are writen in the framework where for several sets the expression data are a vector of
# lists, with each list having a component "data" in which the actual expression data for the set are
# stored.

if (exists("plotHclustColors")) rm(plotHclustColors);
plotHclustColors=function(dendro, color, rowLabels=NULL, cex.rowLabels = 0.9, ...) 
{
  options(stringsAsFactors=FALSE);
  if (length(dendro$order) != dim(color)[[1]] ) 
  { 
    stop("ERROR: length of color vector not compatible with no. of objects in the hierarchical tree");
  } else {
    nSets = dim(color)[[2]];
    C = color[dendro$order, ]; 
    step = 1/dim(color)[[1]];
    ystep = 1/nSets;
    barplot(height=1, col = "white", border=F, space=0, axes=F, ...)
    for (j in 1:nSets)
    {
      ind = (1:dim(C)[1]);
      xl = (ind-1) * step; xr = ind * step; 
      yb = rep(ystep*(j-1), dim(C)[1]); yt = rep(ystep * j, dim(C)[1]);
      rect(xl, yb, xr, yt, col = as.character(C[,j]), border = as.character(C[,j]));
      if (is.null(rowLabels))
      {
        text(as.character(j), pos=2, x=0, y=ystep*(j-0.5), cex=cex.rowLabels, xpd = TRUE);
      } else {
        text(rowLabels[j], pos=2, x=0, y=ystep*(j-0.5), cex=cex.rowLabels, xpd = TRUE);
      }
    }
    for (j in 1:nSets) lines(x=c(0,1), y=c(ystep*j,ystep*j));
  }
}

# ===================================================
#The function moduleEigengenes finds the first principal component (eigengene) in each 
# module defined by the colors of the input vector "color".
# It also reports the variances explained by the first 5 principal components.
# And it yields a measure of module conformity for each gene,
# which is highly correlated to the within module connectivity.
# The theoretical underpinnings are described in Horvath, Dong, Yip (2005)
# http://www.genetics.ucla.edu/labs/horvath/ModuleConformity/
# This requires the R library impute
# The output is 
# 1) a data frame of module eigengenes (MEs), 

if(exists("moduleEigengenes")) rm(moduleEigengenes);

moduleEigengenes=function(expr, color, impute = TRUE, nPC = 1, align = "along average",
                          verbose = 0, indent = 0)
{
  spaces = indentSpaces(indent);
  
  if (is.null(expr))
  {  
    stop("moduleEigengenes: Error: expr is NULL. ");
  }
  if (is.null(color))
  {  
    print("moduleEigengenes: Error: color is NULL. ");
    stop()
  }
  alignRecognizedValues =  c("", "along average");
  if (!is.element(align, alignRecognizedValues)) {
    printFlush(paste("ModulePrincipalComponents: Error:",
                     "parameter align has an unrecognised value:", 
                     align, "; Recognized values are ", alignRecognizedValues));
    stop()
  }
  
  maxVarExplained = 5;
  nVarExplained = min(nPC, maxVarExplained);
  modlevels=levels(factor(color))
  PrinComps = data.frame(matrix(NA,nrow=dim(expr)[[1]], ncol= length(modlevels))) 
  averExpr = data.frame(matrix(NA,nrow=dim(expr)[[1]], ncol= length(modlevels))) 
  varexplained= data.frame(matrix(NA, nrow= nVarExplained, ncol= length(modlevels)))
  names(PrinComps)=paste("ME",modlevels,sep="")
  names(averExpr)=paste("AE",modlevels,sep="")
  for(i in c(1:length(modlevels)) )
  {
    if (verbose>0) 
      printFlush(paste(spaces, "moduleEigengenes : Working on ME for module", modlevels[i]));
    modulename = modlevels[i]
    restrict1 = as.character(color)== modulename
    datModule = t(expr[, restrict1])
    if (impute)
    {
      saved.seed = .Random.seed;
      datModule=impute.knn(as.matrix(datModule))
      .Random.seed = saved.seed;
    }
    datModule=t(scale(t(datModule)));
    n = dim(datModule)[1]; p = dim(datModule)[2];
    svd1=svd(datModule, nu = min(n, p, nPC), nv = min(n, p, nPC));
    #mtitle=paste("MEs of ", modulename," module", sep="");
    varexplained[,i]= (svd1$d[1:min(n,p,nVarExplained)])^2/sum(svd1$d^2)
    # this is the first principal component
    PrinComps[,i] = svd1$v[,1]
    if (align == "along average")
    {
      if (verbose>4) printFlush(paste(spaces,
                                      " .. aligning module eigengene with average expression."))
      averExpr[,i] = apply(datModule, 2, mean);  #datModule is transposed...
      if (cor(averExpr[,i], PrinComps[,i])<0) PrinComps[,i] = -PrinComps[,i]
    }
  }
  
  list(eigengenes = PrinComps, averageExpr = averExpr, varExplained = varexplained, nPC = nPC)
}

#-------------------------------------------------------------------------------------
#
#  ModulePrincipalComponents
#
#-------------------------------------------------------------------------------------
# Has been superseded by moduleEigengenes above.

# ===================================================
# This function collects garbage

if (exists("collectGarbage")) rm(collectGarbage);
collectGarbage=function(){while (gc()[2,4] != gc()[2,4] | gc()[1,4] != gc()[1,4]){}}

#--------------------------------------------------------------------------------------
#
# orderMEs
#
#--------------------------------------------------------------------------------------
#
# performs hierarchical clustering on MEs and returns the order suitable for plotting.

orderMEs = function(MEs, greyLast = TRUE, greyName = "MEgrey", orderBy = 1, order = NULL, 
                    useSets = NULL)
{
  if (!is.null(useSets)) 
    if (is.na(match(orderBy, useSets))) orderBy = useSets[1];
  
  if (is.null(order))
  {
    printFlush(paste("orderMEs: order not given, calculating using given set", orderBy));
    corPC = cor(MEs[[orderBy]]$data, use="p")
    disPC = as.dist(1-corPC);
    clust = hclust(disPC, method = "average");
    order = clust$order;
  } 
  
  if (length(order)!=dim(MEs[[orderBy]]$data)[2])
    stop("orderMEs: given MEs and order have incompatible dimensions.");
  
  if (greyLast)
  {
    printFlush("orderMEs:: Putting grey module last");
    ind.grey = 0;
    PCNames = names(MEs[[orderBy]]$data);  
    for (i in 1:length(PCNames))
    {
      if (PCNames[order[i]]==greyName) {order.grey = i; ind.grey = order[i]; }
    }
    if (ind.grey==0)
    {
      print(paste("orderMEs:: Error: The grey ME name", greyName, 
                  "was not found among the names of the given MEs:", PCNames));
      stop();
    }
    
    if (order.grey<length(order))
    { 
      for (i in order.grey:(length(order)-1))
      {
        order[i] = order[i+1];
      }
    }
    order[length(order)] = ind.grey;
  }  
  nSets = length(MEs);
  orderedMEs = MEs;
  if (is.null(useSets)) useSets = c(1:nSets);
  for (subset in useSets) 
  {
    orderedMEs[[subset]]$data = MEs[[subset]]$data[,order];
    names(orderedMEs[[subset]]$data) = names(MEs[[subset]]$data)[order];
    if (!is.null(MEs[[subset]]$averageExpr))
    {
      orderedMEs[[subset]]$averageExpr = MEs[[subset]]$averageExpr[, order]
      names(orderedMEs[[subset]]$averageExpr) = names(MEs[[subset]]$data)[order];
    }
  }
  orderedMEs;
}

#---------------------------------------------------------------------------------------------
#
# consensusOrderMEs
#
#---------------------------------------------------------------------------------------------
# Orders MEs by the dendrogram of their consensus dissimilarity.

consensusOrderMEs = function(MEs, useAbs = FALSE, useSets = NULL, greyLast = TRUE, 
                             greyName = "MEgrey", method = "consensus")
{
  Diss = consensusMEDissimilarity(MEs, useAbs = useAbs, useSets = useSets, method = method);
  h = hclust(as.dist(Diss));
  order = h$order;
  orderMEs(MEs, greyLast = greyLast, greyName = greyName, order = order, useSets = useSets);
} 

#---------------------------------------------------------------------------------------------
#
# consensusMEDissimilarity
#
#---------------------------------------------------------------------------------------------
# This function calcualtes a consensus dissimilarity (i.e., correlation) among sets of MEs (more generally,
# any sets of vectors). 
# CAUTION: when not using absolute value, the minimum similarity will favor the large negative values!

consensusMEDissimilarity = function(MEs, useAbs = FALSE, useSets = NULL, method = "consensus")
{
  methods = c("consensus", "majority");
  m = charmatch(method, methods);
  if (is.na(m))
    stop("Unrecognized method given. Recognized values are", paste(methods, collapse =", "));
  
  nSets = length(MEs);
  MEDiss = vector(mode="list", length = nSets);
  if (is.null(useSets)) useSets = c(1:nSets);
  for (set in useSets)
  {
    if (useAbs)
    {
      diss = 1-abs(cor(MEs[[set]]$data, use="p"));
    } else
    {
      diss = 1-cor(MEs[[set]]$data, use="p");
    }
    MEDiss[[set]] = list(Diss = diss);
  }
  
  for (set in useSets)
    if (set==useSets[1])
    {
      ConsDiss = MEDiss[[set]]$Diss;
    } else {
      if (m==1) {
        ConsDiss = pmax(ConsDiss, MEDiss[[set]]$Diss);
      } else {
        ConsDiss = ConsDiss + MEDiss[[set]]$Diss;
      }
    }
  
  if (m==2) ConsDiss = ConsDiss/nSets;
  
  ConsDiss = as.data.frame(ConsDiss);
  names(ConsDiss) = names(MEs[[useSets[1]]]$data);
  rownames(ConsDiss) = names(MEs[[useSets[1]]]$data);
  
  ConsDiss;
}

#======================================================================================================
# ColorHandler.R
#======================================================================================================

# A set of global variables and functions that should help handling color names for some 400+ modules.
# A vector called .GlobalStandardColors is defined that holds color names with first few entries 
# being the well-known and -loved colors. The rest is randomly chosen from the color names of R,
# excluding grey colors.

#---------------------------------------------------------------------------------------------------------
#
# .GlobalStandardColors 
#
#---------------------------------------------------------------------------------------------------------
# This code forms a vector of color names in which the first entries are given by BaseColors and the rest
# is "randomly" chosen from the rest of R color names that do not contain "grey" nor "gray".

BaseColors = c("turquoise","blue","brown","yellow","green","red","black","pink","magenta",
               "purple","greenyellow","tan","salmon","cyan", "midnightblue", "lightcyan",
               "grey60", "lightgreen", "lightyellow", "royalblue", "darkred", "darkgreen",
               "darkturquoise", "darkgrey",
               "orange", "darkorange", "white", "skyblue", "saddlebrown", "steelblue", 
               "paleturquoise", "violet", "darkolivegreen", "darkmagenta", "white" );

RColors = colors()[-grep("grey", colors())];
RColors = RColors[-grep("gray", RColors)];
InBase = match(BaseColors, RColors);
ExtraColors = RColors[-c(InBase[!is.na(InBase)])];
nExtras = length(ExtraColors);

# Here is the vector of colors that should be used by all functions:

.GlobalStandardColors = c(BaseColors, ExtraColors[rank(sin(13*c(1:nExtras) +sin(13*c(1:nExtras))) )] );

standardColors = function(n = NULL)
{
  if (is.null(n)) return(.GlobalStandardColors);
  if ((n>0) && (n<=length(.GlobalStandardColors))) 
  {
    return(.GlobalStandardColors[c(1:n)]);
  } else {
    stop("Invalid number of standard colors requested.");
  }
}

rm(BaseColors, RColors, ExtraColors, nExtras, InBase);

#---------------------------------------------------------------------------------------------------------
#
# normalizeLabels
#
#---------------------------------------------------------------------------------------------------------
# "Normalizes" numerical labels such that the largest group is labeled 1, the next largest 2 etc.
# If KeepZero == TRUE, label zero is preserved.

normalizeLabels = function(labels, keepZero = TRUE)
{
  if (keepZero)
  {
    NonZero = (labels!=0);
  }
  else
  {
    NonZero = rep(TRUE, length(labels));
  }
  f = as.numeric(factor(labels[NonZero]));
  t = table(labels[NonZero]);
  # print(t)
  r = rank(-as.vector(t));
  norm_labs = rep(0, times = length(labels));
  norm_labs[NonZero] = r[f];
  norm_labs;
}

#---------------------------------------------------------------------------------------------------------
#
# labels2colors
#
#---------------------------------------------------------------------------------------------------------
# This function converts integer numerical labels labels into color names in the order either given by
# colorSeq,
# or (if colorSeq==NULL) by standardColors(). If GreyIsZero == TRUE, labels 0 will be assigned
# the color grey; otherwise presence of labels below 1 will trigger an error.

labels2colors = function(labels, zeroIsGrey = TRUE, colorSeq = NULL)
{
  if (is.null(colorSeq)) colorSeq = standardColors();
  if (zeroIsGrey) minLabel = 0 else minLabel = 1
  if (sum( (labels>=minLabel) & (labels <= length(colorSeq)) )!= length(labels))
    stop(paste("Input error: something's wrong with labels. Either they are not a numeric vector,",
               "or some values are below", minLabel, 
               "or some values are above the maximum number of colors", length(colorSeq)));
  colors = rep("grey", length(labels));
  colors[labels!=0] = colorSeq[labels[labels!=0]];
  colors;
}

#========================================================================================
#
# MergeCloseModules
#
#========================================================================================

#---------------------------------------------------------------------------------
#
# moduleNumber
#
#---------------------------------------------------------------------------------
# Similar to modulecolor2 above, but returns numbers instead of colors, which is oftentimes more useful.
# 0 means unassigned.
# Return value is a simple vector, not a factor.
# Caution: the module numbers are neither sorted nor sequential, the only guarranteed fact is that grey
# probes are labeled by 0 and all probes belonging to the same module have the same number.

moduleNumber = function(dendro, cutHeight = 0.9, minSize = 50)
{
  Branches = cutree(dendro, h = cutHeight);
  NOnBranches = table(Branches);
  TrueBranch = NOnBranches >= minSize;
  Branches[!TrueBranch[Branches]] = 0;
  
  Branches;
}

#--------------------------------------------------------------------------------------
#
# fixDataStructure
#
#--------------------------------------------------------------------------------------
# Check input data: if they are not a vector of lists, put them into the form of a vector of lists.

fixDataStructure = function(data, verbose = 0, indent = 0)
{
  spaces = indentSpaces(indent);
  if ((class(data)!="list") || (class(data[[1]])!="list"))
  {
    if (verbose>0)
      printFlush(paste(spaces, 
                       "fixDataStructure: data is not a vector of lists: converting it into one."));
    x = data;
    data = vector(mode = "list", length = 1);
    data[[1]] = list(data = x);
    rm(x);
  }
  data;
}

#-------------------------------------------------------------------------------------------
#
# checkSets
#
#-------------------------------------------------------------------------------------------
# Checks sets for consistency and returns some diagnostics.

checkSets = function(data, checkStructure = FALSE)
{
  nSets = length(data);
  if (nSets<=0) stop("No data given.");
  structureOK = TRUE;
  if ((class(data)!="list") || (class(data[[1]])!="list"))
  {
    if (checkStructure)
    {
      structureOK = FALSE;
      nGenes = 0; nSamples = 0;
    } else {
      stop("data does not appear to have the correct format. Consider using fixDataStructure",
           "or setting checkStructure = TRUE when calling this function.");
    }
  } else {
    nSamples = vector(length = nSets);
    nGenes = dim(data[[1]]$data)[2];
    for (set in 1:nSets)
    {
      if (nGenes!=dim(data[[set]]$data)[2])
        stop(paste("Incompatible number of genes in set 1 and", set));
      nSamples[set] = dim(data[[set]]$data)[1];
    }
  }
  
  list(nSets = nSets, nGenes = nGenes, nSamples = nSamples, structureOK = structureOK);
}


#--------------------------------------------------------------------------------------
#
# multiSetMEs
#
#--------------------------------------------------------------------------------------

multiSetMEs = function(exprData, colors, universalColors = NULL, useSets = NULL, nPC = 1, 
                       align = "along average", 
                       verbose = 1, indent = 0)
{
  spaces = indentSpaces(indent);
  
  nSets = length(exprData);
  setsize = checkSets(exprData);
  nGenes = setsize$nGenes;
  nSamples = setsize$nSamples;
  
  if (verbose>0) printFlush(paste(spaces,"multiSetMEs: Looking for module MEs."));
  MEs = vector(mode="list", length=nSets);
  if (is.null(useSets)) useSets = c(1:nSets);
  for (set in useSets) {
    if (verbose>0) printFlush(paste(spaces,"  Working on set", as.character(set), "...")); 
    if (is.null(universalColors)) {
      setColors = colors[,set];
    } else {
      setColors = universalColors; 
    }
    setMEs = moduleEigengenes(expr = exprData[[set]]$data,
                              color = setColors, nPC = nPC, align = align, verbose = verbose-1,
                              indent = indent+1);
    MEs[[set]] = list(data = setMEs$eigengenes, averageExpr = setMEs$averageExpr, 
                      varExplained = setMEs$varExplained, nPC = setMEs$nPC);
  }
  MEs;
}

#---------------------------------------------------------------------------------------------
#
# MergeCloseModules
#
#---------------------------------------------------------------------------------------------
# This function merges modules whose MEs fall on one branch of a hierarchical clustering tree

mergeCloseModules = function(exprData, colors, cutHeight = 0.2, MEs = NULL, useAbs = TRUE, 
                             iterate = TRUE,
                             relabel = FALSE, colorSeq = NULL, getNewMEs = TRUE,
                             useSets = NULL, checkDataFormat = TRUE, 
                             verbose = 1, indent=0)
{
  
  MEsInSingleFrame = FALSE;
  spaces = indentSpaces(indent);
  
  if (verbose>0) printFlush(paste(spaces, 
                                  "MergeCloseModules: Merging modules whose distance is less than", cutHeight));
  
  if (!checkSets(exprData, checkStructure = TRUE)$structureOK)
  {
    if (checkDataFormat)
    {
      exprData = fixDataStructure(exprData);
      MEsInSingleFrame = TRUE;
    } else {
      stop("Given exprData appear to be misformatted.");
    }
  }
  
  setsize = checkSets(exprData);
  nSets = setsize$nSets;
  
  if (!is.null(MEs))
  {
    checkMEs = checkSets(MEs, checkStructure = TRUE);
    if (checkMEs$structureOK)
    {
      if (nSets!=checkMEs$nSets)
        stop("Input error: numbers of sets in exprData and MEs differ.")
      for (set in 1:nSets)
      {
        if (checkME$nSamples[set]!=setsize$nSamples[set])
          stop(paste("Number of samples in MEs is incompatible with subset length for set", set));
      }
    } else {
      if (MEsInSingleFrame)
      {
        MEs = fixDataStructure(MEs);
        checkMEs = checkSets(MEs);
      } else {
        stop("MEs do not have the appropriate structure (same as exprData). ");
      }
    }
  }
  
  if (setsize$nGenes!=length(colors))
    stop("Number of genes in exprData is different from the length of original colors. They must equal.");
  
  if ((cutHeight <0) | (cutHeight>(1+as.integer(useAbs)))) 
    stop(paste("Given cutHeight is out of sensible range between 0 and", 1+as.integer(useAbs) ));
  
  done = FALSE; iteration = 1;
  while (!done)
  {
    if (is.null(MEs)) 
    {
      MEs = multiSetMEs(exprData, colors = NULL, universalColors = colors,
                        useSets = useSets, 
                        verbose = verbose-1, indent = indent+1);
      MEs = consensusOrderMEs(MEs, useAbs = useAbs, useSets = useSets, greyLast = TRUE);
      collectGarbage();
    } else if (nlevels(as.factor(colors))!=checkMEs$nGenes)
    {
      if ((iteration==1) & (verbose>0)) printFlush(paste(spaces, "  Number of given module colors", 
                                                         "does not match number of given MEs => recalculating the MEs."))
      MEs = multiSetMEs(exprData, colors = NULL, universalColors = colors,
                        useSets = useSets, 
                        verbose = verbose-1, indent = indent+1);
      MEs = consensusOrderMEs(MEs, useAbs = useAbs, useSets = useSets, greyLast = TRUE);
      collectGarbage();
    }
    if (iteration==1) oldMEs = MEs;
    
    # Cluster the found module eigengenes and merge ones that are too close to one another _in both sets_.
    
    MEDiss = vector(mode="list", length = nSets);
    if (is.null(useSets)) useSets = c(1:nSets);
    for (set in useSets)
    {
      useMEs = c(1:dim(MEs[[set]]$data)[2])[names(MEs[[set]]$data)!="MEgrey"];
      if (useAbs)
      {
        diss = 1-abs(cor(MEs[[set]]$data[, useMEs], use = "p"));
      } else {
        diss = 1-cor(MEs[[set]]$data[, useMEs], use = "p");
      }
      MEDiss[[set]] = list(Diss = diss);
    }
    
    for (set in useSets)
      if (set==useSets[1])
      {
        ConsDiss = MEDiss[[set]]$Diss;
      } else {
        ConsDiss = pmax(ConsDiss, MEDiss[[set]]$Diss);
      }
    
    nOldMods = dim(ConsDiss)[1]; 
    Tree = hclust(as.dist(ConsDiss), method = "average");
    if (iteration==1) oldTree = Tree;
    TreeBranches = as.factor(moduleNumber(dendro = Tree, cutHeight = cutHeight, minSize = 1));
    UniqueBranches = levels(TreeBranches);
    nBranches = nlevels(TreeBranches)
    NumberOnBranch = table(TreeBranches);
    MergedColors = colors;
    
    # Merge modules on the same branch
    
    for (branch in 1:nBranches) if (NumberOnBranch[branch]>1)
    {
      ModulesOnThisBranch = names(TreeBranches)[TreeBranches==UniqueBranches[branch]];
      ColorsOnThisBranch = substring(ModulesOnThisBranch, 3);
      if (verbose>3) printFlush(paste(spaces, "  Merging original colors", paste(ColorsOnThisBranch, 
                                                                                 collapse=", ")));
      for (color in 2:length(ColorsOnThisBranch))
        MergedColors[MergedColors==ColorsOnThisBranch[color]] = ColorsOnThisBranch[1];
    }
    
    nNewMods = nlevels(as.factor(MergedColors));
    
    if (nNewMods<nOldMods)
    {
      RawModuleColors = levels(as.factor(MergedColors));
      if (relabel) 
      {
        # relabel the merged colors to the usual order based on the number of genes in each module
        if (is.null(colorSeq)) colorSeq = standardColors();
        
        nGenesInModule = rep(0, nNewMods);
        for (mod in 1:nNewMods) nGenesInModule[mod] = sum(MergedColors==RawModuleColors[mod]);
        
        SortedRawModuleColors = RawModuleColors[order(-nGenesInModule)]
        
        # Change the color names to the standard sequence, but leave grey grey (that's why rank in general
        # does not equal color)
        MergedNewColors = MergedColors;
        if (verbose>3) print(paste(spaces, "   Changing original colors:"));
        rank = 0;
        for (color in 1:length(SortedRawModuleColors)) if (SortedRawModuleColors[color]!="grey")
        {
          rank = rank + 1;
          if (verbose>3) print(paste(spaces, "      ", SortedRawModuleColors[color], 
                                     "to ", colorSeq[rank]));
          MergedNewColors[MergedColors==SortedRawModuleColors[color]] = colorSeq[rank];
        }
      } else {
        MergedNewColors = MergedColors; 
      }
      if (iterate) 
      {
        colors = MergedNewColors;
      } else {
        done = TRUE;
      }
    } else {
      done = TRUE;
    }
    
    iteration = iteration+1;
    
  } 
  if (getNewMEs)
  {
    if (nNewMods<nOldMods)
    {
      if (verbose>0) printFlush(paste(spaces, "  Calculating new MEs..."));
      NewMEs = multiSetMEs(exprData, colors = NULL, universalColors = MergedNewColors,
                           useSets = useSets, 
                           verbose = verbose-1, indent = indent+1);
      newMEs = consensusOrderMEs(NewMEs, useAbs = useAbs, useSets = useSets, greyLast = TRUE);
    } else {
      newMEs = MEs;
    }
  } else {
    newMEs = NULL;
  }
  
  if (MEsInSingleFrame) 
  {
    newMEs = newMEs[[1]]$data;
    oldMEs = oldMEs[[1]]$data;
  }
  
  list(colors = MergedNewColors, dendro = Tree, oldDendro = oldTree, cutHeight = cutHeight, 
       oldMEs = oldMEs, newMEs = newMEs);
}


setwd(WorkingDirectory);
