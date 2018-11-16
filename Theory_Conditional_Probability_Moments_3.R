#---------------------------------------R API -----------------------------------------------------------
library(xtable);library(vars);library(tseries);library(urca);library(fUnitRoots);library(tsDyn)
library(tseriesChaos);library(fNonlinear);library(stats);library(Fractal);library(fracdiff)
library(fSeries);library(stringr);library(stringi);library(corrplot);library(ggplot2)
#---------------------------------------------------Data----------------------------------------------------

Gene.Expression<-read.delim("dataPlusScores_all5.txt")
Gene.Description<-as.array(Gene.Expression[,2])
indx<-grep('UBC',Gene.Description)
indx.2<-grep('UBE2D2',Gene.Description)
Gene.Description[grep('UBC',Gene.Description)]
Gene.List<-strsplit(as.character(Gene.Description)," ")

UBC.Gene.Name<-Gene.List[[indx[2]]][2]
UBC.Gene.Description<-Gene.Description[indx[2]]
UBE2D2.Gene.Name<-Gene.List[[indx.2[2]]][2]
UBE2D2.Gene.Description<-Gene.Description[indx.2[2]]

UBC<-t(Gene.Expression[indx[2],53:100])
UBE2D2<-t(Gene.Expression[indx.2[2],53:100])
#--------------------------------------------------Modeling and Analysis-----------------------------------
Gene.Study<-as.matrix(cbind(UBC,UBE2D2))
colnames(Gene.Study)<-c('UBC','UBE2D2')
Gene.Study.Trans<-as.matrix(cbind(diff(UBC),diff(UBE2D2)))
M <- cor(Gene.Study)
Gene.Network.Causality<-gene.causality(Gene.Study,1)
Gene.Network.Impulses<-gene.impulse.response(Gene.Study)

Gene.Network.Causality.df<-data.frame()
Gene.Network.Causality.df<-rbind(c(Gene.Network.Causality$Granger.hypothesis,
                                 Gene.Network.Causality$Granger.statistic,
                                 Gene.Network.Causality$Granger.pvalue),
                                 c(Gene.Network.Causality$Instant.hypothesis,
                                   Gene.Network.Causality$Instant.statistic,
                                   Gene.Network.Causality$Instant.pvalue))
colnames(Gene.Network.Causality.df)<-c("Hypothesis","F-Test","P-values")
#---------------------------------------------------Tables------------------------------------------------
Table.1<-xtable(Gene.Network.Causality.df)
#---------------------------------Figures for Classroom-------------------------------------------------
Figure.1<-plot(UBC,type="l",lty=1,col=1,ylab="x values")
lines(UBE2D2,lty=2,col=2)
grid()
legend("bottomright", col = 1:2, 
       lty = 1:4, cex=0.75,
       legend = c("UBC","UBE2D2"))
Figure.2<-corrplot(M, order="hclust", addrect=3,col=1:5,t1.cex=0.8)
title(main="Gene Correlation Plot")
Figure.3<-boxplot(Gene.Study)
Figure.4<-plot(Gene.Network.Impulses$irfs)
#---------------------------Function Library for Student Modification------------------------------------------
gene.causality<-function(x,lagOrder)
{
    var.1<-VAR(x, lag.max = lagOrder, type="const")
    roots(var.1)
    causal<-causality(var.1)
    return(list(roots=roots(var.1),
              Granger.hypothesis=causal$Granger$method,
              Granger.statistic=causal$Granger$statistic,
              Granger.pvalue=causal$Granger$p.value,
              Instant.hypothesis=causal$Instant$method,
              Instant.statistic=causal$Instant$statistic,
              Instant.pvalue=causal$Instant$p.value))
}
gene.impulse.response<-function(x)
{
  var.1<-VAR(x, 1, type="const")
  irf.1<-irf(var.1)
  return(list(irfs=irf.1))
}
