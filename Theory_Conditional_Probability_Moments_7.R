library(xtable)
library(bayou)
library(Matrix)
library(BCE)
#---------------------------------------Data--------------------------------------------------------
Measurements<-c(0,0.25,0.5,0.75,1,1.25,1.5,1.75,2)
Species.seq.1<-sample(Measurements,6)
Species.seq.2<-sample(Measurements,6)
Species.seq.3<-sample(Measurements,6)

Sample.seq.1<-sample(Measurements,6)
Sample.seq.2<-sample(Measurements,6)
Sample.seq.3<-sample(Measurements,6)

Species.sequence.group.matrix<-matrix(c(Species.seq.1,
                                          Species.seq.2,
                                          Species.seq.3),nrow=3,ncol=6)

Sample.sequence.group.matrix<-matrix(c(Sample.seq.1,
                                         Sample.seq.2,
                                         Sample.seq.3),nrow=3,ncol=6)

Sample.sequence.group.matrix.df<-as.data.frame(Sample.sequence.group.matrix)
Species.sequence.group.matrix.df<-as.data.frame(Species.sequence.group.matrix)

colnames(Sample.sequence.group.matrix.df)<-c("Metric 1","Metric 2","Metric 3","Metric 4","Metric 5","Metric 6")
colnames(Species.sequence.group.matrix.df)<-c("Metric 1","Metric 2","Metric 3","Metric 4","Metric 5","Metric 6")

rownames(Sample.sequence.group.matrix.df)<-c("Sample 1","Sample 2","Sample 3")
rownames(Species.sequence.group.matrix.df)<-c("Species 1","Species 2","Species 3")

#--------------------------------------Modeling--------------------------------------------------

model.1 <- BCE(Species.sequence.group.matrix.df,
               Sample.sequence.group.matrix.df,
                 relsdRat=.2,relsdDat=.2,
                 iter=1000,outputlength=5000,
                 jmpX=.01,
                 jmpRat=.01)

model.2 <- BCE(Species.sequence.group.matrix.df,
               Sample.sequence.group.matrix.df,
                 relsdRat=.2,relsdDat=.2,
                 iter=1000,outputlength=5000,
                 jmpX=.1,
                 jmpRat=.002)

model.3 <- BCE(Species.sequence.group.matrix.df,
               Sample.sequence.group.matrix.df,
                relsdRat=.2,relsdDat=.2,
                iter=1000,outputlength=5000,
                jmpX=.02,
                jmpRat=.2*(.2*Species.sequence.group.matrix.df))

model.4 <- BCE(Species.sequence.group.matrix.df,
               Sample.sequence.group.matrix.df,
               relsdRat=.2,relsdDat=.2,
               iter=1000,outputlength=5000,
               jmpX=.02,
               jmpRat=.8*(.2*Species.sequence.group.matrix.df))

model.5<-tlsce(Species.sequence.group.matrix.df,
               Sample.sequence.group.matrix.df)

model.6<-bce1(Species.sequence.group.matrix.df,
              Sample.sequence.group.matrix.df,niter=500,
              initX=matrix(1/ncol(Species.sequence.group.matrix.df),
                           ncol(Species.sequence.group.matrix.df),
                           ncol(Sample.sequence.group.matrix.df)),jmpA=.01,jmpX=.01)

Summary.model.1 <-summary(model.1)
Summary.model.2 <-summary(model.2)
Summary.model.3 <-summary(model.3)
Summary.model.4 <-summary(model.4)
Summary.model.5 <-summary(model.5)
Summary.model.6 <-summary(model.6)

#----------------------------------Analysis------------------------------------------------------------

Analysis.1.df<-data.frame()
Analysis.1.df<-rbind(c("")
                     )
colnames(Analysis.1.df)<-c("")
rownames(Analysis.1.df)<-c("")


#---------------------------------------Tables------------------------------------------------------

Table.1<-xtable(Analysis.1.df)

# --------------------------------------Figures----------------------------------------------------

Figure.1
palette(rainbow(12, s = 0.6, v = 0.75))
mp     <- apply(Species.sequence.group.matrix.df,MARGIN=2,max)
mp2    <- apply(Sample.sequence.group.matrix.df,MARGIN=2,max)
pstars <- rbind(t(t(Species.sequence.group.matrix.df)/mp),t(t(Sample.sequence.group.matrix.df)/mp2))
stars(pstars, len = 0.9, key.loc = c(7.2, -2),scale=FALSE,
      ncol=5,ylim=c(0,3),main = "Input: species + samples", 
      draw.segments = TRUE, flip.labels=FALSE)

Figure.2<-plot(model.1)
Figure.3<-plot(model.2)
Figure.4<-plot(model.3)
Figure.5<-plot(model.4)

Figure.5A<-plot(model.6)
Figure.6
xlim <- range(c(Summary.model.1$lbX,Summary.model.1$ubX))
dotchart(x=t(Summary.model.1$meanX),xlim=xlim,main="Taxonomic composition",sub="using bce",pch=16)
nr <- nrow(Summary.model.1$meanX)
nc <- ncol(Summary.model.1$meanX)
for (i in 1:nr) 
{ip <-(nr-i)*(nc+2)+1
cc <- ip : (ip+nc-1)
segments(t(Summary.model.1$lbX[i,]),cc,t(Summary.model.1$ubX[i,]),cc)
}

par(mfrow = c(2,2), mar = c(0.5,1,0.5,0.5))
Figure.7<-pairs(model.1,sample=1,main="Sample 1")
Figure.8<-pairs(model.2,sample=1,main="Sample 1")
Figure.9<-pairs(model.3,sample=1,main="Sample 1")
Figure.10<-pairs(model.4,sample=1,main="Sample 1")
Figure.11<-pairs(model.6,sample=1,main="Sample 1")
#--------------------------------------References--------------------------------------------------

Reference.1<-c("Van den Meersche, K., K. Soetaert and J.J. Middelburg (2008)",
               "A Bayesian compositional estimator for microbial taxonomy based on biomarkers",
               "Limnology and Oceanography Methods 6, 190-199")


#-------------------------------------Function Library---------------------------------------------


