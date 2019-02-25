#----------------R Code To Modify in the Classroom Lecture with Students-------------------
#------------------------------------------R API -------------------------------------------
library(wmtsa);library(wavelets);library(wavethresh);library(waveslim);library(wavemulcor)
library(PearsonDS);library(xtable);library(psych);library(adwave);library(biwavelet);library(BootWPTOS)
library(liftLRD);library(mwaved);library(mvLSW);library(unbalhaar)
#---------------------------------------Wavelet Theory-------------------------------------

#---------------------------------------Data----------------------------------------------
pearson.N<-512
pIIIpars <- list(shape=1, location=1, scale=1)
error.pearson.3<-rpearsonIII(pearson.N,params=pIIIpars)
experimental.data.1<-window(error.pearson.3, end=256)

cusps.1 <- function(x) -0.5*abs(x-1)^0.5 -0.5* abs(x-2)^0.5 + 0.5*x + 3.5
cusps.2 <- function(x) -1*abs(x-1)^0.75 -0.5* abs(x-8)^0.75 + 1*x + 3     

x <- seq(0, 10, length=512)
y.1 <- splus2R::signalSeries(cusps.1(x), x)
y.2 <- splus2R::signalSeries(cusps.2(x), x)
y.3<- makeHeaviSine(pearson.N)
#---------------------------------------Transforms---------------------------------------
x.cwt.1 <- wavCWT( as.ts(experimental.data.1), wavelet="gaussian1")
x.cwt.2 <- wavCWT( as.ts(experimental.data.1), wavelet="gaussian2")
x.cwt.3 <- wavCWT( as.ts(experimental.data.1), wavelet="Haar")
x.cwt.4 <- wavCWT( as.ts(experimental.data.1), wavelet="morlet")
#-------------------------------------Trees----------------------------------------------
W.tree.1 <- wavCWTTree(x.cwt.1)
W.tree.2 <- wavCWTTree(x.cwt.2)
W.tree.3 <- wavCWTTree(x.cwt.3)
W.tree.4 <- wavCWTTree(x.cwt.4)
#-----------------------------------Holder Exponent---------------------------------------
holder.1 <- holderSpectrum(W.tree.1)
holder.2 <- holderSpectrum(W.tree.2)
holder.3 <- holderSpectrum(W.tree.3)
holder.4 <- holderSpectrum(W.tree.4)
#------------------------------------Tables----------------------------------------------

distribution.moments.df<-data.frame()
distribution.moments.df<-rbind(c(empMoments(experimental.data.1)),
                               c(empMoments(y.1)),
                               c(empMoments(y.2)))

Table.1<-xtable(distribution.moments.df)

distribution.Holder.df<-data.frame()
distribution.Holder.df<-rbind(c(empMoments(holder.1$exponent)),
                                c(empMoments(holder.2$exponent)),
                                c(empMoments(holder.3$exponent)),
                                c(empMoments(holder.4$exponent)))
colnames(distribution.Holder.df)<-c("mean","variance","skewness","kurtosis")
Table.2<-xtable(distribution.Holder.df)
#-----------------Figures for Classroom Presentation-----------------------------------------------
Figure.1<-plot(W.tree.1$`3`$extrema, type="l", lty=1, 
                 col="black",
                 xlab="Observations",
                 ylab="Extrema")
lines(W.tree.2$`3`$extrema,lty=2, col="green")
lines(W.tree.3$`3`$extrema,lty=3, col="blue")
lines(W.tree.4$`3`$extrema,lty=4, col="red")
rug(W.tree.1$`3`$extrema, side=4, col="black")
rug(W.tree.2$`3`$extrema, side=4, col="green")
rug(W.tree.3$`3`$extrema, side=4, col="blue")
rug(W.tree.4$`3`$extrema, side=4, col="red")
legend("topright",
       c("Branch 5_1=Gaussian 1",
         "Branch 5_2=Gaussian 2",
         "Branch 5_3=Haar",
         "Branch 5_4=Morlet"),
       inset = .01,
       col=c("black","green","blue","red"),
       lwd=2,
       cex=0.8)

op <- par(mfrow = c(3, 3),
          mar=c(1,1,1,1))
Figure.2<-plot(W.tree.1, xlab="(a)")
plot(W.tree.2, xlab="(b)")
plot(W.tree.3, xlab="(c)")
plot(W.tree.4, xlab="(d)")
par(op)

Figure.2<-densityBy(as.data.frame(as.data.frame(y.1)),xlab="CUSP Expression 1", ylab="Simulated",col="green")
Figure.3<-plot(holder.1$exponent, type="l",
                 ylim=c(-4,1),
                 lty=1,
                 xlab="Observations", 
                 ylab="Exponent",
                 col="black")
lines(holder.2$exponent,lty=2, col="green")
lines(holder.3$exponent,lty=3, col="blue")
lines(holder.4$exponent,lty=4, col="red")
rug( holder.1$exponent, side=4, col="black")
rug( holder.2$exponent, side=4, col="green")
rug( holder.3$exponent, side=4, col="blue")
rug( holder.4$exponent, side=4, col="red")
legend("bottomright",
       c("Holder_1=Gaussian 1",
         "Holder_2=Gaussian 2",
         "Holder_3=Haar",
         "Holder_4=Morlet"),
       inset = .01,
       col=c("black","green","blue","red"),
       lwd=2,
       cex=0.8)


Figure.4<-plot(experimental.data.1, type="l",
              lty=1,
              ylim<-c(-1,5),
              xlab="Observations", 
              ylab="value",
              col="black")
lines(y.1,lty=2, col="green")
lines(y.2,lty=3, col="blue")
lines(y.3,lty=4,col="red")
rug(experimental.data.1, side=4, col="black")
rug(y.1, side=4, col="green")
rug(y.2, side=4, col="blue")
rug(y.3, side=4, col="red")
legend("topleft",
       c("Pearson III",
         "CUSP 1",
         "CUSP 2",
         "HeaviSine"),
       inset = .01,
       col=c("black","green","blue","red"),
       lwd=2,
       cex=0.8)
