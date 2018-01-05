library(wmtsa)
library(wavelets)
library(wavethresh)
library(waveslim)
library(wavemulcor)
library(PearsonDS)
library(xtable)
library(psych)
library(adwave)
library(biwavelet)


#---------------------------------------Wavelet Theory-------------------------------------

#---------------------------------------Data----------------------------------------------
pearson.N<-512
pIIIpars.1 <- list(shape=1, location=1, scale=1)
pIIIpars.2 <- list(shape=1, location=2, scale=1)
pIIIpars.3 <- list(shape=1, location=3, scale=1)

error.pearson.3.1<-rpearsonIII(pearson.N,params=pIIIpars.1)
error.pearson.3.2<-rpearsonIII(pearson.N,params=pIIIpars.2)
error.pearson.3.3<-rpearsonIII(pearson.N,params=pIIIpars.3)

experimental.data.1<-window(error.pearson.3.1, end=256)
experimental.data.2<-window(error.pearson.3.2, end=256)
experimental.data.3<-window(error.pearson.3.3, end=256)

experimental.data.1<-cbind(1:256, experimental.data.1)
experimental.data.2<-cbind(1:256, experimental.data.2)
experimental.data.3<-cbind(1:256, experimental.data.3)


#---------------------------------------Compute wavelet spectra-----------------------
wt.t1 <- wt(experimental.data.1)
wt.t2 <- wt(experimental.data.2)
wt.t3 <- wt(experimental.data.3)
w.arr <- array(dim = c(3, NROW(wt.t1$wave), 
                       NCOL(wt.t1$wave)))
w.arr[1, , ] <- wt.t1$wave
w.arr[2, , ] <- wt.t2$wave
w.arr[3, , ] <- wt.t3$wave

#--------------------------------------Multiresolution Decomposition------------------

z.decomposition.1<-wavMRDSum(as.vector(experimental.data.1), levels=3:6)


#----------------------------------------Partial wavelet coherence of 1 and 2---------
pwtc.1.2 <- pwtc(experimental.data.1,
                 experimental.data.2, 
                 experimental.data.3,
                 nrands = 0)
#----------------------------------------Partial wavelet coherence of 1 and 3----------
pwtc.1.3 <- pwtc(experimental.data.1,
                 experimental.data.3, 
                 experimental.data.2,
                 nrands = 0)

#----------------------------------------Compute dissimilarity--------------------------
dissimilarity.index<-c(wdist(wt.t1$wave, wt.t2$wave),wdist(wt.t1$wave, wt.t3$wave))
w.arr.dis <- wclust(w.arr)

#----------------------------------------Compute Cross-wavelet--------------------------
xwt.1.2 <- xwt(experimental.data.1, experimental.data.2)

#----------------------------------------Tables-----------------------------------------
index.dissimilarity.df<-data.frame()
index.dissimilarity.df<-cbind(dissimilarity.index)
colnames(index.dissimilarity.df)<-c("Index")
Table.1<-xtable(index.dissimilarity.df)

#----------------------------------Figures-----------------------------------------------

Figure.1<-plot(experimental.data.1, type="l", lty=1, 
               col="black",
               xlab="Observations",
               ylab="Value")
lines(experimental.data.2,lty=2, col="green")
rug(experimental.data.1, side=4, col="black")
rug(experimental.data.1, side=4, col="green")
legend("topright",
       c("Pearson III A",
         "Pearson III B"),
       inset = .01,
       col=c("black","green"),
       lwd=2,
       cex=0.8)

par(mfrow = c(2,1), oma = c(4, 0, 0, 1),
    mar = c(1, 4, 4, 5), mgp = c(1.5, 0.5, 0))

Figure.2<-plot(pwtc.1.2, xlab = "", plot.cb = TRUE,
     main = "Partial wavelet coherence of 1 and 2 | 3")

plot(pwtc.1.3, plot.cb = TRUE,
     main = "Partial wavelet coherence of 1 and 3 | 2")

par(oma = c(0, 0, 0, 1), mar = c(5, 4, 4, 5) + 0.1)
Figure.3<-plot(wt.t1, 
                 plot.cb = TRUE, 
                 plot.phase = TRUE)
Figure.4<-plot(wt.t2, 
               plot.cb = TRUE, 
               plot.phase = TRUE)
Figure.5<-plot(wt.t3, 
               plot.cb = TRUE, 
               plot.phase = TRUE)

Figure.6<-plot(xwt.1.2, plot.cb = TRUE, plot.phase = TRUE,
                 main = "Plot cross-wavelet and phase difference (arrows)")

Figure.7<-plot(hclust(w.arr.dis$dist.mat, method = "ward.D"),
                 sub = "", main = "", ylab = "Dissimilarity", hang = -1)

Figure.8<-ifultools::stackPlot(x=as.vector(experimental.data.1),
                                 y=list(experimental.data.1=x, "D3+D4+D5"=z.decomposition.1))
