#------------------------------------------------R API ------------------------------------------------
library(SpectralMap);library(PearsonDS);library(moonsun);library(magicaxis);library(astro)
library(astrodatR);library(cosmoFns);library(astrolibR);library(sphereplot);library(stellaR)
library(TED);library(spectral);library(BayesSpec);library(beyondWhittle);library(bootspecdens)
library(pavo);library(ftsspec);library(frequencyConnectedness);library(EMD);library(colorSpec)
library(bspec);library(quantspec);library(spectral.methods);library(abd);library(xtable)

library(TSA);library(tsDyn);library(DescTools)
#-------------------------------------Generate Data--------------------------------------------------------------------

#--------------------------------------Noisy signal with amplitude modulation-----------------------------------------
x <- seq(0,1, length.out=1000)
p.1<--1
p.2<-0
p.3<-1
q=3
sine.data.1 <- x^(p.1)*sin((2^(q))*pi*x)
sine.data.2 <- x^(p.2)*sin((2^(q))*pi*x)
sine.data.3 <- x^(p.3)*sin((2^(q))*pi*x)

#-------------------------------------sunspot----------------------------------------------------------------------------
#data("sunspot.classic")
#data("ss08")
#data("ss08.1850")

experimental.data<-sunspots
sunspot.N<-length(experimental.data)
experimental.data.1<-window(experimental.data)
experimental.data.train<-window(experimental.data)
experimental.data.test<-window(experimental.data)

#-------------------------------------Descriptive Statistics-------------------------------------------------------------
distribution.moments.1.df<-data.frame()
distribution.moments.1.df<-rbind(c("Sunspot","X1",empMoments(experimental.data)),
                                 c("Sunspot","Subset",empMoments(experimental.data.1)),
                                 c("Sunspot","Train",empMoments(experimental.data.train)),
                                 c("Sunspot","Test",empMoments(experimental.data.test)))
colnames(distribution.moments.1.df)<-c("Data","Variable","Moment 1", "Moment 2","Moment 3","Moment 4")

#-----------------------------------------------------Stationarity Test--------------------------------------------------

test.1<-pp.test(experimental.data.train )
test.2<-kpss.test(experimental.data.train )
test.3<-adf.test(experimental.data.train ,alternative="stationary")
test.4<-adf.test(experimental.data.train ,alternative="explosive")


stationarity.test.df<-data.frame()
stationarity.test.df<-rbind(c(test.1$statistic,test.1$p.value,
                              test.2$statistic,test.2$p.value,
                              test.3$statistic,test.3$p.value,
                              test.4$statistic,test.4$p.value))
colnames(stationarity.test.df)<-c("Phillips-Perron","P Value","KPSS","P Value",
                                  "ADF1","P Value","ADF2","P Value")

#---------------------------------------------------Nonlinear Test--------------------------------------------------------
bds.test<-bds.test(experimental.data.train)
K.test<-Keenan.test(abs(experimental.data.train))
mcleod.test<-McLeod.Li.test(y=experimental.data.train)

nonlinear.test.df<-data.frame()
nonlinear.test.df<-rbind(c(K.test$test.stat,K.test$p.value))
colnames(nonlinear.test.df)<-c("K","P Value")

#------------------------------------Fourier Transform--------------------------
x <- seq(1,length(experimental.data))
FT <- spec.fft(experimental.data, x)

#------------------------------------AdaptSPec----------------------------------
model1 <- adaptspec(nloop = 10, nwarmup = 10,nexp_max = 2, x = experimental.data)
str(model1)
summary(model1$nexp_curr)
#------------------------------------Tests-----------------------------------------------
w<-24
sunspot.noise.test.1<-noiseTests(experimental.data,w,'white',parallel=FALSE)
sunspot.noise.test.2<-noiseTests(experimental.data,w,'red',parallel=FALSE)

data.measure.df<-data.frame()
data.measures.df<-cbind(measures(experimental.data))

#-----------------------------------Fit an AR(p) model to sunspot data
ar.order <- 2
mcmc <- gibbs_AR(data=experimental.data, Ntotal=500, burnin=100, ar.order=ar.order)
ar.order <-2^2
mcmc.1 <- gibbs_AR(data=experimental.data, Ntotal=500, burnin=100, ar.order=ar.order)

N <- length(mcmc$psd.median)
pdgrm <- (abs(fft(experimental.data))^2 / (2*pi*length(experimental.data)))[1:N]

#----------------------------------------Tables-----------------------------------------
Table.1<-xtable(data.measures.df)
Table.2<-xtable(distribution.moments.1.df)
Table.3<-xtable(stationarity.test.df)
Table.4<-xtable(nonlinear.test.df)
#----------------------------------Figures for Classroom Presentation-----------------------------------------------

Figure.A<-plot.ts(sine.data.1, col="blue",lty=1, ylim=c(-5,5))
lines(sine.data.2,col="green",lty=2)
lines(sine.data.3,col="red",lty=3)
legend(x="topright", 
       legend=c("Sine 1","Sine 2","Sine3"), 
       lty=c(1,2,3), col=c("blue", "green","red"))


Figure.1<-plot.ts(experimental.data, col="blue",lty=1)
lines(experimental.data.1,col="green",lty=2)
lines(experimental.data.train,col="red",lty=3)
lines(experimental.data.test,col="black",lty=4)
legend(x="topright", 
       legend=c("Sunspot 1","Sunspot Subset","Sunspot Train","Sunspot Test"), 
       lty=c(1,2,3,4), col=c("blue", "green","red","black"))

par(mfrow = c(2, 1))
Figure.2<-plot(x, experimental.data, type = "l", main = "Sunspots",col="blue")
Figure.2A<-plot(
  FT,
  ylab = "Amplitude",
  xlab = "Frequency",
  type = "l",
  main = "Spectrum", ylim=c(0,15), xlim=c(0,0.05),col="green")

Figure.3<-plot(sunspot.noise.test.1)

Figure.4<-plot(model1$nexp_curr, typ="l")

Figure.5<-plot.ts(log(pdgrm[-1]), col="blue", 
        main=paste0("AR Model Sunspot 1 on logarithmic scale for p=", ar.order))
lines(log(mcmc$psd.median[-1]),lty=2, col="green")
lines(log(mcmc$psd.p05[-1]),lty=3,col="red")
lines(log(mcmc$psd.p95[-1]),lty=3,col="red")
lines(log(mcmc$psd.u05[-1]),lty=4,col="black")
lines(log(mcmc$psd.u95[-1]),lty=4,col="black")
legend(x="topright", 
       legend=c("periodogram", "pointwise median", 
                              "pointwise CI (5,95)", "uniform CI(5,95)"), 
       lty=c(1,2,3,3,4,4), col=c("blue", "green", "red", "red","black","black"))
