library(moonsun);library(magicaxis);library(astro);library(astrodatR);library(cosmoFns);library(astrolibR);library(sphereplot)
library(stellaR);library(PearsonDS);library(readr);library(HistogramTools);library(DescTools)
#----------------------------------HYG star database archive Data-----------------------------------------
data(starcat)
data(bright)
#------------------Read the Data Set in the Classroom for Student Mathematical Notes---------------------
hygdata.v3 <- read_csv("HYG-Database-master/hygdata_v3.csv")
View(hygdata.v3)

summary(na.omit(hygdata.v3$ci))
(kernels <- eval(formals(density.default)$kernel))
h<-hist(na.omit(hygdata.v3$ci), plot=FALSE)

#---------------------------------Analysis---------------------------------------------------------------
hygdata.distribution.MLE<-pearsonFitML(na.omit(hygdata.v3$ci))

moments<-c(mean=mean(na.omit(hygdata.v3$ci)),
           variance=var(na.omit(hygdata.v3$ci)),
           skewness=Skew(na.omit(hygdata.v3$ci)),
           kurtosis=Kurt(na.omit(hygdata.v3$ci)))
hygdata.distribution.M<-pearsonFitM(moments=moments)

hygdata.v3.B<-hygdata.v3[grep("B",hygdata.v3$spect),]
hygdata.v3.O<-hygdata.v3[grep("O",hygdata.v3$spect),]

summary(na.omit(hygdata.v3.B$ci))
summary(na.omit(hygdata.v3.O$ci))

h.B<-hist(na.omit(hygdata.v3.B$ci), plot=FALSE)
h.O<-hist(na.omit(hygdata.v3.O$ci), plot=FALSE)

#intersect.dist(h.B, h.O)
#kl.divergence(h.B, h.O)
#jeffrey.divergence(h.B, h.O)

#-------------------------------Planet Pathways------------------------------------------
angle(bright)
angle(planets())
j=jd(2017,1,1,length = 365)
j.1=jd(1600,1,1,length=365)

ephem.mercury <- mercury(jd(2017,1,1,length = 365))
ephem.venus <- venus(jd(2017,1,1,length = 365))
ephem.mars<-mars(jd(2017,1,1,length = 365))

#------------------Tables to Be Done By Students--------------------------------------------------


#-------------Figures for Presentation in Classroom-------------------------------------------------

par(mfrow=c(1,2))
Figure.1<-plot(h, xlab="(B-V) Color Index",main="Histogram for All Spectral Types")
Figure.2<-plot(HistToEcdf(h), main="CDF") 

par(mfrow=c(1,1))
Figure.3<-plot(density(na.omit(hygdata.v3$ci)),col="black", xlim=c(-1,3), main="Smooth Kernel Density Estimation of the Color Index")
lines(density(na.omit(hygdata.v3$ci),kernel="epanechnikov"), col = "blue")
lines(density(na.omit(hygdata.v3$ci),kernel="rectangular"), col = "red")
lines(density(na.omit(hygdata.v3$ci),kernel="triangular"), col = "green")
lines(density(na.omit(hygdata.v3$ci),kernel="biweight"), col = "orange")
lines(density(na.omit(hygdata.v3$ci),kernel="cosine"), col = "purple")
lines(density(na.omit(hygdata.v3$ci),kernel="optcosine"), col = "yellow")
legend(1.9,.4, legend = kernels, col = seq(kernels),
       lty = 1, cex = .8, y.intersp = 1)
par(mfrow=c(1,2))
Figure.4<-plot(h.B, xlab="Color Index",main="Histogram for B Spectral Types")
Figure.5<-plot(h.O, xlab="Color Index",main="Histogram for O Spectral Types")

Figure.6<-plot(angle(mercury(j),venus(j)),lty=1,col=1)
lines(angle(venus(j),mars(j)),lty=2,col=2)
lines(angle(mercury(j),mars(j)),lty=3,col=3)
lines(angle(mercury(j.1),venus(j.1)),lty=4,col=4)

Figure.7<-plot(angle(sun(j),mercury(j)),lty=1,col=1)

par(mfrow=c(2,2))
Figure.8<-track(sun(jd(2017,1,1,length = 365)),col.track = "green", lwd = 1, mag = 5, starcat = starcat, bright = bright)
Figure.9<-track(ephem.mercury, col.track = "green", lwd = 1, mag = 5, starcat = starcat, bright = bright)
Figure.10<-track(ephem.venus, col.track = "blue", lwd = 2, mag = 5, starcat = starcat, bright = bright)
Figure.11<-track(ephem.mars, col.track = "green", lwd = 3, mag = 5, starcat = starcat, bright = bright)

par(mfrow=c(2,2))
Figure.12<-plot(bright)
Figure.13<-plot(as.ecc(bright))
Figure.14<-plot(as.hoc(bright))
Figure.15<-plot(as.lt(rst(bright)))

#------------------Reference To be Added by Students-------------------------------------------

#------------------Function Library To be Added by Students------------------------------------



