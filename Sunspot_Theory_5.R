library(rgl)
library(xtable)
library(readr)
library(spdep)
library(spatstat)
library(scatterplot3d)
library(gstat)
#---------------------------------------Data from Reference 1-------------------------------------------------

Shapley.galaxy <- read.table("Shapley_galaxy.dat", header=T)

Shapley.galaxy.Mean<-Shapley.galaxy[Mean(Shapley.galaxy$R.A.)& Mean(Shapley.galaxy$Dec.) ,]
Shapley.galaxy.Median<-Shapley.galaxy[Median(Shapley.galaxy$R.A.)& Median(Shapley.galaxy$Dec.) ,]
Shapley.galaxy.Mode<-Shapley.galaxy[Mode(Shapley.galaxy$R.A.)& Mode(Shapley.galaxy$Dec.) ,]
Shapley.galaxy.Cluster<-Shapley.galaxy[Mean(Shapley.galaxy$R.A.)& 
                                      Mean(Shapley.galaxy$Dec.) & 
                                        Mean(Shapley.galaxy$Vel),]

summary(Shapley.galaxy)
rgl.open()
rgl.points(scale(Shapley.galaxy[,1]), scale(Shapley.galaxy [,2]), scale(Shapley.galaxy [,4]))
rgl.bbox()

#------------------------------------Spatial Dependence Analysis-----------------------------

Shapley.galaxy.Mean.mat <- as.matrix(Shapley.galaxy.Mean[,1:2])
Shapley.galaxy.Mean.nn.shap.lo <- knearneigh(Shapley.galaxy.Mean.mat, k=1) 
str(Shapley.galaxy.Mean.nn.shap.lo)
Shapley.galaxy.Mean.nb<- knn2nb(Shapley.galaxy.Mean.nn.shap.lo)
Shapley.galaxy.Mean.nb.wt<- nb2listw(Shapley.galaxy.Mean.nb,style='B')
summary(Shapley.galaxy.Mean.nb.wt)

#-----------------------------------Spatial Statistics---------------------------------------

Shapley.galaxy.Mean.win <- owin(range(Shapley.galaxy.Mean[,1]), range(Shapley.galaxy.Mean[,2]))
centroid.owin(Shapley.galaxy.Mean.win) ; area.owin(Shapley.galaxy.Mean.win)
Shapley.galaxy.Mean.ppp <- as.ppp(Shapley.galaxy.Mean[,c(1,2,4)], Shapley.galaxy.Mean.win) # planar point pattern
summary(Shapley.galaxy.Mean.ppp)

#-----------------------------------Moran and Geary C Statistics-------------------------------

Shapley.galaxy.Mean.Moran<-moran(Shapley.galaxy.Mean.mat[,1], Shapley.galaxy.Mean.nb.wt, n=length(Shapley.galaxy.Mean.nb), 
      S0=Szero(Shapley.galaxy.Mean.nb.wt))
Shapley.galaxy.Mean.Moran.test<-moran.test(Shapley.galaxy.Mean.mat[,1], Shapley.galaxy.Mean.nb.wt)
Shapley.galaxy.Mean.Moran.mc<-moran.mc(Shapley.galaxy.Mean.mat[,1], Shapley.galaxy.Mean.nb.wt, nsim=100)
Shapley.galaxy.Mean.Geary<-geary(Shapley.galaxy.Mean.mat[,1], Shapley.galaxy.Mean.nb.wt, n=length(Shapley.galaxy.Mean.nb), 
      n1=length(Shapley.galaxy.Mean.nb)-1, S0=Szero(Shapley.galaxy.Mean.nb.wt))
Shapley.galaxy.Mean.Geary.test<-geary.test(Shapley.galaxy.Mean.mat[,1], Shapley.galaxy.Mean.nb.wt)

Shapley.galaxy.Mean.Statistics.MoranNGeary.df<-data.frame()
Shapley.galaxy.Mean.Statistics.MoranNGeary.df<-rbind(c(Shapley.galaxy.Mean.Moran.test$statistic,Shapley.galaxy.Mean.Moran.test$p.value),
                                                     c(Shapley.galaxy.Mean.Geary.test$statistic,Shapley.galaxy.Mean.Geary.test$p.value))
colnames(Shapley.galaxy.Mean.Statistics.MoranNGeary.df)<-c("Statistic","P Value")
rownames(Shapley.galaxy.Mean.Statistics.MoranNGeary.df)<-c("Moran","Geary")

#-----------------------------------Variogram Analysis------------------------------------------------

Shapley.galaxy.Mean.variog <- variogram(Vel~1, locations=~R.A.+Dec., data=Shapley.galaxy.Mean)
Shapley.galaxy.Mean.variog.mod1 <- vgm(7e+07, "Gau", 3.0,2e+07)
Shapley.galaxy.Mean.variog.fit <- fit.variogram(Shapley.galaxy.Mean.variog,Shapley.galaxy.Mean.variog.mod1) 
Shapley.galaxy.Mean.variog.fit

#-----------------------------------K, L and J functions for the Mean Density Region-------------------

Shapley.galaxy.Mean.K <- Kest(Shapley.galaxy.Mean.ppp, correction='isotropic')
Shapley.galaxy.Mean.K.bias <- Kest(Shapley.galaxy.Mean.ppp, correction='none')
Shapley.galaxy.Mean.L <- Lest(Shapley.galaxy.Mean.ppp, correction='isotropic')
Shapley.galaxy.Mean.L.bias <- Lest(Shapley.galaxy.Mean.ppp, correction='none')
Shapley.galaxy.Mean.J<-Jest(Shapley.galaxy.Mean.ppp)

#--------------------------------------Tables------------------------------------------------

Table.1<-xtable(Shapley.galaxy.Mean.Statistics.MoranNGeary.df)

#-------------------------------------Figures------------------------------------------------
Figure.0<-plot(Shapley.galaxy.Mean[,1], Shapley.galaxy.Mean[,2], cex=(scale(Shapley.galaxy.Mean[,4])+1.5)/2, 
     xlab='Right Ascension (degrees)', ylab='Declination (degrees)')
Figure.0.1<-scatterplot3d(Shapley.galaxy.Mean[,c(1,2,4)], pch=20, cex.symbols=0.7, type='p', angl=40, 
              zlim=c(0,50000))
rgl.snapshot('Figure1.png')
rgl.close()
Figure.2<-plot.nb(Shapley.galaxy.Mean.nb, Shapley.galaxy.Mean[,1:2])
Figure.3<-plot(density(Shapley.galaxy.Mean.ppp,0.3), col=topo.colors(20), main='', xlab='R.A.', 
     ylab='Dec.')
Figure.4<-plot(Shapley.galaxy.Mean.ppp, lwd=2, add=T)
opa <- par(mfrow=c(2,2),mar=c(2,2,2,2))
Figure.5<-plot.fv(Shapley.galaxy.Mean.K, lwd=2, col='black',  main='', xlab='r (degrees)', legend=F)
Figure.6<-plot.fv(Shapley.galaxy.Mean.K.bias, add=T, lty=3, lwd=2, col='black', legend=F)
Figure.7<-plot(Shapley.galaxy.Mean.L$r, (Shapley.galaxy.Mean.L$iso - Shapley.galaxy.Mean.L$r), lwd=2, col='black',
     main='',xlab='r (degrees)', ylab='L*(r)', ty='l', ylim=c(-0.2,0.2))
lines(Shapley.galaxy.Mean.L$r, (Shapley.galaxy.Mean.L$theo - Shapley.galaxy.Mean.L$r), lwd=2, lty=2)
lines(Shapley.galaxy.Mean.L$r, (Shapley.galaxy.Mean.L.bias$un - Shapley.galaxy.Mean.L$r), lwd=2, lty=3)
par(opa)
Figure.8<-plot(Shapley.galaxy.Mean.variog,model <- Shapley.galaxy.Mean.variog.fit, col='black', pch=20, 
     xlab='Distance (degree)', ylab="Semivariance (km/s*km/s)", lwd=2)
#------------------------------------Reference-----------------------------------------------

Reference.1<-c("Eric D. Feigelson & G. Jogesh Babu",
               "Modern Statistical Methods for Astronomy With R Applications", 
               "Cambridge University Press (2012)") 


#------------------------------------Function Library----------------------------------------