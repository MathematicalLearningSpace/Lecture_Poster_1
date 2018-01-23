library(xtable)
library(quantreg)
library(quantregForest)
library(readr)
library(readxl)
library(nlstools)
library(car)
library(MASS)
library(rgl)
library(jsonlite)
library(PearsonDS)
library(spatstat)
library(splines)
#---------------------------------Data From Reference 1-------------------------------------

Data.6D<-read.table("Data_6D.dat",header=T,fill=T)
Data.6D.name<-c("Data_6D")

Data.6D.names<-Data.6D[,1] # name of math object
Feature.1<-Data.6D[,2] # x coordinate
Feature.2<-Data.6D[,3] # y coordinate
Feature.3<-Data.6D[,4] # z coordinate
Feature.4<-Data.6D[,5] # time
Feature.5<-Data.6D[,6] # color

#--------------------------------Analysis------------------------------------
summary(Data.6D)
Data.6D.Statistics.df<-data.frame()
Data.6D.correlation<-cor(Data.6D[,2:5], method='kendall')

Data.6D.pattern.1<-pp3(Data.6D[,2],Data.6D[,3],Data.6D[,4],box3(c(0,1)))
Data.6D.pattern.1.d<-pairdist(Data.6D.pattern.1,period=TRUE)
Data.6D.pattern.1.summary<-summary(Data.6D.pattern.1)

Data.6D.pattern.2.window.polygon<-ppp(Data.6D[,2],Data.6D[,3], poly=list(x=c(0,10,0), y=c(0,0,10)))
Data.6D.pattern.2.window.polygon.fit <- ppm(Data.6D.pattern.2.window.polygon, ~ polynom(x,y,2), Poisson())
Data.6D.pattern.2.window.polygon.fit.1<-ppm(Data.6D.pattern.2.window.polygon ~ bs(x,df=5))
Data.6D.pattern.2.window.polygon.fit.2<-ppm(Data.6D.pattern.2.window.polygon ~ polynom(x,2))

Data.6D.pattern.2.window.polygon.fit.summary.df<-data.frame()
Data.6D.pattern.2.window.polygon.fit.summary.df<-rbind(Data.6D.pattern.2.window.polygon.fit$coef)

Data.6D.pattern.2.window.polygon.fit.auc<-auc(Data.6D.pattern.2.window.polygon.fit)
Data.6D.pattern.2.window.polygon.fit.auc.1<-auc(Data.6D.pattern.2.window.polygon.fit.1)
Data.6D.pattern.2.window.polygon.fit.auc.2<-auc(Data.6D.pattern.2.window.polygon.fit.2)

Data.6D.pattern.2.window.polygon.fit.lambda <- predict(Data.6D.pattern.2.window.polygon.fit, type="trend")
Data.6D.pattern.2.window.polygon.fit.lambda.2 <- predict(Data.6D.pattern.2.window.polygon.fit,
                                                         locations=Data.6D.pattern.2.window.polygon,
                                                         type="trend")
Data.6D.pattern.2.window.polygon.Ki.0 <- Kinhom(Data.6D.pattern.2.window.polygon, 
                                                Data.6D.pattern.2.window.polygon.fit.lambda )
Data.6D.pattern.2.window.polygon.Ki.1 <- Kinhom(Data.6D.pattern.2.window.polygon, sigma=0.1)
lamfun <- function(x,y) { 50 + 100 * x }
Y <- rpoispp(lamfun, 150, owin())
Data.6D.pattern.2.window.polygon.Ki.2 <- Kinhom(Y, lamfun)
#-------------------------------Tables--------------------------------------------------------

Table.1<-xtable(Data.6D.Statistics.df)
Table.2<-xtable(Data.6D.pattern.2.window.polygon.fit.summary.df)

#----------------------------------Figures----------------------------------------------------

Figure.1<-dotchart(Feature.1, labels=Data.6D.names, cex=0.9, xlab=expression(Feature.1~(units)))
Figure.2<-plot(Feature.1, ylim=c(0,8), xlab=Data.6D.name, ylab=expression(Feature.1~(units)),pch=20)
Figure.3<-boxplot(Data.6D[,2:3], varwidth=T, notch=T, xlab=Data.6D.name, ylab="Feature.1 and Feature.2",
                  pars=list(boxwex=0.3,boxlwd=1.5,whisklwd=1.5,staplelwd=1.5,outlwd=1.5,
                            font=2))
Figure.4<-qqnorm(Data.6D[,4], main='', xlab='Feature.3', ylab='Value', cex.lab=1.3, cex.axis=1.3, pch=20)
par(mfrow=c(2,3))
histograms(Data.6D[,2:5])

open3d() 
Figure.6<-plot3d(Data.6D[,2:5])

snapshot3d(file='Figure6.png')

Figure.7<-plot(roc(Data.6D.pattern.2.window.polygon.fit),
               main="ROC of Model 1: Spatial Trend and Polynomial")
lines(roc(Data.6D.pattern.2.window.polygon.fit.1),col="blue")
lines(roc(Data.6D.pattern.2.window.polygon.fit.2),col="green")
legend("topright", legend =c("Spatial Trend","B spline","Polynomial"), col =c("black","blue","green"),
       lty = 1:3, cex = .8, y.intersp = 1)

par(mfrow=c(2,2))
Figure.8<-plot(Data.6D.pattern.2.window.polygon.Ki.0,main="predict intensity at all locations")
Figure.9<-plot(Data.6D.pattern.2.window.polygon.Ki.1,main="intensity function estimated by heavy smoothing")
Figure.10<-plot(Data.6D.pattern.2.window.polygon.Ki.2, main="Inhomogeneous K function")
Figure.11<-plot(density(Data.6D.pattern.2.window.polygon))
#----------------------------------References-------------------------------------------------
Reference.1<-c("",
               "", 
               "") 

#----------------------------------Function Library-------------------------------------------
histograms <- function(x) { 
  par(mfrow=n2mfrow(dim(x)[2]))  
  for (i in 1:dim(x)[2]) { name=names(x)[i]                
  hist(x[,i], main='', breaks='FD', ylab='', xlab=name) }}
