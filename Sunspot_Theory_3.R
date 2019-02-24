#-----------------------------R Code To Modify in the Classroom Lecture with Students-----------------------
#---------------------------------------R API ------------------------------------------------------------
library(xtable);library(quantreg);library(quantregForest);library(readr);library(readxl);library(nlstools);library(car);library(MASS)
#------------------------------------Data From Reference 1----------------------------------------------------
NGC4472 <- read.table("NGC4472_profile.dat", header=T) #NGC 4472 elliptical galaxy surface brightness profile Reference.1
NGC4406 <- read.table("NGC4406_profile.dat", header=T) #NGC 4406 Reference.1
NGC4551 <- read.table("NGC4551_profile.dat", header=T) #NGC 4551 Reference.1
#-----------------------------------Model 1 From Reference 1--------------------------------------------------

NGC4472.fit.nls <-  nls(surf_mag ~ -2.5*log10(I.e * 10^(-(0.868*n-0.142)*((radius/r.e)^{1/n}-1))) + 26, 
                    data=NGC4472, 
                    start=list(I.e=20.,r.e=120.,n=4.), model=T, trace=T)

summary(NGC4472.fit.nls)
NGC4472.fit.nls.coef.df<-as.data.frame(coef(NGC4472.fit.nls))
deviance(NGC4472.fit.nls)
logLik(NGC4472.fit.nls)

NGC4472.boot <- nlsBoot(NGC4472.fit.nls)
summary(NGC4472.boot)

#--------------------------------Model 2 and 3----------------------------------------------------

boxcox(surf_mag~log(radius), data =NGC4406,lambda = seq(-0.25, 0.25, length = 10))
NGC4406.fit <- lm(surf_mag ~ boxCoxVariable(log(radius)),  NGC4406)
summary(NGC4406.fit)
NGC4551.fit <- lm(surf_mag ~ boxCoxVariable(radius),  NGC4551)
summary(NGC4551.fit)
#----------------------------------Model 1 Analysis-----------------------------------------

qqnorm(residuals(NGC4472.fit.nls) / summary(NGC4472.fit.nls)$sigma) 
abline(a=0,b=1)
shapiro.test(residuals(NGC4472.fit.nls) / summary(NGC4472.fit.nls)$sigma)
shapiro.test(residuals(NGC4406.fit) / summary(NGC4406.fit)$sigma)
shapiro.test(residuals(NGC4551.fit) / summary(NGC4551.fit)$sigma)
#----------------------------------Tables----------------------------------------------------

Table.1<-xdata(NGC4472.fit.nls.coef.df)

#-----------------Figures for Presentation in the Classroom----------------------------------------------------

Figure.1<-plot(NGC4472$radius, NGC4472$surf_mag,lty=1,col=1,xlab="r  (arcsec)", 
               ylab=expression(mu ~~ (mag/sq.arcsec)), ylim=c(16,28), 
               cex.lab=1.5, cex.axis=1.5)
lines(NGC4406$radius, NGC4406$surf_mag,lty=2,col=2)
lines(NGC4551$radius, NGC4551$surf_mag,lty=3,col=3)
grid()
legend("bottomright", col = 1:3, 
       lty = 1:3, cex=0.75,
       legend = c("NGC4472","NGC4406","NGC4551"))

Figure.2<-plot(NGC4472.fit.nls$model$radius, NGC4472.fit.nls$model$surf_mag, pch=20, 
     xlab="r  (arcsec)", ylab=expression(mu ~~ (mag/sq.arcsec)), ylim=c(16,28), 
     cex.lab=1.5, cex.axis=1.5,col=1)
lines(NGC4472.fit.nls$model$radius, fitted(NGC4472.fit.nls),col=2)

Figure.3<-plot(NGC4472.fit.nls$model$radius,residuals(NGC4472.fit.nls), xlab="r (arcsec)", 
     ylab="Residuals", pch=20, cex.lab=1.5, cex.axis=1.5)
lines(supsmu(NGC4472.fit.nls$model$radius, residuals(NGC4472.fit.nls), span=0.05), lwd=2)

Figure.4<-curve(dnorm(x,m=5.95, sd=0.10)*58/5.95, xlim=c(5.6,6.4), ylim=c(0,50))
hist(NGC4472.boot$coefboot[,3], breaks=50, add=T)

Figure.5<-plot(NGC4406.fit$model$`boxCoxVariable(log(radius))`, NGC4406.fit$model$surf_mag, pch=20, 
xlab="transformed r  (arcsec)", ylab=expression(mu ~~ (mag/sq.arcsec)), ylim=c(16,24), 
cex.lab=1.5, cex.axis=1.5,col=1)
lines(NGC4406.fit$model$`boxCoxVariable(log(radius))`, NGC4406.fit$fitted.values,col=2)

Figure.6<-plot(NGC4551.fit$model$`boxCoxVariable(radius)`, NGC4551.fit$model$surf_mag, pch=20, 
               xlab="transformed r  (arcsec)", ylab=expression(mu ~~ (mag/sq.arcsec)), ylim=c(16,24), 
               cex.lab=1.5, cex.axis=1.5,col=1)
lines(NGC4551.fit$model$`boxCoxVariable(radius)`, NGC4551.fit$fitted.values,col=2)

#------------------Additional References to be Added by Students-------------------------------------------------
Reference.1<-c("Eric D. Feigelson & G. Jogesh Babu",
               "Modern Statistical Methods for Astronomy With R Applications", 
               "Cambridge University Press (2012)") 

#----------------------------------Function Library-------------------------------------------

