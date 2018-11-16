#---------------------------------------------------R API --------------------------------------------------------
library(SpectralMap); library(PearsonDS);library(moonsun);library(magicaxis);library(astro);library(astrodatR);library(cosmoFns)
library(astrolibR);library(sphereplot);library(stellaR)
#-------------------------------------------------------Generate Two Datasets---------------------------------------
pearson.N <- 100
#-------------------------------------------------------Attractor  I------------------------------------------------
k.1<-2^2 + 2^1
k.2<-2^2
#-------------------------------------------------------Attractor II------------------------------------------------
k.1<-2^3+2^2 + 2^1
k.2<-2^3
#-------------------------------------------------------Noise Specification-----------------------------------------
pIIpars.1 <- list(a=2, location=1, scale=2) 
pIIpars.2 <- list(a=2, location=1, scale=4) 
error.pearson.2.1<-rpearsonII(pearson.N,params=pIIpars.1)
error.pearson.2.2<-rpearsonII(pearson.N,params=pIIpars.2)
#-------------------------------------------------------Equation System---------------------------------------------
K<-k.2/k.1
theta <- 2*pi*seq(from=0, to=1-1/n, by=1/n)
r.1=(1 + cos(k.1*theta)/k.2)
r.2 = (1 + cos(k.1*theta)/k.2) + error.pearson.2.1
r.3 = (1 + cos(k.1*theta)/k.2) + error.pearson.2.2

#-------------------------------------------------------X is a circle------------------------------------------------
X1 = cos(theta)
X2 = sin(theta)
X<-data.frame(X1,X2)

#-------------------------------------------------------Y is a hexagon-----------------------------------------------
Y1.1 = r.1*cos(theta);Y2.1 = r.1*sin(theta)
Y1.2 = r.2*cos(theta);Y2.2 = r.2*sin(theta)
Y1.3 = r.3*cos(theta);Y2.3 = r.3*sin(theta)

Y.1<-data.frame(Y1.1,Y2.1)
Y.2<-data.frame(Y1.2,Y2.2)
Y.3<-data.frame(Y1.3,Y2.3)

Data.1<-list(X,Y.1)
Data.2<-list(X,Y.2)
Data.3<-list(X,Y.3)
#------------------------------------------------------Create the spectral map from X to Y
par(mfrow = c(2,2))
example2.1<-SpectralMap(Data.1, epsilon=0.1, range=1:2, Plot2D=TRUE, Plot3D=FALSE)
example2.2<-SpectralMap(Data.2, epsilon=0.1, range=1:2, Plot2D=TRUE, Plot3D=FALSE)
example2.3<-SpectralMap(Data.3, epsilon=0.1, range=1:2, Plot2D=TRUE, Plot3D=FALSE)
#------------------------Tables-------------------------------------------------------
#------------------------Figures for Classroom Presentation------------------------------------------------------
par(mfrow = c(2,2))

Figure.1<-plot(r.1, type="l")
Figure.1A<-plot(r.2, type="l")
Figure.1B<-plot(r.3, type="l")

Figure.2<-plot(X)
Figure.3<-plot(Y.1)
Figure.4<-plot(Y.2)
Figure.5<-plot(Y.3)
