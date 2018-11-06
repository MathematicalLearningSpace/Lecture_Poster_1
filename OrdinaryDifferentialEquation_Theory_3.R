#-----------------RPI 
library(deSolve)
library(RecTran)
library(rootSolve)
library(Sim.DiffProc)
library(bvpSolve)
library(odeintr)
library(scaRabee)
library(yuima)
library(Sim.DiffProc)
library(fptdApprox)
library(rpgm)
library(yuima)
#--------------------Bioconductor libraries go here for preprocessing the gene and protein interaction networks
#
library(Biostrings)
library(seqinr)
library(igraph)
library(Matrix)
library(rgl)
#
#-------------------Combinatorics
#
library(combinat)
library(GeomComb)

#Admin and presentation
library(ggplot2)
library(xtable)

#---------------------------------------------------------Data------------------------------------------------------------------

Pearson.N<-100
#--------------------------------------------------------Brownian motion
set.seed(1234)
X.1 <- BM(M = Pearson.N)
#-------------------------------------------------------Brownian bridge
X.2 <- BB(M =Pearson.N)
#-------------------------------------------------------Geometric Brownian motion
X.3 <- GBM(M = Pearson.N)
#-------------------------------------------------------Arithmetic Brownian motion
X.4<- ABM(M = Pearson.N)

#---------------------------------------------------------Models and Expression-------------------------------------------------------------
#f(t,w(t)) = int(exp(w(t) - 0.5*t) * dw(s)) with t in [0,1]
f <- expression( exp(w-0.5*t) )
mod1 <- st.int(expr=f,
               type="ito",
               M=50,
               lower=0,
               upper=1)
summary(mod1)

#---------------------------------------------------------Parameter Values------------------------------------
a1<-1
a2<-1
a3<-1
d1<-1
d2<-1
d3<-1
parameters.1<-c(a1,a2,a3)
parameters.2<-c(d1,d2,d3)

#--------------------------------------------------------Drift------------------------------------------------
fx <- expression(a1*(-d1-x)*y , 
                 a2*(d2-y)*x , 
                 a3*(d3-z)*y)
#--------------------------------------------------------Diffusion---------------------------------------------
gx <- rep(expression(0.5),3)
mod3d.1 <- snssde3d(drift=fx,
                    diffusion=gx,
                    M=500)

#------------------------------------------------------- boundary specification--------------------------------
beta1<-1
beta2<-1
gamma1<-1
St.1 <- expression(beta1+beta2*t^gamma1)
mod3d.1.boundary <- fptsde3d(mod3d.1,boundary=St.1)

#--------------------------------------------------------Bridge Estimation------------------------------------------
M1<-10^1
M2<-10^2
M3<-10^3
#--------------------------------------------------------Models------------------------------------------------------
x0<-c(0,-1,0.5)
y<-c(0,-2,0.5)

bridge.model.1<- bridgesde3d(x0=c(0,-1,0.5),
                    y=c(0,-2,0.5),
                    drift=fx,
                    diffusion=gx,
                    M=M2)
#------------------------------------Marginal density of X(t-t0)|X(t0) = 0, X(T) = 0 at time t = 0.75-----------------
s=0.75
bridge.denM.1 <- dsde3d(bridge.model.1,pdf="M",at =s)
#--------------------------------------------------Joint density------------------------------------------------------ 
bridge.denJ.1 <- dsde3d(bridge.model.1,pdf="J",at=0.75)



parameters.df<-data.frame()
parameters.df<-cbind(parameters.1,parameters.2,x0,y0)

#--------------------------------------------Tables-------------------------------------------------------------------

Table.1<-xtable(parameters.df)
#--------------------------------------------Figures-----------------------------------------------------------------

Figure.1<-plot(X.1,plot.type="single")
lines(as.vector(time(X.2)),rowMeans(X.2),col="red")
lines(as.vector(time(X.3)),rowMeans(X.3),col="green")
lines(as.vector(time(X.4)),rowMeans(X.4),col="blue")
legend("topleft",c("Brownian Motion","Brownian Bridge",
         "Geometric Brownian Motion", paste("Arithmetic"," ","Brownian Motion")),
       inset = .01,
       col=c(4,2,5),
       lwd=2,
       cex=0.8)

Figure.2<-plot(bridge.model.1)
lines(time(bridge.model.1),apply(bridge.model.1$X,1,mean),col=2,lwd=2)
lines(time(bridge.model.1),apply(bridge.model.1$X,1,bconfint,level=0.95)[1,],col=4,lwd=2)
lines(time(bridge.model.1),apply(bridge.model.1$X,1,bconfint,level=0.95)[2,],col=4,lwd=2)
legend("topleft",
       c("mean path",
         paste("bound of", 95," percent confidence")),
       inset = .01,col=c(2,4),lwd=2,cex=0.8)

Figure.3<-plot3D(mod3d.1,display = "persp",main="No-Boundary")


Figure.4<-plot3D(bridge.model.1,display = "persp",main="3D Bridge SDE's")

Figure.5<- plot(ts.union(bridge.model.1$X[,1],
                         bridge.model.1$Y[,1],
                         bridge.model.1$Z[,1]),
                col=1:3,lty=3,
                plot.type="single",
                type="l",ylab= "",xlab="time",axes=T)
Figure.6<-plot(bridge.denM.1, main="Marginal Density")
Figure.7<-plot(bridge.denJ.1,display="rgl")
#--------------------------------------------Function Library--------------------------------------------------------
