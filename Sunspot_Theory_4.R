#-----------------------------R Code To Modify in the Classroom Lecture with Students-----------------------
#----------------------------------R API -------------------------------------------------------------------
library(moonsun);library(magicaxis);library(astro);library(astrodatR);library(cosmoFns);library(astrolibR);library(sphereplot)
library(stellaR);library(PearsonDS);library(readr);library(HistogramTools);library(DescTools);library(spatstat);library(G2Sd)
library(plot3D)
#----------------------------------HYG star database archive Data-----------------------------------------
data(starcat)
data(bright)
hygdata.v3 <- read_csv("HYG-Database-master/hygdata_v3.csv")
stellar.mag.1<-as.data.frame(cbind(c(hygdata.v3$spect),
                                   c(hygdata.v3$x),
                                   c(hygdata.v3$y),
                                   c(hygdata.v3$z),
                                   c(hygdata.v3$mag),
                                   c(hygdata.v3$ci)   
))
colnames(stellar.mag.1)<-c("Spectral","x","y","z","Magnitude","ColorIndex")

box.1<-box3(c(0,max(as.numeric(stellar.mag.1$x))),
            c(0,max(as.numeric(stellar.mag.1$y))),
            c(0,max(as.numeric(stellar.mag.1$z))))
diameter(box.1) 
volume(box.1)
sidelengths(box.1)
shortside(box.1)
hd <- shortside(box.1)/2
eroded.volumes(box.1, seq(0,hd, length=10))
X.stellar.mag<- pp3(as.numeric(stellar.mag.1$x), 
                    as.numeric(stellar.mag.1$y),
                    as.numeric(stellar.mag.1$z), 
                    box.1)
#---------------------------------Analysis---------------------------------------------------------------
summary(X.stellar.mag)
X.stellar.mag.O<-stellar.mag.1[grep("O",stellar.mag.1$Spectral),]
Mean(as.numeric(X.stellar.mag.O$Magnitude))
MeanCI(as.numeric(X.stellar.mag.O$Magnitude), sides="left")
RunsTest(as.numeric(X.stellar.mag.O$Magnitude) > mean(as.numeric(X.stellar.mag.O$Magnitude)))

Median(as.numeric(X.stellar.mag.O$Magnitude))
MedianCI(as.numeric(X.stellar.mag.O$Magnitude))
RunsTest(as.numeric(X.stellar.mag.O$Magnitude) > median(as.numeric(X.stellar.mag.O$Magnitude)))

X.stellar.mag.O.table.freq<-xtable(Freq(as.numeric(X.stellar.mag.O$Magnitude), ord="desc"))

X.stellar.mag.O.table.correlation <- cor(as.numeric(X.stellar.mag.O$Magnitude),
                         as.numeric(X.stellar.mag.O$ColorIndex), 
                         use="pairwise.complete.obs")
X.stellar.mag.O.cor.test.1<-cor.test(as.numeric(X.stellar.mag.O$Magnitude), 
                     as.numeric(X.stellar.mag.O$ColorIndex),
                     method = "kendall", 
                     alternative = "greater")
X.stellar.mag.O.cor.test.2<-cor.test(as.numeric(X.stellar.mag.O$Magnitude), 
                     as.numeric(X.stellar.mag.O$ColorIndex),
                     method = "pearson", 
                     alternative = "greater")
X.stellar.mag.O.cor.test.3<-cor.test(as.numeric(X.stellar.mag.O$Magnitude), 
                     as.numeric(X.stellar.mag.O$ColorIndex),
                     method = "spearman", 
                     alternative = "greater")
ShapiroFranciaTest(as.numeric(X.stellar.mag.O$Magnitude))
ShapiroFranciaTest(as.numeric(X.stellar.mag.O$ColorIndex))
#---------------------------------Pattern Analysis--------------------------------------
X.stellar.mag.O.1<- pp3(as.numeric(X.stellar.mag.O$x), 
                        as.numeric(X.stellar.mag.O$y),
                        as.numeric(X.stellar.mag.O$z), 
                        box.1)
summary(X.stellar.mag.O.1)

X.stellar.mag.O.1.x_c <- cut(X.stellar.mag.O.1$data$x, 50)
X.stellar.mag.O.1.y_c <- cut(X.stellar.mag.O.1$data$y, 50)
X.stellar.mag.O.1.z <- table(X.stellar.mag.O.1.x_c, X.stellar.mag.O.1.y_c)

window.1=ellipse(a=50000,b=50000,centre=c(500,500),phi=pi/1)
plot(window.1)
X.stellar.mag.O.1.stats<-allstats(ppp(X.stellar.mag.O.1$x,X.stellar.mag.O.1$y,window=window.1))
X.stellar.mag.O.1.pattern.1<-ppp(X.stellar.mag.O.1$x,X.stellar.mag.O.1$y,window=window.1)

#--------------------------------Tables--------------------------------------------------

#------------------Figures for Scientific Visualization-------------------------------------------------
Figure.1<-plot( (as.numeric(X.stellar.mag.O$Magnitude) < median(as.numeric(X.stellar.mag.O$Magnitude))) - 0.5, 
                type="s", 
                ylim=c(-1,1),ylab="X < median(x)-0.5" )
opa <- par(mfrow=c(2,2),mar=c(2,2,2,2))

Figure.2<-PlotViolin(as.numeric(X.stellar.mag.O$Magnitude), col = "purple",main="A")

Figure.3<-PlotViolin(as.numeric(X.stellar.mag.O$x), col = "green",na.rm=TRUE,main="B")

Figure.4<-PlotViolin(as.numeric(X.stellar.mag.O$y), col = "red",na.rm=TRUE,main="C")

Figure.5<-PlotViolin(as.numeric(X.stellar.mag.O$z), col = "blue",na.rm=TRUE,main="D")
#---------------------Figure Groups----------------------------------
par(opa)
testpos <- seq(0, 1.98*pi, length=270)

Figure.6<-PlotPolar(as.numeric(X.stellar.mag.O$Magnitude), 
          testpos, type="l", 
          main="Stellar Magnitude O Spectral Type", col="blue")

Figure.7<-PolarGrid(ntheta=9, col="grey", lty="solid", lblradians=TRUE)
#---------------------Figure Groups----------------------------------
opa <- par(mfrow=c(2,2),mar=c(2,2,2,2))
Figure.8<-hist(X.stellar.mag.O.1$data$x)
Figure.9<-hist(X.stellar.mag.O.1$data$y)
Figure.10<-hist(X.stellar.mag.O.1$data$z)
Figure.11<-hist3D(z=z, border="black")
#-----------------------Figure Groups-----------------------------
par(opa)
opa <- par(mfrow=c(2,2),mar=c(1,1,1,1))
Figure.12<-plot(X.stellar.mag.O.1, theta=45, phi=45,main="Perspective 1",
     type="p",
     box.back=list(col="pink"),
     box.front=list(col="blue", lwd=1),cex=0.75, pch=1)
Figure.13<-plot(X.stellar.mag.O.1, theta=90, phi=45,main="Perspective 2",
     type="p",
     box.back=list(col="pink"),
     box.front=list(col="blue", lwd=1),cex=0.25, pch=1)
Figure.14<-plot(X.stellar.mag.O.1, theta=0, phi=90,main="Perspective 3",
     type="p",
     box.back=list(col="pink"),
     box.front=list(col="blue", lwd=1),cex=0.25, pch=1)
Figure.15<-plot(X.stellar.mag.O.1, theta=25, phi=15,main="Perspective 4",
     type="p",
     box.back=list(col="pink"),
     box.front=list(col="blue", lwd=1),cex=0.25, pch=1)
#------------------------------Figure Groups------------------------
par(opa)
opa <- par(mfrow=c(2,2),mar=c(2,2,2,2))
Figure.16<-plot(X.stellar.mag.O.1.stats$`K function`, col=1,main="K Function")
Figure.17<-plot(X.stellar.mag.O.1.stats$`F function`, col=2,main="F Function")
Figure.18<-plot(X.stellar.mag.O.1.stats$`G function`, col=3,main="G Function")
Figure.19<-plot(X.stellar.mag.O.1.stats$`J function`, col=4,main="J Function")

#----------------------------Reference-------------------------------------------



#----------------------------Function Library------------------------------------



