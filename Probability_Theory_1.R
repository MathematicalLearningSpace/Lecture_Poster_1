#-----------------------------------R API -----------------------------------------------------------
library(PearsonDS); library(fitdistrplus); library(xtable)
#-------------------------------------Probability Theory: Pearson Distributions----------------------------------------------
pearson.N<-512
#-------------------------------------Moment Parameters for the Distributions-------------------------------------------------
p0pars <- list(mean=1, sd=1)
pIpars <- list(a=1, b=1, location=1, scale=1) 
pIIpars <- list(a=1, location=1, scale=1) 
pIIIpars <- list(shape=1, location=1, scale=1)
pIVpars <- list(m=1, nu=1, location=1, scale=1)
pVpars <- list(shape=1, location=1, scale=1) 
pVIpars <- list(a=1, b=1, location=1, scale=1) 
pVIIpars <- list(df=10, location=1, scale=1)

#-------------------------------------Generate Random variables from the Distributions----------------------------------------
error.pearson.0<-rpearson0(pearson.N,params=p0pars)
error.pearson.1<-rpearsonI(pearson.N,params=pIpars)
error.pearson.2<-rpearsonII(pearson.N,params=pIIpars)
error.pearson.3<-rpearsonIII(pearson.N,params=pIIIpars)
error.pearson.4<-rpearsonIV(pearson.N,params=pIVpars)
error.pearson.5<-rpearsonV(pearson.N,params=pVpars)
error.pearson.6<-rpearsonVI(pearson.N,params=pVIpars)
error.pearson.7<-rpearsonVII(pearson.N,params=pVIIpars)

#----------------Artificial Data Generation and Windows for Data Segmentation----------------------------------------------------
sample.pct<-1
K<-sample.pct*pearson.N
experimental.data.1<-window(error.pearson.1, end=256)
experimental.data.2<-window(error.pearson.2, end=256)
experimental.data.3<-window(error.pearson.3,end=256)

#------------------------------------Descriptive Statistics ---------------------------------------------------------------------
distribution.moments.1.df<-data.frame()
distribution.moments.1.df<-rbind(c("X1","1",empMoments.mutation(experimental.data.1)),
                                 c("X2","2",empMoments.mutation(experimental.data.2)),
                                 c("X3","3",empMoments.mutation(experimental.data.3))
)
colnames(distribution.moments.1.df)<-c("Variable","Type", "Moment 1", "Moment 2","Moment 3","Moment 4","Theta")
Table.1<-xtable(distribution.moments.1.df)

#------------------------------------Fit Pearson Family of Distributions-------------------------------------------------------
pearson.MLE.1 <- pearsonFitML(experimental.data.1)
pearson.MLE.2 <- pearsonFitML(experimental.data.2)
pearson.MLE.3 <- pearsonFitML(experimental.data.3)

pearson.MLE.dist<-c(pearson.MLE.1$type,
                    pearson.MLE.2$type,
                    pearson.MLE.3$type)

pearson.MOM.1 <-pearsonFitM(moments=empMoments(experimental.data.1))
pearson.MOM.2 <-pearsonFitM(moments=empMoments(experimental.data.2))
pearson.MOM.3 <-pearsonFitM(moments=empMoments(experimental.data.3))

pearson.MOM.dist<-c(pearson.MOM.1$type,
                    pearson.MOM.2$type,
                    pearson.MOM.3$type)

Data.1.CDF.Categories.1.df<-data.frame()
Data.1.CDF.Categories.1.df<-cbind(pearson.MLE.dist,
      pearson.MOM.dist)
row.names(Data.1.CDF.Categories.1.df)<-c("x1","x2","x3")
colnames(Data.1.CDF.Categories.1.df)<-c("MLE","Method Of Moments")
Table.2<-xtable(Data.1.CDF.Categories.1.df)
#-------------------------------------Figures for Presentation in the Classroom --------------------------------------
Figure.1<-plot(experimental.data.1, type="l", 
               lty=1,col="black", 
               ylab="Simulated Value", 
               xlab="Temporal Position")
lines(experimental.data.2, lty=2,col="green")
lines(experimental.data.3, lty=3,col="blue")
rug(experimental.data.1, side=4, col="black")
rug(experimental.data.2, side=4, col="green")
rug(experimental.data.3, side=4, col="blue")
legend("bottomleft",
       c("Type 1",
         "Type 2",
         "Type 3"),
       inset = .01,
       col=c("black","green","blue","red"),
       lwd=2,
       cex=0.5)
#-----------------------------------Function Library for Classroom Modification------------------------------------------------------------------------
empMoments.mutation<-function (x) 
{
  n <- length(x)
  moment.1 <- mean(x)
  moment.2  <- var(x) * (n - 1)/n
  if (moment.2 > 0) 
    moment.3 <- sum((x - moment.1)^3/sqrt(moment.2)^3)/n
  else moment.3 <- 0
  if (moment.2 > 0) 
    moment.4 <- sum((x - moment.1)^4/moment.2^2)/n
  else moment.4 <- 0
  if(moment.4> 0)
  theta.1<-moment.3/moment.4 
  else theta.1<-0
  c(average = moment.1, 
    variance = moment.2, 
    skewness = moment.3, 
    kurtosis = moment.4, 
    theta=theta.1)
}



