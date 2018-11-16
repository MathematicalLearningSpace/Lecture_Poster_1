#--------------------------------R API ------------------------------------
library(xtable); library(fitdistrplus); library(PearsonDS);library(lmom);library(pracma)
library(ExtDist)
#-------------------Data--------------------------------------------------
#-------------------Geometric Stable Distributions-------------------------
sequence<-seq(1,100,1)
Algorithm.1.GS.1<-Algorithm.1.GS(sequence,100,0.25)
Algorithm.1.GS.2<-Algorithm.1.GS(sequence,100,0.5)
Algorithm.1.GS.3<-Algorithm.1.GS(sequence,100,0.75)
Algorithm.1.GS.4<-Algorithm.1.GS(sequence,100,0.95)


Algorithm.2.GS.strict.1<-Algorithm.2.GS.strict(100,0.25,0.5,1)
Algorithm.2.GS.strict.2<-Algorithm.2.GS.strict(100,0.5,0.5,1)
Algorithm.2.GS.strict.3<-Algorithm.2.GS.strict(100,0.75,0.5,1)
Algorithm.2.GS.strict.4<-Algorithm.2.GS.strict(100,0.95,0.5,1)


Algorithm.3.GS.Mittag_Leffler.1<-Algorithm.3.GS.Mittag_Leffler(100,1,0.25,1)
Algorithm.3.GS.Mittag_Leffler.2<-Algorithm.3.GS.Mittag_Leffler(100,1,0.5,1)
Algorithm.3.GS.Mittag_Leffler.3<-Algorithm.3.GS.Mittag_Leffler(100,1,0.75,1)
Algorithm.3.GS.Mittag_Leffler.4<-Algorithm.3.GS.Mittag_Leffler(100,1,0.95,1)


Algorithm.4.GS.beta.1<-Algorithm.4.GS.beta(100,1,0.25)
Algorithm.4.GS.beta.2<-Algorithm.4.GS.beta(100,1,0.5)
Algorithm.4.GS.beta.3<-Algorithm.4.GS.beta(100,1,0.75)
Algorithm.4.GS.beta.4<-Algorithm.4.GS.beta(100,1,0.95)


Moments.Algorithms.df<-data.frame()
Moments.Algorithms.df<-rbind(c("0.25",empMoments.mutation(Algorithm.4.GS.beta.1)),
                                   c("0.5",empMoments.mutation(Algorithm.4.GS.beta.2)),
                                   c("0.75",empMoments.mutation(Algorithm.4.GS.beta.3)),
                                   c("0.95",empMoments.mutation(Algorithm.4.GS.beta.4))
)
colnames(Moments.Algorithms.df)<-c("Index Rho","Moment 1","Moment 2","Moment 3","Moment 4","Theta")

#-------------------Tables------------------------------------------------
Table.1<-xtable(Moments.Algorithms.df)
#------Figures for Presentation in the Classroom------------------------------------------------

Figure.1<-plot(Algorithm.1.GS.1, type="l", 
               lty=1,col="black", 
               ylab="Simulated Value", 
               xlab="Temporal Position")
lines(Algorithm.1.GS.2, lty=2,col="green")
lines(Algorithm.1.GS.3, lty=3,col="blue")
lines(Algorithm.1.GS.4, lty=4,col="red")

rug(Algorithm.1.GS.1, side=4, col="black")
rug(Algorithm.1.GS.2, side=4, col="green")
rug(Algorithm.1.GS.3, side=4, col="blue")
rug(Algorithm.1.GS.4, side=4, col="red")

legend("bottomleft",
       c("rho=0.25",
         "rho=0.5",
         "rho=0.75",
         "rho=0.95"),
       inset = .01,
       col=c("black","green","blue","red"),
       lwd=2,lty=1:4,
       cex=0.5)

Figure.2<-plot(Algorithm.2.GS.strict.1, type="l", 
               lty=1,col="black", 
               ylab="Simulated Value", 
               xlab="Temporal Position")
lines(Algorithm.2.GS.strict.2, lty=2,col="green")
lines(Algorithm.2.GS.strict.3, lty=3,col="blue")
lines(Algorithm.2.GS.strict.4, lty=4,col="red")

rug(Algorithm.2.GS.strict.1, side=4, col="black")
rug(Algorithm.2.GS.strict.2, side=4, col="green")
rug(Algorithm.2.GS.strict.3, side=4, col="blue")
rug(Algorithm.2.GS.strict.4, side=4, col="red")

legend("bottomleft",
       c("alpha=0.25",
         "alpha=0.5",
         "alpha=0.75",
         "alpha=0.95"),
       inset = .01,
       col=c("black","green","blue","red"),
       lwd=2,lty=1:4,
       cex=0.5)

Figure.3<-plot(Algorithm.3.GS.Mittag_Leffler.1, type="l", 
               lty=1,col="black", 
               ylab="Simulated Value", 
               xlab="Temporal Position")
lines(Algorithm.3.GS.Mittag_Leffler.2, lty=2,col="green")
lines(Algorithm.3.GS.Mittag_Leffler.3, lty=3,col="blue")
lines(Algorithm.3.GS.Mittag_Leffler.4, lty=4,col="red")

rug(Algorithm.3.GS.Mittag_Leffler.1, side=4, col="black")
rug(Algorithm.3.GS.Mittag_Leffler.2, side=4, col="green")
rug(Algorithm.3.GS.Mittag_Leffler.3, side=4, col="blue")
rug(Algorithm.3.GS.Mittag_Leffler.4, side=4, col="red")

legend("bottomleft",
       c("lambda=0.25",
         "lambda=0.5",
         "lambda=0.75",
         "lambda=0.95"),
       inset = .01,
       col=c("black","green","blue","red"),
       lwd=2,lty=1:4,
       cex=0.5)

Figure.4<-plot(Algorithm.4.GS.beta.1, type="l", 
               lty=1,col="black", 
               ylab="Simulated Value", 
               xlab="Temporal Position")
lines(Algorithm.4.GS.beta.2, lty=2,col="green")
lines(Algorithm.4.GS.beta.3, lty=3,col="blue")
lines(Algorithm.4.GS.beta.4, lty=4,col="red")

rug(Algorithm.4.GS.beta.1, side=4, col="black")
rug(Algorithm.4.GS.beta.2, side=4, col="green")
rug(Algorithm.4.GS.beta.3, side=4, col="blue")
rug(Algorithm.4.GS.beta.4, side=4, col="red")

legend("bottomleft",
       c("beta=0.25",
         "beta=0.5",
         "beta=0.75",
         "beta=0.95"),
       inset = .01,
       col=c("black","green","blue","red"),
       lwd=2,lty=1:4,
       cex=0.5)

#------Reference Pattern for Classroom Presentation-----------------------------------------------
Reference.1<-c("Tomasz J. Kozubowski",
               "Computer simulation of geometric stable distributions",
               "Journal of Computational and Applied Mathematics 116 (2000) 221-229")

#---------------Function Library-----------------------------------------
Algorithm.1.GS<-function(x,N, rho)
{
  y<-NULL
  for(i in 1:N)
  {
    y[i]<-sin(pi * rho)/rho*pi*(x[i]^2+2*x[i]*cos(pi*rho)+1)
  }
  return(y)
}

sequence<-runif(100)
test.Algorithm.1.GS<-Algorithm.1.GS(sequence,100,0.5)
test.Algorithm.1.GS

Algorithm.2.GS.strict<-function(N,alpha,lambda,tau)
{
  y<-NULL;Z<-NULL;W<-NULL; U<-NULL
  if(abs(tau)<= min(1,(2/alpha)-1))
     {
  p<-(1+tau)/2
  #Z<-quaexp(runif(N), c(0,1))
  Z<-rExp(N, scale=1)
  for(i in 1:N)
  {
  U[i]<-runif(1)
  if(U[i] < p)
  {
    rho<-alpha*p
    I<-1
  }
  else{
    rho<-alpha*(1-p)
    I<-1
  }
  if (rho==1)
  {
    W<-1
  }
  else
  {
    V<-runif(1)
    W[i]<-sin(pi*rho)*cot(pi*rho*V)-cos(pi*rho)
  }
  
  y[i]<-I*Z[i]*(lambda*W[i])^(1/alpha)
  }
  }
  return(y)
}

test.Algorithm.2.GS.strict<-Algorithm.2.GS.strict(100,0.5,0.5,0.5)
test.Algorithm.2.GS.strict

Algorithm.3.GS.Mittag_Leffler<-function(N,alpha,lambda,rho)
{
  y<-NULL;Z<-NULL;V<-NULL;W<-NULL
  #Z<-quaexp(runif(N), c(0,1))
  Z<-rExp(N, scale=1)
  V<-runif(N)
  for(i in 1:N)
  {
    W[i]<-sin(pi*rho)*cot(pi*rho*V[i])-cos(pi*rho)
    y[i]<-Z[i]*(lambda*W[i])^(1/alpha)
  }
  return(y)
}

test.Algorithm.3.GS.Mittag_Leffler<-Algorithm.3.GS.Mittag_Leffler(100,0.5,0.5,0.5)
test.Algorithm.3.GS.Mittag_Leffler

Algorithm.4.GS.beta<-function(N,alpha,beta)
{
  y<-NULL
  #Z<-quaexp(runif(N), c(0,1))
  Z<-rExp(N, scale=1)
  V<-runif(N,-pi/2,pi/2)
  for(i in 1:N)
  {
  if(alpha==1)
  {
    y[i]<-(2/pi)*((pi/2) + beta *V[i])*tan(V[i])-beta*log(Z[i]*cos(V[i])/(pi/2+beta*V[i]))
  }
  else{
    gamma<-atan(beta*tan(pi*alpha/2))/alpha
    delta<-(1+beta^2*tan((pi*alpha))/2)^(1/(2*alpha))
    y[i]<-delta*sin(alpha*(V[i]+gamma))/(cos(V[i])^(1/alpha))*
      (cos(V[i]-alpha*(V[i]+gamma))/Z[i])^((1-alpha)/alpha)
  }
  }
  return(y)
}

test.Algorithm.4.GS.beta<-Algorithm.4.GS.beta(100,0.75,0.5)
test.Algorithm.4.GS.beta
#--------------------------------Function Library for Classroom Modification by Students-----------

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
  
