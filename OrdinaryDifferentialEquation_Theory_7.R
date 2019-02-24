#-----------------------------R Code To Modify in the Classroom Lecture with Students-----------------------
#----------------------------------R API --------------------------------------------------
library(xtable); library(deSolve); library(igraph); library(Matrix); library(visNetwork);library(corrplot);library(plot3D);library(scatterplot3d)
library(rgl)
#------------------------------------Data--------------------------------------------------
parms.names<-c("forward transition rate",
               "backward transition rate",
               "attachment rate",
               "detachment rate")
parms<-c(l.0=1,g.0=0.5,a.1=0.2,b.1=0.2,
         l.1=1.25,g.1=0.5,a.2=0.2,b.2=0.2,
         l.2=1,g.2=0.5,a.3=0.2,b.3=0.2,
         l.3=1.25,g.3=0.5,a.4=0.2,b.4=0.2,
         l.4=1.25,g.4=0.5,a.5=0.2,b.5=0.2,
         l.5=1.25,g.5=0.5,a.6=0.2,b.6=0.2,
         l.6=1.25,g.6=0.5)
param.functionals<-function(l,g,t)
{
  parm.l<-NULL
  parm.g<-NULL
  parm.a<-NULL
  parm.b<-NULL
  for(i in 1:6)
  {
    parm.l[i]<-1+sin(pi*t/4)
    parm.g[i]<-sin(pi*t/2)
    parm.a[i]<-sin(pi*t/2)
    parm.b[i]<-1+sin(pi*t/4)
  }
  parms.periodic<-c(l.0=l,g.0=g,a.1=parm.a[1],b.1=parm.b[1],
           l.1=parm.l[1],g.1=parm.g[1],a.2=parm.a[2],b.2=parm.b[2],
           l.2=parm.l[2],g.2=parm.g[2],a.3=parm.a[3],b.3=parm.b[3],
           l.3=parm.l[3],g.3=parm.g[3],a.4=parm.a[4],b.4=parm.b[4],
           l.4=parm.l[4],g.4=parm.g[4],a.5=parm.a[5],b.5=parm.b[5],
           l.5=parm.l[5],g.5=parm.g[5],a.6=parm.a[6],b.6=parm.b[6],
           l.6=parm.l[6],g.6=parm.g[6])
  return(parms.periodic)
}

#--------------------Nonlinear Tridiagonal mRNA translation Ribosome Model----------------------------

Ribosome.model.1 <- function(t, x, parms, input)  {
  with(as.list(c(parms, x)), {
    import <- input(t)
    dx1 <- l.0*(1-x1)+g.1*x2*(1-x1)+b.1*(1-x1)-l.1*x1*(1-x2)-g.0*x1-a.1*x1
    dx2 <- l.1*x1*(1-x2)+g.2*x3*(1-x2)+b.2*(1-x2)-l.2*x2*(1-x3)-g.1*x2*(1-x1)-a.2*x2     
    dx3 <- l.2*x2*(1-x3)+g.3*x4*(1-x3)+b.3*(1-x3)-l.3*x3*(1-x4)-g.2*x3*(1-x2)-a.3*x3
    dx4 <- l.3*x3*(1-x4)+g.4*x5*(1-x4)+b.4*(1-x4)-l.4*x4*(1-x5)-g.3*x4*(1-x3)-a.4*x4
    dx5 <- l.4*x4*(1-x5)+g.5*x6*(1-x5)+b.5*(1-x5)-l.5*x5*(1-x6)-g.4*x5*(1-x4)-a.5*x5
    dx6 <- l.5*x5*(1-x6)+g.6*(1-x6)+b.6*(1-x6)-l.6*x6-g.5*x6*(1-x5)-a.6*x6
    res <- c(dx1, dx2, dx3,dx4,dx5,dx6)
    list(res=res)
  })
}

#------------------Periodic function for x2------------------------------------------
Ribosome.model.2 <- function(t, x, parms, input)  {
  with(as.list(c(parms, x)), {
    import <- input(t)
    dx1 <- l.0*(1-x1)+g.1*x2*(1-x1)+b.1*(1-x1)-l.1*x1*(1-x2)-g.0*x1-a.1*x1
    dx2 <- import*x1*(1-x2)+g.2*x3*(1-x2)+b.2*(1-x2)-l.2*x2*(1-x3)-g.1*x2*(1-x1)-a.2*x2     
    dx3 <- l.2*x2*(1-x3)+g.3*x4*(1-x3)+b.3*(1-x3)-l.3*x3*(1-x4)-g.2*x3*(1-x2)-a.3*x3
    dx4 <- l.3*x3*(1-x4)+g.4*x5*(1-x4)+b.4*(1-x4)-l.4*x4*(1-x5)-g.3*x4*(1-x3)-a.4*x4
    dx5 <- l.4*x4*(1-x5)+g.5*x6*(1-x5)+b.5*(1-x5)-l.5*x5*(1-x6)-g.4*x5*(1-x4)-a.5*x5
    dx6 <- l.5*x5*(1-x6)+g.6*(1-x6)+b.6*(1-x6)-l.6*x6-g.5*x6*(1-x5)-a.6*x6
    res <- c(dx1, dx2, dx3,dx4,dx5,dx6)
    list(res=res)
  })
}
#----------------External signal with rectangle impulse and periodic function-----------------------
#
sequence <- seq(0, 10, 0.1)
signal <- data.frame(times = sequence,import = rep(0, length(sequence)))
signal.2<-2^(-1)+2^(-1)*sin((2^3)*pi*sequence/length(sequence))
signal.2.df <- data.frame(times = sequence,import = signal.2)
interval.impulse.time<-50
interval.impulse.value<-0.5
signal$import[signal$times >= interval.impulse.time & signal$times <= interval.impulse.time+1] <-interval.impulse.value 
sigimp <- approxfun(signal$times, signal$import, rule = 2)
sigimp.2 <- approxfun(signal.2.df$times, signal.2.df$import, rule = 1)

#-----------------Specification of the Spline functional for the Periodic function----------

periodic.function.1<- sin((sequence-0.5)*pi)
periodic.function.2<-2^(-1)+2^(-1)*sin((2^3)*pi*sequence/length(sequence))

periodic.function.1.spline<- splinefun(sequence, periodic.function.1)
periodic.function.2.spline<-splinefun(sequence, periodic.function.2)

#---------------------------------System Solution-------------------------------------------
xstart <- c(x1 = 1, x2 = 1, x3 = 1,x4=1,x5=1,x6=1)

Ribosome.system.solution.1 <- ode(y = xstart, times = sequence,
                             func = Ribosome.model.1, 
                             parms = parms, 
                             input = sigimp)
Ribosome.system.solution.1


Ribosome.system.solution.2 <- ode(y = xstart, times = sequence,
                                  func = Ribosome.model.2, 
                                  parms = parms, 
                                  input = periodic.function.1.spline)
Ribosome.system.solution.2


Ribosome.system.solution.3 <- ode(y = xstart, times = sequence,
                                  func = Ribosome.model.2, 
                                  parms = param.functionals(1.25,0.5,c(1,1,1,1,1,1)), 
                                  input = periodic.function.2.spline)
Ribosome.system.solution.3


parms.df<-as.data.frame(parms)

Rate.output<-(parms[25]+parms[23])*Ribosome.system.solution.1[,7] -(parms[22]+parms[24])*(1-Ribosome.system.solution.1[,7])
#------------------------------------Tables-----------------------------------------------
parameter.table.df<-t(as.data.frame(parms))
parameter.table.df<-cbind(parameter.table.df, t(as.data.frame(xstart)))
row.names(parameter.table.df)<-c("Parameter Value 1")

Table.1<-xtable(parameter.table.df)

distribution.moments.1.df<-data.frame()
distribution.moments.1.df<-rbind(c("ES_1","X1",empMoments.mutation(Ribosome.system.solution.1[,2])),
                                 c("ES_1","X2",empMoments.mutation(Ribosome.system.solution.1[,3])),
                                 c("ES_1","X4",empMoments.mutation(Ribosome.system.solution.1[,4])),
                                 c("ES_1","X5",empMoments.mutation(Ribosome.system.solution.1[,5])),
                                 c("ES_1","X6",empMoments.mutation(Ribosome.system.solution.1[,6])))
colnames(distribution.moments.1.df)<-c("Equation","Variable","Moment 1", "Moment 2","Moment 3","Moment 4","theta")

Table.2<-xtable(distribution.moments.1.df)

#--------------------------Figures for Classroom Presentation--------------------------------------------------------

Figure.1<-plot(periodic.function.1)
Figure.1A<-plot(periodic.function.2)


Figure.1<-plot(Ribosome.system.solution.2[,7], type="l", 
               lty=1,col="black", ylim=c(0,1),
               ylab="Simulated Value", 
               xlab="Temporal Position")
lines(Ribosome.system.solution.2[,2], lty=2,col="green")
lines(Ribosome.system.solution.2[,3], lty=3,col="red")
lines(Ribosome.system.solution.2[,4], lty=4,col="blue")
lines(Ribosome.system.solution.2[,5], lty=5,col="yellow")
lines(Ribosome.system.solution.2[,6], lty=6,col="purple")
rug(Ribosome.system.solution.2[,7], side=4, col="black")
rug(Ribosome.system.solution.2[,2], side=4, col="green")
rug(Ribosome.system.solution.2[,3], side=4, col="red")
rug(Ribosome.system.solution.2[,4], side=4, col="blue")
rug(Ribosome.system.solution.2[,5], side=4, col="yellow")
rug(Ribosome.system.solution.2[,6], side=4, col="purple")
legend("bottomleft",
       c("x6",
         "x1",
         "x2",
         "x3",
         "x4",
         "x5"),
       inset = .01,
       col=c("black","green","red","blue","yellow","purple"),
       lwd=2,lty=1:6,
       cex=0.5)

Figure.2<-hist(Ribosome.system.solution.2, mfrow = c(2, 3))
mtext(outer = TRUE, side = 3, "Species Histograms", cex = 1.5)

Figure.2A<-plot(Rate.output, type="l", 
               lty=1,col="black", 
               ylab="Simulated Output", 
               xlab="Temporal Position")
rug(Rate.output, side=4, col="black")
legend("bottomleft",
       c("Output"),
       inset = .01,
       col=c("black"),
       lwd=2,
       cex=0.5)

Figure.3<-plot(Ribosome.system.solution.1[,2],Ribosome.system.solution.1[,7])

cc <- colors()
Figure.4<-scatterplot3d(Ribosome.system.solution.2[,7], 
                        Ribosome.system.solution.2[,2], 
                        Ribosome.system.solution.2[,4], highlight.3d = TRUE, col.axis = "blue",
                        col.grid = "lightblue", main = "Variable Solutions", pch = 20,angle=120,
                        xlab = "x(t)", 
                        ylab = "y(t)", zlab = "z(t)")

rgl.viewpoint(0, 20)
Figure.5<-plot3d(Ribosome.system.solution.2[,7], 
                 Ribosome.system.solution.2[,2], 
                 Ribosome.system.solution.2[,4], 
                 pch = 1, cex = 1, xlab = "x(t)", 
                 ylab = "y(t)", zlab = "z(t)")


#----------------------------------References---------------------------------------------

Reference.1<-c("Yoram Zarai, Michael Margaliot, Tamir Tuller (2017)",
               "A deterministic mathematical model for bidirectional excluded flow with Langmuir kinetics",
               "PLos ONE 12(8): e0182178. https:// doi.org/10.1371/journal.pone.0182178")

#----------------------------------Function Library---------------------------------------
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
