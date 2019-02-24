#-----------------------------R Code To Modify in the Classroom Lecture with Students-----------------------
#----------------------------------R API -------------------------------------------------------------------
library(NMOF);library(xtable)
#------------------------------------Data----------------------------------------------------------
N=100
y<-mRN(N,1)
X<-mRN(N,2)
h <- (N+2) %/% 2
data.simulated<-list(y=y,X=X,h=h)
#-----------------------------------Particle Swarm Optimization------------------------------------
PSO.options<-list(min = rep(-5, 2), max = rep( 5, 2),
                  c1 = 1.0, c2 = 2.0,
                  iner = 0.7, initV = 1, maxV = 3, 
                  nP = 100, nG = 30, loopOF = FALSE)
PSO.options.1<-list(min = rep(-5, 2), max = rep( 5, 2),
                  c1 = 1.0, c2 = 2.0,
                  iner = 0.7, initV = 1, maxV = 6, 
                  nP = 100, nG = 30, loopOF = FALSE)
PSO.options.2<-list(min = rep(-5, 2), max = rep( 5, 2),
                  c1 = 1.0, c2 = 2.0,
                  iner = 0.7, initV = 1, maxV = 9, 
                  nP = 100, nG = 30, loopOF = FALSE)
PSO.options.3<-list(min = rep(-5, 2), max = rep( 5, 2),
                  c1 = 1.0, c2 = 2.0,
                  iner = 0.7, initV = 1, maxV = 12, 
                  nP = 100, nG = 30, loopOF = FALSE)

PSO.algorithm.solution<-PSopt(OF=objective.function, algo=PSO.options,data=data.simulated)
PSO.algorithm.solution.1<-PSopt(OF=objective.function, algo=PSO.options.1,data=data.simulated)
PSO.algorithm.solution.2<-PSopt(OF=objective.function, algo=PSO.options.2,data=data.simulated)
PSO.algorithm.solution.3<-PSopt(OF=objective.function, algo=PSO.options.3,data=data.simulated)


PSO.algorithm.solution.accuracy.error <- sort((y - X %*% as.matrix(PSO.algorithm.solution$xbest))^2)[h]
PSO.algorithm.solution.accuracy.error.1 <- sort((y - X %*% as.matrix(PSO.algorithm.solution.1$xbest))^2)[h]
PSO.algorithm.solution.accuracy.error.2 <- sort((y - X %*% as.matrix(PSO.algorithm.solution.2$xbest))^2)[h]
PSO.algorithm.solution.accuracy.error.3 <- sort((y - X %*% as.matrix(PSO.algorithm.solution.3$xbest))^2)[h]


PSO.solution.df<-data.frame()
PSO.solution.df<-rbind(c("3",PSO.algorithm.solution.accuracy.error,PSO.algorithm.solution$xbest),
                       c("6",PSO.algorithm.solution.accuracy.error.1,PSO.algorithm.solution.1$xbest),
                       c("9",PSO.algorithm.solution.accuracy.error.2,PSO.algorithm.solution.2$xbest),
                       c("12",PSO.algorithm.solution.accuracy.error.3,PSO.algorithm.solution.3$xbest))
colnames(PSO.solution.df)<-c("Velocity","Error","Parameter 1","Parameter 2")

#------------------------------------Tables---------------------------------------------------------

Table.1<-xtable(PSO.solution.df)

#------------------------------------Figures--------------------------------------------------------

par(mfrow = c(1,1), mar = c(5,5,5,5))
Figure.1<-plot(seq(1,N,1),PSO.algorithm.solution$popF,xlab='Population',ylab='Objective Function Population Values',
               main='',col='black',type="l")
lines(PSO.algorithm.solution.1$popF,lty=2,col="red")
lines(PSO.algorithm.solution.2$popF,lty=3,col="green")
lines(PSO.algorithm.solution.3$popF,lty=4,col="blue")
legend("topright", legend=c("v=3","v=6","v=9","v=12")
       , bty = "n",lwd=2, 
       cex=0.75, lty=1:4,col=c("black","red","green","blue"))
Figure.1.caption<-c("Population Values by Velocity")

#------------------------------------Function Library-----------------------------------------------

mRN<-function(m,n) array(rnorm(m*N),dim=c(m,n))

objective.function<-function(param,data.simulated)
{
  X<-data.simulated$X
  y<-data.simulated$y
  accuracy<-as.vector(y)-X%*%param
  accuracy<-accuracy*accuracy
  accuracy<-apply(accuracy,2,sort,partial=data.simulated$h)
  return(accuracy[h,])
}




