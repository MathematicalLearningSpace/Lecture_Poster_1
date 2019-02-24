library(deSolve);library(ReacTran);library(rootSolve);library(fda);library(xtable);library(CHNOSZ);library(changepoint)

#------Solution of Delayed Differential Equation System for Classroom-------------
y.init<-0.1
times <- seq(from = 0, to = 100, by = 0.1) 
solution.1.mg.test<- dede(y = y.init, 
                          times = times, 
                          func =function.mg.test, 
                          parms = c(0.2,9,0.1), 
                          tau = 5)
solution.2.mg.test<- dede(y = y.init, 
                          times = times, 
                          func = function.mg.test, 
                          parms = c(0.2,10,0.1), 
                          tau = 10)
solution.3.mg.test<- dede(y = y.init, 
                          times = times, 
                          func = function.mg.test, 
                          parms = c(0.2,8,0.1), 
                          tau = 10)
solution.4.mg.test<- dede(y = y.init, 
                          times = times, 
                          func = function.mg.test, 
                          parms = c(0.2,6,0.1), 
                          tau = 10)
#----------------------Changepoint Analysis-----------------------------------------
y<-c(solution.1.mg.test[,-1],solution.2.mg.test[,-1])
ansmean=cpt.mean(y)
ansvar=cpt.var(y)
ansmeanvar=cpt.meanvar(y)
#----------------------Tables for Classroom--------------------------------------------------------
#----------------------Figures for Classroom------------------------------------------------------
Figure.1<-plot(solution.1.mg.test, 
               lwd = 2, 
               main = "MG Solution Trajectories", 
               ylab = "y", 
               mfrow = c(2, 2), which = 1) 
lines(solution.2.mg.test, type = "l", lwd = 2, xlab = "y",lty=2,col="red")
lines(solution.3.mg.test, type = "l", lwd = 2, xlab = "y",lty=3,col="green")
lines(solution.4.mg.test, type = "l", lwd = 2, xlab = "y",lty=4,col="blue")
grid()
legend("bottomright", col = c("black", "red", "green","blue"), 
       lty = 1:4, cex=0.75,
       legend = c("MG-9","MG-10", "MG-8", "MG-6"))

Figure.2<-plot(solution.1.mg.test[,-1], type = "l", lwd = 2, xlab = "y") 
Figure.3<-plot(solution.2.mg.test, lwd = 2, main = "tau=20", ylab = "y", 
               mfrow = NULL, which = 1) 
Figure.4<-plot(solution.2.mg.test[,-1], type = "l", lwd = 2, xlab = "y")

Figure.5<-plot(solution.1.mg.test[,-1], type = "l", lwd = 2, xlab = "y", lty=1,col="black", xlim=c(0,1.5),ylim=c(0,1.5))
lines(solution.2.mg.test[,-1], type = "l", lwd = 2, xlab = "y",lty=2,col="red")
lines(solution.3.mg.test[,-1], type = "l", lwd = 2, xlab = "y",lty=3,col="green")
lines(solution.4.mg.test[,-1], type = "l", lwd = 2, xlab = "y",lty=4,col="blue")
grid()
legend("bottomleft", col = c("black", "red", "green","blue"), 
       lty = 1:4, cex=0.75,
       legend = c("MG-9","MG-10", "MG-8", "MG-6"))

par(mfrow = c(2,1))
Figure.6<-plot(ansmean,cpt.col='blue')
Figure.7<-plot(ansvar)
Figure.8<-plot(ansmeanvar,cpt.width=3)
#--------------------------------Function Library for Classroom---------------------------------------------

#---------------------------------Example of Mackey-Glass Equation--------------------------------------------------
#-------------Leon Glass and Michael Mackey (2010) Mackey-Glass equation. Scholarpedia, 5(3):6908.
function.mg.test <-function(t,y,parms, tau)
{
  tlag<-t-tau
  if(tlag <=0)
    ylag<-0.5
  else
    ylag <-lagvalue(tlag)
  dy<-parms[1]*ylag * 1/(1+ylag^parms[2])-parms[3]*y
  list(dy=dy,ylag=ylag)
}


