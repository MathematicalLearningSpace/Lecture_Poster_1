#-----------------------------R Code To Modify in the Classroom Lecture with Students-----------------------
#-----------------------------R API ------------------------------------------------------------------------
#-------------------------------------Ordinary Differential Equations-----------
library(deSolve);library(ReacTran);library(rootSolve);library(fda);library(phaseR)
#---------------------------TeX Formatting----------------------------------
library(xtable)
library(tseriesChaos)
#----------------------------Scientific  Visualization-----------------------
library(corrplot)
library(plot3D)
library(scatterplot3d)
library(rgl)
#---------------Reusable Ordinary Differential Equation Model------------------
#-------------- Three Variable Specification------------------------------------

BCC.model.1 <- function(t, x, parms, input)  {
  with(as.list(c(parms, x)), {
    import <- input(t)
    dx1 <- a1*x1*x2 + a2*x3 + a3*import
    dx2 <- a4*x2*x1  - a5*x2*x3          
    dx3 <- a6*x2*x3  - a7*x1*x3-a8*x3            
    res <- c(dx1, dx2, dx3)
    list(res)
  })
}
BCC.model.2 <- function(t, x, parms)  {
  with(as.list(c(parms, x)), {
    dx1 <- a1*x1*x2 + a2*x3 + a3
    dx2 <- a4*x2*x1  - a5*x2*x3          
    dx3 <- a6*x2*x3  - a7*x3            
    res <- c(dx1, dx2, dx3)
    list(res)
  })
}
#----------------Parameter Value Examples-------------------------------
parms <- c(a1 = 0.0001, 
           a2 = 0.1, 
           a3 = 0.0, 
           a4 = 0.1, 
           a5 = 0.1, 
           a6 = 0.1,
           a7= 0.1,
           a8= 0.5)
## vector of timesteps
sequence <- seq(0, 10, 0.1)
#--------------------External signal with rectangle impulse for noise generation-------------
signal <- data.frame(times = sequence,
                     import = rep(0, length(sequence)))
interval.impulse.time<-50
interval.impulse.value<-0.5
signal$import[signal$times >= interval.impulse.time & signal$times <= interval.impulse.time+1] <-interval.impulse.value 

sigimp <- approxfun(signal$times, signal$import, rule = 2)
#
#----------------------Start values for steady state
xstart <- c(x1 = 1, x2 = 1, x3 = 1)
#---------------------------------------Estimation-----------------------------------------------------
BCC.system.solution.1 <- ode(y = xstart, times = sequence,
                              func = BCC.model.1, 
                              parms = parms, 
                              input = sigimp)

#------------------------Qualitative Analysis in the Classroom-------------------------------------------
#------------------------Example of Classroom Tables---------------------------------------------------------

parameter.table.df<-t(as.data.frame(parms))
parameter.table.df<-cbind(parameter.table.df, t(as.data.frame(xstart)))
row.names(parameter.table.df)<-c("Parameter Value 1")
Table.1<-xtable(parameter.table.df)
#-------------Example of Empirical Moment Estimation from Function Library---------------------
distribution.moments.1.df<-data.frame()
distribution.moments.1.df<-rbind(c("ES_1","X1",empMoments.mutation(BCC.system.solution.1[,2])),
                                 c("ES_1","X2",empMoments.mutation(BCC.system.solution.1[,3])),
                                 c("ES_1","X3",empMoments.mutation(BCC.system.solution.1[,4])))
colnames(distribution.moments.1.df)<-c("Equation","Variable","Moment 1", "Moment 2","Moment 3","Moment 4","theta")
Table.2<-xtable(distribution.moments.1.df)
#-----------------------Examples of Classroom Figures--------------------------------------------------------
Figure.1<-plot(BCC.system.solution.1[,2], type="l", 
                lty=1,col="black", 
                ylab="Simulated Value", 
                xlab="Temporal Position")
lines(BCC.system.solution.1[,3], lty=2,col="green")
lines(BCC.system.solution.1[,4], lty=3,col="blue")
rug(BCC.system.solution.1[,2], side=4, col="black")
rug(BCC.system.solution.1[,3], side=4, col="green")
rug(BCC.system.solution.1[,4], side=4, col="blue")
legend("bottomleft",
       c("x1",
         "x2",
         "x3"),
       inset = .01,
       col=c("black","green","blue"),
       lwd=2,
       cex=0.5)
Figure.2<-hist(BCC.system.solution.1, mfrow = c(2, 3))
mtext(outer = TRUE, side = 3, "Species Histograms", cex = 1.5)
cc <- colors()
Figure.3<-scatterplot3d(BCC.system.solution.1[,2], 
              BCC.system.solution.1[,3], 
              BCC.system.solution.1[,4], highlight.3d = TRUE, col.axis = "blue",
              col.grid = "lightblue", main = "Species Trajectories", pch = 20,angle=120,
              xlab = "x(t)", 
              ylab = "y(t)", zlab = "z(t)")
rgl.viewpoint(0, 20)
Figure.4<-plot3d(BCC.system.solution.1[,2], 
       BCC.system.solution.1[,3], 
       BCC.system.solution.1[,4], 
       pch = 1, cex = 1, xlab = "x(t)", 
       ylab = "y(t)", zlab = "z(t)")


rgl.snapshot("BCC_Model_1.png")
#--------------------------------------Function Library-----------------------------------------------
