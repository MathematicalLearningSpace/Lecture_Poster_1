#----------------------------------R API ----------------------------------------------
library(readr);library(xtable);library(deSolve);library(ReacTran);library(rootSolve);library(fda);library(phaseR)
library(xtable);library(tseriesChaos);library(corrplot);library(plot3D);library(scatterplot3d)
library(tm);library(topicmodels);library(wordcloud);library(rgl);library(mRMRe)
#----------------------------Data---------------------------------------------------
article.files <- list.files(patt='publications_*.*csv$')
article.files
print(article.files[1])
#------------------------------publications_MYC.csv Example For Classroom ---------------------------------
z.df<-data.frame()
z.df<-read_csv(article.files[1])
View(z.df)
print(z.df$Title)
#--------ODE Model Pattern for MYC without Equation and Parameter Specification To be Completed in Classroom-----------------
Gene.ROI<-c('MYc', 'Miz1','p15','CKS1','Skp2','p27','CyclinD1','CDK4','CDK2','RB', 'E2F')

model.1 <- function(t, x, parms)  {
  with(as.list(c(parms, x)), {
    
    dx1 <-a11*x1 + a12*x1*x2 + runif(1,0,1)
    dx2 <-a21*x2 + a22*x2*x3      
    dx3 <-a31*x3 + a32*x3*x4
    
    dx4<-a41*x4 + a42*x4*x5
    dx5<-a51*x5 + a52*x5*x6
    dx6<-a61*x6 + a62*x6*x7
    
    dx7<-a71*x7 + a72*x7*x8
    dx8<-a81*x8 + a82*x8*x9
    dx9<-a91*x9 + a92*x9*x10
    
    dx10<-a101*x10 + a102*x10*x11
    dx11<-a111*x11 + a112*x11*x12
    dx12<-a121*x12 
    
    res <- c(dx1, dx2, dx3, dx4,dx5,dx6,dx7,dx8, dx9, dx10, dx11, dx12)
    list(res)
  })
}
model.2 <- function(t, x, parms)  {
  with(as.list(c(x)), {
    
    dx1 <-parms[1]*x1 + parms[2]*x1*x2 + runif(1,0,1)
    dx2 <-parms[3]*x2 + parms[4]*x2*x3      
    dx3 <-parms[5]*x3 + parms[6]*x3*x4
    
    dx4<-parms[7]*x4 + parms[8]*x4*x5
    dx5<-parms[9]*x5 + parms[10]*x5*x6
    dx6<-parms[11]*x6 + parms[12]*x6*x7
    
    dx7<-parms[13]*x7 + parms[14]*x7*x8
    dx8<-parms[15]*x8 + parms[16]*x8*x9
    dx9<-parms[17]*x9 + parms[18]*x9*x10
    
    dx10<-parms[19]*x10 + parms[20]*x10*x11
    dx11<-parms[21]*x11 + parms[22]*x11*x12
    dx12<-parms[23]*x12 
    
    res <- c(dx1, dx2, dx3, dx4,dx5,dx6,dx7,dx8, dx9, dx10, dx11, dx12)
    list(res)
  })
}
#------------------------------------- Parameter Values--------------------------
parms <- c(a11 = 0.01, a12 = -0.1, 
           a21 = 0.0, a22 = -0.1, 
           a31 = -0.1, a32 = 0.1, 
           a41= 0.1, a42= -0.5,
           a51 = 0.1, a52 = 0.1, 
           a61= 0.1, a62= -0.5,
           a71 = 0.1, a72 =- 0.25, 
           a81= -0.1, a82= -0.5,
           a91 = 0.1, a92 = 0.1, 
           a101= -0.1, a102=0.5,
           a111= 0.1, a112= 0.5,
           a121= -0.1)

parms.sim.unif<-NULL
parms.sim.norm<-NULL
parms.sim.gamma<-NULL
for(j in 1:length(parms))
{
  parms.sim.unif[j]<-runif(1,0,1)
  parms.sim.norm[j]<-rnorm(1,0,0.5)
  parms.sim.gamma[j]<-rgamma(1,2)
}
Table.Parameters.df<-data.frame()
Table.Parameters.df<-cbind(parms,parms.sim.unif,parms.sim.norm,parms.sim.gamma)
colnames(Table.Parameters.df)<-c("Point","Uniform",'Normal',"Gamma")

## vector of timesteps
sequence <- seq(0, 10, 0.1)
#---------------------Start values for steady state-----------------------------------
xstart <- c(x1 = 0.5, x2 = 1, x3 =0.5,x4=1,x5=1,
            x6=1,x7=1,x8=0.5,x9=1,x10=1,x11=1,x12=1)
#------------------------Estimation-----------------------------------------------------
system.solution.1 <- ode(y = xstart, times = sequence,
                             func = model.1, 
                             parms = parms)
system.solution.2<- ode(y = xstart, times = sequence,
                            func = model.2, 
                            parms = parms.sim.unif)
#---------------------------Tables----------------------------------------------------
Table.1<-xtable(z.df$Title)
Table.2<-xtable(Table.Parameters.df)
Table.2
#----------------------------Figures---------------------------------------------------
par(mfcol = c(1, 2))
Figure.1<-plot(system.solution.1[,2], type="l", 
               lty=1,col="black", ylim=c(0,10),
               ylab="Simulated Value", 
               xlab="Time")
for(i in 3:nrow(system.solution.1))
{
lines(system.solution.1[,i], lty=i)
}
legend("topright",
       c("x1","x2","x3","x4","x5","x6",
         "x7","x8","x9","x10","x11","x12"),
       inset = .01,
       lwd=2,
       cex=0.5)

Figure.2<-plot(system.solution.2[,2], type="l", 
               lty=1,col="black", ylim=c(0,10),
               ylab="Simulated Value", 
               xlab="Time")
for(i in 3:nrow(system.solution.2))
{
  lines(system.solution.2[,i], lty=i)
}
legend("topright",
       c("x1","x2","x3","x4","x5","x6",
         "x7","x8","x9","x10","x11","x12"),
       inset = .01,
       lwd=2,
       cex=0.5)

par(mfcol = c(1, 1))
Figure.2<-hist(system.solution.1, mfrow = c(3, 3))
mtext(outer = TRUE, side = 3, "Species Histograms", cex = 1.5)
Figure.3<-hist(system.solution.2, mfrow = c(3, 3))
mtext(outer = TRUE, side = 3, "Species Histograms", cex = 1.5)

par(mfcol = c(1, 1))
cc <- colors()
Figure.4<-scatterplot3d(system.solution.1[,2], 
                        system.solution.1[,3], 
                        system.solution.1[,4], highlight.3d = TRUE, col.axis = "blue",
                        col.grid = "lightblue", main = "Species Trajectories", pch = 20,angle=120,
                        xlab = "x1(t)", 
                        ylab = "x2(t)", zlab = "x3(t)")
rgl.viewpoint(0, 20)
Figure.5<-plot3d(system.solution.2[,2], 
                 system.solution.2[,3], 
                 system.solution.2[,4], 
                 pch = 1, cex = 1, xlab = "x1(t)", 
                 ylab = "x2(t)", zlab = "x3(t)")

#----------------------------LinkOut References------------------------------------------------
References.1<-c("MYC drives overexpression of telomerase RNA (hTR/TERC) in prostate cancer")
#---------------------------Function Library-------------------------------------------


