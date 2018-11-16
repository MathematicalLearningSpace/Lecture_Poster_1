#----------------------------------------R API ----------------------------------------------------------
library(DiffusionRimp);library(Langevin);library(xtable);library(igraph);library(costat);library(HMM);library(markovchain)
#-------------------------------------Data----------------------------------------------------------------
Data.Captions<-c("The number of bins in the data is 20",
                 "The number of steps is a sequence from 1 to 5",
                 "The sampling frequency is 10 to the power of 3",
                 "The sequence t is from 1 to the number of observations e^7")
data.bins <- 20
data.steps <- c(1:5)
data.sf <- 10^3
N<-exp(7)
t <- 1:N
#-------------Temporal Data Generation for Cubic Drift and Quadratic Diffusion Polynomials----------------
Temporal.Sequence.Caption<-c("Several temporal sequences are generated.",
                             "Sequence 1 has coefficients for both the drift and diffusion polynomials",
                             "Sequence 2 has alpha_1=-1 for drift polynomial and 
                              beta_0=1 for diffusion polynomial")
set.seed(1111)
temporal.sequence.1 <- timeseries1D(N=N, sf=1, dt=0.01)
# drift D^1 = -x and constant diffusion D^2 = 1
temporal.sequence.2 <- timeseries1D(N=N, 
                                    d11=-1, 
                                    d20=1, 
                                    sf=data.sf)

#------------------------------------Test of Two Conditions-------------------------------------------

temporal.sequence.1.Property.1<-verifyMarkovProperty(temporal.sequence.1, verbose = TRUE)
wilcox.test(temporal.sequence.1[2:length(temporal.sequence.1)],diff(temporal.sequence.1),paired = TRUE, alternative = "less")

temporal.sequence.1.Property.2<-assessStationarity(temporal.sequence.1, 1, verbose = TRUE)

temporal.sequence.1.Property.2.xews <- ewspec(temporal.sequence.1[1:2^6], smooth.dev=var)$S
temporal.sequence.1.Property.2A<-TOSts(temporal.sequence.1.Property.2.xews)


temporal.sequence.1.Property.1.df<-data.frame()
temporal.sequence.1.Property.1.df<-rbind(c(temporal.sequence.1.Property.1$statistic,
                                           temporal.sequence.1.Property.1$dof,
                                           temporal.sequence.1.Property.1$p.value))
temporal.sequence.1.Property.2.df<-data.frame()
temporal.sequence.1.Property.2.df<-rbind(c(temporal.sequence.1.Property.2$statistic,
                                           temporal.sequence.1.Property.2$dof,
                                           temporal.sequence.1.Property.2$p.value))
temporal.sequence.1.Property.3.df<-data.frame()
temporal.sequence.1.Property.3.df<-rbind(c(temporal.sequence.1.Property.2A))

#--------------------Model Estimation of Drift and Diffusion Vectors-----------------------------------

Model.Estimation.Caption<-c("Three different models are estimated based on both 1D and 2D temporal sequences",
                            "Model 1 fits Langevin model to Temporal Sequence 1 with bins=40, steps=1:3 and sampling frequency of 10^3",
                            "Model 2 fits Langevin model to Temporal Sequence 2 with bins=40, steps=1:3")

model.estimation.1D.1<-Langevin1D(temporal.sequence.1, data.bins, data.steps, data.sf, reqThreads=2)
model.estimation.1D.2<-Langevin1D(temporal.sequence.2, data.bins, data.steps)

summary(model.estimation.1D.1)
summary(model.estimation.1D.2)

print(model.estimation.1D.1, digits = max(3, getOption("digits") - 3))

#------------------------------Estimate Model Coefficients-----------------------------------

Model.Estimation..Coefficients.Caption<-c("The drift coefficients for Model 2 are estimated by linear regression 
                                          from mean bin from estimation of Model 2",
                                          "The diffusion coefficients for Model 2 are estimated by linear regression 
                                          from mean bin from estimation of Model 2",
                                          "Temporal Sequence 4 is generated based on both estimated drift and diffusion coefficient models")

model.estimation.1D.2.drift <- coef(lm(model.estimation.1D.2$D1 ~ model.estimation.1D.2$mean_bin + I(model.estimation.1D.2$mean_bin^2) + 
                   I(model.estimation.1D.2$mean_bin^3),weights = 1/model.estimation.1D.2$eD1))
model.estimation.1D.2.diffusion <- coef(lm(model.estimation.1D.2$D2 ~ model.estimation.1D.2$mean_bin + 
                   I(model.estimation.1D.2$mean_bin^2), weights =1/model.estimation.1D.2$eD2))

#------------------------------Test Third Condition: Ratio of D4 to D2^2-----------------------------------
model.estimation.1D.2$D4 < model.estimation.1D.2$D2^2
summary(model.estimation.1D.2)

#------------------------------Generate Temporal sequence from estimated coefficients----------------------

temporal.sequence.4<-timeseries1D(N =N, 
                                  d10 = model.estimation.1D.2.drift[1], 
                                  d11 = model.estimation.1D.2.drift[2],
                                  d12 = model.estimation.1D.2.drift[3], 
                                  d13 = model.estimation.1D.2.drift[4], 
                                  d20 = model.estimation.1D.2.diffusion[1], 
                                  d21 = model.estimation.1D.2.diffusion[2],
                                  d22 = model.estimation.1D.2.diffusion[3], 
                                  sf = data.sf)

Model.estimation.parameters.df<-data.frame()
Model.estimation.parameters.df<-rbind(c(model.estimation.1D.2.drift,model.estimation.1D.2.diffusion))

#---------------------------------------Tables-------------------------------------------------------------

Table.1.Caption<-c("Table 1 shows values of the estimated coefficients fore each of the models",
                   "fit a cubic function to the estimated drift coefficient and a quadratic function to the diffusion coefficient")
Table.1<-xtable(Model.estimation.parameters.df)

Table.2.Caption<-c("Table 2 shows values of statistic (the chi - square statistic)",
                   "dof (degrees of freedom)",
                   "corresponding p value")
Table.2<-xtable(temporal.sequence.1.Property.1.df)
Table.3.Caption<-c("Table 3 shows value statistic (the chi - square statistic)",
                   "dof (degrees of freedom)",
                   "corresponding p value")
Table.3<-xtable(temporal.sequence.1.Property.2.df)
Table.4.Caption<-c("Table 4 shows the value of the statistic")
Table.4<-xtable(temporal.sequence.1.Property.3.df)

#---------------------------------------Figures------------------------------------------------------------

par(mfrow = c(1,2), mar = c(1,1,1,1))
Figure.Group.1.Caption<-c("The following set of 4 figures show the relationships of",
                          "Figure A shows Temporal Sequence 1",
                          "Figure B shows Temporal Sequence 2",
                          "")
Figure.1<-plot(t, temporal.sequence.1, type="l",xlab='phi',ylab='psi',
                         main='',col='Blue')
legend("bottomright", legend=c(c(paste("mean:", mean(temporal.sequence.1), "\n var:", var(temporal.sequence.1))))
       , bty = "n",lwd=2, 
       cex=0.75, lty=1:1)
Figure.2<-plot(t, temporal.sequence.2,  type="l", xlab='phi',ylab='psi',
               main='',col='Blue')
legend("bottomright", legend=c(c(paste("mean:", mean(temporal.sequence.2), "\n var:", var(temporal.sequence.2))))
       , bty = "n",lwd=2, 
       cex=0.75, lty=1:1)

par(mfrow = c(2,2), mar = c(2,2,2,2))
x<-rnorm(100)
Figure.Group.2.Caption<-c("One-dimensional Langevin Approach:",
                          "(a1) drift coefficient, D(1)(x) =x^3+x^2+x+1.",
                          "(b1) diffusion coefficient,D(2)(x)=x^2+x+1.",
                          "(a2) drift coefficient, D(1)(x) =-x.",
                          "(b2) diffusion coefficient,D(2)(x)=1.")
  
Figure.1<-plot(x, cubic.drift.equation(x,1,1,1,1), xlab='phi',ylab='psi',
               main='',col='Blue')
legend("topright", legend=c("")
       , bty = "n",lwd=2, 
       cex=0.75, lty=1:1)
Figure.2<-plot(x, quadratic.diffusion.equation(x,1,0,1), xlab='phi',ylab='psi',
               main='',col='Green')
legend("topright", legend=c("")
       , bty = "n",lwd=2, 
       cex=0.75, lty=1:1)
Figure.3<-plot(x, cubic.drift.equation(x,0,0,-1,0), xlab='phi',ylab='psi',
               main='',col='Blue')
legend("topright", legend=c("")
       , bty = "n",lwd=2, 
       cex=0.75, lty=1:1)
Figure.4<-plot(x, quadratic.diffusion.equation(x,0,0,1), xlab='phi',ylab='psi',
               main='',col='Green')
legend("topright", legend=c("")
       , bty = "n",lwd=2, 
       cex=0.75, lty=1:1)

#---------------------------------------References-------------------------------

Reference.1<-c("Kleinhans D (2012).", 
               "Estimation of Drift and Diffusion Functions from Time Series Data: A Maximum Likelihood Framework.",
               "Physical Review E, 85 (2), 026705. doi:10.1103/physreve.85.026705.")

Reference.2<-c("Scholz T., Raischel F., Lopes V, Lehle B, Wachter M., Peinke J., Lind P. (2015).", 
               "Parameter-free Resolution of the Superposition of Stochastic Signals.",
               "submitted URL http://arxiv.org/abs/1510.07285.")
               
Reference.3<-c("Friedrich R, Peinke J (1997)", 
               "Description of a Turbulent Cascade by a Fokker-Planck Equation.", 
               "Physical Review Letters, 78, 863-866.doi:10.1103/PhysRevLett.78.863.")

Reference.4<-c("Friedrich R, Peinke J, Sahimi M, Reza Rahimi Tabar M (2011).", 
               "Approaching Complexity by Stochastic Methods: From Biological Systems to Turbulence.",
               "Physics Reports, 506 (5), 87-162. ISSN 0370-1573. doi:10.1016/j.physrep.2011.05.003.")

#---------------------------------------Function Library------------------------

Function.Caption<-c("Equation 1 shows a cubic polynomial.",
                    "Equation 2 shows a quadratic polynomial")
cubic.drift.equation<-function(x,alpha_3,alpha_2,alpha_1,alpha_0)
{
  
  y<-alpha_3*x^3+alpha_2*x^2+alpha_1*x+alpha_0
  return(y)
}
quadratic.diffusion.equation<-function(x,beta_2,beta_1,beta_0)
{
  y<-beta_2*x^2+beta_1*x+beta_0
  return(y)
}
