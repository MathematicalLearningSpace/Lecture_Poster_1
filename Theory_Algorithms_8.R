#-----------------------------R Code To Modify in the Classroom Lecture with Students-----------------------
#----------------------------------R API -------------------------------------------------------------------
library(xtable);library(smfsb);library(VGAM);library(nlstools);library(datasets);library(seqRFLP);library(stats);library(graphics);library(org.Hs.eg.db)
#-------------------------------Data-----------------------------------------

enz<-org.Hs.egENZYME
mapped_genes <- mappedkeys(enz)
enz.genes <- as.list(enz[mapped_genes])
EC.Level<-c("EC 1 oxidoreductases","EC 2 transferases", "EC 3 hydrolases","EC 4 lyases","EC 5 isomerases","EC 6 ligases")

data(enzdata)
data(spnModels)
Puromycin
Puromycin.caption<-c("reaction velocity versus substrate concentration in an enzymatic reaction with untreated cells or cells treated with Puromycin.",
                     "substrate concentrations (ppm)",
                     "instantaneous reaction rates (counts/min/min)",
                     "levels treated untreated")

#----------------------------Michaelis-Menten Model------------------------

PurTrt.treated <- Puromycin[ Puromycin$state == "treated", ]

SSmicmen(PurTrt.treated$conc, 200, 0.05)  # response only
Vm <- 200; K <- 0.05
SSmicmen(PurTrt.treated$conc, Vm, K)      # response and gradient
print(getInitial(rate ~ SSmicmen(conc, Vm, K), data = PurTrt.treated), digits = 3)

model.1 <- nls(rate ~ SSmicmen(conc, Vm, K), data = PurTrt.treated)
model.2 <- nls(rate ~ SSmicmen(conc, Vm, K), data = Puromycin,
           subset = state == "treated")
model.summary.1<-summary(model.1)
model.summary.2<-summary(model.2)

model.3 <- nls(rate ~ Vm * conc/(K + conc), data = Puromycin,
           subset = state == "treated",
           start = c(Vm = 200, K = 0.05))
model.4 <- nls(rate ~ Vm * conc/(K + conc), data = Puromycin,
           subset = state == "untreated",
           start = c(Vm = 160, K = 0.05))
model.5 <- nls(rate ~ conc/(K + conc), data = Puromycin,
           subset = state == "treated", start = c(K = 0.05),
           algorithm = "plinear")

model.summary.3<-summary(model.3)
model.summary.4<-summary(model.4)
model.summary.5<-summary(model.5)

model.summary.3$coefficients
model.summary.4$coefficients
model.summary.5$coefficients

#--------------Stochastic Petri Net Model-----------------------------------
parms.E.1 <-c(k_forward = 0.00166, k_reverse  = 0.01831564, k_cat = 0.1) 
parms.E.2 <-c(k_forward = 0.00166, k_reverse = 0.01831564, k_cat = 0.1) 
parms.E.3 <-c(k_forward = 0.00166, k_reverse = 0.01831564, k_cat = 0.1) 
parms.E.4 <-c(k_forward = 0.00166, k_reverse = 0.01831564, k_cat = 0.1) 

params.E.1.caption<-c("rate of formation of ES",
              "rate of decomposition of ES",
              "maximum number of enzymatic reactions catalysed per second.")

Michaelis.constant<-(parms.E.1[2]+parms.E.1[3])/parms.E.1[1] # Ribonuclease 0.0079
Kinetic.efficiency<-parms.E.1[3]/Michaelis.constant
Dissocation.constant<-parms.E.1[2]/parms.E.1[1] #change with temperature, pH and salt concentration for ligand-protein interaction

stepMM.kinetic<-StepCLE(MM)
MM$h<-Michaelis.Menten.enzyme.kinetic.modified
Michaelis.Menten.enzyme.kinetic.simulation.1<-simTs(MM$M,0,20,0.1,stepMM.kinetic,parms=parms.E.1)

parms.D.1<-c(alpha = 0.00166, beta = 0.2)

stepDim.kinetic<-StepCLE(Dimer)
Dimer$h<-Dimerisation.kinetic.modified
Dimerisation.enzyme.kinetic.simulation.1<-simTs(Dimer$M,0,20,0.1,stepDim.kinetic,parms=parms.D.1)

#--------------Simulate the Drift and Diffusion relationship-----------------

stepProc = StepSDE(drift,diffusion)
sequence<-seq(0,20,0.1)

simulation.1 = simTs(c(x=5,y=0.1,z=0.1),0,20,0.1,stepProc,parms=c(delta=1,lambda=1,alpha=1,mu=0.1,sigma=0.1,kappa=1))
simulation.2 = simTs(c(x=5,y=0.1,z=0.1),0,20,0.1,stepProc,parms=c(delta=1,lambda=1,alpha=1,mu=0.1,sigma=0.1,kappa=1.25))
simulation.3 = simTs(c(x=5,y=0.1,z=0.1),0,20,0.1,stepProc,parms=c(delta=1,lambda=1,alpha=1,mu=0.1,sigma=0.1,kappa=1.5))
simulation.4 = simTs(c(x=5,y=0.1,z=0.1),0,20,0.1,stepProc,parms=c(delta=1,lambda=1,alpha=1,mu=0.1,sigma=0.1,kappa=1.75))
simulation.5 = simTs(c(x=5,y=0.1,z=0.1),0,20,0.1,stepProc,parms=c(delta=1,lambda=1,alpha=1,mu=0.1,sigma=0.1,kappa=2))

Summary.Statistics.df<-data.frame()
Summary.Statistics.df<-rbind(c(summary(simulation.1)),
                             c(summary(simulation.2)),
                             c(summary(simulation.3)),
                             c(summary(simulation.4)),
                             c(summary(simulation.5)))
colnames(Summary.Statistics.df)<-c("X:Min","X:1stQ","X:Median","X:Mean","X:3rdQ","X:Max",
                                   "Y:Min","Y:1stQ","Y:Median","Y:Mean","Y:3rdQ","Y:Max",
                                   "Z:Min","Z:1stQ","Z:Median","Z:Mean","Z:3rdQ","Z:Max")
rownames(Summary.Statistics.df)<-c("Simulation 1",
                                   "Simulation 2",
                                   "Simulation 3",
                                   "Simulation 4",
                                   "Simulation 5")

simulation.sample.1 = simSample(25,c(x=5,y=0.1,z=0.1),0,20,stepProc,parms=c(delta=1,lambda=1,alpha=1,mu=0.1,sigma=0.1,kappa=1))

#------------------------------Tables---------------------------------------

Table.1<-xtable(Summary.Statistics.df)

#-----------------------------Figures--------------------------------------

par(mfrow = c(2,1), mar = c(5,5,5,5))
Figure.1<-plot(sequence,simulation.1[,2],xlab='',ylab='',
               main='',col='black',type="l")
lines(simulation.2[,2],lty=2,col="red")
#lines(simulation.3[,2],lty=3,col="green")
lines(simulation.4[,2],lty=4,col="blue")
#lines(simulation.5[,2],lty=4,col="purple")
legend("bottomright", legend=c("k=1","k=1.25","k=1.5","k=1.75","k=2")
       , bty = "n",lwd=2, 
       cex=0.75, lty=1:5,col=c("black","red","green","blue","purple"))
Figure.1A<-plot(sequence,simulation.1[,1],xlab='',ylab='',
               main='',col='black',type="l")
lines(simulation.2[,1],lty=2,col="red")
#lines(simulation.3[,1],lty=3,col="green")
#lines(simulation.4[,1],lty=4,col="blue")
#lines(simulation.5[,1],lty=4,col="purple")
legend("bottomright", legend=c("k=1","k=1.25","k=1.5","k=1.75","k=2")
       , bty = "n",lwd=2, 
       cex=0.75, lty=1:5,col=c("black","red","green","blue","purple"))
Figure.1.caption<-c("(a) x2 species and (b) x1 species for Drift and Diffusion Parameters of 
                    a range of kappa from (1,2) with 
                    x1=5,y1=0.1 and lambda=1,alpha=1,mu=0.1,and sigma=0.1")

Figure.2<-plot(sequence,simulation.1[,1],xlab='',ylab='',
               main='',col='black',type="l",ylim=c(0,10))
lines(simulation.1[,2],lty=2,col="red")
lines(simulation.1[,3],lty=3,col="green")
#lines(simulation.4[,2],lty=4,col="blue")
#lines(simulation.5[,2],lty=4,col="purple")
legend("bottomright", legend=c("Simulation 1:X","Simulation 1:Y","Simulation 1:Z")
       , bty = "n",lwd=2, 
       cex=0.75, lty=1:5,col=c("black","red","green","blue","purple"))

par(mfrow = c(2,2), mar = c(5,5,5,5))
Figure.3<-hist(simulation.sample.1[,"x"])
Figure.4<-hist(simulation.sample.1[,"y"])
Figure.5<-hist(simulation.sample.1[,"z"])

Figure.6<-plot(sequence,Michaelis.Menten.enzyme.kinetic.simulation.1[,1],xlab='',ylab='',
               main='',col='black',type="l")
lines(Michaelis.Menten.enzyme.kinetic.simulation.1[,2],lty=2,col="red")
lines(Michaelis.Menten.enzyme.kinetic.simulation.1[,3],lty=3,col="green")
lines(Michaelis.Menten.enzyme.kinetic.simulation.1[,4],lty=4,col="blue")
legend("bottomright", legend=c("E","S","ES","P")
       , bty = "n",lwd=2, 
       cex=0.75, lty=1:5,col=c("black","red","green","blue"))

Figure.7<-plot(sequence,Dimerisation.enzyme.kinetic.simulation.1[,1],xlab='',ylab='',
               main='',col='black',type="l",ylim=c(0,300))
lines(Dimerisation.enzyme.kinetic.simulation.1[,2],lty=2,col="red")
legend("bottomright", legend=c("x1","x2")
       , bty = "n",lwd=2, 
       cex=0.75, lty=1:5,col=c("black","red"))


Figure.9<-plot(rate ~ conc, data = Puromycin, las = 1,
               xlab = "Substrate concentration (ppm)",
               ylab = "Reaction velocity (counts/min/min)",
               pch = as.integer(Puromycin$state),
               col = as.integer(Puromycin$state),
               main = "Puromycin data and fitted Michaelis-Menten curves")
conc <- seq(0, 1.2, length.out = 101)
lines(conc, predict(model.3, list(conc = conc)), lty = 1, col = 1)
lines(conc, predict(model.4, list(conc = conc)), lty = 2, col = 2)
lines(conc, predict(model.5, list(conc = conc)), lty = 2, col = 2)
legend(0.8, 120, levels(Puromycin$state),
       col = 1:3, lty = 1:3, pch = 1:3)

#----------------------------References----------------------------------

Reference.1<-c("Handbook of Mathematical Functions",
                  "https://dlmf.nist.gov/")

Reference.2<-c("Wikipedia: Enzyme Kinetics",
              "https://en.wikipedia.org/wiki/Enzyme_kinetics")

Reference.3<-c("Treloar, M. A. (1974)", 
               "Effects of Puromycin on Galactosyltransferase in Golgi Membranes", 
               "M.Sc. Thesis, U. of Toronto.")
Reference.4<-c("Roberts, R.J., Vincze, T., Posfai, J., Macelis, D. (2010)", 
               "REBASE-a database for DNA restriction and modification: enzymes, genes and genomes.", 
               "Nucl. Acids Res. 38: D234-D236. http://rebase.neb.com"
               )

#----------------------------Function Library-----------------------------
Dimerisation.kinetic.modified<-function(x,t,parms.D)
{
  with(as.list(c(x,parms.D)), {
    return(c(1/2*(alpha * x1^2 - alpha*x1), 
             beta * x2))
  })
}

Michaelis.Menten.enzyme.kinetic.modified<-function (x, t, parms.E) 
{
  with(as.list(c(x, parms.E)), {
    return(c(k_forward * S * E, 
             k_reverse * SE, 
             k_cat * SE))
  })
}

drift<-function(x,t,parms)
{
  with(as.list(c(x,parms)),{
    c( lambda - x*y ,
       alpha*(mu-y),
       delta-z)
  })
}
diffusion <- function(x,t,parms)
{
  with(as.list(c(x,parms)),{
    matrix(c( sqrt(lambda + x*y^kappa) , 0, 0,
              0, sigma*sqrt(y^kappa),0,
              0,0,sqrt(z)),ncol=3,nrow=3,byrow=TRUE)
  })
}
