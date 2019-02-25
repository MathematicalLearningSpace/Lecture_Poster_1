#-----------------------------R Code To Modify in the Classroom Lecture with Students-----------------------
#---------------------------------------------------R API --------------------------------------------------------
library(combinat);library(permutations);library(permute);library(xtable);library(igraph);library(Matrix);library(vegan);library(zoo)
library(adegenet);library(ActiveDriver);library(VennDiagram);library(venneuler);library(backShift);library(fields);library(ggplot2)
library(dnet);library(dcGor);library(ggm)
#------------------------------------------------------Data--------------------------------------------------------------

#----------------------------------------------------Cyclic Causal Graphs----------------------------------------------
CCG.Description<-c("Estimate connectivity matrix of hidden variables 
directed graph (cyclic,noncyclic graph) from linear system of observations from different shift interventions.")

CCG.Observations.N.Multiple<-3
CCG.Observations.N<-10^(CCG.Observations.N.Multiple)
CCG.variables.N<-6
CCG.Environments.N<-10
CCG.useCov <- TRUE 
#---------------------------------------------------Threshold Parameters-----------------------------------------------
CCG.thres<-0.75
CCG.thres.1<-0.9
CCG.thres.pe<-0.25
CCG.thres.pe.1<-0.25
CCG.thres.pe.2<-0.15
CCG.thres.pe.3<-0.05
#----------------------------------------------------------------------------------------------------------------------
CCG.providedA <- FALSE
CCG.cyclic<-TRUE
CCG.expNumNeigh<-0.1*CCG.variables.N
#---------------------------------------------------Coefficients-------------------------------------------------------
CCG.minCoef<-0.3
CCG.maxCoef<-0.8
#---------------------------------------------------Generate Graph-----------------------------------------------------
CCG.generate.result<-generateA(CCG.variables.N,
                               CCG.expNumNeigh,
                               CCG.minCoef,
                               CCG.maxCoef,
                               CCG.cyclic)

CCG.generate.result.A<-CCG.generate.result$A
cat("A has a cycle of size", CCG.generate.result$sizeCycle, "\n") 

#---------------------------------------------------Simulation of the Cyclic Causal Graph------------------------------
CCG.simulateObs<-TRUE
CCG.hidden<-FALSE
CCG.knownInterventions<-FALSE
CCG.fracVarInt<-0.5
CCG.intMult<-1.5
CCG.noiseMult<-1
CCG.nonGauss<-FALSE

CCG.simulation.result<- simulateInterventions(CCG.Observations.N, 
                                              CCG.variables.N,
                                              CCG.generate.result.A, 
                                              CCG.Environments.N, 
                                              CCG.intMult, 
                                              CCG.noiseMult, 
                                              CCG.nonGauss, 
                                              CCG.hidden, 
                                              CCG.knownInterventions, 
                                              CCG.fracVarInt, 
                                              CCG.simulateObs)


CCG.simulation.result.X <- CCG.simulation.result$X
CCG.simulation.result.Environment <- CCG.simulation.result$environment
CCG.simulation.result.baseInd <- CCG.simulation.result$configs$indexObservationalData

CCG.simulation.result.backshift.result<- backShift(CCG.simulation.result.X ,
                                                     CCG.simulation.result.Environment, 
                                                     covariance=CCG.useCov, 
                                                     ev=CCG.simulation.result.Environment, 
                                                     threshold=CCG.thres, 
                                                     baseSettingEnv = CCG.simulation.result.baseInd, 
                                                     tolerance = 1e-6, 
                                                     verbose = FALSE)
CCG.simulation.result.backshift.result.1<- backShift(CCG.simulation.result.X ,
                                                   CCG.simulation.result.Environment, 
                                                   covariance=CCG.useCov, 
                                                   ev=CCG.simulation.result.Environment, 
                                                   threshold=CCG.thres.1, 
                                                   baseSettingEnv = CCG.simulation.result.baseInd, 
                                                   tolerance = 1e-6, 
                                                   verbose = FALSE)

CCG.simulation.result.backshift.result.Ahat <- CCG.simulation.result.backshift.result$Ahat
CCG.simulation.result.backshift.result.1.Ahat <- CCG.simulation.result.backshift.result.1$Ahat

CCG.simulation.result.backshift.result.Ahat.structure <- CCG.simulation.result.backshift.result.1$AhatAdjacency
CCG.simulation.result.backshift.result.1.Ahat.structure <- CCG.simulation.result.backshift.result.1$AhatAdjacency


metrics.Thresholded.A.df <- data.frame()
metrics.Thresholded.A.df <-rbind(c(metricsThreshold(CCG.generate.result.A, 
                                                         CCG.simulation.result.backshift.result.Ahat, 
                                                         thres = CCG.thres.pe)),
                                 c(metricsThreshold(CCG.generate.result.A, 
                                                    CCG.simulation.result.backshift.result.Ahat, 
                                                    thres = CCG.thres.pe.1)),
                                 c(metricsThreshold(CCG.generate.result.A, 
                                                    CCG.simulation.result.backshift.result.Ahat, 
                                                    thres = CCG.thres.pe.2)),
                                 c(metricsThreshold(CCG.generate.result.A, 
                                                    CCG.simulation.result.backshift.result.Ahat, 
                                                    thres = CCG.thres.pe.3)))

metrics.Stability.Selection.df <- data.frame()
metrics.Stability.Selection.df<-rbind(c(metricsThreshold(CCG.generate.result.A,CCG.simulation.result.backshift.result.Ahat.structure, 
                                                          thres = 0)),
                                      c(metricsThreshold(CCG.generate.result.A,CCG.simulation.result.backshift.result.1.Ahat.structure, 
                                                         thres = 0)))

#------------------------------------------------------Tables-----------------------------------------------------------
Table.1<-xtable(metrics.Thresholded.A.df)
Table.2<-xtable(metrics.Stability.Selection.df)
#----------------------------------------------------Figures------------------------------------------------------------

Figure.1<-plot(CCG.simulation.result$environment)
Figure.2<-plot(CCG.simulation.result$X)
Figure.3<-plotGraphEdgeAttr(estimate =  CCG.generate.result.A, 
                            plotStabSelec = FALSE, 
                            labels = colnames(CCG.generate.result.A), 
                            thres.point = 0, 
                            thres.stab = CCG.thres, 
                            main = "True graph")

Figure.4<-plotGraphEdgeAttr(estimate = CCG.simulation.result.backshift.result.Ahat, 
                            plotStabSelec = FALSE, 
                            labels = colnames(CCG.generate.result.A), 
                            thres.point = CCG.thres.pe, thres.stab = CCG.thres, 
                            main = paste("Point estimate thresholded at", CCG.thres.pe))

Figure.5<-plotGraphEdgeAttr(estimate = CCG.simulation.result.backshift.result.Ahat.structure , 
                            plotStabSelec = TRUE, 
                            labels = colnames(CCG.generate.result.A), thres.point = CCG.thres.pe, 
                            edgeWeights = CCG.simulation.result.backshift.result.Ahat, thres.stab = CCG.thres, 
                            main = "Stability selection result")

Figure.6<-plotInterventionVars(CCG.simulation.result.backshift.result$varianceEnv, 
                               CCG.simulation.result$interventionVar)

#------------------------------------Figure Group ----------------------------------
par(mfrow = c(3,3), mar = c(0.5,1,0.5,0.5))
for(i in 1:CCG.Environments.N){
  plotDiagonalization(estConnectivity = CCG.simulation.result.backshift.result$Ahat, 
                      X = CCG.simulation.result.X, 
                      env =CCG.simulation.result.Environment, 
                      whichEnv = i)
}

#---------------------------------------------------References---------------------------------------------------------
citation("backshift")

Additional.Reference.1<-c("","","'")

#---------------------------------------------------Function Library---------------------------------------------------
#-------------Function Template Library for Classroom Presentation and Modification---------------------
f.1<-function(X)
 {
  Z<-""
  a<-1
  W<-runif(length(X),0,1)
  for(i in 1:length(X))
  {  
	Z<-stringr::str_c(Z,X[i])
	W[i]<-a*W[i]
  }
  output<-list()
  output$X<-X
  output$a<-a
  output$Z<-Z
  output$W<-W
  return(output)
 } 
test.f.1<-f.1(letters)
test.f.1



