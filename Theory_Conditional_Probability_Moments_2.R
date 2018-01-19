library(markovchain)
library(HMM)
library(adaptMCMC)
library(mcmc)
library(MCMCpack)
library(combinat)
library(coda)
library(rbenchmark)
library(stringi)
library(stringr)
library(Matrix)
library(corrplot)
library(xtable)
library(coda)

#----------------------------------------------------------Data-----------------------------------------------------------------
libraries.R<-c('markovchain','HMM','adaptMCMC','mcmc','MCMCpack','combinat',
               'coda','rbenchmark','stringi','stringr','Matrix',
               "corrplot","xtable","coda")
#----------------------------------------------------------Project Data---------------------------------------------------------

project.description.file <- file("Project_Description.txt", "r", 
                                 blocking = FALSE)
project.note.file<-file("Project_Notes.txt", "r", 
                        blocking = FALSE)

project.Figures <- readPNG(system.file("img", stri_paste("Figure.",i,".png"), package="png"), info=TRUE)
project.Figures.annotation<-project.Figures[1]$text

project.description<-readLines(project.description.file)
project.notes<-readLines(project.note.file)
close(project.description.file)
close(project.note.file)
edit(project.description)
edit(project.notes)

#-----------------------------------------------------Mathematical Biology Data-------------------------------------------------
AA.names<-c('A','R','N','D','C','Q','E','G','H','I','L','K','M','F',
            'P','S','T','W','Y','V')
AA.order<-c('CYS', 'MET', 'PHE', 'ILE', 'LEU', 'VAL', 'TRP', 'TYR', 'ALA', 
                   'GLY', 'THR', 'SER', 'ASN', 'GLN', 'ASP', 'GLU', 'HIS', 'ARG', 'LYS', 'PRO')

AA.environment<-c('HIS','ILe','LYS','MET','PHE','THR','TRP','VAL')

Study.Cystein<-c("CYS","CIR","CME","CMT","CSD","CSO","CSW","CSX","CYM","CYX")

Informatics.Gene.Human<-c('Disorders','Domains','Drugs','Expression','Function',
                          'Genomics', 'Localization','Orthologs','Paralogs','Pathways',
                          'Products','Proteins','Publications','Sources','Summaries',
                          'Transcripts','Variants')
Informatics.Gene.Cards.Human<-c('GeneCards','MalaCard','LifeMap Discovery','PathCards','TGex',
                                'VarElect','GeneAnalytics','GeneALaCart',
                                'GenesLikeMe')
#----------------------------------------------------Markov Model Data and Parameters------------------------------------------

States.names<-c('S_1','S_2','S_3','S_4')
States.Description<-c('Description A','Description B','Description C','Description D')
Alphabet<-c('A','B','C','D')
Alphabet.symbols<-c('+','-','*','/')
Alphabet.description<-c('Description A','Description B','Description C','Description D')

sequence<-c('')
sequence.temporal<-c('A','B','A','B','C','A','A','A','B','C')
sequence.observations<-c('A','B','B','B','C','D','C','A','A','D')
sequence.observations.1<-c('A','B','B','B','C','D','C','A','A','D')
sequence.observations.2<-c('A','B','B','B','C','D','C','A','A','D')
sequence.observations.3<-c('A','B','B','B','C','D','C','A','A','D')

sequence.matrix <- createSequenceMatrix(sequence.observations, sanitize = FALSE)
sequence.matrix.sanitize<-createSequenceMatrix(sequence.observations, sanitize = TRUE)
sequence.matrix.temporal <- createSequenceMatrix(sequence.temporal, sanitize = FALSE)
sequence.matrix.temporal.sanitize<-createSequenceMatrix(sequence.temporal, sanitize = TRUE)
#----------------------------------------------------------Generate Matrices---------------------------------------------------
nrows<-3; ncolumns<-3
Model.1.matrix.adjacency<-diag(4)
Model.1.matrix.transition<-matrix(1,nrows,ncolumns)
Model.1.matrix.emission<- matrix(1,nrows,ncolumns)

Model.1.matrix.transition.A<-matrix(data=Model.1.matrix.adjacency,
                                    byrow=TRUE,nrow=4,dimnames=list(Alphabet,Alphabet))
Model.1.matrix.emission.A<-matrix(data=Model.1.matrix.adjacency,
                                    byrow=TRUE,nrow=4,dimnames=list(States.names,States.names))

Model.1.matrix.generator<-transition2Generator(P=Model.1.matrix.transition.A, t = 1, method = "logarithm")
States.names.scale<-c(1,1,1,1)

Model.1.matrix.hyperparameters.inference<-inferHyperparam(transMatr = Model.1.matrix.transition.A, scale = States.names.scale)
Model.1.matrix.hyperparameters.inference.dataset<-inferHyperparam(data=sequence.observations)

Model.1.hyperMatrix<-Model.1.matrix.adjacency
Model.1.predProb <- predictiveDistribution(sequence.observations[1:5], sequence.observations[6:10], hyperparam =Model.1.hyperMatrix )

#----------------------------------------------------------Model Development---------------------------------------------------
Model.1.markov.model<-initHMM(States.names,Alphabet,startProbs=NULL,
                              transProbs=Model.1.matrix.transition.A,
                              emissionProbs=Model.1.matrix.emission.A)
Model.1.markov.model.posterior.Probabilities<-posterior(Model.1.markov.model,sequence.observations)
Model.1.markov.model.logforward.Probabilities<-forward(Model.1.markov.model,sequence.observations)
Model.1.markov.model.backward.Probabilities<-backward(Model.1.markov.model,sequence.observations)
Model.1.markov.model.estimate.BW<-baumWelch(Model.1.markov.model,sequence.observations,
                                            maxIterations = 1000,delta=0.0001)
Model.1.markov.model.estimate.BW
Model.1.markov.model.estimate.BW$States$transProbs
Model.1.markov.model.estimate.Viterbi<-viterbi(Model.1.markov.model,sequence.observations)
Model.1.markov.model.probabilities<-exp(Model.1.markov.model.logforward.Probabilities)

#-----------------------------------------------------------Dirichlert Prior-----------------------------------------------------
Model.1.distribution.Dirichlet.prior.A<-priorDistribution(Model.1.matrix.transition.A,hyperparam=Model.1.matrix.hyperparameters)
#-----------------------------------------------------------Markov Chains----------------------------------------------------------

Model.1.markov.model.A<-new("markovchain", transitionMatrix = Model.1.matrix.transition.A,states=States.names)
Model.1.markov.model.continuous.A<-new("ctmc", transitionMatrix = Model.1.matrix.transition.A,states=States.names)

nbr.sequence.states<-100
Model.1.sequence.states <- markovchainSequence(n = nbr.sequence.states, 
                                       markovchain = Model.1.markov.model.A, 
                                       t0 = States.names[1])
Model.1.sequence.states.sample <- markovchainSequence(n = nbr.sequence.states, 
                                              markovchain = Model.1.markov.model.A, 
                                              t0 = sample(States.names[1],1))    
nbrboot<-5
mc.confidence.interval<-c(0.90,0.95,0.99)

mc.estimate.MLE.1 <- markovchainFit(data = sequence.observations,name="mc.estimate.MLE.1",
                                    sanitize = FALSE,byrow = TRUE,
                                    confidencelevel = mc.confidence.interval[1],
                                    possibleStates = States.names)
mc.estimate.MAP.1 <- markovchainFit(data = sequence.observations, method = "map",
                                    sanitize = FALSE,byrow = TRUE,
                                    name="mc.estimate.MAP.1",
                                    confidencelevel = mc.confidence.interval[1],
                                    possibleStates = States.names)
mc.estimate.LAP.1 <- markovchainFit(data = sequence.observations, method = "laplace",
                                    sanitize = FALSE,byrow = TRUE,
                                    laplacian = 0, name="mc.estimate.LAP.1",
                                    possibleStates = States.names)

mc.estimate.BSP.1 <- markovchainFit(data = sequence.observations.1, method = "bootstrap", 
                                    nboot = nbrboot, sanitize = FALSE,byrow = TRUE,
                                    name = "mc.estimate.BSP.1",
                                    possibleStates = States.names)
mc.estimate.BSP.2 <- markovchainFit(data = sequence.observations.2, method = "bootstrap", 
                                    nboot = nbrboot,parallel = TRUE, 
                                    sanitize = FALSE,byrow = TRUE,
                                    name = "mc.estimate.BSP.2",
                                    possibleStates = States.names)
mc.estimate.BSP.3 <- markovchainFit(data = sequence.observations.3, method = "bootstrap", 
                                    sanitize = FALSE,byrow = TRUE,
                                    nboot = nbrboot, name = "mc.estimate.BSP.3",
                                    possibleStates = States.names)

mc.list.estimate.nondistinct <- markovchainListFit(data=sequence.matrix.temporal, 
                                                   name = "mc.list.estimate.nondistinct")
mc.list.estimate.distinct <- markovchainListFit(data=sequence.matrix.temporal,byrow = TRUE,
                                                name = "mc.list.estimate.distinct")
mc.order<-2
mc.estimate.order.1<-fitHigherOrder(sequence.observations, order = mc.order)
mc.estimate.order.2<-seq2freqProb(sequence.observations)
mc.estimate.order.3<-seq2matHigh(sequence.observations, mc.order)

sequence.probabilties.df<-data.frame()
sequence.probabilities.df<-cbind(mc.estimate.order.2)

States.distribution <- c(0.9, 0.1)
nbr.samples<-c(5,'Inf')
mc.continuous.estimate<-rctmc(n =nbr.samples[2], ctmc = markov.model.continuous.A, T = 1)
mc.continuous.estimate.A<-rctmc(n = nbr.samples[2], ctmc = markov.model.continuous.A, 
                                initDist = States.distribution,T = 1)

mm.estimate.BW<-baumWelch(Model.1.markov.model,sequence.observations,nbr.sequence.states,maxIterations = 1000,delta=0.0001)
mm.estimate.VT<-viterbiTraining(Model.1.markov.model,sequence.observations,nbr.sequence.states)
#---------------------------------------------------------Model Analysis---------------------------------------------------------

sequence.property<-verifyMarkovProperty(sequence.observations)
sequence.property.order<-assessOrder(sequence.observations)
sequence.property.stationarity<-assessStationarity(sequence.observations, 1)
sequence.property.steadyStates<-steadyStates(Model.1.markov.model.A)
sequence.property.communication<-communicatingClasses(Model.1.markov.model.A)
sequence.property.recurrent<-recurrentClasses(Model.1.markov.model.A)
sequence.property.states.absorbing<-absorbingStates(Model.1.markov.model.A)
sequence.property.states.transient<-transientStates(Model.1.markov.model.A)
sequence.property.canonic<-canonicForm(Model.1.markov.model.A)
sequence.property.irreducible<-is.irreducible(Model.1.markov.model.A)
sequence.property.accessible<-is.accessible(Model.1.markov.model.A, States.names[1], States.names[length(States.names)-1])
sequence.property.period<-period(Model.1.markov.model.A)

sequence.posterior<-posterior(Model.1.markov.model,sequence.observations)
sequence.logForwardProbabilities = forward(Model.1.markov.model,sequence.observations)
sequence.logBackwardProbabilities = backward(Model.1.markov.model,sequence.observations)

sequence.simulation<-simHMM(Model.1.markov.model,length(States.names))

sequence.state.probabilties.df<-data.frame()
sequence.State.probabilities.df<-cbind(sequence.posterior)

sequence.simulation.df<-data.frame()
sequence.simulation.df<-rbind(sequence.observations,sequence.simulation)


Model.1.markov.model.A.ConditionalDistribution<-conditionalDistribution(Model.1.markov.model.A, States.names[1])
nrows.distribution<-20
Model.1.markov.model.A.Passage<-firstPassage(Model.1.markov.model.A, States.names[1], nrows.distribution)

multi.ci.MLE<-multinomialConfidenceIntervals(mc.estimate.MLE.1$estimate@transitionMatrix, sequence.matrix,mc.confidence.interval[1])
multi.ci.MAP<-multinomialConfidenceIntervals(mc.estimate.MAP.1$estimate@transitionMatrix, sequence.matrix,mc.confidence.interval[1])
multi.ci.LAP<-multinomialConfidenceIntervals(mc.estimate.LAP.1$estimate@transitionMatrix, sequence.matrix ,mc.confidence.interval[1])
multi.ci.BSP.1<-multinomialConfidenceIntervals(mc.estimate.BSP.1$estimate@transitionMatrix, sequence.matrix,mc.confidence.interval[1])
multi.ci.BSP.2<-multinomialConfidenceIntervals(mc.estimate.BSP.2$estimate@transitionMatrix, sequence.matrix,mc.confidence.interval[1])
multi.ci.BSP.3<-multinomialConfidenceIntervals(mc.estimate.BSP.3$estimate@transitionMatrix, sequence.matrix,mc.confidence.interval[1])


mm.estimate.probable.path<-viterbi(Model.1.markov.model,sequence.observations)
sequence.path.probabilties.df<-data.frame()
sequence.path.probabilities.df<-cbind(mm.estimate.probable.path)

#----------------------------------------------------------Model Diagnostics------------------------------------------------------------------

#---------------------------------------------------------Tables----------------------------------------------------------------

Table.1<-xtable(sequence.State.probabilities.df)
Table.2<-xtable(sequence.path.probabilities.df)
Table.3<-xtable(sequence.probabilties.df)
Table.4<-xtable(sequence.simulation.df)
Table.5<-print(mc.estimate.BSP.3$estimate)
#---------------------------------------------------------Figures---------------------------------------------------------------
Figure.1<-corrplot(multi.ci.MLE$upperEndpointMatrix)

Figure.2<-plot(sequence.posterior[1,],lty=1,col=1)
title(main='')
lines(sequence.posterior[2,],lty=2,col=2)
lines(sequence.posterior[3,],lty=3,col=3)
lines(sequence.posterior[4,],lty=4,col=4)
legend("topright", legend=c("A","B","C","D")
       , bty = "n",lwd=2, 
       cex=0.75, col=1:4, text.col=1:4, lty=1:4)

#Save Figures
Figure.1.Anotation<-c('')
writePNG(project.Figures,Figure.1.Anotation,asp=1.25)

#---------------------------------------------------------Function Library------------------------------------------------------



