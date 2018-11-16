#--------------------------------R API ------------------------------------------------------------------------
library(dnet);library(dcGor);library(ggm);library(abn);library(catnet);library(GADAG)
library(gRbase);library(HEMDAG);library(pcalg);library(unifDAG);library(igraph);library(xtable)
library(MXM);library(hash)
#-------------------------------------------------Data and Parameters------------------------------------------
DAG.model.random <- rdag(50, 4, 0.2)
DAG.model.random.X <- DAG.model.random$x
DAG.model.random.G<- DAG.model.random$G 
#-------------------------------------------------Topological Sorting------------------------------------------
DAG.model.random.G[DAG.model.random.G== 2] <- 1
DAG.model.random.G[DAG.model.random.G== 3] <- 0
topological_sort(DAG.model.random.G)

DAG.model.random.con <- pc.con(DAG.model.random.X)
DAG.model.random.or <- pc.or(DAG.model.random.con)

DAG.model.random.path.undir.1<-undir.path(DAG.model.random.G, 3, 4)
DAG.model.random.path.undir.2<-undir.path(DAG.model.random.G, 1, 3)

test.data<-list(X=t(matrix(c(rnorm(100),rnorm(100),rnorm(100)),3)),
                G=matrix(c(1,0,0,
                    0,1,1,
                    1,1,1),nrow=3,ncol=3))
k<-2
DAG.model.generation.N<-10^k
DAG.model.tolerance.Shannon<-10^(-6)
DAG.model.population.size<-5*ncol(test.data$X)
DAG.model.maximum.evaluation<-DAG.model.generation.N*DAG.model.population.size

target<-test.data$X[,3]
target.1 <- 3 * test.data$X[,1] + 2 * test.data$X[, 1] + 3 * test.data$X[, 1] 

dataset<-test.data$X[,-3]

DAG.model.auc<-auc(test.data$X[,1],test.data$X[,2],roc=TRUE)
DAG.model.auc.1<-auc(test.data$X[,1],test.data$X[,2],roc=TRUE)

#-----------------------------------------------Independence Tests---------------------------------------------

DAG.model.Test.Independence.Fisher<- testIndFisher(target, dataset, xIndex = 1, csIndex = 2)
DAG.model.Test.Independence.Spearman <- testIndSpearman(target, dataset, xIndex = 1, csIndex = 2)
DAG.model.Test.Independence.Fisher.Permutation <- permFisher(target, dataset, xIndex = 1, csIndex = 2, R = 999)

DAG.model.Test.Independence.df<-data.frame()
DAG.model.Test.Independence.df<-rbind(c(DAG.model.Test.Independence.Fisher$stat,DAG.model.Test.Independence.Fisher$pvalue),
                                      c(DAG.model.Test.Independence.Spearman$stat,DAG.model.Test.Independence.Spearman$pvalue),
                                      c(DAG.model.Test.Independence.Fisher.Permutation$stat,DAG.model.Test.Independence.Fisher.Permutation$pvalue))
colnames(DAG.model.Test.Independence.df)<-c("Statistic","P Value")
rownames(DAG.model.Test.Independence.df)<-c("Fisher","Spearman","Fisher Permutation")

DAG.model.algorithm.SES.Fisher<-SES(target, 
                                    dataset, 
                                    max_k = 3, 
                                    threshold = 0.05, 
                                    test = "testIndFisher")
DAG.model.algorithm.MMPC.Fisher<- MMPC(target, 
                                       dataset, 
                                       max_k = 3, 
                                       threshold = 0.05, 
                                       test="testIndFisher")

summary(DAG.model.algorithm.SES.Fisher)
summary(DAG.model.algorithm.MMPC.Fisher)

#-------------------------------------Permutation Matrix-------------------------------------------------------------
variable.N<-ncol(test.data$G)
variable.N.permutations <- sample(variable.N) # permutation
variable.N.permutations.matrix <- matrix(0,variable.N,variable.N)
variable.N.permutations.matrix[variable.N*0:(variable.N-1) + variable.N.permutations] <- 1
T <- matrix(rnorm(variable.N),variable.N,variable.N)
T[upper.tri(T,diag=TRUE)] <- 0

#-------------------------------------Estimation --------------------------------------------------------------------
return.level <- 1
DAG.optimal.results<-GADAG_Run(X=test.data$X, 
                       lambda=0.1,
                       return.level=return.level,
                       GADAG.control=list(n.gen=DAG.model.generation.N, 
                                          tol.Shannon=DAG.model.tolerance.Shannon,
                                          pop.size = DAG.model.population.size, 
                                          max.eval=DAG.model.maximum.evaluation))

DAG.optimal.results.best<-DAG.optimal.results$G.best
#-------------------------------------Analysis of Genetic Algorithm----------------------------------------------------
plot.evol<-TRUE
plot.graph<-TRUE
DAG.optimal.results.analysis<-GADAG_Analyze(DAG.optimal.results,
                                            G=test.data$G,
                                            X=test.data$X,
                                            plot.control=list(plot.evol=plot.evol,
                                                              plot.graph=plot.graph))

DAG.model.Fitness <- fitness(P=variable.N.permutations.matrix, X=test.data$X, T=T, lambda=0.1)
print(DAG.model.Fitness) 

DAG.optimal.results.df<-data.frame()
DAG.optimal.results.df<-rbind(c(DAG.model.generation.N,
                                DAG.model.population.size,
                                DAG.optimal.results$f.best))
colnames(DAG.optimal.results.df)<-c("Generations","Population","F best")

DAG.model.path.mmpc <- mmpc.path(target, 
                                 dataset, max_ks = NULL, 
                                 thresholds = NULL, 
                                 test = NULL, user_test = NULL, 
                                 robust = FALSE, ncores = 1)

#--------------------------------------Regularisation Analysis---------------------------------------------------
DAG.model.regularisation.ridge.1 <- ridge.reg(target, dataset, lambda = 0, B = 1, newdata = NULL)
DAG.model.regularisation.ridge.2 <- ridge.reg(target, dataset, lambda = 0.5, B = 1, newdata = NULL)
DAG.model.regularisation.ridge.3 <- ridge.reg(target, dataset, lambda = 0.5, B = 100, newdata = NULL) 

DAG.model.regularisation.df<-data.frame()
DAG.model.regularisation.df<-rbind(c(DAG.model.regularisation.ridge.1$beta),
                                      c(DAG.model.regularisation.ridge.2$beta),
                                      c(DAG.model.regularisation.ridge.3$beta))
colnames(DAG.model.regularisation.df)<-c("Beta-1","Beta-2")
rownames(DAG.model.regularisation.df)<-c("Ridge-1","Ridge-2","Ridge-3")

#------------------------------------------------Tables-------------------------------------------------------
Table.1<-xtable(DAG.optimal.results.df)
Table.2<-xtable(DAG.model.Test.Independence.df)
Table.3<-xtable(DAG.model.regularisation.df)

#------------------------Figures for Presentation in Classroom------------------------------------------------------

Figure.1<-plotnetwork(DAG.model.random.G) 
par(mfrow = c(2,2), mar = c(0.5,1,0.5,0.5))
Figure.2<-plot(DAG.model.algorithm.SES.Fisher, mode = "all")
Figure.3<-plot(DAG.model.algorithm.MMPC.Fisher,mode="all")

#-----------------References for Student Completion in Classroom----------------------------------------------------
Reference.1<-c("",
               "",
               "")
#------------------------------------------------Function Library----------------------------------------------
