
library(backShift); library(BiDag);library(CIEE);library(dagbag);library(dagR);library(GADAG);library(GA);library(HEMDAG)
library(unifDAG);library(Matrix);library(xtable);library(PearsonDS);library(igraph);library(qgraph)
#---------------------------------------Data------------------------------------------------------------

pearson.N<-100
#-------------------------------------Moment Parameters for the Distributions-------------------------------------------------

p0pars <- list(mean=1, sd=1)
pIpars <- list(a=1, b=1, location=1, scale=1) 
pIIpars <- list(a=1, location=1, scale=1) 
pIIIpars <- list(shape=1, location=1, scale=1)
pIVpars <- list(m=1, nu=1, location=1, scale=1)
pVpars <- list(shape=1, location=1, scale=1) 
pVIpars <- list(a=1, b=1, location=1, scale=1) 
pVIIpars <- list(df=10, location=1, scale=1)

#-------------------------------------Generate Random variables from the Distributions----------------------------------------

error.pearson.0<-rpearson0(pearson.N,params=p0pars)
error.pearson.1<-rpearsonI(pearson.N,params=pIpars)
error.pearson.2<-rpearsonII(pearson.N,params=pIIpars)
error.pearson.3<-rpearsonIII(pearson.N,params=pIIIpars)
error.pearson.4<-rpearsonIV(pearson.N,params=pIVpars)
error.pearson.5<-rpearsonV(pearson.N,params=pVpars)
error.pearson.6<-rpearsonVI(pearson.N,params=pVIpars)
error.pearson.7<-rpearsonVII(pearson.N,params=pVIIpars)


design.matrix.Pearson<-matrix(error.pearson.7,nrow=600,ncol=6)

g <- make_empty_graph(n = 6) %>%
  add_edges(c(1,2, 1,3, 2,4, 4,5,4,6,5,6)) %>%
  set_edge_attr("color", value = "blue") %>%
  set_edge_attr("weight",value=c(0.5,0.3,0.1,0.2,0.7,0.3))
 
V(g)[[]]
E(g)[[]]

is_dag(g)

design.matrix<-as_adjacency_matrix(g)

experimental.data<-list(A=design.matrix.Pearson,G=design.matrix)

experimental.data$A

n.gen <- 1e10
tol.Shannon <- 1e-10
pop.size <- 10*ncol(experimental.data$G)
max.eval <- n.gen * pop.size 
p.xo<-0.5
p.mut<-0.05
A_results<- GADAG_Run(X=experimental.data$A, lambda=0.1,GADAG.control=list(n.gen=n.gen,
                                                                           pop.size=pop.size,
                                                                           max.eval=max.eval,
                                                                           tol.Shannon=tol.Shannon,
                                                                           p.xo=p.xo,
                                                                           p.mut=p.mut))
plot.graph<-TRUE
plot.evol<-TRUE
plot.png<-FALSE
A_analysis<-GADAG_Analyze(A_results, G=experimental.data$G, X=experimental.data$A,  
                          plot.control = list(plot.graph= plot.graph, 
                                              plot.evol = plot.evol)
                          )

A.results.table.df<-data.frame()
A.results.table.df<-cbind(A_results$f.best)
colnames(A.results.table.df)<-c("Best Fitness")

A.analysis.table.df<-data.frame()
A.analysis.table.df<-cbind(A_analysis)
colnames(A.analysis.table.df)<-c("A")

GA.parameter.table.df<-data.frame()
GA.parameter.table.df<-rbind(c("n.gen","maximal number of population generations","(>0)"),
                              c("pop.size","initial population size for the genetic algorithm", "(>0)"),
                              c("max.eval" , "overall maximal number of calls of the evaluation function" ,"(>0)"),
                              c("tol.Shannon" ,"threshold for the Shannon entropy" , "(>0)"),
                              c("p.xo", "crossover probability of the genetic algorithm" ,"(between 0 and 1)"),
                              c("p.mut" , "mutation probability of the genetic algorithm"  , "(between 0 and 1)"))
colnames(GA.parameter.table.df)<-c("Parameter","Description","Interval")


#--------------------------------------Tables-----------------------------------------------------------

Table.1<-xtable(GA.parameter.table.df)
Table.2<-xtable(A.results.table.df)
Table.3<-xtable(A.analysis.table.df)

#--------------------------------------Figures----------------------------------------------------------

Figure.1<-plot(g,layout=layout_in_circle,
               vertex.color="Blue",vertex.size=10,
               vertex.label.dist=2,edge.arrow.size=0.5,
               main='Experimental Graph')

#-------------------------------------Function Library--------------------------------------------------
