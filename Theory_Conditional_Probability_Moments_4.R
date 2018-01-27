library(dnet)
library(dcGor)
library(ggm)
library(abn)
library(catnet)
library(GADAG)
library(gRbase)
library(HEMDAG)
library(pcalg)
library(unifDAG)
library(igraph)
library(xtable)
library(readr)
#-------------------------------------Data-------------------------------------------------------------------

File.categories<-c("VertexID","Name","Description","VertexID_1","VertexID_2","Direction", "EvidenceType", "Weight","Comment")

Conditional.Moments.Concept.Map <- read_csv("data/Conditional_Moments_Concept_Map.txt")
View(Conditional.Moments.Concept.Map)

#-------------------------------------Graph Specification----------------------------------------------------

UG.1<-UG(~(X + Y + Z)^2)
UG.1.matrix<-as.matrix(UG.1)
UG.1.matrix.cycle<-cycleMatrix(UG.1)
UG.1.Cycles.Fundamental<-fundCycles(UG.1)
UG.1.Path<-findPath(bfsearch(UG.1)$tree, st=1, en=3)

UG.2<-UG(~(W+X*Y + Z)^2)
UG.2.matrix<-as.matrix(UG.2)
UG.3<-UG(~(V+W*X-Y*Z)^2)
UG.3.matrix<-as.matrix(UG.3)

UG.4<-UG(~(U +V+ W + X + Y + Z)^2)
UG.4.matrix<-as.matrix(UG.4)
UG.4.matrix.cycle<-cycleMatrix(UG.4)
UG.4.Cycles.Fundamental<-fundCycles(UG.4)
UG.4.Path<-findPath(bfsearch(UG.1)$tree, st=1, en=6)


UG.formulas.df<-data.frame()
UG.formulas.df<-rbind(c("~(X+Y+Z)^2",UG.1.Path))

#---------------------------------------Mixed Graphs---------------------------------------------------------

Graph.mixed.1 <- makeMG(dg = DG(Y~X, Z~W, W~Z), 
                        ug = UG(~ X*Q), 
                        bg = UG(~ Y*X + X*Q + Q*W + Y*Z))

Graph.mixed.2 <- makeMG(ug=UG(~y0*y1), 
                        dg=DAG(y4~y2, y2~y1), 
                        bg=UG(~y2*y3+y3*y4)) 


Graph.mixed.3 = makeMG(dg = DAG(W ~ X, Q ~ Z), 
                       bg = UG(~ X*Y + Y*Z + W*Q)) 

Graph.mixed.4 = makeMG(dg = DAG(W ~ X, Q ~ Z), 
                       bg = UG(~ X*Y + Y*Z + W*Q)) 


DAG.1.V<- matrix(c(0,1,0,0,0,0,0,0, 
               0,0,1,0,0,0,0,0, 
               0,0,0,0,0,0,0,0, 
               0,0,1,0,1,0,1,0, 
               0,0,0,0,0,1,0,0, 
               0,0,0,0,0,0,0,0, 
               0,0,0,0,0,1,0,0, 
               0,0,0,0,0,1,1,0),8,8,byrow=TRUE)

DAG.1.E <- c("a",1,2,
             "a",2,3,
             "a",4,3, 
             "a",4,5,
             "a",4,7,
             "a",5,6,
             "a",7,6,
             "a",8,6,
             "a",8,7)

DAG.1 <- DAG(y ~ x+z, z~u+v)
DAG.1.E <- edgematrix(DAG.1)
DAG.1.Matrix.1<-adjMatrix(DAG.1.E)
eigen(DAG.1.Matrix.1)

#--------------------------------Fit DAG Model-----------------------------------------------------------------
DAG.model.N<-10^2
DAG.model <- DAG(y ~ x-u, 
                 x ~ z, 
                 u ~ z)
DAG.model.S <- structure(c(2.93, -1.7, 0.76, -0.06,
                           -1.7, 1.64, -0.78, 0.1,
                           0.76, -0.78, 1.66, -0.78,
                           -0.06, 0.1, -0.78, 0.81), 
                           .Dim = c(4,4),
                          .Dimnames = list(c("y", "x", "z", "u"), c("y", "x", "z", "u")))

DAG.model.Fit<-fitDag(DAG.model, DAG.model.S, DAG.model.N)
DAG.model.Fit.Latent<-fitDagLatent(DAG.model, DAG.model.S, n=DAG.model.N, latent="u", seed=4564)
DAG.model.ancestral.graph<-AG(DAG.model,M="u")
DAG.model.Fit.2<-fitAncestralGraph(DAG.model.ancestral.graph, DAG.model.S, n = DAG.model.N)

#--------------------------------Analysis of the Conditional Models--------------------------------------------

DAG.model.max<-Max(DAG.model)
DAG.model.essential<-essentialGraph(DAG.model)
DAG.model.Induced.Covariance<-inducedCovGraph(DAG.model, sel=c("y"), cond=NULL)
DAG.model.Induced.Concentration<-inducedConGraph(DAG.model, sel=c("y"), cond="x")
DAG.model.Induced.Covariance.Overall<-inducedCovGraph(DAG.model)
DAG.model.Concentration.Overall<-inducedConGraph(DAG.model)

cc = list(c("y"))
DAG.model.Induced.Regression<-inducedRegGraph(DAG.model, sel=c("y"), cond=c("x"))
DAG.model.Induced.Chain.LWF<-inducedChainGraph(DAG.model,cc=cc,type="LWF")
DAG.model.Induced.Chain.AMP<-inducedChainGraph(DAG.model,cc=cc,type="AMP")
DAG.model.Induced.Chain.MRG<-inducedChainGraph(DAG.model,cc=cc,type="MRG")

#--------------------------------Tables----------------------------------------------------------
Table.1<-xtable(UG.1.matrix)
Table.2<-xtable(UG.formulas.df)

#--------------------------------Figures----------------------------------------------------------
lo<-layout_in_circle(graph.adjacency(UG.1.matrix))
par(mfrow = c(2,2), mar = c(0.5,1,0.5,0.5))
Figure.1A<-plot(graph.adjacency(UG.1.matrix))
Figure.1B<-plot(graph.adjacency(UG.2.matrix))
Figure.1C<-plot(graph.adjacency(UG.3.matrix))
Figure.1D<-plot(graph.adjacency(UG.4.matrix))

Figure.2<-plotGraph(DAG.model)

Figure.3A<-plotGraph(Graph.mixed.1)
Figure.3B<-plotGraph(Graph.mixed.2)
Figure.3C<-plotGraph(Graph.mixed.3)
Figure.3D<-plotGraph(Graph.mixed.4)

par(mfrow = c(2,2), mar = c(0.5,1,0.5,0.5))
Figure.4A<-plot(graph.adjacency(as.matrix(Graph.mixed.1)))
Figure.4B<-plot(graph.adjacency(as.matrix(Graph.mixed.2)))
Figure.4C<-plot(graph.adjacency(as.matrix(Graph.mixed.3)))
Figure.4D<-plot(graph.adjacency(as.matrix(Graph.mixed.4)))


par(mfrow = c(2,2), mar = c(0.5,1,0.5,0.5))
Figure.5A<-plot(graph.adjacency(as.matrix(DAG.model.Fit$Ahat)))
Figure.5B<-plot(graph.adjacency(as.matrix(DAG.model.Fit$Shat)))
Figure.5C<-plot(graph.adjacency(as.matrix(DAG.model.Fit.2$Ahat)))
Figure.5D<-plot(graph.adjacency(as.matrix(DAG.model.Fit.2$Shat)))

par(mfrow = c(2,2), mar = c(0.5,1,0.5,0.5))
Figure.6A<-plot(graph.adjacency(as.matrix(DAG.model.max)))
Figure.6B<-plot(graph.adjacency(as.matrix(DAG.model.essential)))
Figure.6C<-plot(graph.adjacency(as.matrix(DAG.model.Fit$Shat)))
Figure.6D<-plot(graph.adjacency(as.matrix(DAG.model.Fit$Ahat)))

par(mfrow = c(2,2), mar = c(0.5,1,0.5,0.5))
Figure.7A<-plot(graph.adjacency(as.matrix(DAG.model.Induced.Covariance)))
Figure.7B<-plot(graph.adjacency(as.matrix(DAG.model.Induced.Concentration)))
Figure.7C<-plot(graph.adjacency(as.matrix(DAG.model.Induced.Covariance.Overall)))
Figure.7D<-plot(graph.adjacency(as.matrix(DAG.model.Induced.Concentration.Overall)))

par(mfrow = c(2,2), mar = c(0.5,1,0.5,0.5))
Figure.8A<-plot(graph.adjacency(as.matrix(DAG.model.Induced.Regression)))
Figure.8B<-plot(graph.adjacency(as.matrix(DAG.model.Induced.Chain.LWF)))
Figure.8C<-plot(graph.adjacency(as.matrix(DAG.model.Induced.Chain.AMP)))
Figure.8D<-plot(graph.adjacency(as.matrix(DAG.model.Induced.Chain.MRG)))

#------------------------------------------------------References------------------------------------------------------

Reference.1<-c("","","")

#-----------------------------------------------------Function Library------------------------------------------------