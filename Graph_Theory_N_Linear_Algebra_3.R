#-----------------------------------R Code to be Updated in the Classroom Lecture-------------------
library(Matrix);library(matrixcalc);library(microbenchmark);library(igraph);library(xtable);library(Hmisc);library(devtools);library(XML)
library(inline);library(utils);library(png);library(qlcMatrix);library(qgraph);library(rcrossref);library(pracma)
library(PearsonDS);library(RandomFields)
#------------------------------------------------------Graph Design----------------------------------------------------------
g.ABCD<- make_full_graph(4) %du% make_full_graph(5) %du% make_full_graph(5) %du% make_full_graph(5)
g.ABCD <- add_edges(g.ABCD, c(1,6, 
                              1,11, 
                              6,11, 
                              6,15))
g.1 <- make_full_graph(4) %du% make_full_graph(5) %du% make_full_graph(3)
g.2 <- make_full_graph(3) %du% make_full_graph(5) %du% make_full_graph(5)
g.3 <- make_full_graph(4) %du% make_full_graph(8) %du% make_full_graph(5)
g.4 <- make_full_graph(2) %du% make_full_graph(5) %du% make_full_graph(5)

g.1 <- add_edges(g.1, c(1,6, 1,11, 6, 11))
g.2 <- add_edges(g.2, c(1,6, 1,11, 6, 11))
g.3 <- add_edges(g.3, c(1,6, 1,11, 6, 11))
g.4<-add_edges(g.4, c(1,6, 1,11, 6, 11))
#------------------------------------------------------------Assign Weights---------------------------------------------------------

E(g.ABCD)$weight<-runif(ecount(g.ABCD))
E(g.1)$weight<-runif(ecount(g.1),0,1)
E(g.2)$weight<-runif(ecount(g.2),0,1)
E(g.3)$weight<-runif(ecount(g.3),-3,3)
E(g.4)$weight<-runif(ecount(g.4),0,1)


#------------------------------------------------------------Topological Measures spectra--------------------------------------------
graph.alpha<-alpha_centrality(g.ABCD)
graph.power<-power_centrality(g.ABCD)
graph.hub_score<-hub_score(g.ABCD)$vector
graph.diameter<-diameter(g.ABCD)

graph.topology.df<-data.frame()
graph.topology.df<-cbind(graph.alpha,
                         graph.power,
                         graph.hub_score,
                         graph.diameter)

colnames(graph.topology.df)<-c("Alpa","Power","Hub Score","Diameter")

centralization.1<-centr_eigen(g.ABCD)$centralization
centralization.2<-centr_eigen(g.1)$centralization
centralization.3<-centr_eigen(g.2)$centralization
centralization.4<-centr_eigen(g.3)$centralization
centralization.5<-centr_eigen(g.4)$centralization

graph.centralization.df<-data.frame()
graph.centralization.df<-cbind(centralization.1,
                         centralization.2,
                         centralization.3,
                         centralization.4,
                         centralization.5)
colnames(graph.centralization.df)<-c("1","2","3","4","5")

#--------------------------------------------------------------Clustering Leading Eignevalues-----------------------------------

cle.g.ABCD<-cluster_leading_eigen(g.ABCD)
cle.g.1<-cluster_leading_eigen(g.1)
cle.g.2<-cluster_leading_eigen(g.2)
cle.g.3<-cluster_leading_eigen(g.3)
cle.g.4<-cluster_leading_eigen(g.4)
#-------------------------------------------------------------Walking the Graph------------------------------------------------------
cw.g<-cluster_walktrap(g.ABCD)
eb.g<-edge_betweenness(g.ABCD)
cfg.g<-cluster_fast_greedy(g.ABCD)
cle.g<-cluster_leading_eigen(g.ABCD)

graph.cluster.df<-data.frame()
graph.cluster.df<-cbind(cw.g$modularity,
                        eb.g,
                        cfg.g$modularity,
                        cle.g$modularity)
colnames(graph.cluster.df)<-c("CW","EB","CFG","CLE")

#--------------------------------------------------------------Table--------------------------------------------------------------------------------

Table.1<-xtable(graph.topology.df)
Table.2<-xtable(graph.centralization.df)
Table.3<-xtable(graph.cluster.df)


#-------------------------------------------------------------Figures-------------------------------------------------------------------------------

par(mfrow = c(1,1), mar = c(0.25,1,0.25,0.25))
Figure.1<-plot(g.ABCD,vertex.color="Green",
     vertex.size=10,vertex.label.dist=2.5,
     edge.arrow.size=0.5,main="")
#--------------------------------------------Figure Group-------------------------------
par(mfrow = c(2,2), mar = c(2,2,2,2))
Figure.2<-plot(diversity(g.ABCD), type="l",ylim=c(0,1), xlab="Node Index",
     main="(A): Node Diversity by Graph Type")
lines(diversity(g.1), col="blue")
lines(diversity(g.2), col="green")
lines(diversity(g.3), col="red")
lines(diversity(g.4), col="orange")
legend("bottom", 
       lwd = c(2, 1), 
       lty = c(1, 4),
       col=c("black","blue","green","red","orange"),
       bty = "n",
       legend = c("A", "B", "C","D"))

Figure.3<-plot(hub_score(g.ABCD)$vector, type="l",ylim=c(0,1), xlab="Node Index",
     main="(B): Hub Score by Graph Type")
lines(hub_score(g.1)$vector, col="blue")
lines(hub_score(g.2)$vector, col="green")
lines(hub_score(g.3)$vector, col="red")
lines(hub_score(g.4)$vector, col="orange")
legend("bottom", 
       lwd = c(2, 1), 
       lty = c(1, 4),
       col=c("black","blue","green","red","orange"),
       bty = "n",
       legend = c("A", "B", "C","D"))

Figure.4<-plot(alpha_centrality(g.ABCD), type="l",ylim=c(-2,2), xlab="Node Index",
     main="(C): Alpha Centrality by Graph Type")
lines(alpha_centrality(g.1), col="blue")
lines(alpha_centrality(g.2), col="green")
lines(alpha_centrality(g.3), col="red")
lines(alpha_centrality(g.4), col="orange")
legend("bottom", 
       lwd = c(2, 1), 
       lty = c(1, 4),
       col=c("black","blue","green","red","orange"),
       bty = "n",
       legend = c("A", "B", "C","D"))

Figure.5<-plot(authority_score(g.ABCD)$vector, type="l",ylim=c(0,1), xlab="Node Index",
     main="(D):Authority Score by Graph Type")
lines(authority_score(g.1)$vector, col="blue")
lines(authority_score(g.2)$vector, col="green")
lines(authority_score(g.3)$vector, col="red")
lines(authority_score(g.4)$vector, col="orange")
legend("bottom", 
       lwd = c(2, 1), 
       lty = c(1, 4),
       col=c("black","blue","green","red","orange"),
       bty = "n",
       legend = c("A", "B", "C","D"))

par(mfrow = c(3,3), mar = c(2,2,2,2))
Figure.6<-plot_dendrogram(cle.g)
plot_dendrogram(cle.g.1)
plot_dendrogram(cle.g.2)
plot_dendrogram(cle.g.3)
plot_dendrogram(cle.g.4)

#-------------------------------------------------------------Function Library----------------------------------------------------------------------



