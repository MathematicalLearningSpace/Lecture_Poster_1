#----------------------R Code To Modify in the Classroom Lecture with Students-----------------------
library(xtable);library(igraph);library(Matrix);library(qgraph);library(qlcMatrix);library(rcrossref);library(PearsonDS);library(catnet)
#-------------------------------Graph Theory N Linear Algebra----------------------------------------
#-----------------------Synthetic Data Generation for Classroom Lectures-----------------------------
data.N<-100
set.seed(10)
x1 <- rnorm(data.N,0,1)
x2 <- x1 + rnorm(data.N,0,0.2)
x3 <- x1 + x2 + rnorm(data.N,0,0.2)
x4 <- rnorm(data.N,0,1)
x5 <- x4 + rnorm(data.N,0,0.4)
x6 <- x4 + rnorm(data.N,0,0.4)
x7 <- x1 + x5 + rnorm(data.N,0,0.1)

experimental.Data.1 <- cbind(x1,x2,x3,x4,x5,x6,x7)

#--------------Examples of Graph Motifs for Topological Metrics-----------------
graph.name<-c('Barbell','Balaban','Book', 'Centipede', 'Complete BiPartite','Bull',
              'Cage','Cayley','Cocktail Party','Coxeter','Chvatal',
              'Crown','Cubic Symmetric', 'Cycle','Desargues','Doyle',
              'Folkman','Foster','Franklin','Gear','Gray','Groetzsch',
              'Heawood','Helm','Hoffman-Singleton','Hypercube','Konigsberg',
              'Levi','McGee','Mobius-Kantor','Pan','Pappus','Path','Petersen',
              'Prism','Tetrahedral','Octahedral','Cubical','Icosahedral','Dodecahedral',
              'Robertson','Star','Sun','Sunlet','Tutte-Coxeter','Web','Wheel','Wong')

graph.name.Petersen<-make_graph('Petersen')
graph.Petersen.adjM<-as_adjacency_matrix(graph.name.Petersen,sparse=FALSE)

graph.name.Petersen.generalized.adjM<-matrix(
  c(0,0,0,0,0,0,
    1,0,0,0,0,0,
    1,1,0,0,0,0,
    1,1,1,0,0,0,
    1,1,1,1,0,0,
    1,1,1,1,1,0),6,6)
graph.Petersen.adjM<-as_adjacency_matrix(graph.name.Petersen.generalized.adjM,sparse=FALSE)
graph.name.Petersen.generalized<-graph_from_adjacency_matrix(graph.name.Petersen.generalized.adjM)
g<-graph.name.Petersen

#------------------------------------------Analysis----------------------------------------------------------------
#----------------------------To Be Completed in Class By Students--------------------------------------------------
#-----------------------------------------Distance----------------------------------------------------------------
Graph.Topology.Metrics.Distance.df<-data.frame()
graph.mean_distance<-mean_distance(g)
graph.diameter<-diameter(g)
Graph.Topology.Metrics.Distance.df<-cbind(graph.mean_distance,graph.diameter)
Graph.Topology.Metrics.Distance.df.Names<-c('Mean Distance',"Diameter")
colnames(Graph.Topology.Metrics.Distance.df)<-Graph.Topology.Metrics.Distance.df.Names

#----------------------------------------Connectivity------------------------------------------------------------
Graph.Topology.Metrics.Connectivity.df<-data.frame()
graph.eccentricigraph.degree.all<-degree(g, v = V(g), mode = c("all"),loops = TRUE, normalized = FALSE)
graph.degree.out<-degree(g, v = V(g), mode = c("out"),loops = TRUE, normalized = FALSE)
graph.degree.in<-degree(g, v = V(g), mode = c("in"),loops = TRUE, normalized = FALSE)
graph.degree.total<-degree(g, v = V(g), mode = c("total"),loops = TRUE, normalized = FALSE)
graph.degree_distribution<-degree_distribution(g)
graph.edge_conn<-edge_connectivity(g)
graph.vertex<-vertex_connectivity(g)
graph.adhesion<-adhesion(g)
graph.cohesion<-cohesion(g)

Graph.Topology.Metrics.Connectivity.df<-rbind(graph.eccentricity,
                                              graph.degree.out,
                                              graph.degree.in,
                                              graph.degree.total,
                                              graph.degree_distribution,
                                              graph.edge_conn,
                                              graph.vertex,
                                              graph.adhesion,
                                              graph.cohesion
                                              )
Graph.Topology.Metrics.Connectivity.df.Names<-c('Eccentricity',
                                                "Degree.out",
                                                "Degree.in",
                                                "Degree.total",
                                                "Degree_distribution",
                                                "Edge_conn",
                                                "Vertex",
                                                "Adhesion",
                                                "Cohesion")
colnames(Graph.Topology.Metrics.Connectivity.df)<-Graph.Topology.Metrics.Connectivity.df.Names

#---------------------------------------Spectra-----------------------------------------------------------------
Graph.Topology.Metrics.Spectra.df<-data.frame()
#graph.alpha<-alpha_centrality(g)
#graph.power<-power_centrality(g)
graph.eigen<-eigen_centrality(g)$vector
graph.hub_score<-hub_score(g)$vector
graph.authority<-authority_score(g)$vector
graph.automorphism<-automorphisms(g)
graph.components <- decompose(g, min.vertices=2)
Graph.Topology.Metrics.Spectra.df<-cbind(graph.eigen,
                                         graph.hub_score,
                                         graph.authority,
                                         graph.automorphism$group_size,
                                         graph.components)
Graph.Topology.Metrics.Spectra.df.Names<-c('Eigen Centraility',
                                           "Hub Score",
                                           "Authority",
                                           "Group Size",
                                           "Graph Components")
colnames(Graph.Topology.Metrics.Spectra.df)<-Graph.Topology.Metrics.Spectra.df.Names
#---------------------------------------Miscellaneous-----------------------------------------------------------------
Graph.Topology.Metrics.Miscellaneous.df<-data.frame()
Gc.maxParents<-3
Gc.numCategories<-2
cnet <- cnNew(
  nodes =  c("a", "b", "c"),
  #biological Processes and Gene Function in development
  cats = list(c("1","2"), c("1","2"), c("1","2")), 
  parents = list(NULL, c(1), c(1,2)), 
  probs = list(	c(0.2,0.8), 
                list(c(0.6,0.4),c(0.4,0.6)), 
                list(list(c(0.3,0.7),c(0.7,0.3)),list(c(0.9,0.1),c(0.1,0.9))))
)
Gc <- cnRandomCatnet(numnodes=length(V(g)), maxParents=Gc.maxParents, numCategories=Gc.numCategories)
Graph.complexity<-cnComplexity(Gc)
Graph.entropy<-cnEntropy(Gc)

Graph.Topology.Metrics.Miscellaneous.df<-cbind(graph.complexity,graph.entropy)
Graph.Topology.names<-c('Complexity','Entropy')
names(Graph.Topology.Metrics.Miscellaneous.df)<-Graph.Topology.names
#------------------------------------------Correlation Analysis--------------------------------------------------

CorMat<-cor_auto(experimental.Data.1)
# Compute graph with tuning = 0 (BIC):
BICgraph <- qgraph(CorMat, graph = "glasso", sampleSize = nrow(experimental.Data.1),
                   tuning = 0, layout = "spring", title = "BIC", details = TRUE)

# Compute graph with tuning = 0.5 (EBIC)
EBICgraph <- qgraph(CorMat, graph = "glasso", sampleSize = nrow(experimental.Data.1),
                    tuning = 0.5, layout = "spring", title = "BIC", details = TRUE)

# Compare centrality and clustering:
centralityPlot(list(BIC = BICgraph, EBIC = EBICgraph))
clusteringPlot(list(BIC = BICgraph, EBIC = EBICgraph))

# Compute centrality and clustering:
centrality_auto(BICgraph)
clustcoef_auto(BICgraph)

smallworldness(BICgraph)

#------------------------------------------Tables-----------------------------------------------------------------
Table.1<-xtable(Graph.Topology.Metrics.Distance.df)                                        
Table.2<-xtable(Graph.Topology.Metrics.Connectivity.df) 
Table.3<-xtable(Graph.Topology.Metrics.Spectra.df) 
Table.4<-xtable(Graph.Topology.Metrics.Miscellaneous.df) 

#-----------------------------------------Figures-----------------------------------------------------------------

Figure.1<-plot(graph.name.Petersen, 
     layout=layout_in_circle,
     vertex.color="Blue",vertex.size=10,
     vertex.label.dist=2,edge.arrow.size=0.5,
     main='Petersen Graph')

Figure.2<-plot(graph.name.Petersen.generalized, 
     layout=layout_in_circle,
     vertex.color="Blue",vertex.size=10,
     vertex.label.dist=2,edge.arrow.size=0.5,
     main='Generalized Petersen Graph')

colors <- c("red", "red", "blue", "blue", "white",'green')
Figure.3<-qgraph(cor(experimental.Data.1,method="pearson")
           ,layout="spring"
           ,label.cex=0.9
           ,labels=colnames(df)
           ,label.scale=F
           ,details=T
           ,edge.labels=F
           ,doNotPlot=T
           ,alpha=0.05
           ,minimum='sig'
           ,sampleSize=100,
           colors=colors)


#----------------------------------------Function Library---------------------------------------------------------

