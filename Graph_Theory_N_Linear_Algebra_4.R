#-----------------------------R Code To Modify in the Classroom Lecture with Students-----------------------
#---------------------------------------------------R API --------------------------------------------------------
library(xtable);library(igraph);library(Matrix);library(readr);library(readxl);library(stringr);library(stringi)
#------------------------------------------------Data from Gene Cards and String DB-----------------------------------------------------------
UBC.Enrichment.Components<-read.delim("UBC/enrichment.Component.tsv", header=FALSE, comment.char="#")
UBC.Enrichment.KEGG<-read.delim("UBC/enrichment.KEGG.tsv", header=FALSE, comment.char="#")
UBC.Network.Process <- read.delim("UBC/enrichment.Process.tsv", header=FALSE, comment.char="#")
UBC.Network.Process$V2
UBC.Network.regulation.positive<-UBC.Network.Process$V2[grep("positive regulation",UBC.Network.Process$V2)]
UBC.Network.regulation.negative<-UBC.Network.Process$V2[grep("negative regulation",UBC.Network.Process$V2)]
UBC.Network.regulation.ratio<-length(grep("regulation",UBC.Network.Process$V2))/length(UBC.Network.Process$V2)

UBC.Network.Process.Biological.df<-data.frame()
UBC.Network.Process.DNA.Damage<-UBC.Network.Process$V6[grep('DNA damage', UBC.Network.Process$V2, ignore.case=TRUE)]
UBC.Network.Process$V6
#------------------------------------------------Graph---------------------------------------------------------
nodes<-str_split(UBC.Network.Process.DNA.Damage[3],",")
UBC.network.study<-nodes[[1]]
UBC.vertex.n<-length(nodes[[1]])
UBC.graph<-graph_from_edgelist(cbind(nodes[[1]],nodes[[1]]))
UBC.graph.full<-make_full_graph(UBC.vertex.n)

v<-V(UBC.graph.full)
e<-E(UBC.graph.full)
e.n<-length(E(UBC.graph.full))
e.weight<-E(UBC.graph.full)$weight

vertex_attr(UBC.graph.full)<-list(name=unique(nodes[[1]]))
vertex_attr(UBC.graph.full,"label")<-V(UBC.graph.full)$name

#-----------------------------------------------Transformations on Graphs------------------------------------
prob.low<-1/UBC.vertex.n
prob.middle<-3/UBC.vertex.n
prob.high<-5/UBC.vertex.n
edge.prob<-c(prob.low, prob.middle, prob.high)

UBC.graph.2<-UBC.graph.full%>%delete.edges(seq(1, length(E(UBC.graph.full)), by = 2))
UBC.graph.3<-UBC.graph.2%>%set_edge_attr("weight", value = 1:length(E(UBC.graph.2)))%>%
  set_edge_attr("color", value = "green")

UBC.network.adjacency<-as_adjacency_matrix(UBC.graph.3,attr="weight")
UBC.network.adjacency.eigenvalues<-eigen(UBC.network.adjacency)
UBC.network.adjacency.eigenvalues$values
#------------------------------------------------Tables--------------------------------------------------------

#-----------------------------------------------Figures--------------------------------------------------------

lo<-layout_in_circle(UBC.graph.full)
Figure.1<-plot(UBC.graph.full,layout=lo, vertex.color="Blue",vertex.size=20,
     vertex.label.dist=5,edge.arrow.size=0.5,main='UBC Network')
lo<-layout_in_circle(UBC.graph.full)
Figure.1.caption<-stri_join(UBC.Network.Process.DNA.Damage[3]," ",UBC.Enrichment.Components$V2[1])

Figure.2<-plot(UBC.graph.3,layout=lo, vertex.color="Blue",vertex.size=20,
               vertex.label.dist=5,edge.arrow.size=0.5,main='UBC Network',edge.label = E(UBC.graph.3)$weight)

#-----------------------------Function Library for the Classroom-----------------------------------------------
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



