library(xtable)
library(visNetwork)
library(igraph)
#-----------------------------------------Data-------------------------------------------
Petersen.graph <- graph.famous("Petersen")
Petersen.graph.translation<-toVisNetworkData(Petersen.graph )
#-----------------------------------------Poster Blocks---------------------------------
nodes.df <- data.frame(id = 1:16, 
                    label = paste("Poster Block", 1:16),        
                    group = c("GrA", "GrB"),                                     
                    value = 1:16,                                             
                    shape = c("square", "triangle", "box", "circle", "dot", "star",
                              "ellipse","database"),                   
                    color = c("darkblue","orange", "blue", "purple"),
                    shadow = c(FALSE, TRUE, FALSE, TRUE))                 
edges.df <- data.frame(
  from = c(1,1,1,1,1,1,1,
           2),
  to = c(2,3,4,5,6,7,8,
         9),
  label = paste("E-On", 1:8),
  length=seq(1:8)
)
math.lectures <- visNetwork(nodes.df, edges.df, width = "100%")
visSave(lectures, file = "math_lectures.html")
visSave(Petersen.graph.translation, file = "Petersen_graph.html")
#-----------------------------------------Tables-----------------------------------------

Table.1<-xtable(nodes.df)
Table.2<-xtable(edges.df)

#-----------------------------------------Figures---------------------------------------
Figure.1<-visNetwork(nodes.df, edges.df, width = "100%") %>% 
  visInteraction(navigationButtons = TRUE)%>%
  visOptions(manipulation = TRUE)%>%
  visEdges(arrows = 'to')

Figure.2<-visNetwork(Petersen.graph.translation$nodes, Petersen.graph.translation$edges, width = "100%") %>% 
  visInteraction(navigationButtons = TRUE)%>%
  visOptions(manipulation = TRUE)%>%
  visEdges(arrows = 'to')

#----------------------------------------Function Library------------------------------