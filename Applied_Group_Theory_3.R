library(combinat)
library(permutations)
library(permute)
library(xtable)
library(igraph)
library(Matrix)
library(vegan)
library(zoo)
library(adegenet)
library(ActiveDriver)
library(VennDiagram)
library(venneuler)
#------------------------------------------------------Data--------------------------------------------------------------

Species.Population.Category.n <- c(100,20,10)
Species.population.matrix.probability <- matrix(c(.3,.1,.5,
                                                  .1,.1,.2,
                                                  .6,.8,.3),3)
Probability.Matrix <- as.matrix(data.frame(A=c(0.1, 0.2, 0.4, 0.2, 0.1),
                                           B=c(0, 0.2, 0, 0.8, 0),
                                           C=c(0, 0, 0.3, 0, 0.7)))

VE.Probability.Matrix.1<-venneuler(Probability.Matrix)
VE.Probability.Matrix.2<-venneuler(Probability.Matrix > 0)
#----------------------------------------------------Combinatorics------------------------------------------------------

Species.population.matrix.Counts<-rmultinomial(Species.Population.Category.n,Species.population.matrix.probability)

#----------------------------------------------------Symmetric Groups, Graphs, and Matrices----------------------------------------------------
S4<-allperms(4)
S4.cycle<-as.cycle(allperms(4))
S4.matrix<-as.matrix(S4)
S4.matrix.1<-S4.matrix[,1]
S4.matrix.2<-S4.matrix[,2]
S4.matrix.3<-S4.matrix[,3]
S4.matrix.4<-S4.matrix[,4]

A4.Tetrahedron<-S4[is.even(S4)]
A4.Tetrahedron.matrix<-as.matrix(A4.Tetrahedron)
A4.Tetrahedron.matrix.word<-word(A4.Tetrahedron.matrix)

Coxeter.graph<-make_graph("Coxeter")
Coxeter.matrix<-as_adjacency_matrix(Coxeter.graph)
Coxeter.matrix.eigenvalues<-eigen(Coxeter.matrix)
Coxeter.matrix.eigenvalues$values

Octahedral.graph<-make_graph("Octahedral")
Octahedral.matrix<-as_adjacency_matrix(Octahedral.graph)
Octahedral.matrix.eigenvalues<-eigen(Octahedral.matrix)
Octahedral.matrix.eigenvalues$values

Species.population.matrix.probability.graph<-graph.adjacency(Species.population.matrix.probability)

#------------------------------------------------------Tables-----------------------------------------------------------
Table.1<-xtable(Species.population.matrix.probability)
Table.2<-xtable(Species.population.matrix.Counts)
#----------------------------------------------------Figures------------------------------------------------------------

Figure.1<-plot(Species.population.matrix.probability.graph, main="Species Probability Matrix")
Figure.2<-plot(Coxeter.graph, main="Coxeter Graph")
Figure.3<-plot(Octahedral.graph, main="Octahedral Graph")
Figure.4<-plot(VE.Probability.Matrix.1)
Figure.5<-plot(VE.Probability.Matrix.2)

Figure.6<-plot(S4.matrix.1, type="l", xlab="Placement",
               main="S4 Permutation Sequences by column")
lines(S4.matrix.2, col="blue")
lines(S4.matrix.3, col="green")
lines(S4.matrix.4, col="red")
legend("topright", 
       lwd = c(2, 1), 
       lty = c(1, 4),
       col=c("black","blue","green","red"),
       bty = "n",
       legend = c("A", "B", "C","D"))

#---------------------------------------------------References---------------------------------------------------------

Reference.1<-c("Reingold, E.M., Nievergelt, J., Deo, N. (1977)", 
               "Combinatorial Algorithms: Theory and Practice.", 
               "NJ: Prentice-Hall. pg. 170.")

#---------------------------------------------------Function Library---------------------------------------------------




