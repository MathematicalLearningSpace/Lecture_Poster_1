
library(xtable);library(igraph);library(Matrix);library(qgraph);library(qlcMatrix);library(rcrossref);library(PearsonDS);library(catnet);library(CCA)
#-------------------------------------------------Data Matrix------------------------------------------------------
A.matrix<-matrix(c(0.5,0.3,0.4,0.7,0.7,1,
            0.7,0.5,1,0.8,0.9,0.3,
            0.6,0,0.5,0.5,0.6,0.2,
            0.3,0.2,0.5,0.5,0.3,0.5,
            0.3,0.1,0.4,0.7,0.5,0.3,
            0,0.7,0.8,0.5,0.7,0.5),nrow=6,ncol=6)

#------------------------------------------------EigenValue Analysis with Skewness---------------------------------

A.SyK<-symmpart(A.matrix) + skewpart(A.matrix) 
A.Sy<-symmpart(A.matrix)
A.SK<-skewpart(A.matrix)
A.Sy.Eigen<-eigen(A.Sy)
A.SK.Eigen<-eigen(A.SK)

A.matrix.rho.max<-eigen(A.matrix)$values[1]
A.matrix.S<-symmpart(t(A.matrix))
A.matrix.skew<-skewpart(t(A.matrix)) 

A.delta<-eigen(A.matrix.S)
A.delta.max<-max(A.delta$values)

A.singular<-eigen(A.matrix.skew)
A.singular.max<-A.singular$values[1]

w<-matrix(1,6,1)
#------------------------Matrix Operations----------------------------------
A.Decomp<-w%*%t(w)
A.Var<-(t(w)%*%t(A.SK)%*%A.SK%*%w)/t(w)%*%w
A.alpha<-A.Sy+(1-2*alpha)*A.SK
A.alpha.variance<-(1-2*alpha)^2*A.Var
#
#----------------------------------------------------------Holder-------------------------------------------------------
#
A.Holder<-holder(A.matrix,1,1,1)
A.spectralRadius<-eigen(A.matrix)$values
A.Sy.spectralRadius<-max(A.Sy.Eigen$values)
A.SK.spectralRadius<-max(A.SK.Eigen$values)

#----------------------------------------------------------Perron Vectors-----------------------------------------------
holder<-function(A,a,b,p)
{
  Hp<-(((a^p)*A+(b^p)*t(A))/2)^(1/p)
  return(Hp)
}

alpha<-1/3
A.alpha.Perron.root<-eigen((1-alpha)*A.matrix+alpha*t(A.matrix))
A.alpha.Perron.root.max<-A.alpha.Perron.root$values[1]
A.alpha.Perron.root.vectors.left<-A.alpha.Perron.root$vectors[,1]
A.alpha.Perron.root.vectors.right<-A.alpha.Perron.root$vectors[,5]
A.alpha.Perron.root.Prime<--2*((t(A.alpha.Perron.root.vectors.left)%*%A.matrix.skew%*%A.alpha.Perron.root.vectors.right)/
                                   (t(A.alpha.Perron.root.vectors.left)%*%A.alpha.Perron.root.vectors.right))

Perron.alpha.df<-data.frame()
Perron.alpha.df<-cbind(alpha,min(Re(A.alpha.Perron.root$values)),
                       max(Re(A.alpha.Perron.root.max)))
colnames(Perron.alpha.df)<-c("alpha","Min Eigen","Max Root")


#---------------------------------------------------------Correlation Structure--------------------------
graph.sample.1<-make_graph("Petersen")
graph.adjM.1<-as_adjacency_matrix(graph.sample.1,sparse=FALSE)
graph.adjM.1.eigen<-eigen(graph.adjM.1)
G.Petersen.S<-symmpart(t(graph.adjM.1))
G.Perersen.K<-skewpart(t(graph.adjM.1))

graph.sample.2<-make_ring(10)
graph.adjM.2<-as_adjacency_matrix(graph.sample.2,sparse=FALSE)
#----------------------------------------Second Moment Analysis-------------------------------------------
matrix.correlation<-cor(c(as.matrix(graph.adjM.1)), c(as.matrix(graph.adjM.2)))
matrix.correlation.test<-cor.test(c(as.matrix(graph.adjM.1)), c(as.matrix(graph.adjM.2)))
cancor(as.matrix(graph.adjM.1), as.matrix(graph.adjM.2))

results.correlation.table.header<-c('Graphs','Correlation','Test')
results.correlation.table.row<-c('Petersen-Ring',
                                 matrix.correlation,
                                 matrix.correlation.test$statistic,
                                 matrix.correlation.test$pvalue)
results.correlation.table<-rbind(results.correlation.table.header,
                                 results.correlation.table.row)


#---------------------------------------------Tables------------------------------------------------------

Table.1<-xtable(Perron.alpha.df)
Table.2<-xtable(results.correlation.table)

#---------------------------------------------Figures-----------------------------------------------------

Figure.1<-plot(alpha.range,as.numeric(A.alpha.Perron.root.max),
     type='l',xlab="alpha",ylab="Perron Root",
     main="Perron Root and Alpha")
lines(alpha.range,A)

Figure.2<-plot(graph.sample, 
     vertex.color="Blue",
     vertex.size=20,
     vertex.label.dist=4,
     edge.arrow.size=0.5,
     main='Petersen Graph')

#--------------------------------------------Function Library---------------------------------------------
f.1<-function(X)
{
A.alpha.Perron.root.max<-NULL
alpha.range<-seq(0,0.5,by=0.01)
A<-NULL
for(i in 1:length(alpha.range))
{
  A.alpha.Perron.root<-eigen((1-alpha.range[i])*X+alpha.range[i]*t(X))
  A.alpha.Perron.root.max[i]<-A.alpha.Perron.root$values[1]
  A[i]<-(1-alpha.range[i])+alpha.range[i]
}
output<-list()
output$X<-X
output$A.alpha.Perron.root<-A.alpha.Perron.root
output$A.alpha.Perron.root.max<-A.alpha.Perron.root.max
output$A<-A
return(output)
}
test.f.1<-f.1(A.matrix)
test.f.1

