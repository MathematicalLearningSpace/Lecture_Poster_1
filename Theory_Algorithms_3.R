library(NMOF)
library(xtable)
library(recommenderlab)
library(corrplot)
#------------------------------------Data----------------------------------------------------------------

#-----------------------------------Curated Data Set for Recommendations-----------------------------------------------------
v1 <- rnorm(20)
v2 <- runif(50)
v3 <- rbinom(100, size = 1, prob = 0.4)
#----------------------------------Correlation matrix----------------------------------------------------
matrix.correlation <- array(0.5, dim = c(3, 3))
diag(matrix.correlation ) <- 1
cor.Spearman.1<-cor(resampleC(a = v1, b = v2, v3, size = 100, cormat = matrix.correlation),
    method = "spearman")
cor.spearman.eigenvalues.1<-eigen(cor.Spearman.1, only.values = TRUE)
cor.spearman.eigenvalues.df<-data.frame()
cor.spearman.eigenvalues.df<-rbind(c(cor.spearman.eigenvalues.1$values))
colnames(cor.spearman.eigenvalues.df)<-c("E1","E2","E3")

#------------------------------------A: Differential Evolution---------------------------------------------
Objective.Function.1 <- f
algo <- list(nP = 50L,          ### population size
             nG = 300L,         ### number of generations
             F = 0.6,          ### step size
             CR = 0.9,          ### prob of crossover
             min = c(-10, -10),  ### range for initial population
             max = c( 10,  10))
DE.sol <- DEopt(OF = Objective.Function.1, algo = algo)
DE.sol.sd<-sd(DE.sol$popF)
DE.sol$xbest
DE.sol$OFvalue
#------------------------------------B: Genetic Algorithms-------------------------------------------------
#---------------------------------------String Match-------------------------------------------------------
size <- 30L 
Objective.Function.2 <- function(x, y) sum(x != y)
#y <- rnorm(size,0,2) > 0.5
y <- rbinom(size, size = 1, prob = 0.5)
x <- runif(size) > 0.5

similarity.ratio<-Objective.Function.2(x,y)/size

algo.1 <- list(nB = size, nP = 20L, nG = 300L, prob = 0.002,
             printBar = TRUE)
GA.sol <- GAopt(OF=Objective.Function.2 , algo = algo.1,y=y)
GA.sol$xbest
GA.sol$OFvalue

cat(as.integer(y), "\n", 
    as.integer(GA.sol$xbest), "\n",
    ifelse(y == GA.sol$xbest , " ", "*"), "\n", 
    sep = "")

#------------------------------------C: Grid Search----------------------------------------------------
f.1 <- function(x) x[1L] + x[2L]^5

levels <- list(a = 1:5, b = 1:5)

GS.sol <- gridSearch(fun = f.1, levels)
GS.sol$minfun
GS.sol$minlevels
GS.sol$values

Algorithm.optimization.df<-data.frame()
Algorithm.optimization.df<-rbind(c("DE",algo$nB,algo$nG,algo$nP,algo$F,algo$CR,algo$min[1],algo$max[1], 
                                   sd(DE.sol$popF),DE.sol$xbest,DE.sol$OFvalue),
                                 c("GA",algo.1$nB,algo.1$nG,algo.1$nP,"",algo.1$prob,"","", 
                                   sd(GA.sol$popF),GA.sol$xbest,GA.sol$OFvalue))
colnames(Algorithm.optimization.df)<-c("Algorithm","nB","nG","nP","F","CR","Min","Max","SD","xbest","Fvalue")

#------------------------------------Tables---------------------------------------------------------

Table.1<-xtable(Algorithm.optimization.df)
Table.2<-xtable(cor.spearman.eigenvalues.df)
#------------------------------------Figures--------------------------------------------------------
Figure.1<-corrplot(cor.Spearman.1)
Figure.2<-ts.plot(DE.sol$Fmat, xlab = "generations", ylab = "OF")
Figure.3<-ts.plot(GA.sol$Fmat)
Figure.4<-ts.plot(GS.sol$values)

#------------------------------------Function Library-----------------------------------------------

f<-function(x)
{
  y<-x[2]
  x<-x[1]
  a<-2^0;b<-(2^2+2^0);c<-2^3;d<-2^0;
  a*sin(b*x) + c*cos(d*y)
}
