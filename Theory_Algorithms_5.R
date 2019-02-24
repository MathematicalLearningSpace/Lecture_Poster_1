#-----------------------------R Code To Modify in the Classroom Lecture with Students-----------------------
#----------------------------------R API -------------------------------------------------------------------
library(NMOF);library(xtable)
#------------------------------------Data----------------------------------------------------------
N<-100
categories<-11
categories.desired<-5

X<-array(rnorm(N*categories),dim=c(N,categories))
X.desired<-sample(1:categories,categories.desired)
X.subset<-X[,X.desired,drop=FALSE]

y<-rowMeans(X.subset)

Data <- list(X = X, y = y, nc = categories, nr = N, n = 1)

#----------------------------------Neighborhoods-------------------------------------------------
random.solution.1 <- make.Random.Solution(categories)
neighborhood.1.solution<-neighborhood.1(random.solution.1,Data)
neighborhood.2.solution<-neighborhood.2(random.solution.1,Data)

random.solution.1.eval<-abs(sum(y - rowMeans(X[ , random.solution.1, drop = FALSE]))) 
neighborhood.1.eval<-abs(sum(y - rowMeans(X[ , neighborhood.1.solution, drop = FALSE])))
neighborhood.2.eval<-abs(sum(y - rowMeans(X[ , neighborhood.2.solution, drop = FALSE])))


#----------------------------------Greedy Search-------------------------------------------------
Search.Options <- list(nS = 500, neighbour = neighborhood.1, x0 = random.solution.1,
                       printBar = FALSE, printDetail = FALSE)

GS.solutions <- LSopt(Objective.function, algo = Search.Options, Data = Data)
GS.solutions.iterations <- min(which(GS.solutions$Fmat[ ,2] == GS.solutions$OFvalue))

#----------------------------------Local Search--------------------------------------------------
Search.Options$neighbour<-neighborhood.2
LS.solutions <- LSopt(Objective.function, algo = Search.Options, Data = Data)
LS.solutions.iterations <- min(which(LS.solutions$Fmat[ ,2] == LS.solutions$OFvalue))

#----------------------------------Threshold Accepting-------------------------------------------
Search.Options$nT <- 10
Search.Options$nS <- ceiling(Search.Options$nS/Search.Options$nT)
Search.Options$q <- 0.99

TA.solutions <- TAopt(Objective.function, algo = Search.Options, Data = Data)
TA.solutions.iterations <- min(which(TA.solutions$Fmat[ ,2] == TA.solutions$OFvalue))


#----------------------------------Stochastic Search Analysis------------------------------------

Stochastic.search.all <- sort(unique(c(which(TA.solutions$xbest),
                     which(LS.solutions$xbest),
                     which(GS.solutions$xbest),
                     X.desired)))

ta <- ls <- greedy <- true <- character(length(Stochastic.search.all))
true[  match(X.desired, Stochastic.search.all)] <- "Z"
greedy[match(which(GS.solutions$xbest),  Stochastic.search.all)] <- "Z"
ls[    match(which(LS.solutions$xbest), Stochastic.search.all)] <- "Z"
ta[    match(which(TA.solutions$xbest), Stochastic.search.all)] <- "Z"

Stochastic.Solutions.df<-data.frame(true = true, 
                                    greedy = greedy, 
                                    LS = ls , 
                                    TA = ta,
                                    row.names=Stochastic.search.all)

Solutions.neighborhood.df<-data.frame()
Solutions.neighborhood.df<-rbind(c(random.solution.1.eval),
                                 c(neighborhood.1.eval),
                                 c(neighborhood.2.eval)
)
colnames(Solutions.neighborhood.df)<-c("ABS(Error")
rownames(Solutions.neighborhood.df)<-c("Random","Neighborhood 1","Neighborhood 2")

#------------------------------------Tables---------------------------------------------------------
Table.1<-xtable(Stochastic.Solutions.df)
Table.2<-xtable(Solutions.neighborhood.df)
#------------------------------------Figures--------------------------------------------------------

par(ylog = TRUE, mar = c(5,5,1,6), las = 1)
plot(TA.solutions$Fmat[seq_len(TA.solutions.iterations) ,2],type = "l", log = "y",
     ylim = c(1e-3,
              max(pretty(c(GS.solutions$Fmat,LS.solutions$Fmat,TA.solutions$Fmat)))),
     xlab = "iterations", ylab = "Objective Function Value", col = grey(0.75))
lines(cummin(TA.solutions$Fmat[seq_len(TA.solutions.iterations), 2]), type = "l")
lines(GS.solutions$Fmat[ seq_len(GS.solutions.iterations),  2], type = "p", col = "blue")
lines(LS.solutions$Fmat[seq_len(LS.solutions.iterations), 2], type = "l", col = "red")
legend(x = "bottomleft",
       legend = c("TA best solution", "TA current solution",
                  "Greedy", "LS current/best solution"),
       lty = c(1,1,0,1),
       col = c("black",grey(0.5),"blue","red"),
       pch = c(NA,NA,21,NA))
axis(4, at = c(GS.solutions$OFvalue, LS.solutions$OFvalue, TA.solutions$OFvalue),
     labels = NULL, las = 1)
lines(x = c(GS.solutions.iterations, par()$usr[2]), y = rep(GS.solutions$OFvalue,2),
      col = "blue", lty = 2)
lines(x = c(TA.solutions.iterations, par()$usr[2]), y = rep(TA.solutions$OFvalue,2),
      col = "black", lty = 3)
lines(x = c(LS.solutions.iterations, par()$usr[2]), y = rep(LS.solutions$OFvalue,2),
      col = "red", lty = 4)

#------------------------------------Function Library-----------------------------------------------
Objective.function <- function(xn, Data)
  abs(sum(Data$y - rowMeans(Data$X[ ,xn, drop = FALSE])))


make.Random.Solution <- function(categories) {
  ii <- sample.int(categories, sample.int(categories, 1))
  random.solution <- logical(categories)
  random.solution[ii] <- TRUE
  random.solution
}

neighborhood.1 <- function(xc, Data) {
  of <- function(x)
    abs(sum(Data$y - rowMeans(Data$X[ ,x, drop = FALSE])))
  xbest <- xc
  Fxbest <- of(xbest)
  for (i in 1L:Data$nc) {
    xn <- xc; p <- i
    xn[p] <- !xn[p]
    if (sum(xn) >= 1) {
      Fxn <- of(xn)
      if (Fxn < Fxbest) {
        xbest <- xn
        Fxbest <- Fxn
      }
    }
  }
  xbest
}

neighborhood.2 <- function(xc, Data) {
  xn <- xc
  p <- sample.int(Data$nc, Data$n)
  xn[p] <- !xn[p]
  if (sum(xn) < 1)
    xn <- xc
  xn
}



