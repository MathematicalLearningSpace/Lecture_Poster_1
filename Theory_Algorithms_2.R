#-----------------------------------------------R API ------------------------------
library(NMOF); library(xtable); library(recommenderlab)
#--------------------------------------------Data-------------------------------------------------------------------------------
#--------------------------------------------Correlation and Rating Matrix-------------------------------------------------------
n <- 100
Obs<-100
Features<-5
matrix.rating<-matrix(rbinom(500, size = 5, prob = 0.4),Obs,Features,
                      dimnames=list(user=paste("u", 1:100, sep=''), 
                                    item=paste("i", 1:5, sep='')))
rho <- 0.75
C <- matrix(rho, 2L, 2L); 
diag(C) <- 1
x <- matrix(rnorm(n * 2L), n, 2L) %*% chol(C)
data <- list(x = x, n = n, nmin = 40L)
#--------------------------------------------Objective Function--------------------------------------------------------------------
Objective.Function <- function(xc, data)
  -abs(cor(data$x[xc, ])[1L, 2L] - cor(data$x[!xc, ])[1L, 2L])
x0 <- runif(n) > 0.5
Objective.Function(x0, data)
Objective.Function(neighbour(x0, data), data)
#--------------------------------------------Local Search-------------------------------------------------------------------------
x0 == neighbour(x0, data)
algo <- list(nS = 3000L,
             neighbour = neighbour,
             x0 = x0,
             printBar = FALSE)
LS.sol <- LSopt(Objective.Function, algo = algo, data=data)
LS.sol.OFvalue<-LS.sol$OFvalue

#--------------------------------------------Threshold Accepting-----------------------------------------------------------------
algo$nT <- 10L
algo$nS <- ceiling(algo$nS/algo$nT)
TA.sol <- TAopt(Objective.Function, algo = algo, data = data)
TA.sol$OFvalue

c1 <- cor(data$x[ TA.sol$xbest, ])[1L, 2L]
c2 <- cor(data$x[!TA.sol$xbest, ])[1L, 2L]

#---------------------------------------------Top N Recommender Algorithm--------------------------------------------------------
r <- as(matrix.rating[1:100,], "realRatingMatrix")
getRatingMatrix(r)
r_m <- normalize(r) 
recommend <- Recommender(r, method = "POPULAR") 
names(getModel(recommend))

recommend.topN<-getModel(recommend)$topN
recommend.bestN.10<- bestN(recommend.predict, n = 10)
recommend.predict<-predict(recommend, r,n=5)
#-----------------------------------Evalutation of Filtering Methods------------------------------------------------------------

scheme <- evaluationScheme(r, method="split", train=0.9, given=3, goodRating=4)
r1 <- Recommender(getData(scheme, "train"), "UBCF") 
r2 <- Recommender(getData(scheme, "train"), "IBCF") 
p1 <- predict(r1, getData(scheme, "known"), type="ratings") 
p2 <- predict(r2, getData(scheme, "known"), type="ratings") 
error <- rbind( UBCF = calcPredictionAccuracy(p1, getData(scheme, "unknown")), 
                IBCF = calcPredictionAccuracy(p2, getData(scheme, "unknown")) ) 

Evalutation.recommender.df<-data.frame(error)
#---------------------------------Top 20-----------------------------------------------------------------------------------------
results <- evaluate(scheme, method="POPULAR", type = "topNList", 
                    n=c(1,3,5,10,15,20))
matrix.cf<-getConfusionMatrix(results)[[1]]
matrix.cf.avg<-avg(results)
#-------------------------------------------Tables-------------------------------------------------------------------------------

Table.1<-xtable(Evalutation.recommender.df)
Table.2<-xtable(matrix.cf)
#-------------------------------------------Figures for Presentation in the Classroom------------------------------------------------------------------------------
Figure.1<-par(mfrow = c(1,3), bty = "n")
plot(data$x,
     xlim = c(-3,3), ylim = c(-3,3),
     main = "all data", col = "darkgreen")
lines(data$x[ TA.sol$xbest, ], type = "p", col = "blue")

Figure.2<-plot(data$x[ TA.sol$xbest, ], col = "blue",
     xlim = c(-3,3), ylim = c(-3,3),
     main = paste("subset 1, corr.", format(c1, digits = 3)))

Figure.3<-plot(data$x[!TA.sol$xbest, ], col = "darkgreen",
     xlim = c(-3,3), ylim = c(-3,3),
     main = paste("subset 2, corr.", format(c2, digits = 3)))

#-------------------------------------------compare LS/TA----------------------------------------------------------------------
par(mfrow = c(2,2), bty = "n")
Figure.4<-plot(x0,type="l")
Figure.5<-plot(LS.sol$Fmat[1:500 ,2L],type="l", ylim=c(-1.5,0.5),
     ylab = "OF", xlab = "TA iterations")
lines(TA.sol$Fmat[1:500 ,2L],type = "l", col = "blue")
legend(x = "topright",legend = c("LS", "TA"),
       lty = 1, lwd = 2,col = c("black", "blue"))

Figure.6<- plot(results, annotate=TRUE)
Figure.7<- plot(results, "prec/rec", annotate=TRUE)
Figure.8<-hist(matrix.rating)

#---------------------------------------------------------------
par(mfrow = c(1,2))
image(r, main = "Raw Ratings") 
image(r_m, main = "Normalized Ratings")

#------------------------------------------Function Library--------------------------------------------------------------------
