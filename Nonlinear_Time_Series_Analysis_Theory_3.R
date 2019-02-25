#-----------------------------R Code To Modify in the Classroom Lecture with Students-----------------------
#---------------------------------------------------R API --------------------------------------------------------
library(xtable);library(boot);library(sampling);library(PearsonDS);library(RandomFields);library(dtwclust);library(dtw);library(TSMining)
library(ggplot2);library(deSolve);library(ReacTran);library(rootSolve);library(fda)

#--------------------Data for the 5th Dimension----------------------------------------------------------------------------
pearson.N<-512
pVpars <- list(shape=1, location=0, scale=0.1) 
error.pearson.5<-rpearsonV(pearson.N,params=pVpars)
experimental.data.1.train<-window(error.pearson.5, end=475)
experimental.data.1.test<-window(error.pearson.5, start=476)

v<-seq(1,475,1)
w<-runif(n =475, min = -1, max = 1)
x <- experimental.data.1.train
y<-rnorm(n=475,0,1)
z<-runif(n = 475, min = -1, max = 1)

ts.1<-as.data.frame(cbind(v,w,x,y,z))

#------------------Differential Equation Model Template-----------------------------

equation.parameter.matrix<-matrix(
  c(0,0,0,0,0,0,
    1,0,0,0,0,0,
    1,1,0,0,0,0,
    1,1,1,0,0,0,
    1,1,1,1,0,0,
    1,1,1,1,1,0),6,6)


parameters.test <- c(a = -8/3, b = -10, c = 28,d=-1,e=1)
variables.initial.test <- c(x1 = 1, x2 = 1, x3 = 1,x4=1,x5=1)
times <- seq(0, 100, by = 0.1)


system.equation.model.test<-function(times, variables.intitial.test, parameters.test)
{
  with(as.list(c(parameters.test, variables.intitial.test)), {
    #---------------Equation Matrix for the Classroom--------------------	  
    dx.dt.1 <- a*x1 + x2*x3  
    dx.dt.2 <- b * (x2-x3)
    dx.dt.3 <- -x1*x2 + c*x2 - x3 
    dx.dt.4<- -d*x4
    dx.dt.5<- e*x5+x1
    #---------------Dependent Variable for Ratio of Changes Vector--------	  
    res <- c(dx.dt.1,dx.dt.2,dx.dt.3,dx.dt.4,dx.dt.5)
    list(res)
  })
}

#---------------------------------Template Model Solution--------------------------------

system.equation.model.test.solution.1 <- ode(y = variables.initial.test , times = times, 
                                             func = system.equation.model.test, 
                                             parms = parameters.test) 
system.equation.model.test.solution.1.df<-as.data.frame(system.equation.model.test.solution.1)

summary(system.equation.model.test.solution.1)


#-------------SAX Representations for Motifs-------------------------------

sax.w<-Func.SAX(x = ts.1$w, w = 5, a = 5, eps = .01, norm = TRUE)
sax.x<-Func.SAX(x = ts.1$x, w = 5, a = 5, eps = .01, norm = TRUE)
sax.x3<-Func.SAX(x = system.equation.model.test.solution.1.df$x3, w = 5, a = 5, eps = .01, norm = TRUE)
sax.z<-Func.SAX(x = ts.1$z, w = 5, a = 5, eps = .01, norm = TRUE)

res.1 <- Func.motif(ts.1$w,global.norm = TRUE, local.norm = FALSE,
                    window.size = 2, overlap = 0, 
                    w = 5, a = 3, mask.size = 3, eps = .01)
res.2 <- Func.motif(ts = ts.1$x, 
                    global.norm = TRUE, 
                    local.norm = FALSE,
                    window.size = 3, 
                    overlap = 0, 
                    w = 5, 
                    a = 3, 
                    mask.size = 3, 
                    eps = .01)
data.vis <- Func.visual.SingleMotif(single.ts=ts.1$w, window.size=5, motif.indices=res.1$Indices)
res.multi.local <- Func.motif.multivariate(motif.list = list(res.1$Indices, 
                                                             res.2$Indices),
                                           window.sizes = c(5,3), 
                                           alpha = .8)
res.multi.global <- Func.motif.multivariate(motif.list = list(res.1$Indices, 
                                                              res.2$Indices),
                                            window.sizes = c(5,3), 
                                            alpha = .8)
data.multi.local <- Func.visual.MultiMotif(data = ts.1, 
                                           multi.motifs = res.multi.local, 
                                           index = 1)

#--------------Hierarchical clustering of time series for Classification--------------------------------
ctrl <- new("dtwclustControl", trace = TRUE)
dtw.cluster.ts.1 <- dtwclust(ts.1, k = 2:4,
                  distance = "L2", centroid = "pam",
                  seed = 3247, control = ctrl)
dtw.cluster.hc.ts.1<- dtwclust(ts.1, type = "hierarchical",
                   k = 4, method = "all",
                   distance = "sbd", control = ctrl)

acf_fun <- function(dat) {
  lapply(dat, function(x) as.numeric(acf(x, lag.max = 50, plot = FALSE)$acf))
}
dtw.fuzzy.ts.1 <- dtwclust(ts.1, type = "fuzzy", k = 4,
               preproc = acf_fun, distance = "L2",
               seed = 123, control = ctrl)
z.1<-rnorm(n=100,0,1)

dtw.cluster.predictions.ts.1.df<-data.frame()
dtw.cluster.predictions.ts.1.df<-rbind(c(predict(dtw.fuzzy.ts.1, newdata = z.1)))
colnames(dtw.cluster.predictions.ts.1.df)<-c("Cluster 1","Cluster 2","Cluster 3","Cluster 4")
rownames(dtw.cluster.predictions.ts.1.df)<-c("Fuzzy")
#-------------------------------------------------------------------------Tables-----------------------------------------------------------------------------

Table.1<-xtable(dtw.cluster.predictions.ts.1.df)

#-------------------------------------------------------------------------Figures----------------------------------------------------------------------------

Figure.1<-plot(ts.1$w, type='l', lty=1,col=1,xlab="Observations",ylab="Interval of Values")
lines(ts.1$x,lty=2,col=2)
lines(ts.1$y,lty=3,col=3)
lines(ts.1$z,lty=4,col=4)
legend("topright", col = 1:4, lty = 1:4, legend = c("Unif", "Pearson 5",'Normal','Unif'))

Figure.2<-ggplot(data = data.vis$data.1) +
  geom_line(aes(x = 1:dim(data.vis$data.1)[1], y = X)) +
  geom_point(aes(x = 1:dim(data.vis$data.1)[1], y = X, color=Y))

Figure.3<-ggplot(data = data.vis$data.2[[3]]) + geom_line(aes(x = Time, y = Value, linetype=Instance))
Figure.4<-ggplot(data = data.multi.local) +
  geom_line(aes(x = T, y = X)) +
  geom_point(aes(x = T, y = X, col=Lab, shape=Lab)) + facet_grid(Facet~.)
Figure.5<-plot(dtw.cluster.ts.1[[3L]])
Figure.6<-plot(dtw.cluster.hc.ts.1[[ which.max(sapply(dtw.cluster.hc.ts.1, randIndex, y = colnames(ts.1))) ]])
Figure.7<-plot(dtw.fuzzy.ts.1, data = ts.1, type = "series")
Figure.8<-plot(dtw.fuzzy.ts.1, data = ts.1, type = "series")
         
         
#-------------Function Template Library for Classroom Presentation and Modification---------------------
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
