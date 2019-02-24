#-----------------------------R Code To Modify in the Classroom Lecture with Students-----------------------
#-----------------------------R API ------------------------------------------------------------------------
library(Matrix);library(igraph);library(rbenchmark);library(xtable);library(stringi)
library(PearsonDS);library(h2o);library(darch);library(deepnet);library(caret)
#---------------------------------------------Data----------------------------------------------------------------------
filesToProcess <- dir(pattern = "file.*\\.txt$")
listOfFiles <- lapply(filesToProcess, function(x) read.table(x, header = TRUE))
#--------------------------------------------Choose the Pearson Data-----------------------------------------
pearson.N<-512
pIIIpars.1 <- list(shape=1, location=1, scale=1)
pIIIpars.2 <- list(shape=2, location=2, scale=2)
pIIIpars.3 <- list(shape=3, location=2, scale=3)
pIIIpars.4 <- list(shape=4, location=2, scale=4)

error.pearson.3.1<-rpearsonIII(pearson.N,params=pIIIpars.1)
error.pearson.3.2<-rpearsonIII(pearson.N,params=pIIIpars.2)
error.pearson.3.3<-rpearsonIII(pearson.N,params=pIIIpars.3)
error.pearson.3.4<-rpearsonIII(pearson.N,params=pIIIpars.4)

#--------------------------------------Define the Window-----------------------------------------------------
experimental.data.1<-window(error.pearson.3.1, end=256)
experimental.data.2<-window(error.pearson.3.2, end=256)
experimental.data.3<-window(error.pearson.3.3, end=256)
experimental.data.4<-window(error.pearson.3.4, end=256)

#--------------------------------------Train and Test Data--------------------------------------------------
experimental.data.1.train<-window(error.pearson.3.1, end=400)
experimental.data.2.train<-window(error.pearson.3.2, end=400)

experimental.data.1.test<-window(error.pearson.3.1, start=401)
experimental.data.2.test<-window(error.pearson.3.2, start=401)

#------------------------------------Conditional Nonlinear relationship in the mean-------------------------
y.actual<-experimental.data.1^2
y.actual.train<-experimental.data.1.train^2
y.actual.test<-experimental.data.1.test^2
y.estimate<-experimental.data.1^2+experimental.data.2
y.estimate.train<-experimental.data.1.train^2+experimental.data.2.train
y.estimate.test<-experimental.data.1.train^2+experimental.data.2.train

#--------------------------------------Translate Data Structure------------------------------------------
Var.1<-c(experimental.data.1,experimental.data.2)
Var.2<-c(experimental.data.3,experimental.data.4)
X<-matrix(c(Var.1,Var.2),nrow=length(experimental.data.1),ncol=2)
k<-0.80
Y<- c(rep(1, 200), rep(0,56))

matrix.X<-as.matrix(experimental.data.1)
matrix.Y<-as.matrix(y.estimate)

matrix.X.train<-as.matrix(experimental.data.1.train)
matrix.Y.train<-as.matrix(y.estimate.train)

matrix.X.test<-as.matrix(experimental.data.1.test)
matrix.Y.test<-as.matrix(y.estimate.test)

#--------------------------------------Conditional Moment Models---------------------------------------

model.NN.1<-nn.train(matrix.X,matrix.Y, initW=NULL, initB=NULL, hidden=c(1), 
         activationfun="sigm", 
         learningrate=0.8, 
         momentum=0.5, 
         learningrate_scale=1, 
         output="sigm", 
         numepochs=3, batchsize=100, hidden_dropout=0, visible_dropout=0)

model.RBM.1<-rbm.train(X, 3, numepochs = 20, cd = 10)

model.dbn.1<-dbn.dnn.train(matrix.X, matrix.Y, hidden=c(1), 
              activationfun="sigm", 
              learningrate=0.8, momentum=0.5, 
              learningrate_scale=1, output="sigm", 
              numepochs=3, batchsize=100, 
              hidden_dropout=0, visible_dropout=0, cd=1)

model.sae.1<-sae.dnn.train(matrix.X, matrix.Y, hidden=c(1), 
              activationfun="sigm", 
              learningrate=0.8, momentum=0.5, 
              learningrate_scale=1, output="sigm", 
              sae_output="linear", numepochs=3, batchsize=100, 
              hidden_dropout=0, visible_dropout=0)

model.darch.1<- darch(matrix.X, Y,
                rbm.numEpochs = 0,
                rbm.batchSize = 100,
                rbm.trainOutputLayer = F,
                layers = c(784,100,10),
                darch.batchSize = 100,
                darch.learnRate = 2,
                darch.retainData = T,
                darch.numEpochs = 20 )
#----------------------------------Grid Search-----------------------------------------------------
tc <- trainControl(method = "boot", number = 2, allowParallel = F,verboseIter = T)
parameters <- data.frame(parameter = c("layers", "bp.learnRate", "darch.unitFunction"),
                         class = c("character", "numeric", "character"),
                         label = c("Network structure", "Learning rate", "unitFunction"))
grid <- function(x, y, len = NULL, search = "grid")
{
  df <- expand.grid(layers = c("c(0,3,0)","c(0,3,3,0)","c(0,3,3,3,0)"),bp.learnRate = c(1,2,5,10))
  df[["darch.unitFunction"]] <- rep(c("c(tanhUnit, softmaxUnit)",
                                      "c(tanhUnit, tanhUnit, softmaxUnit)",
                                      "c(tanhUnit, tanhUnit, tanhUnit, softmaxUnit)"), 4)
  df
}
model.darch.1.caret<-train(Y ~ ., data = matrix.X, tuneLength = 12, trControl = tc,
                    method = darchModelInfo(parameters, grid), preProc = c("center", "scale"),
                    darch.numEpochs = 15, darch.batchSize = 6)


#-------------------------------Model Summaries, Classification and Prediction-------------------------------------------------------------
y.estimate.NN.1.predict <- nn.predict(model.NN.1, experimental.data.3)
y.estimate.sae.1.predict <- nn.predict(model.sae.1, experimental.data.3)

Y.predictions <- predict(model.darch.1, newdata = matrix.X, type = "class")

err <- nn.test(model.NN.1, experimental.data.4, y.estimate,t=0.5)

states.H.1 <- c(0.1, 0.1, 0.1)
states.H.V.1 <- rbm.down(model.RBM.1, states.H.1)
states.H.V.2 <- c(0.2, 0.8)
states.H.2 <- rbm.up(model.RBM.1, states.H.V.2)

#-----------------------------Classification-----------------------------------------------------
cat(paste("Incorrect classifications:", sum(Y.predictions != Y)))

Model.Statistics.Classification.df<-data.frame()
Model.Statistics.Classification <- darchTest(model.darch.1, newdata = matrix.X,targets = Y)

Model.Statistics.Classification.df<-rbind(c(Model.Statistics.Classification$error,
                                            Model.Statistics.Classification$percentIncorrect,
                                            Model.Statistics.Classification$numIncorrect))
colnames(Model.Statistics.Classification.df)<-c("Error","Incorrect.PCT","Incorrect.Num")

Model.Statistics.Prediction.df<-data.frame()
Model.Statistics.Prediction.df<-rbind(c(Metric.RMS(y.estimate.NN.1.predict,y.estimate),
                                      Metric.RMS(y.estimate.sae.1.predict,y.estimate),
                                      Metric.RMS(Y.predictions,y.estimate)))
colnames(Model.Statistics.Prediction.df)<-c("NeuralNetwork.1","SAE.1","DARCH.1")
rownames(Model.Statistics.Prediction.df)<-c("RMS")

#-----------------------------Tables--------------------------------------------------------------------
Table.1<-xtable(Model.Statistics.Prediction.df)
Table.2<-xtable(Model.Statistics.Classification.df)
#-------------Figures for Classroom Presentation--------------------------------------------------------

Figure.1<-plot(model.NN.1$e,lty=1,col=1)
title(main='')
lines(model.sae.1$e,lty=2,col=2)
cols <- c("black","red")
legend("topright", legend=c("NN.e","SAE.e")
       , bty = "n",lwd=2, 
       cex=0.75, col=cols, text.col=cols, lty=1:2)

Figure.2<-plot(model.darch.1)
Figure.3<-plot(model.darch.1, "class")
Figure.4<-plot(model.darch.1, "time")
Figure.5<-plot(model.darch.1, "momentum")
Figure.6<-plot(model.darch.1, "net")

#--------Function Library for Student Modification in the Classroom---------------------------------------------------------
Metric.RMS <- function(x.actual, x.predict) {
  res <- (mean((x.actual-x.predict)^2))^(1/2)
  return(res)
}

