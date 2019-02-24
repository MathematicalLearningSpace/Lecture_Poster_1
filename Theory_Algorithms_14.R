#-----------------------------R Code To Modify in the Classroom Lecture with Students-----------------------
#----------------------------------R API -------------------------------------------------------------------
library(darch);library(e1071);library(klaR);library(ROCR);library(caret);library(recipes);library(deepnet);library(Matrix)
#---------------------Data------------------------------------------
u<-runif(100,0,1)
v<-runif(100,2,3)
w<-runif(100,0,4)
x<-runif(100,1,5)
y<-runif(100,1,3)
z<-runif(100,0,4)
l<-c(rep("A",35),rep("B",15),rep("C",25),rep("D",25))

group.X.df<-as.data.frame(cbind(x,y,z,l))

X <- matrix(c(x,y,z),nrow=100,ncol=3)
l.id <- as.numeric(group.X.df[, 4])
sample.idx <- sample(1:100, 80)

#--------------------Classification Models with Deep Learning--------------------------

model <- darch(l ~ ., data=group.X.df,retainData=F, layers=3,darch.numEpochs = 25)
model.A <- darch(l ~ ., data=group.X.df,retainData=F, layers=4,darch.numEpochs = 25)
model.B <- darch(l ~ ., data=group.X.df,retainData=F, layers=5,darch.numEpochs = 25)
model.C <- darch(l ~ ., data=group.X.df,retainData=F, layers=6,darch.numEpochs = 25)
model.D <- darch(l ~ ., data=group.X.df,retainData=F, layers=7,darch.numEpochs = 25)

#--------------------Other Models-----------------------------------------------------

model.2 <-naiveBayes(group.X.df[,-4], group.X.df[,4])
model.3<- rbm.train(X, 10, numepochs = 20, learningrate = 0.8, 
                    learningrate_scale = 1, momentum = 0.5,cd = 10)
model.4<-dbn.dnn.train(X, l.id, hidden = c(5,5))
model.5<-sae.dnn.train(X, l.id, hidden = c(5,5), 
                       activationfun = "sigm", learningrate = 0.8, 
                       momentum = 0.5, learningrate_scale = 1, 
                       output = "sigm", sae_output = "linear", 
                       numepochs = 3, batchsize = 100, hidden_dropout = 0, 
                       visible_dropout = 0)
model.6<-nn.train(X, l.id, hidden = c(5))

#------------------Model Predictions---------------------------------------------------

model.1.predictions <- predict(model, newdata = group.X.df, type = "class")
table(model.1.predictions)

model.2.predictions<-predict(model.2,newdata = group.X.df, type = "class")
table(model.2.predictions)

classificationStats <- darchTest(model)

model.4.predict<-nn.test(model.4, X,l.id)
model.5.predict<-nn.test(model.5, X,l.id)
model.6.predict<-nn.test(model.6, X,l.id)
#-----------------Multiple Deep Learning Models on a Grid--------------------------------------

tc <- trainControl(method = "boot", 
                   number = 1, 
                   allowParallel = F,
                   verboseIter = F)

parameters <- data.frame(parameter = c("layers", "bp.learnRate"),
                         class = c("character", "numeric"),
                         label = c("Network structure", "Learning rate"))

darch.grid <- function(x, y, len = NULL, search = "grid")
{
  df <- expand.grid(layers = c("c(0,20,0)",
                               "c(0,10,10,0)"),
                    bp.learnRate = c(1,2))
  df
}

darch.models <- train(l~ ., data = group.X.df, trControl = tc,
                    method = darchModelInfo(parameters, darch.grid), 
                    preProc = c("center", "scale"),
                    darch.numEpochs = 2)
darch.models$bestTune
#--------------------Tables----------------------------------------

Table.1<-xtable(darch.models$results)

#--------------------Figures---------------------------------------

Figure.1<-plot(model@stats$trainErrors$raw)
lines(model.A@stats$trainErrors$raw)
lines(model.B@stats$trainErrors$raw)
lines(model.C@stats$trainErrors$raw)
lines(model.D@stats$trainErrors$raw)

#--------------------References------------------------------------

Reference.1<-c("Hinton, G. E., S. Osindero, Y. W. Teh", 
"A fast learning algorithm for deep belief nets", 
"Neural Computation 18(7), S. 1527-1554, DOI: 10.1162/neco.2006.18.7.1527 2006.")


#-------------------Function Library-------------------------------
