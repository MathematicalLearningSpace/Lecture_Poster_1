library(TSA)
library(PearsonDS)
library(tsDyn)
library(DescTools)

#------------------------------------------------------Data Generation----------------------------------------------------
pearson.N <- 100
pIIpars.1 <- list(a=2, location=1, scale=2) 
pIIpars.2 <- list(a=2, location=1, scale=4) 
error.pearson.2.1<-rpearsonII(pearson.N,params=pIIpars.1)
error.pearson.2.2<-rpearsonII(pearson.N,params=pIIpars.2)
x.ts<-ts(error.pearson.2.1,freq=1)
experimental.data.train <- window(x.ts, end = 80)
experimental.data.test<-window(x.ts,start=81)

#-----------------------------------------------------Stationarity Test--------------------------------------------------

test.1<-pp.test(experimental.data.train )
test.2<-kpss.test(experimental.data.train )
test.3<-adf.test(experimental.data.train ,alternative="stationary")
test.4<-adf.test(experimental.data.train ,alternative="explosive")


stationarity.test.df<-data.frame()
stationarity.test.df<-rbind(c(test.1$statistic,test.1$p.value,
                              test.2$statistic,test.2$p.value,
                              test.3$statistic,test.3$p.value,
                              test.4$statistic,test.4$p.value))
colnames(stationarity.test.df)<-c("Phillips-Perron","P Value","KPSS","P Value",
                                  "ADF1","P Value","ADF2","P Value")

#---------------------------------------------------Nonlinear Test--------------------------------------------------------
bds.test<-bds.test(experimental.data.train)
K.test<-Keenan.test(abs(experimental.data.train))
mcleod.test<-McLeod.Li.test(y=experimental.data.train)

nonlinear.test.df<-data.frame()
nonlinear.test.df<-rbind(c(K.test$test.stat,K.test$p.value))
colnames(nonlinear.test.df)<-c("K","P Value")

#-----------------------------------------------------Specification and Estimation of Model----------------------------------------------

model.tar.1<-tar(y=experimental.data.train ,
               p1=4,
               p2=4,
               d=1,
               a=.1,
               b=.9,print=TRUE)
model.tar.2<-tar(y=experimental.data.train ,
               p1=4,
               p2=4,
               d=2,
               a=.1,
               b=.9,print=TRUE)
model.tar.3<-tar(y=experimental.data.train ,
               p1=4,
               p2=4,
               d=3,
               a=.1,
               b=.9,print=TRUE)


model.estimation.df<-data.frame()
row.1<-c(model.tar.1$thd,model.tar.1$thdindex)
row.2<-c(model.tar.2$thd,model.tar.2$thdindex)
row.3<-c(model.tar.3$thd,model.tar.3$thdindex)
model.estimation.df<-rbind(row.1,row.2,row.3)
colnames(model.estimation.df)<-c("Threshold","Threshold Index")
#------------------------------------------------------Prediction--------------------------------------------------------

mod.predict.1<-predict(model.tar.1,n.ahead=20,n.sim=10)
mod.predict.2<-predict(model.tar.2,n.ahead=20,n.sim=10)
mod.predict.3<-predict(model.tar.3,n.ahead=20,n.sim=10)

experimental.data.predict.1=ts(c(experimental.data.train,mod.predict.1$fit),frequency=1,start=1)

tar.forecasts.accuracy.df<-data.frame()
tar.forecasts.accuracy.df<-rbind(c(Model="TAR",p1="4",p2="4",d="1",a=".1",b=".9",MAE=MAE(mod.predict.1$fit,experimental.data.test), 
                                   MSE=MSE(mod.predict.1$fit,experimental.data.test), 
                                   RMSE=RMSE(mod.predict.1$fit,experimental.data.test)),
                                 c(Model="TAR",p1="4",p2="4",d="1",a=".1",b=".9",MAE=MAE(mod.predict.2$fit,experimental.data.test), 
                                   MSE=MSE(mod.predict.2$fit,experimental.data.test), 
                                   RMSE=RMSE(mod.predict.2$fit,experimental.data.test)),
                                 c(Model="TAR",p1="4",p2="4",d="1",a=".1",b=".9",MAE=MAE(mod.predict.3$fit,experimental.data.test), 
                                   MSE=MSE(mod.predict.3$fit,experimental.data.test), 
                                   RMSE=RMSE(mod.predict.3$fit,experimental.data.test))
)

#-----------------------------------------------------Neural Network Prediction Comparision------------------------------

mod.nnet.1<- nnetTs(experimental.data.train, m=2, size=3)
mod.nnet.2<- nnetTs(experimental.data.train, m=6, size=3)
mod.nnet.3<- nnetTs(experimental.data.train, m=8, size=3)

mod.nn.predict.1<-predict(mod.nnet.1, n.ahead=20)
mod.nn.predict.2<-predict(mod.nnet.2, n.ahead=20)
mod.nn.predict.3<-predict(mod.nnet.3, n.ahead=20)

#---------------------------------------------------Tables---------------------------------------------------------------

Table.1<-xtable(tar.forecasts.accuracy.df)
Table.2<-xtable(stationarity.test.df)
Table.3<-xtable(nonlinear.test.df)
#-------------------------------------------------Figures----------------------------------------------------------------

Figure.1<-plot(experimental.data.predict.1,type='n',ylim=range(c(experimental.data.predict.1,mod.predict.1$pred.interval)),
    ylab='Actial Values/Predicted Values',
    xlab=expression(t),col="black")
lines(experimental.data.train,lty=2, col="red")
lines(window(experimental.data.predict.1, start=end(experimental.data.train)[1]+1),col="green", lty=3)
lines(ts(mod.predict.1$pred.interval[2,],
         start=end(experimental.data.train)[1]+1),lty=4, col="blue")
lines(ts(mod.predict.1$pred.interval[1,],
         start=end(experimental.data.train)[1]+1),lty=5,col="orange")
grid()
legend("bottomleft", col = c("black", "red", "green","blue","orange"), 
       lty = 1:3, cex=0.75,
       legend = c("Model 1", "Model 2", "Model 3"))


Figure.2<-plot(model.tar.1$residuals,
               type = "l", col = "black", lwd = 2,
               xlab = "Observations", ylab = "Time", 
               main = "TAR Model Residuals")
lines(model.tar.2$residuals,col="red", lty=2)
lines(model.tar.3$residuals,col="green", lty=3)
grid()
legend("bottomleft", col = c("black", "red", "green"), 
       lty = 1:3, cex=0.75,
       legend = c("Model 1", "Model 2", "Model 3"))

Figure.3<-plot(window(experimental.data.predict.1, start=end(experimental.data.train)[1]+1),lty=3,
               type = "l", col = "black", lwd = 2,
               xlab = "Observations", ylab = "Time", 
               main = "TAR Model Prediction With Neural Networks")
lines(mod.nn.predict.1,col="red", lty=2)
lines(mod.nn.predict.2,col="green", lty=3)
grid()
legend("bottomleft", col = c("black", "red", "green"), 
       lty = 1:3, cex=0.75,
       legend = c("TAR 1", "NN 1", "NN 2"))


