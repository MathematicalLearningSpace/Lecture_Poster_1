library(tseries)
library(costat)
library(locits)
library(wbsts)
library(forecast)
library(tsoutliers)
library(jmotif)
library(TSclust)
library(TSMining)
library(ggplot2)
library(tsDyn)
library(tseriesChaos)
library(yuima)
library(DescTools)

library(xtable)
library(PearsonDS)
library(fitdistrplus)
library(psych)

#------------------------------------------Nonlinear Time Series Analysis--------------------------------------

#------------------------------------------Data-----------------------------------------------------------------
pearson.N<-512
pVpars <- list(shape=1, location=1, scale=0.1) 
error.pearson.5<-rpearsonV(pearson.N,params=pVpars)
experimental.data.1.train<-window(error.pearson.5, end=475)
experimental.data.1.test<-window(error.pearson.5, start=476)

#-----------------------------------------Transformations-----------------------------------------------------

experimental.data.1.train<-diff(experimental.data.1.train,1)
experimental.data.1.test<-diff(experimental.data.1.test,1)

#-----------------------------------------Nonlinear Tests-----------------------------------------------------

setar.test.1<-setarTest(experimental.data.1.train, m=3, thDelay=0:1, nboot=5,trim=0.1, test="1vs")
setar.test.2<-setarTest(experimental.data.1.train, m=7, thDelay=0:1, nboot=5,trim=0.1, test="1vs")
setar.test.3<-setarTest(experimental.data.1.train, m=11, thDelay=0:1, nboot=5,trim=0.1, test="1vs")

setar.test.df<-data.frame()
setar.test.df<-rbind(c(setar.test.1$SSRs,setar.test.1$Ftests,setar.test.1$firstBests,setar.test.1$secBests),
                     c(setar.test.2$SSRs,setar.test.2$Ftests,setar.test.2$firstBests,setar.test.2$secBests),
                     c(setar.test.3$SSRs,setar.test.3$Ftests,setar.test.3$firstBests,setar.test.3$secBests))

colnames(setar.test.df)<-c("AR","TAR(1)","TAR(2)","1vs2","1vs3","2vs3","tau1","tau1A","tau2")
Table.1<-xtable(setar.test.df)

delta.test<-delta(experimental.data.1.train[1:450], m=3, eps=sd(experimental.data.1.train[1:450]))
delta.LinTest<-delta.lin(experimental.data.1.train[1:450], m=3)
terasvirta.test<-terasvirta.test(as.ts(experimental.data.1.train[1:450]),lag = 1,type = c("Chisq","F"))
white.test<-white.test(as.ts(experimental.data.1.train[1:450]),lag = 1, qstar = 2, q = 10, range = 4, 
                       type = c("Chisq","F"))

bbc.test<-BBCTest(experimental.data.1.train[1:450], m=3, test="Wald", grid="minPerc")

nonlinear.test.df<-data.frame()
nonlinear.test.df<-rbind(c(delta.test,delta.LinTest,bbc.test$statistic,
                           white.test$statistic, white.test$p.value,
                           terasvirta.test$statistic,terasvirta.test$p.value))
colnames(nonlinear.test.df)<-c("Delta","Delta Linear","BBC","white","p Value","Terasvirta","p value")
Table.2<-xtable(nonlinear.test.df)

#--------------------------------------------------------Specification of Nonlinear Time Series Models----------------------------------

mod.setar.1 <- setar(experimental.data.1.train, m=2, thDelay=1)
mod.setar.2 <- setar(experimental.data.1.train, m=3, thDelay=1)
mod.setar.3 <- setar(experimental.data.1.train, m=4, thDelay=1)

pred_setar_naive.1<- predict(mod.setar.1, n.ahead=35)
pred_setar_naive.2<- predict(mod.setar.1, n.ahead=35)
pred_setar_naive.3<- predict(mod.setar.1, n.ahead=35)

pred_setar_boot.1<- predict(mod.setar.1, n.ahead=35, type="bootstrap", n.boot=100)
pred_setar_boot.2<- predict(mod.setar.2, n.ahead=35, type="bootstrap", n.boot=100)
pred_setar_boot.3<- predict(mod.setar.3, n.ahead=35, type="bootstrap", n.boot=100)

pred_setar_Bboot.1 <- predict(mod.setar.1, n.ahead=35, type="block-bootstrap", n.boot=100)
pred_setar_Bboot.2 <- predict(mod.setar.2, n.ahead=35, type="block-bootstrap", n.boot=100)
pred_setar_Bboot.3 <- predict(mod.setar.3, n.ahead=35, type="block-bootstrap", n.boot=100)

pred_setar_MC.1 <- predict(mod.setar.1, n.ahead=35, type="MC", n.boot=100)
pred_setar_MC.2 <- predict(mod.setar.2, n.ahead=35, type="MC", n.boot=100)
pred_setar_MC.3 <- predict(mod.setar.3, n.ahead=35, type="MC", n.boot=100)

regime(mod.setar.1)
regime(mod.setar.2)
regime(mod.setar.3)

regime(mod.setar.1, time=FALSE, initVal=FALSE)

setar.forecasts.accuracy.df<-data.frame()
setar.forecasts.accuracy.df<-rbind(c(Model="1",Forecast="Naive",MAE=MAE(pred_setar_naive.1,experimental.data.1.test), 
                                     MSE=MSE(pred_setar_naive.1,experimental.data.1.test), 
                                     RMSE=RMSE(pred_setar_naive.1,experimental.data.1.test)),
                                   c(Model="1",Forecast="Bootstrap",MAE=MAE(pred_setar_boot.1$pred,experimental.data.1.test), 
                                     MSE=MSE(pred_setar_boot.1$pred,experimental.data.1.test), 
                                     RMSE=RMSE(pred_setar_boot.1$pred,experimental.data.1.test)),
                                   c(Model="1",Forecast="Block Bootstrap",MAE=MAE(pred_setar_Bboot.1$pred,experimental.data.1.test), 
                                     MSE=MSE(pred_setar_Bboot.1$pred,experimental.data.1.test), 
                                     RMSE=RMSE(pred_setar_Bboot.1$pred,experimental.data.1.test)),
                                   c(Model="1",Forecast="MC",MAE=MAE(pred_setar_MC.1$pred,experimental.data.1.test), 
                                     MSE=MSE(pred_setar_MC.1$pred,experimental.data.1.test), 
                                     RMSE=RMSE(pred_setar_MC.1$pred,experimental.data.1.test))
                                    )
Table.3<-xtable(setar.forecasts.accuracy.df)


#-----------------------------------------Figure--------------------------------------------------------------
op <- par(mfrow = c(2, 2))
densityBy(as.data.frame(experimental.data.1.train[1:450]),
          xlab="x", 
          ylab="Simulated",col="black")
densityBy(as.data.frame(experimental.data.1.train),
          xlab="x", 
          ylab="Simulated",col="blue")
densityBy(as.data.frame(experimental.data.1.train),
          xlab="x", 
          ylab="Simulated",col="green")
densityBy(as.data.frame(experimental.data.1.train),
          xlab="x", 
          ylab="Simulated",col="red")
par(op) 

Figure.2<-plot(regime(mod.setar.1))
Figure.3<-hist(regime(mod.setar.1))

## Plot to compare results:
plot.new()
pred_range <- range(pred_setar_naive.1, 
                    pred_setar_boot.1$pred, 
                    pred_setar_MC.1$pred, na.rm=TRUE)
Figure.4<-plot(experimental.data.1.test,type="l", 
     ylab="Pearson V",
     main="Pearson V Test Data Set")

Figure.5<-plot(pred_setar_naive.1, lty=2, col=2, main="Comparison of forecasts methods from same SETAR")
lines(pred_setar_boot.1$pred, lty=3, col=3)
lines(pred_setar_Bboot.1$pred, lty=4, col=4)
lines(pred_setar_MC.1$pred, lty=5, col=5)
legLabels <- c("Naive F", 
               "Bootstrap F",
               "Block-Bootstrap F", 
               "MC F")
legend("topright", leg=legLabels, lty=1:5, col=1:5, cex=0.5)