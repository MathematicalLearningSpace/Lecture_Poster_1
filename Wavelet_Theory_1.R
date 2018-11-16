#---------------------------------------R API ---------------------------------------------
library(wmtsa);library(wavelets);library(wavethresh);library(waveslim);library(wavemulcor)
library(PearsonDS);library(xtable)
#---------------------------------------Wavelet Theory-------------------------------------

#---------------------------------------Data----------------------------------------------
pearson.N<-512
pIIIpars <- list(shape=1, location=1, scale=1)
error.pearson.3<-rpearsonIII(pearson.N,params=pIIIpars)
experimental.data.1<-window(error.pearson.3, end=256)
#---------------------------------------Transforms---------------------------------------

x.dwt <- wavDWT(as.ts(experimental.data.1), n.levels = 8)
x.modwt <- wavMODWT( as.ts(experimental.data.1), n.levels = 8)
x.cwt.1 <- wavCWT( as.ts(experimental.data.1), 
                   wavelet="gaussian2")

#---------------------------------------Tables------------------------------------------

wavelet.statistics.df<-data.frame()
wavelet.statistics.moments.df<-data.frame()
wavelet.statistics.moments.df<-rbind(c(empMoments(x.modwt$data$d1)),
                                     c(empMoments(x.modwt$data$d2)),
                                     c(empMoments(x.modwt$data$d3)),
                                     c(empMoments(x.modwt$data$d4)),
                                     c(empMoments(x.modwt$data$d5)),
                                     c(empMoments(x.modwt$data$d6))
)
wavelet.statistics.Hurst.df<-data.frame()
wavelet.statistics.Hurst.df<-cbind(c(HurstK(x.modwt$data$d1),
                                     HurstK(x.modwt$data$d2),
                                     HurstK(x.modwt$data$d3),
                                     HurstK(x.modwt$data$d4),
                                     HurstK(x.modwt$data$d5),
                                     HurstK(x.modwt$data$d6)))
colnames(wavelet.statistics.Hurst.df)<-c("Hurst")

wavelet.statistics.entropy.df<-data.frame()
wavelet.statistics.entropy.df<-cbind(c(Shannon.entropy(x.modwt$data$d1),
                                       Shannon.entropy(x.modwt$data$d2),
                                       Shannon.entropy(x.modwt$data$d3),
                                       Shannon.entropy(x.modwt$data$d4),
                                       Shannon.entropy(x.modwt$data$d5),
                                       Shannon.entropy(x.modwt$data$d6)))
colnames(wavelet.statistics.entropy.df)<-c("Shannon Entropy")

wavelet.statistics.df<-cbind(wavelet.statistics.moments.df,
                             wavelet.statistics.Hurst.df,
                             wavelet.statistics.entropy.df)
Table.1<-xtable(wavelet.statistics.df)


pearson.MLE.1 <- pearsonFitML(x.modwt$data$d1)
pearson.MLE.2 <- pearsonFitML(x.modwt$data$d2)
pearson.MLE.3 <- pearsonFitML(x.modwt$data$d3)
pearson.MLE.4 <- pearsonFitML(x.modwt$data$d4)
pearson.MLE.5 <- pearsonFitML(x.modwt$data$d5)
pearson.MLE.6 <- pearsonFitML(x.modwt$data$d6)
pearson.MLE.7 <- pearsonFitML(x.modwt$data$s3)

pearson.MLE.dist<-c(pearson.MLE.1$type,
                    pearson.MLE.2$type,
                    pearson.MLE.3$type,
                    pearson.MLE.4$type,
                    pearson.MLE.5$type,
                    pearson.MLE.6$type)

#--------------------------------------------------Distribution MOM Identification

pearson.MOM.1<-pearsonFitM(moments=empMoments(x.modwt$data$d1))
pearson.MOM.2<-pearsonFitM(moments=empMoments(x.modwt$data$d2))
pearson.MOM.3<-pearsonFitM(moments=empMoments(x.modwt$data$d3))
pearson.MOM.4 <- pearsonFitM(moments=empMoments(x.modwt$data$d4))
pearson.MOM.5 <- pearsonFitM(moments=empMoments(x.modwt$data$d5))
pearson.MOM.6 <- pearsonFitM(moments=empMoments(x.modwt$data$d6))
pearson.MOM.7<-pearsonFitM(moments=empMoments(x.modwt$data$s3))

pearson.MOM.dist<-c(pearson.MOM.1$type,
                    pearson.MOM.2$type,
                    pearson.MOM.3$type,
                    pearson.MOM.4$type,
                    pearson.MOM.5$type,
                    pearson.MOM.6$type)

distribution.table.df<-data.frame()
distribution.table.df<-cbind(pearson.MLE.dist,
                             pearson.MOM.dist)
Table.2<-xtable(distribution.table.df)

#-----------------MultiResolution Decomposition for Classroom Presentation--------------------
Figure.1<-plot(wavMRD(x.modwt, level = c(1:6)),col=c(1:7))
abline(v = 100, lty = 2)
abline(v = 250, lty = 3, col="Orange")
abline(v = 400, lty = 2)

par(mfrow = c(2, 2))
Figure.2<-plotdist(x.modwt$data$d1, histo = TRUE, demp = TRUE)
plotdist(x.modwt$data$d2, histo = TRUE, demp = TRUE)
plotdist(x.modwt$data$d3, histo = TRUE, demp = TRUE)
plotdist(x.modwt$data$d4, histo = TRUE, demp = TRUE)
plotdist(x.modwt$data$d5, histo = TRUE, demp = TRUE)
plotdist(x.modwt$data$d6, histo = TRUE, demp = TRUE)


