library(stsm)
library(xtable)
library(Deriv)
library(numDeriv)
library(pracma)
library(Matrix)
library(forecast)
library(PearsonDS)

#---------------------------------------------------------Generate Data-------------------------------------------------
pearson.N <- 100
pIIpars.1 <- list(a=2, location=1, scale=2) 
pIIpars.2 <- list(a=2, location=1, scale=4) 
error.pearson.2.1<-rpearsonII(pearson.N,params=pIIpars.1)
error.pearson.2.2<-rpearsonII(pearson.N,params=pIIpars.2)
x.ts<-ts(error.pearson.2.1,freq=4)

#-------------------------------------------------------# Define local level plus seasonal model--------------------
pars <- c("var1" = 2^2, "var2" = 2^3, "var3" = 2^4)
m.1 <- stsm.model(model = "llm+seas", y = x.ts, pars = pars)
mloglik.fd(model = m.1)
# 'logLik' returns the value of the log-likelihood function
logLik(object = m1, domain = "frequency")

model.likelihood.1 <- mloglik.fd.deriv(model = m.1, gradient = TRUE, hessian = TRUE)

model.estimation.df<-data.frame()
model.estimation.df<-cbind(c(model.likelihood.1$pars),
                        c(model.likelihood.1$std.errors),
                         model.likelihood.1$loglik,
                        model.likelihood.1$convergence,
                          model.likelihood.1$iter)
model.likelihood.1$gradient
model.likelihood.1$hessian
model.Hessian<-as.matrix(model.likelihood.1$hessian)

#--------------------------------------------------------Barrier work-------------------------------------------

bar.1 <- barrier.eval(m.1, type = "1", mu = 1)
bar.2 <- barrier.eval(m.1, type = "2", mu = 2)
bar.3 <- barrier.eval(m.1, type = "2", mu = 3)

barrier.design.df<-data.frame()
barrier.design.df<-cbind(bar.1$barrier,
      bar.2$barrier,
      bar.3$barrier)
colnames(barrier.design.df)<-c("A","B","C")

#--------------------------------------------------------Simulated Hessian-------------------------------------
ss <- char2numeric(m.1)
sim.1 <- datagen.stsm(n = 500, model = list(Z = ss$Z, T = ss$T, H = ss$H, Q = ss$Q), 
                      n0 = 20, freq = 4, old.version = TRUE)$data

#--------------------------------------------------------Spectral generating function of the local level plus seasonal model
res <- stsm.sgf(m.1, gradient=TRUE)
res$sgf

#--------------------------------------------------------Adding a Barrier term with different specifications
#
model.likelihood.1<-maxlik.fd.scoring(m = m.1, information = "expected",
                                      control = list(maxit = 100, tol = 0.001))
model.likelihood.2<-maxlik.fd.scoring(m= m.1, information = "expected", barrier = list(type = "1", mu = 2),
                                      control = list(maxit = 100, tol = 0.001))
model.likelihood.3<-maxlik.fd.scoring(m= m.1, information = "expected", barrier = list(type = "2", mu = 10),
                                      control = list(maxit = 100, tol = 0.001))

row.1<-c(model.likelihood.2$pars,model.likelihood.1$std.errors,
                           model.likelihood.1$loglik,
                           model.likelihood.1$convergence,
                           model.likelihood.1$iter)
row.2<-c(model.likelihood.2$pars,model.likelihood.2$std.errors,
       model.likelihood.2$loglik,
       model.likelihood.2$convergence,
       model.likelihood.2$iter)
row.3<-c(model.likelihood.3$pars,model.likelihood.2$std.errors,
        model.likelihood.3$loglik,
        model.likelihood.3$convergence,
        model.likelihood.3$iter)
model.estimation.df<-rbind(row.1,row.2,row.3)
#--------------------------------------------------------Estimated components with 95% confidence bands
comps <- tsSmooth(model.likelihood.1)

#----------------------------------------------------------Plot predictions eight periods ahead
pred.1<- predict(model.likelihood.1, n.ahead = 20, se.fit = TRUE)
pred.2<- predict(model.likelihood.2, n.ahead = 20, se.fit = TRUE)
pred.3<- predict(model.likelihood.3, n.ahead = 20, se.fit = FALSE)

#------------------------------------------------------------Tables--------------------------------------------------------

Table.1<-xtable(model.Hessian)
Table.2<-xtable(model.estimation.df)
Table.3<-xtable(barrier.design.df)
#-----------------------------------------------------------Figures--------------------------------------------------------

Figure.1<-plot(x.ts)

Figure.2<-plot(res$sgf)
title(main = "Spectral Generating Function")

Figure.3<-tsdiag(model.likelihood.1)
Figure.4<-tsdiag(model.likelihood.2)
Figure.5<-tsdiag(model.likelihood.3)

Figure.6<-plot(comps)
title(main = "smoothed trend and seasonal components")

Figure.7<-plot(pred.1$pred,type = "l", col = "black", lwd = 2,
     xlab = "", ylab = "", main = "Structural Model")
lines(pred.2$pred,col="red", lty=2)
lines(pred.3$pred,col="green", lty=3)
grid()
legend("bottomleft", col = c("black", "red", "green"), 
       lty = 1:3, cex=0.75,
       legend = c("No Barrier", "Barrier=1", "Barrier=2"))


