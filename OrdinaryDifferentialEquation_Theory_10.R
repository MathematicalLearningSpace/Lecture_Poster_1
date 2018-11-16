#----------------------------------R API ------------------------------------------------
library(deBInfer); library(xtable); library(actuar)
#--------------------------------------------Data----------------------------------------
r<-rpois(25,1);s<-rpois(25,1);u<-rpois(25,1);v<-rpoisinvgauss(25,2);w<-rpoisinvgauss(25,4)
x<-rpois(25,5);y<-rpois(25,1);z<-rpois(25,1)
cnt<-c(r,s,u,v,w,x,y,z)
t<-c(rep(1,25),rep(2,25),rep(3,25),rep(4,25),rep(5,25),rep(6,25),rep(8,25),rep(9,25))

group.X.df<-as.data.frame(cbind(t,cnt))
#--------------------------------------------Parameters----------------------------------
x1 <- debinfer_par(name = "x1", var.type = "de", fixed = FALSE, value = 2, 
                   prior = "gamma", 
                   hypers = list(shape = 5, rate = 1), prop.var = c(3,4), 
                   samp.type = "rw-unif")
x2 <- debinfer_par(name = "x2", var.type = "de", fixed = FALSE, value = 0.5, 
                   prior = "beta", 
                   hypers = list(shape1 = 1, shape2 = 1), prop.var = 0.01, 
                   samp.type = "ind")
x3 <- debinfer_par(name = "x3", var.type = "de", fixed = FALSE, value = 2, 
                   prior = "gamma", 
                   hypers = list(shape = 1, rate = 1), prop.var = 0.1, 
                   samp.type = "rw")
x4 <- debinfer_par(name = "x4", var.type = "de", fixed = FALSE, value = 1, 
                   prior = "gamma", 
                   hypers = list(shape = 5, rate = 1), prop.var = c(4,5), 
                   samp.type = "rw-unif")
x5 <- debinfer_par(name = "x5", var.type = "de", fixed = FALSE, value = 10, 
                   prior = "gamma", 
                   hypers = list(shape = 1, rate = 0.25), prop.var = 5, 
                   samp.type = "rw")
x6 <- debinfer_par(name = "x6", var.type = "de", fixed = FALSE, value = 3, 
                   prior = "unif", 
                   hypers = list(min = 2, max = 10), prop.var = 0.05, 
                   samp.type = "rw")
x7 <- debinfer_par(name = "x7", var.type = "de", fixed = FALSE, value = 3, 
                   prior = "unif", 
                   hypers = list(min = 2, max = 10), prop.var = 0.05, 
                   samp.type = "rw")
x8 <- debinfer_par(name = "x8", var.type = "de", fixed = FALSE, value = 3, 
                   prior = "unif", 
                   hypers = list(min = 2, max = 10), prop.var = 0.05, 
                   samp.type = "rw")
x9 <- debinfer_par(name = "x9", var.type = "de", fixed = FALSE, value = 3, 
                   prior = "unif", 
                   hypers = list(min = 2, max = 10), prop.var = 0.05, 
                   samp.type = "rw")
x10 <- debinfer_par(name = "x10", var.type = "de", fixed = FALSE, value = 3, 
                   prior = "unif", 
                   hypers = list(min = 2, max = 10), prop.var = 0.05, 
                   samp.type = "rw")
#------------------------------Delayed Differential Equation Specifications---------------------
group.X.dede<-function(t,y,p){
  
  x1 <- p["x1"] 
  x2 <- p["x2"] 
  x3 <- p["x3"]
  x4 <- p["x4"]
  x5 <- p["x5"] 
  x6 <- p["x6"] 
  x7 <- p["x7"]
  x8 <- p["x8"] 
  x9 <- p["x9"] 
  x10 <- p["x10"]
  
  Rs <- Ms <- 0 
  lag1 <- lag2 <- 0
  if (t>x6) { 
    lag1 <- lagvalue(t - x6) 
    Rs <- x1 * x2 * lag1[1] 
    }
  phiZ <- x5 * y[2] 
  
  dy1 <- -(x4 + x1) * y[1] 
  dy2 <- Rs - Ms - x3 * y[2] 
  dy3 <- phiZ - (x4 + x1) * y[3]
  
  if(y[1]<0) 
    dy1<-0 
  if (y[2] < 0) { 
    dy2 <- Rs - Ms 
    dy3 <- -(x4 + x1) * y[3] 
    } 
  if (y[3] < 0) { 
    dy3 <- dy3 + (x4 + x1) * y[3] 
    }
  
  list(c(dy1,dy2,dy3)) 
}

group.X.dede.2<-function(t,y,p){
  
  x1 <- p["x1"] 
  x2 <- p["x2"] 
  x3 <- p["x3"]
  x4 <- p["x4"]
  x5 <- p["x5"] 
  x6 <- p["x6"] 
  x7 <- p["x7"]
  x8 <- p["x8"] 
  x9 <- p["x9"] 
  x10 <- p["x10"]
  
  Parameter.1<-Parameter.2 <- 0 
  lag1 <-lag2<-lag3<-lag4 <- 0
  if (t>x6) { 
    lag1 <- lagvalue(t - x6) 
    Parmeter.1 <- x1 * x2 * lag1[1] 
  }
  if(t > x7)
  {
    lag3<-lagvalue(t-x7)
    Parameter.3 <- x5 * lag3[1]
  }
  
  dy1 <- x1* y[1] 
  dy2 <- x3 * y[2] 
  dy3 <- -(x4 + x1) * y[3]
  dy4 <- x4*y[4]
  dy5 <- x5*y[5]
  dy6 <- x6*y[6]
  dy7 <- x7*y[7]
  dy8 <- x8*y[8]
  dy9 <- x9*y[9]
  dy10 <- x10*y[10]
  #
  #----------------------------------------------Inequalites based on context and value
  #
  if(y[1]<0) 
    dy1<-0 
  if (y[2] < 0) { 
    dy2 <- Parameter.1 - Parameter.2 
    dy3 <- -(x4 + x1) * y[3] 
  } 
  if (y[3] < 0) { 
    dy3 <- dy3 + (x4 + x1) * y[3] 
  }
  if(y[4]<0) 
    dy5<-0
  if(y[5]<0) 
    dy5<-0 
  if(y[6]<0) 
    dy6<-0 
  if(y[7]<0) 
    dy7<-0 
  if(y[8]<0) 
    dy8<-0
  if(y[9]<0) 
    dy9<-0
  if(y[10]<0) 
    dy10<-0 
  
  list(c(dy1,dy2,dy3,dy4,dy5,dy6,dy7,dy8,dy9,dy10)) 
}

Pop.initialize<-100
R <- debinfer_par(name = "R", var.type = "init", fixed = TRUE, value = 0) 
S <- debinfer_par(name = "S", var.type = "init", fixed = TRUE, value = 0) 
T <- debinfer_par(name = "T", var.type = "init", fixed = TRUE, value = 0)
U <- debinfer_par(name = "U", var.type = "init", fixed = TRUE, value = 0) 
V <- debinfer_par(name = "V", var.type = "init", fixed = TRUE, value = 0) 
W <- debinfer_par(name = "W", var.type = "init", fixed = TRUE, value = 0)
X <- debinfer_par(name = "X", var.type = "init", fixed = TRUE, value = Pop.initialize) 
Y <- debinfer_par(name = "Y", var.type = "init", fixed = TRUE, value = 0) 
Z <- debinfer_par(name = "Z", var.type = "init", fixed = TRUE, value = 0)

#------------------------Markov Chain Monte Carlo Inference---------------------------------------
mcmc.pars <- setup_debinfer(x1, x2, x3, x4, x5, x6,x7,x8,x9,x10, 
                            X, Y, Z)
mcmc.pars.2 <- setup_debinfer(x1, x2, x3, x4, x5, x6,x7,x8,x9,x10, 
                            R,S,T,U,V,W,X, Y, Z)
iter <- 50
mcmc.model <- de_mcmc(N = iter, data = group.X.df, 
                    de.model = group.X.dede, obs.model = group.X.obs.model, all.params = mcmc.pars, 
                    Tmax = max(group.X.df$t), data.times = c(0,group.X.df$t), 
                    cnt = 50, plot = FALSE, sizestep = 0.1, solver = "dede", verbose.mcmc = FALSE)

mcmc.model.summary<-summary(mcmc.model)

mcmc.model.2 <- de_mcmc(N = iter, data = group.X.df, 
                    de.model = group.X.dede, obs.model = group.X.obs.model.2, all.params = mcmc.pars, 
                    Tmax = max(group.X.df$t), data.times = c(0,group.X.df$t), 
                    cnt = 50, plot = FALSE, sizestep = 0.1, solver = "dede", verbose.mcmc = FALSE)

mcmc.model.2.summary<-summary(mcmc.model.2)
#----------------------------------Simulation--------------------------------------------
burnin <- 10
mcmc.model.post_traj <- post_sim(mcmc.model, n = 10, 
                      times = seq(0,10,by = 0.1), 
                      burnin = burnin, output = "all", prob = 0.95)
mcmc.model.post_traj.2 <- post_sim(mcmc.model.2, n = 10, 
                      times = seq(0,10,by = 0.1), 
                      burnin = burnin, output = "all", prob = 0.95)
mcmc.model.post_traj$HDI$X
#----------------------------------------Tables--------------------------------------
Table.1<-xtable(mcmc.model.summary$statistics)
Table.2<-xtable(mcmc.model.2.summary$statistics)
Table.3<-xtable(head(mcmc.model.post_traj$HDI$X))
#----------------------------------------Figures for the Classroom------------------------------------
Figure.1<-plot(group.X.df, xlab = "Time (days)", ylab = "Test Data", xlim = c(0,10))
par(mfrow = c(3,4)) 
Figure.2<-plot(mcmc.model, ask = FALSE, auto.layout = FALSE)
pairs(mcmc.model, burnin = burnin, scatter = TRUE, trend = TRUE)
Figure.3<-post_prior_densplot(mcmc.model, burnin = burnin)
par(mfrow = c(2,3), mgp = c(2.2, 0.8, 0))
ylabel = expression(paste(Pr,"(", theta,"|", "Y", ")")) 

Figure.4<-post_prior_densplot(mcmc.model, param = "x1",
                              xlab = expression(theta), 
                              ylab = ylabel, show.obs = FALSE, xlim = c(0,8), 
                              main = "x1") 
legend("topright", legend = c("Posterior","Prior"), lty = 1, col = c("black", "red"))

Figure.5<-post_prior_densplot(mcmc.model, param = "x2",
                              xlab = expression(theta), ylab = ylabel, show.obs = FALSE, 
                              xlim = c(-0.1,1.1),main = "x2")
legend("topright", legend = c("Posterior","Prior"), lty = 1, col = c("black", "red"))

Figure.6<-post_prior_densplot(mcmc.model, param = "x3",
                              xlab = expression(theta), ylab = ylabel, show.obs = FALSE, xlim = c(0,3), 
                              main = "x3") 
legend("topright", legend = c("Posterior","Prior"), lty = 1, col = c("black", "red"))

Figure.7<-post_prior_densplot(mcmc.model, param = "x4",
                              xlab = expression(theta), ylab = ylabel, show.obs = FALSE, xlim = c(0,6), 
                              main = "x4") 
legend("topright", legend = c("Posterior","Prior"), lty = 1, col = c("black", "red"))

Figure.8<-post_prior_densplot(mcmc.model, param = "x5",
                              xlab = expression(theta), ylab = ylabel, show.obs = FALSE, xlim = c(0,50), 
                              ylim = c(0,0.2), main = "x5")
legend("topright", legend = c("Posterior","Prior"), lty = 1, col = c("black", "red"))
Figure.9<-post_prior_densplot(dede_rev, param = "x6",
                              xlab = expression(theta), ylab = ylabel, show.obs = FALSE, xlim = c(1.5,10.5), 
                              main = "x6")
legend("topright", legend = c("Posterior","Prior"), lty = 1, col = c("black", "red"))

Figure.10<-post_prior_densplot(mcmc.model, param = "x7",
                              xlab = expression(theta), ylab = ylabel, show.obs = FALSE, xlim = c(1.5,10.5), 
                              main = "x7")
legend("topright", legend = c("Posterior","Prior"), lty = 1, col = c("black", "red"))

Figure.11<-post_prior_densplot(mcmc.model, param = "x8",
                              xlab = expression(theta), ylab = ylabel, show.obs = FALSE, xlim = c(1.5,10.5), 
                              main = "x8")

Figure.12<-post_prior_densplot(dede_rev, param = "x9",
                              xlab = expression(theta), ylab = ylabel, show.obs = FALSE, xlim = c(1.5,10.5), 
                              main = "x9")
legend("topright", legend = c("Posterior","Prior"), lty = 1, col = c("black", "red"))

Figure.13<-post_prior_densplot(mcmc.model, param = "x10",
                              xlab = expression(theta), ylab = ylabel, show.obs = FALSE, xlim = c(1.5,10.5), 
                              main = "x10")
legend("topright", legend = c("Posterior","Prior"), lty = 1, col = c("black", "red"))
par(mfrow = c(1,3)) 
Figure.14<-plot(mcmc.model.post_traj, plot.type = "medianHDI", auto.layout = FALSE)
Figure.15<-plot(mcmc.model.post_traj, plot.type = "ensemble", col = "#FF000040")

Figure.16<-plot(group.X.df, xlab = "Time (days)", ylab = "Test Data", xlim = c(0,10))

for(i in seq_along(mcmc.model.post_traj$sims)) { 
  mcmc.model.sims <- as.data.frame(mcmc.model.post_traj$sims[i]) 
  lines( mcmc.model.sims[,2] ~  mcmc.model.sims[,1]) 
  lines( mcmc.model.sims[,3] ~  mcmc.model.sims[,1],col = "green") 
  lines( mcmc.model.sims[,4] ~  mcmc.model.sims[,1],col = "red") 
}

#--------------------------------------------References----------------------------------

Reference.1<-c("Boersch-Supan, P.H., Ryan, S.J. and Johnson, L.R. (2016).", 
               "deBInfer: Bayesian inference for dynamical models of biological systems.", 
               "arXiv, 1605.00021.")

Reference.2<-c("Johnson, L.R., Pecquerie, L. and Nisbet, R.M. (2013).", 
               "Bayesian inference for bioenergetic models.", 
               "Ecology, 94, 882-894.")

#-------------------------------------------Function Library------------------------------

group.X.obs.model <- function(data, sim.data, samp) {
  ec <- 0.01 
  l.Z <- 0 
  for(i in unique(data$time)){ 
    try(l.Z <- l.Z +sum(dpois(data$count[data$time == i], 
                                    lambda = (sim.data[,"Z"][sim.data[,"time"] == i] + ec), 
                                    log = TRUE))) } 
  llk <- l.Z 
  return(llk)
}

group.X.obs.model.2<-function(data,sim.data,samp) {
  ec<-0.01
  llk.Z <- 0 
  for(i in unique(data$time)){ 
    try(llk.Z <- llk.Z +sum(dpoisinvgauss(data$count[data$time == i], 1,
                                    lambda = (sim.data[,"Z"][sim.data[,"time"] == i] + ec), 
                                    log = TRUE))) } 
  llk <- llk.Z 
  return(llk)
}

