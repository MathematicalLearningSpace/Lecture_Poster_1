
library(xtable)
library(Matrix)
library(combinat)
library(igraph)

#------------------Data----------------------------------------------------------------
n<-100;h<-1;t<-(n+1)*h
parameters.1<-c(n,h,t)
sequence<-seq(1,t,1)
#-----------------Parameter Table------------------------------------------------------
sigma<-0.1181;rho<-1.131;eta<-20.19;mu<-0.00311;delta<-0.3743;alpha<-1.636;beta<-0.002

sigma_LB<-0.05;rho_LB<-1;eta_LB<-20;mu_LB<-0.001;delta_LB<-0.3;alpha_LB<-1;beta_LB<-0.001;
sigma_UB<-0.2;rho_UB<-2;eta_UB<-25;mu_UB<-0.005;delta_UB<-0.5;alpha_UB<-2;beta_UB<-0.003;

parameter.2<-c(sigma,rho,eta,mu,delta,alpha,beta)

parameter.2.lb<-c(sigma_LB,rho_LB,eta_LB,mu_LB,delta_LB,alpha_LB,beta_LB)
parameter.2.ub<-c(sigma_UB,rho_UB,eta_UB,mu_UB,delta_UB,alpha_UB,beta_UB)

biotic.capacity<-1/beta

Parameter.Table.df<-data.frame()
Parameter.Table.df<-rbind(c("sigma","Source rate of effector cells",sigma,sigma_LB,sigma_UB,"Reference.1"),
                          c("rho","Maximum rate of effector cells proliferation",rho,rho_LB,rho_UB,"Reference.1"),
                          c("eta","Half saturation constant",eta,eta_LB,eta_UB,"Reference.1"),
                          c("mu","Effector cells inactivation rate by mutated cells",mu,mu_LB,mu_UB,"Reference.1"),
                          c("delta","The natural survival rate of effector cells",delta,delta_LB,delta_UB,"Reference.1"),
                          c("alpha","The maximum growth rate of mutated cells",alpha,alpha_LB,alpha_UB,"Reference.1"),
                          c("beta","1/beta is the maximal carry capacity",beta,beta_LB,beta_UB,"Reference.1"))
colnames(Parameter.Table.df)<-c("Parameters","Description", "Values","Lower Bound", "Upper Bound","References")

#-----------------Equation Specification of System---------------------------
f.x<-function(x,y,sigma, rho, eta,mu, delta)
{
  dx.dt<-sigma + (rho*x*y/(eta+y)) - mu*x*y - delta*x
  return(dx.dt)
}
f.y<-function(x,y,alpha,beta)
{
  dy.dt<-alpha*y*(1-beta*y)-x*y
  return(dy.dt)
}

#-----------------Test the Algorithm-------------------------------------
simulation.parameters.1<-c(t,n,0.75,h,1,1)
test.algorithm<-Grunwald.Letnikov.algorithm(simulation.parameters.1,
                                 parameter.2)
test.algorithm

Jacobian.E<-Jacobian.equilibrium(sigma/delta,0,sigma, rho, eta,mu, delta,alpha,beta)

Jacobian.E

#------------------Simulation of the System---------------------------------------
simulation.parameters.1<-c(t,n,0.25,h,1,10)

sigma<-0.1181;rho<-1.131;eta<-20.19;mu<-0.00311;delta<-0.3743;alpha<-1.636;beta<-0.002
parameter.2<-c(sigma,rho,eta,mu,delta,alpha,beta)

simulation.1<-Grunwald.Letnikov.algorithm(simulation.parameters.1,parameter.2)
simulation.1
#----------------Variations in sigma and initial conditions----------------------
simulation.parameters.5<-c(t,n,0.75,h,1,20)
parameter.2.1<-c(0.01,rho,eta,mu,delta,alpha,beta)

simulation.5<-Grunwald.Letnikov.algorithm(simulation.parameters.5,parameter.2.1)

simulation.5

simulation.results.df<-data.frame()
simulation.results.df<-rbind(c(simulation.parameters.1,parameter.2),
                             c(simulation.parameters.5,parameter.2.1)
)
colnames(simulation.results.df)<-c("t","n","q","h","x.initial","y.initial","sigma","rho","eta","mu","delta","alpha","beta")

#-------------------------Tables-----------------------------------------------------------

Table.1<-xtable(Parameter.Table.df)
Table.1.caption<-c("Parameter values")
Table.2<-xtable(simulation.results.df)
Table.2.caption<-c("Simulation Results")

#------------------------Figures----------------------------------------------------------

par(mfrow = c(1,2), mar = c(5,5,5,5))
Figure.1<-plot(sequence,simulation.1$x,xlab='t-time Days',ylab='Effector Cells',
               main='',col='black',type="l")
#lines(simulation.2$x,lty=2,col="red")
#lines(simulation.3$x,lty=3,col="green")
#lines(simulation.4$x,lty=4,col="blue")
legend("topright", legend=c("q=0.25","q=0.5","q=0.75","q=1")
       , bty = "n",lwd=2, 
       cex=0.75, lty=1:4,col=c("black","red","green","blue"))
Figure.1A<-plot(sequence,simulation.1$y,xlab='t-time Days',ylab='Mutant Cells',
                main='',col='black',ylim=c(0,300),type="l")
#lines(simulation.2$y,lty=2,col="red")
#lines(simulation.3$y,lty=3,col="green")
#lines(simulation.4$y,lty=4,col="blue")
legend("bottomright", legend=c("q=0.25","q=0.5","q=0.75","q=1")
       , bty = "n",lwd=2, 
       cex=0.75, lty=1:4,col=c("black","red","green","blue"))
Figure.1.caption<-c("Effect of parameter alpha_1 and alpha_2")


Figure.3<-plot(simulation.1$x, simulation.1$y,type="l", xlab='Effector Cells',ylab='Mutant Cells',
               main='',col='black')
legend("topright", legend=c("q=0.25")
       , bty = "n",lwd=2, 
       cex=0.75, lty=1:1)
Figure.3.caption<-c("Phase Portraits with x and y with (delta, sigma and q)")

par(mfrow = c(1,2), mar = c(5,5,5,5))
Figure.5<-plot(sequence, simulation.5$x, xlab='',ylab='Cells',
               main='',col='Green',type="l")
legend("bottomright", legend=c("")
       , bty = "n",lwd=2, 
       cex=0.75, lty=1:1)
Figure.5A<-plot(sequence, simulation.5$y, xlab='',ylab='Cells',
                main='',col='Green',type="l")
legend("bottomright", legend=c("")
       , bty = "n",lwd=2, 
       cex=0.75, lty=1:1)
Figure.5.caption<-c("(a) X and (b) Y with sigma and initial conditions")

#------------------------------------References----------------------------------------------------

Reference.1<-c("Sadia Arshad1, Dumitru Baleanu Jianfei Huang, Yifa Tang and Maysaa Mohamed Al Qurashi",
               "Dynamical analysis of fractional order model of immunogenic tumors",
               "Advances in Mechanical Engineering 2016, Vol. 8(7) 1-13")
Reference.2<-c("D'Onofrio A.", 
               "A general framework for modeling tumor immune system competition and immunotherapy: mathematical analysis and biomedical inferences.", 
               "Physica D 2005; 208: 220-235.")
Reference.3<-c("Scherer R, Kalla SL, Tang Y, et al.", 
               "The Grunwald-Letnikov method for fractional differential equations.", 
               "Comput Math Appl 2011; 62: 902-917.")

#---------Uncurated Function Library (Work in Progress with Updates)------------------------------------

Jacobian.equilibrium<-function(xbar,ybar,sigma,rho,eta,mu,delta,alpha,beta)
{
  Jacobian<-matrix(c((rho*ybar/(eta+ybar)-mu*ybar -delta),
                     ((rho*xbar*eta)/(eta+ybar)^2-mu*xbar),
                     (-1)*ybar,
                     (alpha -2*alpha*beta*ybar - xbar)),nrow=2,ncol=2)
  P1<-Jacobian[1,1]
  P2<-Jacobian[2,2]
  P3<-Jacobian[1,2]*Jacobian[2,1]
  m1<--(P1+P2)
  m2<-P1*P2 + P3
  
  lambda1 <- 0.5*(-m1+(m1^2-4*m2)^(1/2))
  lambda2 <- 0.5*(-m1-(m1^2-4*m2)^(1/2))
  
  return(list(Jacobian=Jacobian,lamda1=lambda1,lambda2=lambda2))
}

Grunwald.Letnikov.algorithm<-function(parms,parms.2)
{
  x.n<-NULL;y.n<-NULL;w.k<-NULL;sum.x<-NULL;sum.y<-NULL;
  g.correction<-NULL;
  x.n[1]<-parms[5]
  y.n[1]<-parms[6]
  w.k[1]<-parms[3]
  sum.x<-0
  sum.y<-0
  h<-parms[4]
  q<-parms[3]
  #-----------------Step 1
  for(i in 1:n)
  {
    #-------------- Step 2
    w.k[i+1]<-(1-((q+1)/i))*w.k[i]
    for(j in 1:i)
    {
      #-------------Step 3
      sum.x[j+1]<-sum.x[j]+w.k[j]*x.n[i+1-j]
      sum.y[j+1]<-sum.y[j]+w.k[j]*y.n[i+1-j]
    }
    #---------------Step 4
    #g.correction[i]<-(parms[2]^(i+1))/gamma(i+1-q)
    g.correction[i]<-0
    #---------------Step 5
    x.n[i+1]<-(h^q)*f.x(x.n[i],y.n[i],parms.2[1],parms.2[2],parms.2[3],parms.2[4],parms.2[5])+ sum.x[i] + g.correction[i]*x.n[1]
    y.n[i+1]<-(h^q)*f.y(x.n[i],y.n[i],parms.2[6],parms.2[7]) + sum.y[i] +  g.correction[i]*y.n[1]
  }
  
  return(list(x=x.n,
              y=y.n,
              x.sum=sum.x,
              y.sum=sum.y,
              w.k=w.k,
              g.correction=g.correction))
}

