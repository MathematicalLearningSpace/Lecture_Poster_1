#-----------------------------R Code To Modify in the Classroom Lecture with Students-----------------------
#----------------------------------R API -------------------------------------------------------------------
library(xtable);library(rnn);library(fpp2);library(astrodatR);library(moonsun);library(magicaxis);library(astro)
library(astrodatR);library(cosmoFns);library(astrolibR);library(sphereplot);library(stellaR);library(PearsonDS);library(readr)
library(HistogramTools);library(DescTools);library(spatstat);library(G2Sd);library(plot3D);library(rgp);library(SPOT)
#---------------------------------Data------------------------------------------------------------------------
data("Sun_spot_num")
Sunspot.1<-sunspot.month
Sunspot.2<-sunspot.year
Sunspot.3<-sunspotarea
Sunspot.4<-Sun_spot_num

x.1<-seq(0,4*pi,length.out=length(Sunspot.2))
x.2<-seq(0,4*pi,length.out=length(Sunspot.2))
k.1<-2^0
k.2<-2^3
y.1<-sin(k.1*x.1)+cos(k.2*x.2)
trig.data<-data.frame(y=y.1,x1=x.1,x2=x.2)

#--------------------------------Pearson Distribution Classification-------------------------------------------
pearson.MLE.1 <- pearsonFitML(Sunspot.3)

pearson.N=length(Sunspot.3)
pIIpars.1 <- list(a=pearson.MLE.1$a,b=pearson.MLE.1$b, location=pearson.MLE.1$location, scale=pearson.MLE.1$scale) 
sunspot.pearson.1<-rpearsonI(pearson.N,params=pIIpars.1)
#--------------------------------Genetic Programming Models----------------------------------------------------
operators<- functionSet("+", "*", "-")
variable.1 <- inputVariableSet("x")
factory.1 <- constantFactorySet(function() rpearsonI(pearson.N,params=pIIpars.1))
interval.1 <- seq(from = 1, to = pearson.N, by = 1) 
fitnessFunction.1 <- function(f) rmse(f(interval.1), sin(interval.1))
set.seed(1) 
gpResult.1 <- geneticProgramming(functionSet = operators, 
                                inputVariables = variable.1, 
                                constantSet = factory.1, 
                                fitnessFunction = fitnessFunction.1, 
                                stopCondition = makeTimeStopCondition(1 * 60))

bestSolution.1 <- gpResult.1$population[[which.min(gpResult.1$fitnessValues)]]
gpResult.1

#----------------------------Symbolic Regression-------------------------------------------------------------

model.1 <- symbolicRegression(Sunspot.3 ~ interval.1, data = Sunspot.3, 
                                stopCondition = makeTimeStopCondition(1 * 60))

result.1 <- symbolicRegression(y ~ x1 + x2, 
                               data = trig.data, 
                              stopCondition = makeTimeStopCondition(30))

best.Fitness.1 <- min(sapply(result.1$population, result.1$fitnessFunction))

bestModel.1 <- model.1$population[[which.min(model.1$fitnessValues)]] 

#---------------------------Symbolic Representations---------------------------------------------------------

variable.4.boolean <- inputVariableSet( 
  "x1" %::% st("logical"), 
  "x2" %::% st("logical"), 
  "x3" %::% st("logical"),
  "x4" %::% st("logical"))

threshold.1<-0.5
factory.1.function <- function() runif(1) > threshold.1 
Fitness.function.parity <- makeBooleanFitnessFunction(parity.3)
factory.1.boolean <- constantFactorySet( "factory.1.function" %::% (list() %->% st("logical")))
function.set.boolean <- functionSet( "&" %::% (list(st("logical"), st("logical")) %->% st("logical")),
                                   "|" %::% (list(st("logical"), st("logical")) %->% st("logical")), 
                                   "!" %::% (list(st("logical")) %->% st("logical")))

typedGp.Result.1 <- typedGeneticProgramming(Fitness.function.parity, st("logical"),
                                          functionSet = function.set.boolean, 
                                          inputVariables = variable.4.boolean, 
                                          constantSet = factory.1.boolean, 
                                          stopCondition = makeTimeStopCondition(30))

bestFunction.1 <- typedGp.Result.1$population[[which.min(typedGp.Result.1$fitnessValues)]]

Best.functions.df<-data.frame()
Best.functions.df<-rbind(c("!(x4 | x3 & x2)"))
colnames(Best.functions.df)<-c("Equations")

#-------------------------------Tables--------------------------------------------------------

Table.1<-xtable(Best.functions.df)

#----------------------------------Figures----------------------------------------------------
Figure.1<-plot(Sunspot.3,main="",ltp="l",lty=1,col="red")
lines(Sunspot.2,lty=2,col="green")
lines(sunspot.pearson.1,lty=3,col="blue")
legend("topright", legend =c("Annual Mean Sunspot Area","Annual Sunspot","Pearson Type 1"), col =c("red","green","blue"),lty = 1:3, cex = .8, y.intersp = 1)

testpos <- seq(0, 1.98*pi, length=270)
                                
Figure.2<-PlotPolar(as.numeric(Sunspot.3), testpos, type="l", main="Annual Mean Sunspot Area", col="blue")
                                
Figure.3<-plot(y = bestSolution.1(interval.1), x = interval.1, type = "l", lty = 1, xlab = "x", ylab = "y") 
lines(y = sin(interval.1), x = interval.1, lty = 2)

Figure.4<-plot(y = bestModel.1(interval.1), x = interval.1, type = "l", lty = 1, xlab = "x", ylab = "y") 
lines(y = Sunspot.3, x = interval.1, lty = 2)

#----------------------------------References-------------------------------------------------
Reference.1<-c("",
               "", 
               "") 
Reference.2<-citation("rgp")
#----------------------------------Function Library-------------------------------------------
parity <- function(x) { 
  numberOfOnes <- sum(sapply(x, function(bit) if (bit) 1 else 0)) 
  numberOfOnes %% 2 != 0 
}
parity.4 <- function(x1, x2, x3,x4) 
{
  parity(c(x1, x2, x3,x4))
}



