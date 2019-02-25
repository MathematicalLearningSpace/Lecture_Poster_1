#-----------------------------R Code To Modify in the Classroom Lecture with Students-----------------------
#---------------------------------------------------R API --------------------------------------------------------
library(xtable);library(BNPTSclust);library(readr)
#--------------------------------------Data--------------------------------------------
n1<-10^2
n2<-10^1
Q<-matrix(0,n1,n2)
X<-matrix(0,n1,n2)
Y<-matrix(0,n1,n2)
Z<-matrix(0,n1,n2)
for(i in 1:n1)
{
  for(j in 1:n2)
  {
    Q[i,j]<-rnorm(1) + rgamma(1,1)
    X[i,j]<-runif(1,0,10) +rnorm(1)
    Y[i,j]<-rgamma(1,1)+runif(1,0,10)
    Z[i,j]<-cos(rgamma(1,5))
  }
}
Q.df<-as.data.frame(Q)
X.df<-as.data.frame(X)
Y.df<-as.data.frame(Y)
Z.df<-as.data.frame(Z)
#-----------------------------------Models---------------------------------------------

L.1 <- tseriescm(Q.df,maxiter =10,level = FALSE,trend = TRUE,deg=1,seasonality = TRUE,
                 c0eps = 2, c1eps = 1, c0beta = 2, c1beta = 1, 
                 c0alpha = 2,c1alpha = 1, priora = FALSE, pia = 0.5, q0a = 1,
                 q1a = 1, priorb = FALSE, q0b = 1, q1b = 1, a = 0.25, b = 0,
                 indlpml = FALSE) 
L.2 <- tseriescm(X.df,maxiter = 10,level = FALSE,trend = TRUE,deg=1,seasonality = TRUE,
                 c0eps = 2, c1eps = 1, c0beta = 2, c1beta = 1, 
                 c0alpha = 2,c1alpha = 1, priora = FALSE, pia = 0.5, q0a = 1,
                 q1a = 1, priorb = FALSE, q0b = 1, q1b = 1, a = 0.25, b = 0,
                 indlpml = FALSE) 
L.3 <- tseriescm(Y.df,maxiter = 10,level = FALSE,trend = TRUE,deg=1,seasonality = TRUE,
                 c0eps = 2, c1eps = 1, c0beta = 2, c1beta = 1, 
                 c0alpha = 2,c1alpha = 1, priora = FALSE, pia = 0.5, q0a = 1,
                 q1a = 1, priorb = FALSE, q0b = 1, q1b = 1, a = 0.25, b = 0,
                 indlpml = FALSE) 
L.4 <- tseriescm(Z.df,maxiter = 10,level = FALSE,trend = TRUE,deg=1,seasonality = TRUE,
                 c0eps = 2, c1eps = 1, c0beta = 2, c1beta = 1, 
                 c0alpha = 2,c1alpha = 1, priora = FALSE, pia = 0.5, q0a = 1,
                 q1a = 1, priorb = FALSE, q0b = 1, q1b = 1, a = 0.25, b = 0,
                 indlpml = FALSE) 

#-------------------------------Analysis------------------------------------------------

Group.Analysis.df<-data.frame()
Group.Analysis.df<-rbind(c("1",L.1$mstar,L.1$HM),
                         c("2",L.2$mstar,L.2$HM),
                         c("3",L.3$mstar,L.3$HM),
                         c("4",L.4$mstar,L.4$HM)
                         )
colnames(Group.Analysis.df)<-c("Model","Groups","Heterogeneity Measure")
Group.Analysis.df

#--------------------------------------Tables------------------------------------------

Table.1<-xtable(Group.Analysis.df)
Table.1
#--------------------------------------Figures------------------------------------------

Figure.1<-plot(Q.df)
Figure.2<-plot(X.df)
Figure.3<-plot(Y.df)
Figure.4<-plot(Z.df)

Figure.5<-clusterplots(L.1,Q.df)
Figure.6<-diagplots(L.1)

Figure.7<-clusterplots(L.2,X.df)
Figure.8<-diagplots(L.2)

Figure.9<-clusterplots(L.3,Y.df)
Figure.10<-diagplots(L.3)

Figure.11<-clusterplots(L.4,Z.df)
Figure.12<-diagplots(L.4)
#------------------------------------Figure Group-------------------------------------
par(mfcol = c(2, 2))
Figure.13<-hist(L.1$gnstar)
Figure.14<-hist(L.2$gnstar)
Figure.15<-hist(L.3$gnstar)
Figure.16<-hist(L.4$gnstar)

#--------------------------------------References---------------------------------------
Reference.1<-c("Nieto-Barajas, L.E. & Contreras-Crist´an A. (2014)",
"A Bayesian Nonparametric Approach for Time Series Clustering.", 
"Bayesian Analysis, Vol. 9, No. 1 pp. 147-170.")

#-------------------------------------Function Library----------------------------------
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
