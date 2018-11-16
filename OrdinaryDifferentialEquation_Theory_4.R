#----------------------------------R API ----------------------------------------------------------------
library(pracma);library(xtable);library(deSolve);library(ReacTran);library(rootSolve);library(fda);library(igraph)
library(boot);library(sampling);library(PearsonDS);library(KEGGgraph);library(keggorthology);library(KEGGREST)
library(Rgraphviz);library(igraph);library(jsonlite);library(tsDyn);library(tseriesChaos)
library(wavelets);library(waveslim);library(readr)
#---------------------------------------Sample Gene List--------------------------------------------------------------------------------------------------
Gene.Study<-c("LD","ATRp","p53","p53p","Mdm2c","Mdm2cp","Mdm2n","Aktp","PIP3","PTEN","DDB2","p21tot","p21CE","CycEtot","E2F1","Rbp","Bax","CytoC",
"Apaf1","Apops","Casp9","Casp3")

Table.Parameter.p21.df <- as.data.frame(t(read.csv("Parameter_Table_P21.txt", sep="")))
View(Table.Parameter.p21.df)
#---------------------------------------Temporal Resolution-----------------------------------------------------------------------------------
times <- seq(0, 10^2, by = 0.1)
#-------------------------------------Classroom Scientific Article Example----------------------------------------------------------
#-------------------------------------DNA Damage Repair Example Adapted from-------------------------------------------------------------------------------- 
#-------------------------------------Journal Article:---------------------------------------------------------------------------------------
#-------------------------------------Li H1, Zhang XP, Liu F. (2013). Coordination between p21 and DDB2 
#-------------------------------------in the cellular response to UV radiation PLoS One. 
#-------------------------------------8(11):e80111. doi: 10.1371/journal.pone.0080111. eCollection 2013.-------------------------------------
#-------------------------------------https://www.ncbi.nlm.nih.gov/pubmed/24260342-----------------------------------------------------------
ATR<-0.1;Akt<-0.1;PIP2<-0.1;p21<-0.1;CycE<-0.1;kdmdm2n<-0.1;RE<-0.1;Rb<-0.1
alpha.1<-4;alpha.2<-2;alpha.3<-7
variables.initial.test.5 <- c(LD=5,ATRp=0.1,p53=0.1,p53p=0.1,Mdm2c=0.1,Mdm2cp=0.1,Mdm2n=0.1,Aktp=0.1,PIP3=0.1,
                              PTEN=0.1,DDB2=0.1,p21tot=0.1,p21CE=0.1,CycEtot=0.1,E2F1=0.1,Rbp=0.1,
                              Bax=0.1,CytoC=0.1,Apaf1=0.1,Apops=0.1,Casp9=0.1,Casp3=0.1)
variables.initial.test.20 <- c(LD=20,ATRp=0,p53=0,p53p=0,Mdm2c=0,Mdm2cp=0,Mdm2n=0,Aktp=0,PIP3=0,
                               PTEN=0,DDB2=0,p21tot=0,p21CE=0,CycEtot=0,E2F1=0,Rbp=0,Bax=0,CytoC=0,Apaf1=0,Apops=0,Casp9=0,Casp3=0)
parameters.test <- c(krepair<-0.01, 
                     kacatr0<-0.001,kacatr<-2.0,jacatr<-1.0,kdeatr<-1.5,jdeatr<-2.5,jd<-2.0,jatr<-1.0,ATRtot<-5.0,   
                     ksp53<-0.04,kdp530<-0.03,kdp53<-0.2,jdp53<-0.1,kdp53p<-0.01,jdp53p<-0.1,kacp530<-0.2,kdep53<-0.1,   
                     kdmdm2<-0.009,ksmdm20<-0.002,ksmdm2<-0.01,jsmdm2<-1.0,kpmdm2<-4.0,jpmdm2<-0.3,kdpmdm2<-0.3,   
                     jdpmdm2<-0.1,kinmdm2c<-0.06,koutmdm2n<-0.09,kdmdm2n0<-0.05,kpakt<-0.25,kdpakt<-0.1,jpakt<-0.1,   
                     jdpakt<-0.2,Akttot<-1.0,PIPtot<-1.0,kp2<-0.1,kp3<-0.45,jp2<-0.2,jp3<-0.4,   
                     ksPTEN0<-0.001, ksPTEN<-0.05,jsPTEN<-2.0,kdPTEN<-0.01,ksDDB20<-0.01,ksDDB2<-0.6,jsDDB2<-2.0,kdDDB2<-0.15,  
                     ksp210<-0.01,ksp21<-0.4,jsp21<-0.7,kdp210<-0.1,kdp21<-0.2,kscyce0<-0.0005,kscyce<-0.0275,jscyce<-0.2, 
                     kdcyce<-0.005,kasre<-0.5,kdsre<-0.05,kprb<-0.05,jprb<-0.1,kdprb<-0.025,jdprb<-0.1,kasp21ce<-0.5,   
                     kdsp21ce<-0.05,Rbtot<-2.0,E2F1tot<-1.0,ksbax0<-0.01,ksbax<-0.8,jsbax<-2.0,kdbax<-0.1,jbax<-3.0,
                     kaccytoc0<-0.01,kaccytoc<-1.0,kdecytoc<-0.1,ksapaf10<-0.001,ksapaf1<-0.6,jsapaf1<-0.6,kdapaf1<-0.1,  
                     kacapops<-5.0,kdeapops<-0.5,kaccasp90<-0.001,kaccasp9<-3.0,jcasp3<-0.5,kdecasp9<-0.05,kaccasp30<-0.001, 
                     kaccasp3<-0.1,jcasp9<-0.5, kdecasp3<-0.07,CytoCtot<-5.0,Casp9tot<-3.0, Casp3tot<-3.0)
parameter.test.df<-as.data.frame(parameters.test)
rownames(parameter.test.df)<-c("EQ1:krepair",
                               "EQ2:kacatr0","EQ2:kacatr","EQ2:jacatr","EQ2:kdeatr","EQ2:jdeatr","EQ2:jd",
                               "EQ3:jatr","EQ2:ATRtot",   
                               "EQ3:ksp53","EQ3:kdp530","EQ3:kdp53","EQ3:jdp53","EQ3:kdp53p","EQ3:jdp53p",
                               "EQ3:kacp530","EQ3:kdep53","EQ3:kdmdm2",
                               "EQ4:ksmdm20","EQ4:ksmdm2","EQ4:jsmdm2","EQ4:kpmdm2","EQ4:jpmdm2","EQ4:kdpmdm2",   
                               "EQ4:jdpmdm2","EQ4:kinmdm2c","EQ4:koutmdm2n","EQ4:kdmdm2n0",
                               "EQ5:kpakt","EQ5:kdpakt","EQ5:jpakt","EQ5:jdpakt","EQ5:Akttot","EQ5:PIPtot",
                               "EQ6:kp2","EQ6:kp3","EQ6:jp2","EQ6:jp3",   
                               "EQ7:ksPTEN0","EQ7:ksPTEN","EQ7:jsPTEN","EQ7:kdPTEN","EQ7:ksDDB20",
                               "EQ7:ksDDB2","EQ7:jsDDB2","EQ7:kdDDB2",  
                               "EQ8:ksp210","EQ8:ksp21","EQ8:jsp21","EQ8:kdp210","EQ8:kdp21","EQ8:kscyce0","EQ8:kscyce","EQ8:jscyce", 
                               "EQ9:kdcyce","EQ9:kasre","EQ9:kdsre","EQ9:kprb","EQ9:jprb","EQ9:kdprb","EQ9:jdprb","EQ9:kasp21ce",   
                               "EQ9:kdsp21ce","EQ9:Rbtot",
                               "EQ10:E2F1tot",
                               "EQ11:ksbax0","EQ11:ksbax","EQ11:jsbax","EQ11:kdbax","EQ11:jbax",
                               "EQ12:kaccytoc0","EQ12:kaccytoc","EQ12:kdecytoc","EQ12:ksapaf10","EQ12:ksapaf1","EQ12:jsapaf1","EQ12:kdapaf1",  
                               "EQ12:kacapops","EQ12:kdeapops",
                               "EQ13:kaccasp90","EQ13:kaccasp9","EQ13:jcasp3","EQ13:kdecasp9","EQ13:kaccasp30", 
                               "EQ13:kaccasp3","EQ13:jcasp9", "EQ13:kdecasp3","EQ13:CytoCtot","EQ13:Casp9tot", "EQ13:Casp3tot")
#-------------------------------------------------Classroom Model A-----------------------------------------------------------
system.equation.model.test<-function(times, 
                                     variables.intitial.test, 
                                     parameters.test)
{
  with(as.list(c(parameters.test, variables.intitial.test)), 
       {
#----------------------------------------------EQ1:LD Model------------------------------------------------------------
          d.1.LD.dt.1 = - krepair*H(LD) 
#----------------------------------------------EQ2:ATR Model-----------------------------------------------------------
          d.1.ATRp.dt.1 = (kacatr0 + kacatr*(LD/(LD + jd))*ATRp)*(ATR/(ATR + jacatr)) - kdeatr* (ATRp/(ATRp + jdeatr)) 
          ATR = ATRtot - ATRp 
#----------------------------------------------EQ3:p53 Model---------------------------------------------------------------
         kacp53 = kacp530*(ATRp/(ATRp + jatr)) 
         d.1.p53.dt.1 = ksp53 - kdp530*p53 - kdp53*Mdm2n * (p53/(p53 + jdp53)) - kacp53*p53 + kdep53*p53p 
         d.1.p53p.dt.1 = kacp53*p53 - kdep53*p53p - kdp53p*Mdm2n *(p53p/(p53p + jdp53p)) 
#----------------------------------------------EQ4:MDM2 Model--------------------------------------------------------------
         d.1.Mdm2c.dt.1 = ksmdm20 + ksmdm2*(p53p^(alpha.1)/(p53p^(alpha.1) + jsmdm2^(alpha.1))) 
         - kdmdm2*Mdm2c + kdpmdm2*(Mdm2cp/(Mdm2cp + jdpmdm2))
         - kpmdm2*(Aktp)*(Mdm2c/(Mdm2c + jpmdm2)) 
         d.1.Mdm2cp.dt.1 = kpmdm2*Aktp*(Mdm2c/ (Mdm2c + jpmdm2)) - kdpmdm2*(Mdm2cp/(Mdm2cp + jdpmdm2))
         - kinmdm2c*Mdm2cp + koutmdm2n*Mdm2n - kdmdm2*Mdm2cp 
         d.1.Mdm2n.dt.1 = kinmdm2c*Mdm2cp - koutmdm2n*Mdm2n - (kdmdm2 + kdmdm2n)*Mdm2n 
         kdmdm2n = kdmdm2n0*(ATRp/(ATRp + jatr))
#--------------------------------------------EQ5:AKT Model----------------------------------------------------------------
         d.1.Aktp.dt.1= kpakt*PIP3*(Akt/(Akt + jpakt)) - kdpakt*(Aktp/(Aktp + jdpakt))
         Akt = Akttot - Aktp 
#--------------------------------------------EQ6:PIP3 Model----------------------------------------------------------------
         d.1.PIP3.dt.1 = kp2*(PIP2/(PIP2 + jp2)) - kp3*PTEN*(PIP3/(PIP3 + jp3)) 
         PIP2= PIPtot - PIP3 
#--------------------------------------------EQ7:Pten Model---------------------------------------------------------------
         d.1.PTEN.dt.1 = ksPTEN0 + ksPTEN*(p53p^(alpha.1)/(p53p^(alpha.1) + jsPTEN^(alpha.1))) - kdPTEN*PTEN 
         d.1.DDB2.dt.1 = ksDDB20 + ksDDB2*(p53p^(alpha.1)/(p53p^(alpha.1) + jsDDB2^(alpha.1))) - kdDDB2*DDB2 
#--------------------------------------------EQ8:P21 model----------------------------------------------------------------
         d.1.p21tot.dt.1 = ksp210 + ksp21*(p53p^(alpha.1)/(p53p^(alpha.1) + jsp21^(alpha.1))) 
         - kdp210*p21tot - kdp21*DDB2*p21tot 
         d.1.p21CE.dt.1 = kasp21ce*p21*CycE - kdsp21ce*p21CE 
         p21 = p21tot - p21CE 
#------------------------------------------EQ9:Cyc Model------------------------------------------------------------------
         d.1.CycEtot.dt.1 = kscyce0 + kscyce*(E2F1^(alpha.2)/(E2F1^(alpha.2) + jscyce^(alpha.2))) - kdcyce*CycEtot 
         CycE = CycEtot-p21CE
         d.1.E2F1.dt.1 = - kasre*Rb*E2F1 + kdsre*RE 
#------------------------------------------EQ10:RE-rbp-rb model-----------------------------------------------------------
         RE = E2F1tot - E2F1 
         d.1.Rbp.dt.1 = kprb*CycE*(Rb/(Rb + jprb)) - kdprb*(Rbp/(Rbp + jdprb)) 
         Rb = Rbtot - Rbp - RE 
#-------------------------------------------EQ11:BAX model------------------------------------------------------------------
         d.1.Bax.dt.1 = ksbax0 + ksbax*(p53p^(alpha.1)/(p53p^(alpha.1) + jsbax^(alpha.1))) - kdbax*Bax 
#-------------------------------------------EQ12:Cytoc-Apaf1-Apops Model---------------------------------------------------
         d.1.CytoC.dt.1 = kaccytoc0 + kaccytoc*(Bax^(alpha.1)/(Bax^(alpha.1) + jbax^(alpha.1)))*(CytoCtot - CytoC) - kdecytoc*CytoC 
         d.1.Apaf1.dt.1 = ksapaf10 + ksapaf1*(E2F1^(alpha.2)/(E2F1^(alpha.2) + jsapaf1^(alpha.2))) - kdapaf1*Apaf1 
         d.1.Apops.dt.1 = kacapops*((CytoC - (alpha.3)*Apops)*(Apaf1 - (alpha.3)*Apops))^(alpha.3) - kdeapops*Apops 
#-------------------------------------------EQ13:Casp Model---------------------------------------------------------------- 
         d.1.Casp9.dt.1 = (kaccasp90 + kaccasp9*(Casp3^(alpha.1)/(Casp3^(alpha.1) + jcasp3^(alpha.1))))*(Casp9tot - Casp9) - kdecasp9*Casp9
         d.1.Casp3.dt.1 = (kaccasp30 + kaccasp3*(Casp9^(alpha.1)/(Casp9^(alpha.1) + jcasp9^(alpha.1))))*(Casp3tot - Casp3) - kdecasp3*Casp3
         
         res <- c(d.1.LD.dt.1,
                  d.1.ATRp.dt.1,
                  d.1.p53.dt.1,
                  d.1.p53p.dt.1,
                  d.1.Mdm2c.dt.1,
                  d.1.Mdm2cp.dt.1,
                  d.1.Mdm2n.dt.1,
                  d.1.Aktp.dt.1,
                  d.1.PIP3.dt.1,
                  d.1.PTEN.dt.1,
                  d.1.DDB2.dt.1,
                  d.1.p21tot.dt.1,
                  d.1.p21CE.dt.1,
                  d.1.CycEtot.dt.1,
                  d.1.E2F1.dt.1,
                  d.1.Rbp.dt.1,
                  d.1.Bax.dt.1,
                  d.1.CytoC.dt.1,
                  d.1.Apaf1.dt.1,
                  d.1.Apops.dt.1,
                  d.1.Casp9.dt.1,
                  d.1.Casp3.dt.1)
         list(res)
       })
}
#-----------------------------------------------Solution of Equation System------------------------------------------
system.equation.model.test.solution.1 <- ode(y = variables.initial.test.5 , 
                                             times = times, 
                                             func = system.equation.model.test, 
                                             parms = parameters.test) 
system.1.summary<-summary(system.equation.model.test.solution.1.df)
system.equation.model.test.solution.1.df<-as.data.frame(system.equation.model.test.solution.1)
system.equation.model.df<-data.frame(system.1.summary)
# speed after 0.5, 1, 1.5, 2 seconds
system.summary.time.P53<-cbind(c(0.5,1,1.5,2), 
                               -deval(system.equation.model.test.solution.1.df$time, 
                                      system.equation.model.test.solution.1.df$p53, 
                                      c(0.5, 1, 1.5, 2)))

#--------------------------------------------Analysis------------------------------------------------------------------

hurst.p53<-hurstexp(system.equation.model.test.solution.1.df$p53)
hurst.Mdm2cp<-hurstexp(system.equation.model.test.solution.1.df$Mdm2cp)

entropy.p53<-approx_entropy(system.equation.model.test.solution.1.df$p53, edim = 2)
entropy.Mdm2cp<-approx_entropy(system.equation.model.test.solution.1.df$Mdm2cp, edim = 2)

System.Analysis.Artifacts.df<-data.frame()
System.Analysis.Artifacts.df<-rbind(c(hurst.p53,entropy.p53),
                                    c(hurst.Mdm2cp,entropy.Mdm2cp))
rownames(System.Analysis.Artifacts.df)<-c("p53","Mdm2cp")
colnames(System.Analysis.Artifacts.df)<-c("Hs","Hrs","He","Hal","ht","Entropy")

#------------------------------------------Classification-------------------------------------------------------------
path.classification.df<-data.frame()
path.classification.df<-rbind(c(classify.orbit(times,system.equation.model.test.solution.1.df$p53,
                                               "p53",
                                               min(system.equation.model.test.solution.1.df$p53),max(system.equation.model.test.solution.1.df$p53))))

rownames(path.classification.df)<-c("p53")
#--------------------------------------Tables-------------------------------------------------------------------------
Table.1<-xtable(parameters.test.df)
Table.2<-xtable(system.1.summary)
Table.3<-xtable(System.Analysis.Artifacts.df)
#--------------------------------------Figures-------------------------------------------------------------------------
Figure.1<-plot(seq(0,10,0.1),f(seq(0,10,0.1),1,4),lty=1,col="black")
lines(f(seq(0,10,0.1),0.5,4),lty=2,col="red")
lines(f(seq(0,10,0.1),1.5,3),lty=3,col="blue")
lines(f(seq(0,10,0.1),1.5,3),lty=4,col="green")
Figure.2<-plot(system.equation.model.test.solution.1.df[,2:7])
Figure.3<-plot(times,system.equation.model.test.solution.1.df$p53p)
Figure.3.caption<-c("Temporal evolution of the levels of DNA damage, ATRp, p53p, 
p21tot, Bax, and Casp3 (from top to bottom) at  (A) or 20 (B)","'DNA damage', 'ATRp', 'p53p', 'p21tot', 'Bax', and 'Casp3'")
Figure.4<-hist(system.equation.model.test.solution.1[,2:14], 
               col = grey(seq(0, 1, by = 0.1)), 
               mfrow = c(2, 3))
mtext(outer = TRUE, side = 3, "Species Histograms", cex = 1.5)
require(car)
Figure.5<-scatterplot.matrix(system.equation.model.test.solution.1[,10:13])

#-------------------------------------Function Library Examples For Classroom-------------------------------------------------------
H<-function(x)
{
  result<-0
  if( x >=0){result<-1}else{result<-0}
  return(result)
}

f<-function(x,a,k)
{
  result<-0
  result<-x^k/(x^k+a^k)
  return(result)
}

classify.orbit<-function(t,x,geneName,a,b)
{
  j=0;
  category<-c("A","B","C","D")
for(i in 1:length(x))
{
  if(geneName=="p53")
  {
    if(x[i]>= a & x[i] <= b)
    {
      j=j+1
    }
  }
  
}
  intervalRatio<-j/length(x)
  if(intervalRatio < 0.25)
  {
    category="A"
  }
  else if (intervalRatio > 0.25 & intervalRatio < 0.5 ) category="B"
  else if (intervalRatio > 0.5 & intervalRatio < 0.75)  category="C"
  else if (intervalRatio > 0.75)  category="D"
  
  return(list(ratio=intervalRatio,
              classify=category))
}

#-----------------------------------Learning Space-----------------------------------------------------------------------



 
