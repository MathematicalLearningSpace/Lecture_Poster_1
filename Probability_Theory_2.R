#-----------------------------R Code To Modify in the Classroom Lecture with Students-----------------------
#------------------------------------R API ----------------------------------------------------------------
library(PearsonDS); library(copula);library(VineCopula);library(vines);library(xtable)
#------------------------------------Data Generation-------------------------------------------------------
copula.families<-c("0 = independence copula",
"1 = Gaussian copula",
"2 = Student t copula (t-copula)" ,
"3 = Clayton copula", 
"4 = Gumbel copula",
"5 = Frank copula" ,
"6 = Joe copula" ,
"7 = BB1 copula" ,
"8 = BB6 copula" ,
"9 = BB7 copula" ,
"10 = BB8 copula" ,
"13 = rotated Clayton copula (180 degrees; "survival Clayton")" ,
"14 = rotated Gumbel copula (180 degrees; "survival Gumbel")", 
"16 = rotated Joe copula (180 degrees; "survival Joe")",
"17 = rotated BB1 copula (180 degrees; "survival BB1")",
"18 = rotated BB6 copula (180 degrees; "survival BB6")",
"19 = rotated BB7 copula (180 degrees; "survival BB7")",
"20 = rotated BB8 copula (180 degrees; "survival BB8")",
"23 = rotated Clayton copula (90 degrees)",
"24 = rotated Gumbel copula (90 degrees)" ,
"26 = rotated Joe copula (90 degrees)",
"27 = rotated BB1 copula (90 degrees)" ,
"28 = rotated BB6 copula (90 degrees)" ,
"29 = rotated BB7 copula (90 degrees)" ,
"30 = rotated BB8 copula (90 degrees)" ,
"33 = rotated Clayton copula (270 degrees)",
"34 = rotated Gumbel copula (270 degrees)" ,
"36 = rotated Joe copula (270 degrees)" ,
"37 = rotated BB1 copula (270 degrees)" ,
"38 = rotated BB6 copula (270 degrees)" ,
"39 = rotated BB7 copula (270 degrees)" ,
"40 = rotated BB8 copula (270 degrees)" ,
"104 = Tawn type 1 copula" ,
"114 = rotated Tawn type 1 copula (180 degrees)" ,
"124 = rotated Tawn type 1 copula (90 degrees)" ,
"134 = rotated Tawn type 1 copula (270 degrees)" ,
"204 = Tawn type 2 copula",
"214 = rotated Tawn type 2 copula (180 degrees)", 
"224 = rotated Tawn type 2 copula (90 degrees)" ,
"234 = rotated Tawn type 2 copula (270 degrees)" )

#------------------------------------Data--------------------------------------------------------------------------

##-----------------------------------Simulate from a bivariate Student-t copula-------------------------------------
data.N<-100
experimental.data.distribution.t<- BiCop(family = 2, par = -0.7, par2 = 4)
# ----------------------------------rotated Tawn T2 copula with parameters------------------------------------------
copulaFromFamilyIndex(224, -2, 0.5)
experimental.data.1<- BiCopSim(data.N, experimental.data.distribution.t)
experimental.data.2<- matrix(runif(5 * 100), ncol = 5, nrow = 100)
colnames(experimental.data.2) <- c("A", "B", "C", "D", "E")
#----------------------------------- Compute the empirical Kendall's taus----------------------------------------
M.1<-TauMatrix(experimental.data.1)
M.2<-TauMatrix(experimental.data.2)
#------------------------------------Compute conditional distributions--------------------------------------------
hist(experimental.data.1[, 1])
hist(experimental.data.1[, 2]) 

# simulate (U2 | U1 = 0.5)
experimental.data.3 <- BiCopCondSim(data.N, cond.val = 0.5, cond.var = 1, experimental.data.distribution.t)
# simulate (U1 | U2 = 0.5)
experimental.data.4 <- BiCopCondSim(data.N, cond.val = 0.5, cond.var = 2, experimental.data.distribution.t)
#----------------------------------------Fit the Model---------------------------------------------------------
fit <- BiCopSelect(experimental.data.1[, 1], experimental.data.1[, 2])
summary(fit)
#round(BiCopPDF(experimental.data.1[, 1], experimental.data.1[, 2], fit), 3)
model.pdf.df<-data.frame()
model.pdf.df<-cbind(fit$familyname,
                    fit$par,
                    fit$logLik,
                    fit$AIC,
                    fit$p.value.indeptest,
                    fit$emptau,
                    fit$beta)
colnames(model.pdf.df)<-c("Family","Parameter","LogLik","AIC","P value","Tau","beta")
#-----------------------------------Tables-----------------------------------------------------------------
Table.1<-xtable(model.pdf.df)

#-------------Figures for Classroom Presentation---------------------------------------------------------------
Figure.1<-matplot(experimental.data.1, type="l", log = "y",main = "",xlab = "")
Figure.2<-matplot(experimental.data.2, type="l", log = "y",main = "",xlab = "")
Figure.3<-hist(experimental.data.3)
Figure.4<-hist(experimental.data.4)
Figure.5<-plot(experimental.data.distribution.t)  
Figure.6<-contour(experimental.data.distribution.t) 
Figure.7<-contour(experimental.data.distribution.t, margins = "unif") 

#-----------------------------------Function Library-------------------------------------------------------
selectCopula <- function (vine, j, i, x, y) {
  data <- cbind(x, y)
  fit <- fitCopula(normalCopula(), data, method = "itau")
  fit@copula
}
