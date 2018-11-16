#-------------------------------------R API ------------------------------------------------------------------
library(pracma); library(xtable)
#--------------------------------------------Data-------------------------------------------------------------

parm.1<-c(1,-1)
parm.2<-c(1,-1,0.5)
f.1<- function(t, y)
  as.matrix(c(parm.1[1]*y[1]*y[2],
              parm.1[2]*y[2]))
f.2 <- function(t, y)
  as.matrix(c(parm.2[1]*y[2]*y[3], 
              parm.2[2]*y[1]*y[3], 
              parm.2[3]*y[1]*y[2]))
y0.1 <- as.matrix(c(0,1))
y0.2 <- as.matrix(c(0, 1, 1))
t0 <- 0; tf <- 20
#-------------------------------------------Computation------------------------------------------------------

solution.1<- ode23(f.1, t0, tf, y0.1, rtol=1e-5, atol=1e-10)
solution.2<- ode23.mutation(f.2, t0, 0.95, tf,y0.2, rtol=1e-5, atol=1e-10)

model.solution.df<-data.frame()
model.solution.df<-cbind(solution.2$step,solution.2$h,solution.2$threshold)
colnames(model.solution.df)<-c("Step Size","h","Threshold")

#--------------------------------------------Tables-----------------------------------------------------------

Table.1<-xtable(model.solution.df)

#--------------------------------------------Figures for Classroom Presentation----------------------------------------------------------

Figure.1<-matplot(solution.1$t, solution.1$y, type = "l", lty = 1, lwd = c(2, 1),
        col = c("darkred", "darkblue"),
        xlab = "Time [min]", ylab= "",
        main = "Two-Objects")
grid()
Figure.2<-matplot(solution.2$t, solution.2$y, type = "l", lty = 1, lwd = c(2, 1, 1),
                  col = c("darkred", "darkblue", "darkgreen"),
                  xlab = "Time [min]", ylab= "",
                  main = "Three-Objects")
grid()

#------Function Library pracma for student modification in class R Studio-----------------------------

ode23.mutation<-function (f, t0, stepSize,tfinal, y0, ..., rtol = 0.001, atol = 1e-06) 
{
  stopifnot(is.numeric(y0), is.numeric(t0), length(t0) == 1, 
            is.numeric(tfinal), length(tfinal) == 1)
  #-----------------Validation----------------------------------------------------
  if (is.vector(y0)) {
    y0 <- as.matrix(y0)
  }
  else if (is.matrix(y0)) {
    if (ncol(y0) != 1) 
      stop("Argument 'y0' must be a vector or single column matrix.")
  }
  #------------------------------------------------------------
  fun <- match.fun(f)
  f <- function(t, y) fun(t, y, ...)
  if (length(f(t0, y0)) != length(y0)) 
    stop("Argument 'f' does not describe a system of equations.")
  eps <- .Machine$double.eps
  realmin <- 1e-100
  tdir <- sign(tfinal - t0)
  threshold <- atol/rtol
  hmax <- abs(0.1 * (tfinal - t0))
  t <- t0
  tout <- t
  y <- y0
  yout <- t(y)
  s1 <- f(t, y)
  r <- max(abs(s1/max(abs(y), threshold))) + realmin
  h <- tdir * stepSize * rtol^(1/3)/r
  while (t != tfinal) {
    hmin <- 16 * eps * abs(t)
    if (abs(h) > hmax) {
      h <- tdir * hmax
    }
    else if (abs(h) < hmin) {
      h <- tdir * hmin
    }
    if (1.1 * abs(h) >= abs(tfinal - t)) 
      h <- tfinal - t
    s2 <- f(t + h/2, y + h/2 * s1)
    s3 <- f(t + 3 * h/4, y + 3 * h/4 * s2)
    tnew <- t + h
    ynew <- y + h * (2 * s1 + 3 * s2 + 4 * s3)/9
    s4 <- f(tnew, ynew)
    e <- h * (-5 * s1 + 6 * s2 + 8 * s3 - 9 * s4)/72
    err <- max(abs(e/max(max(abs(y), abs(ynew)), threshold))) + 
      realmin
    #-----------------------------Check the Error---------------------------------------
    if (err <= rtol) {
      t <- tnew
      y <- ynew
      tout <- c(tout, t)
      yout <- rbind(yout, t(y))
      s1 <- s4
    }
    
    h <- h * min(5, stepSize * (rtol/err)^(1/3))
    if (abs(h) <= hmin) {
      warning("Step size too small.")
      t <- tfinal
    }
  }
  
  return(list(t = c(tout), 
              y = yout,
              step=stepSize,
              h=h,
              threshold=threshold))
}
