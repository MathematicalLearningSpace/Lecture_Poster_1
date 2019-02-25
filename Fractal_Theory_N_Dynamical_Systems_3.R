#-----------------------------R Code To Modify in the Classroom Lecture with Students-----------------------
library(fractal);library(xtable);library(fracdiff);library(pracma);library(RandomFields)

#-------------------------------------------------Data-----------------------------------------------------------
N<-8
z <- fractalcurve(4, which = "hilbert")
z.2 <- fractalcurve(N, which="molecule")
z.3 <- fractalcurve(N, which="dragon")
#-------------------------------------------------Hurst Exponent------------------------------------------------
hurst.1<-hurstexp(f1(z$x, z$y,1,1))
hurst.2<-hurstexp(f1(z$x, z$y,2,2))
#-------------------------------------------------Transformations------------------------------------------------
f <- function(x) x * cos(0.1*exp(x)) * sin(0.1*pi*exp(x))
F <- function(t) c(t, f(t))
L <- arclength(F, 0, 5, tol = 1e-12, nmax = 25)
print(L$length, digits = 16)

t<-seq(0,3*pi)
f.arclength<-arclength(f, 0, 3*pi)
#------------------------------------------------Sine------------------------------------------------------------
alpha<-1;beta<-1
f <- function(t) c(t, alpha*sin(beta*t))
sine.length<-arclength(f, 0, pi)          
#------------------------------------------------Ellipse--------------------------------------------------------
a <- 1.0; b <- 0.5
ellipse.ratio<-a/b
f <- function(t) c(a*cos(t), b*sin(t))
L <- arclength(f, 0, 2*pi, tol = 1e-10)   
#------------------------------------------------Elliptic integral of the second kind---------------------------
ellipticity <- sqrt(1 - b^2/a^2)  # ellipticity
L <- 4 * a * ellipke(ellipticity^2)$e               
#------------------------------------------------Fermat Spiral--------------------------------------------------

f <- function(t) 0.25 * sqrt(t) * c(cos(t), sin(t))
t1 <- 0; t2 <- 6*pi
a  <- 0; 
b  <- arclength(f, t1, t2)$length

fParam <- function(w) {
  fct <- function(u) arclength(f, a, u)$length - w
  urt <- uniroot(fct, c(a, 6*pi))
  output<-list()
  output$w<-w
  output$fct<-fct
  output$Uniroot<-urt$root
  return(output)
}

ts <- linspace(0, 6*pi, 250)

#------------------------------------------------Fit the Curve---------------------------------------------------
N <- 50
k1<-1
k2<-0.05
u <- linspace(0, pi, N)
x <- cos(u) + k1 * randn(1, N)
y <- sin(u) + k2 * randn(1, N)
n <- 8
cfit1 <- curvefit(u, x, y, n)
hurst.3<-hurstexp(f1(cfit1$xp, cfit1$yp,2,2))
hurst.4<-hurstexp(f1(cfit1$xp, cfit1$yp,1,1))

hurst.transformations.df<-data.frame()
hurst.transformations.df<-rbind(c("f1","Hilbert","1","1",hurst.1),
                                c("f1","Hilbert","1","1",hurst.2),
                                c("f1","FitCurve","2","2",hurst.3),
                                c("f1","FitCurve","1","1",hurst.4))
colnames(hurst.transformations.df)<-c("Function","Fractal","P1","P2","Hs","Hrs","He","Hal","Ht")

#------------------------------------------------Archimedian spiral-----------------------------------------------
n <- 8
u <- linspace(0, 3*pi, 50)
a <- 1.0
k1<-1.0
x <- as.matrix(a*u*cos(u)+k1* randn(1, N))
y <- as.matrix(a*u*sin(u))

#------------------------------------------------Tables----------------------------------------------------------

Table.1<-xtable(hurst.transformations.df)
#------------------------------------------------Figures--------------------------------------------------------
opar <- par(mfrow=c(2,2), mar=c(2,2,1,1))
#---------------------------Figure Group-----------------------------------
Figure.1<-plot(f1(z$x, z$y,1,1.25), type = 'l', col = "black", lwd = 2, lty=1,
     ylim = c(-1, 6), main = "Hilbert curve Transformations")
lines(f2(z$x, z$y,2,2), col = "red", lwd = 2, lty=2)
lines(f3(z$x, z$y,1.5,2.5), col = "green", lwd = 2, lty=3)
lines(f4(z$x, z$y,1,1), col = "blue", lwd = 2, lty=4)
grid()
legend("bottomright", col = c("black", "red", "green","blue"), 
       lty = 1:4, cex=0.75,
       legend = c("f1","f2", "f3", "f4"))
#---------------------------Figure 2------------------------------------
Figure.2<-plot(z.2$x, z.2$y, type='l', col="darkblue")
title("Hexagon Molecule Curve")
#---------------------------Figure 3-----------------------------------
Figure.3<-plot(z.3$x, z.3$y, type='l', col="darkblue")
title("Dragon Curve")
#---------------------------Figure 4-----------------------------------
Figure.4<-plot(matrix(f(ts), ncol=2), type='l', col="blue", 
     asp=1, xlab="", ylab = "",
     main = "Fermat's Spiral", sub="20 subparts of equal length")

for (i in seq(0.05, 0.95, by=0.05)) {
  v <- fParam(i*b); fv <- f(v)
  points(fv[1], f(v)[2], col="darkred", pch=20)
} 

#---------------------------Figure Group--------------------------------------------
opar <- par(mfrow=c(2,2), mar=c(2,2,1,1))
Figure.5<-plot(x, y, col = "green", pch = 19, asp = 1)
xp <- cfit1$xp; yp <- cfit1$yp
lines(xp, yp, col="blue")
grid()
Figure.6<-plot(x, y, type = "p", pch = 19, col = "green", asp = 1)
lines(x, y, col = "green", lwd = 3)
cfit <- curvefit(u, x, y, n)
px <- c(cfit$px); py <- c(cfit$py)
v <- linspace(0, 3*pi, 200)
xs <- polyval(px, v)
ys <- polyval(py, v)
lines(xs, ys, col = "navy")
Figure.7<-plot(f1(cfit1$xp, cfit1$yp,2,2))
grid()
Figure.8<-plot(f1(cfit1$xp, cfit1$yp,1,1))
grid()
#------------------------------------------------Function Library-----------------------------------------------
f1 <- function(x, y,p1,p2) x^p1 + y^p2
f2 <- function(x, y,p1,p2) x^p1 - y^p2
f3 <- function(x, y,p1,p2) x^p1 * y^p2
f4 <- function(x, y,p1,p2) x^p1 / y^p2

f <- function(t) c(sin(2*t), cos(t), t)

Fractal.Algorithm<-function (n,which = c("hilbert", "dragon", "molecule")) 
{
  curve <- match.arg(which)
  if (curve == "hilbert") {
    a <- 1 + (0+1i)
    b <- 1 - (0+1i)
    z <- 0
    for (k in 1:n) {
      w <- (0+1i) * Conj(z)
      z <- c(w - a, z - b, z + a, b - w)/2
    }
  }
  else if (curve == "dragon") {
    a <- (1 + (0+1i))/2
    b <- (1 - (0+1i))/2
    c <- sqrt(1/2 + (0+0i))
    z <- c(1 - c, c)
    for (k in 1:n) {
      w <- rev(z)
      z <- c(a * z, 1 - b * w)
    }
  }
  else if (curve == "molecule") {
    a <- (1 + sqrt(-3 + (0+0i)))/2
    b <- (1 - sqrt(-3 + (0+0i)))/2
    c <- c(1, a, -b, -1, -a, b)
    u <- 0
    for (k in 1:n) {
      u <- c(u + 1, -u, u - 1)
    }
    u <- c(u, 1 - u, 2 + u, 3 - u, 4 + u, 5 - u)
    u <- mod(u, 6)
    z <- cumsum(c[u + 1])
    z <- c(0, z/2^n)
  }
  return(list(x = Re(z), y = Im(z)))
  
  #----------------------------------------------Learning Space-------------------------------------------
  
  f <- function(xy) {
    x <- xy[1]; 
    y <- xy[2]
    3*(1-x)^2 * exp(-(x^2) - (y+1)^2) -
      10*(x/5 - x^3 - y^5) * exp(-x^2 - y^2) -
      1/3 * exp(-(x+1)^2 - y^2)
  }
  ezcontour(f, col = "navy")
  ezcontour(f, filled = TRUE)
  ezmesh(f)
  ezmesh(f, col="lightblue", theta = -15, phi = 30)
  output<-list()
  output$xy<-xy
  return(output)
}
