Ovenrun <- c(1,1,1,2,2,2,3,3,4,4,4); Ovenrun
x1 <- c(-1,-1,-1,1,1,1,0,0,0,0,0); x1
x2 <- c(-1,1,0,-1,1,0,-1,1,0,0,0); x2
y <- c(2.7,2.5,2.7,2.9,1.3,2.2,3.7,2.9,2.9,2.8,2.9); y
x1sq=x1^2
x2sq=x2^2
cake=cbind(Ovenrun,x1,x2,y,x1sq,x2sq)
cake=data.frame(cake)
library(rsm)
mod <- rsm(y ~ SO(x1, x2),  data=cake)
summary(mod)
ridge <- steepest(mod, dist=seq(0, 1.5, by=.1), descent=FALSE) # maximum
ridge
# Contour Plot
contour(mod, x1~x2)

summary(mod)

library(daewr)
data(cake)
library(lme4)
mmod <- lmer(y ~ x1 + x2 + x1:x2 + x1sq + x2sq + (1|Ovenrun), data = cake)
mmod
summary(mmod)


par (mfrow=c(1,1)) # temp x click contour plot of nonlinear analysis
x1 <- seq(-1, 1, length=30); x1
x2 <- seq(-1, 1, length=30); x2
y <- matrix(rep(0, 900), nrow=30); y
for (i in 1:30)    {
  for (j in 1:30) {
    y[i,j] <- 3.1312-0.25*x1[j]-0.433*x2[i]-0.684*x1[j]*x1[j]-0.096*x2[i]*x2[i]-0.35*x1[j]*x2[i] 
  }
}
library(graphics)
contour(x1, x2, y, xlab="time", ylab="temp")

