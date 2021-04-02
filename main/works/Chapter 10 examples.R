# Example 1 p. 389 Figure 10.3
library(rsm)
rotd <- ccd(3, n0 = c(4,2), alpha = "rotatable", randomize = FALSE)
rotdm <- rotd[ , 3:5]
library(Vdgraph)
Vdgraph(rotdm)

# Figure 10.5 p. 392
# create design with rsm
library(rsm)
bbd3 <- bbd(3,randomize=FALSE,n0=3)
library(Vdgraph)
Vdgraph(bbd3[ , 3:5])

# Figure 10.6 p. 392
library(Vdgraph)
FDSPlot(bbd3[ , 3:5], mod=2)

# Figure 10.8 p. 394
library(Vdgraph)
data(SCDH2)
Vdgraph(SCDH2)

# Figure 10.9 p. 395
# first get the ccd2
library(rsm)
ccd2<- ccd(2, n0 = c(3,2), alpha = "rotatable", randomize = FALSE)[ ,3:4]
Compare2FDS(SCDH2,ccd2,"SCD","CCD",mod=2)

# Figure 10.10 p. 396
library(Vdgraph)
data(D310)
Vdgraph(D310)

# Example 2 p. 396-397
library(rsm)
ccd.pick(k=3)

# Example 3 p. 397
library(rsm)
ccd.up<-ccd(y~x1+x2+x3,n0=c(4,2),alpha="rotatable",randomize=FALSE)
Vdgraph(ccd.up[ , 3:5])

# Example 4 p. 398
ccd.up<-ccd(y~x1+x2+x3,n0=c(4,2),alpha="rotatable",
            coding=list(x1~(Temp-150)/10,x2~(Press-50)/5, x3~(Rate-4/1)),randomize=FALSE)
head(ccd.up)

# Example 5 p. 398
library(rsm)
Treb<-bbd(y ~ x1 + x2 + x3, randomize = FALSE, n0 = 3,
          coding=list(x1~(A-6)/2, x2~(B-15)/5, x3~(C-2.5)/.5))
head(Treb)

# Example 6 p. 399
library(Vdgraph)
data(D310)
D310

# Example 7 p. 400
des <- transform(D310, Temp=10*x1+150, Press=5*x2+50, Rate = x3 + 4)
des[sample(1:10) ,4:6]

# Example 8 p. 401
library(daewr)
data(qsar)

# Example 9 p.403 and Figure 10.13 p. 404
library(AlgDesign)
desgn1<-optFederov(~quad(.),data=qsar,nTrials=15,center=TRUE, 
                   criterion="D",nRepeats=40)
desgn2<-optFederov(~quad(.),data=qsar,nTrials=15,center=TRUE, 
                   criterion="I",nRepeats=40)
Compare2FDS(desgn1$design, desgn2$design, "D-optimal", "I-optimal", mod=2)

# Example 10 p. 404
desgn1$design

# Example 11 p. 406
k1 <- .15; k2 <- .72; gamma0 <- 2.65; t0 <- 0.41
x <- c(seq(1:25))
dfdk1 <- c(rep(0, 25))
dfdk2 <- c(rep(0, 25))
dfdgamma0 <- c(rep(0, 25))
dfdt0 <- c(rep(0, 25))
for (i in 1:25) {
  dfdk1[i] <- -1 * gamma0 * exp(-k1 * (x[i] - t0)) *(x[i] - t0)
  dfdk2[i] <-gamma0 * exp(-k2 * (x[i] - t0)) * (x[i] - t0)
  dfdgamma0[i] <- exp(-k1 * (x[i] - t0)) - exp( -k2 *                                                  ( x[i] - t0))
  dfdt0[i] <- gamma0 * exp(-k1 * (x[i] - t0)) * k1 - gamma0 * 
    exp(-k2 * (x[i] - t0)) * k2; }
grid <- data.frame(x, dfdk1, dfdk2, dfdgamma0, dfdt0); grid
library(matrixcalc)
X <- as.matrix(grid); I <- t(X) %*% X; I
is.singular.matrix(I, tol = 1e-08)
det(I)

# Example 12 p. 407
library(AlgDesign)
desgn2<-optFederov(~ +dfdk1+dfdk2+dfdgamma0+dfdt0,data=grid,nTrials=5,center=TRUE, criterion="D",nRepeats=20)
desgn2

# Example 13 p. 408
library(daewr)
data(cement)
head(cement)

# Example 14 p. 408
library(rsm)
grout.lin <- rsm(y ~ SO(x1, x2, x3),data = cement, subset = (Block == 1))
anova(grout.lin)

# Example 15 p. 409
library(daewr)
data(Treb)
head(Treb)

# Example 16 p. 410
library(rsm)
treb.quad <- rsm(y ~ SO(x1, x2, x3), data = Treb)
summary(treb.quad)

# Example 17 p. 411 
library(rsm)
library(daewr)
data(cement)
grout.quad <- rsm(y ~ Block + SO(x1,x2,x3), data = cement)
summary(grout.quad)

# Example 18 p. 412 - 413
library(daewr)
data(Tet)
mod.nln1 <-nls(Conc ~ gamma0 * (exp( -k1 * (Time - t0)) - exp( -k2 * (Time - t0))), 
               data = Tet, 
               start = list(gamma0 = 10, k1 = .12, k2 = .5, t0 = .5),
              trace = TRUE); mod.nln1
summary(mod.nln1)

# Example 20 p. 413-414 Figure 10.15
library(daewr)
data(Treb)
library(rsm)
treb.quad <- rsm(y ~ SO(x1, x2, x3), data = Treb)
par (mfrow=c(2,2))
contour(treb.quad, ~ x1+x2+x3 )

# Example 21 p. 414-415 Figure 10.16
par (mfrow=c(2,2))
persp(treb.quad, ~ x1+x2+x3, zlab="Distance", contours=list(z="bottom") )

# Example 22 p. 415-416, Figure 10.17
par (mfrow=c(1,2))
contour(treb.quad, x1~x3, at=list(x2=1))
persp(treb.quad, x1~x3, at=list(x2=1),zlab="Distance", contours=list(z="bottom"))

# Example 23 p. 418
# slice at stationary point instead of coding value
par (mfrow=c(2,2))
contour(treb.quad, ~ x1+x2+x3, at = xs(treb.quad) )

# Example 24 p. 419
ridge<-steepest(treb.quad, dist=seq(0, 1.412, by=.1), descent=FALSE)
ridge

# Figure 10.20 p. 420
par (mfrow=c(2,1))
leg.txt<-c("A","B","C")
plot(ridge$dist,ridge$yhat, type="l",xlab="radius",ylab="Max. Predicted")
plot(ridge$dist,seq(1,22,by=1.5), type="n", xlab="radius", ylab="Factors")
lines(ridge$dist,ridge$A,lty=1)
lines(ridge$dist,ridge$B,lty=2)
lines(ridge$dist,ridge$C,lty=3)
legend(1.2,19,leg.txt,lty=c(1,2,3))

#Example 25 p. 421
start<-c(20, 80, 0.5); start
prod <- function(x) {
  temp <- x[1]
  click <- x[2]
  hard <- x[3]
  f <- exp( -temp/15 * (click - 60)) + exp( -k * (click - 60))
}
ui <- matrix(c(1, -1, 0, 0, 0, 0, 0, 0, 1, -1, 0, 0, 0, 0, 0, 0, 1, -1), 6, 3); ui
ci <- c(-10, -50, 40, -120, 0.1, -0.9); ci
constrOptim(start, prod, NULL, ui, ci)


# Figure 10.21 p. 422
par (mfrow=c(1,1))
time<-seq(0,25,length=52)
Temp<-seq(375,425,length=52)
p<-matrix(rep(0,2704),nrow=52)
for (i in 1:52)    {
  for (j in 1:52) {
    k1<-.523*exp(-9847*((1/Temp[j]-1/400)))
    k2<-.2*exp(-12327*((1/Temp[j]-1/400)))
    p[i,j]<--132*(exp(-k1*time[i])-exp(-k2*time[i]))*k1/(k1-k2)
  }
}
contour(time,Temp,p,xlab="Time, hrs",ylab="Temp, Deg K")

# Example 26 p. 424
# maximize tsr for cost <=10 cents
# Note: since no R packages allow nonlinear equality constraints
# I included this constraint in the tsr function by adding the 
# 100*maxumum of (cost -10) and zero 
start<-c(100,.6,40)
tsrcost<-function(x) {
  Temp<-x[1]
  Ratio<-x[2]
  AOPPM<-x[3]
  tsrcost<- -(-226+3.375*Temp+86.5*Ratio+2.646*AOPPM-.0128*Temp*Temp-17.5*Ratio*Ratio-.0121*AOPPM*AOPPM
              +           -.3857*Ratio*Temp-.0126*AOPPM*Temp-.0333*AOPPM*Ratio)+
    +           100*max((.8313+1.27*Ratio+.37*Ratio*AOPPM-10),0)
}
ui<-matrix(c(1,-1,0,0,0,0,0,0,1,-1,0,0,0,0,0,0,1,-1),6,3)
ci<-c(70,-140,.5,-1.5,5,-65)
constrOptim(start,tsrcost,NULL,ui,ci)
# note this minimizes tsrcost, subject to ui%*%x-ci>=0

# Example 27 p. 426
x<-c(0,0,0)
saPred<-function(x) 125.4106 -8.1233*x[1]+17.0266*x[2]+.4277*x[3]+2.4184*x[1]*x[2]
-8.4376*x[1]*x[3]+9.0134*x[2]*x[3]+33.88054*x[1]^2+14.81976*x[2]^2+13.07001*x[3]^2

pvPred<-function(x) .661354-.1963*x[1]-.02016*x[2]-.00291*x[3]+.02399*x[1]*x[2]
+.010327*x[1]*x[3]-.0374*x[2]*x[3]+.15126*x[1]^2+.118423*x[2]^2+.0679*x[3]^2

dpPred<-function(x) 39.35608+3.19547*x[1]+.21729*x[2]-1.46979*x[3]+.58873*x[1]*x[2]
-.62136*x[1]*x[3]-1.53234*x[2]*x[3]+.41413*x[1]^2-2.39408*x[2]^2-2.36399*x[3]^2

library(desirability)
saD<-dMax(100, 217)
pvD<-dMax(0.6, 1.3)
dpD<-dTarget(38, 40, 42)
overallD<-dOverall(saD, pvD, dpD)

# Code on web page referred to on p. 426
rsmOpt <- function(x, dObject, space = "square")
{
  sa <- saPred(x)
  pv <- pvPred(x)
  dp <- dpPred(x)
  
  out <- predict(dObject, data.frame(sa = sa, pv = pv, dp = dp))
  
  if(space == "circular")
  {
    if(sqrt(sum(x^2)) > 1.0) out <- 0
  } else if(space == "square") if(any(abs(x) > 1.0)) out <- 0
  out
}


searchGrid <- expand.grid(Mixtime = seq(-1.0, 1.0, length = 5),
                          filTtime = seq(-1.0, 1.0, length = 5),
                          PackMth = seq(-1.0, 1.0, length = 5))
for(i in 1:dim(searchGrid)[1])
{
  tmp <- optim(as.vector(searchGrid[i,]),
               rsmOpt,
               dObject = overallD,
               space = "square",
               control = list(fnscale = -1))
  if(i == 1)
  {
    best <- tmp
  } else {
    if(tmp$value > best$value) best <- tmp
  }
}

predOutcomes <-c(saPred(c(-0.2924993,-1,-1)), pvPred(c(-0.2924993,-1,-1)), dpPred(c(-0.2924993,-1,-1)))
print(predOutcomes)

# Example 28  p. 427
#Create Blocked bbd
library(rsm)
# four-factor design
bbd(4,n0=1,randomize=FALSE)
# five-factor design
bbd(5,n0=3,randomize=FALSE)
#Create Blocked CCD designs
# two-factor design
ccd(2,n0=c(3,3),alpha="orthogonal",randomize=FALSE)
# three-factor design
ccd(3, n0=2,alpha="orthogonal",randomize=FALSE, blocks=Block~( x1*x2*x3))
# four-factor design
ccd(4, n0=2, alpha="orthogonal",randomize=FALSE, blocks=Block~(x1*x2*x3*x4))
# five-factor design
ccd(4, generators=(x5~x1*x2*x3*x4), n0=c(6,1),alpha="orthogonal",randomize=FALSE)

# Example 29 p. 428
library(AlgDesign)
fact<-gen.factorial(levels=2,nVars=3)
fact<-rbind(fact,fact)
center<-data.frame(matrix(rep(c(0,0,0),6),ncol=3))
star<-data.frame(rbind(diag(3),-diag(3)))
cand<-rbind(fact,center,star)
bdesign<-optBlock(~quad(.),cand,blocksizes=c(4,4,4,4,4,4,4),criterion="Dpc",nRepeats=1000)
bdesign

# Example 30 p. 429
library(daewr)
data(qsar) 
library(AlgDesign)
desgn1<-optBlock(~quad(.),qsar,blocksizes=c(4,4,4,4,4),criterion="Dpc",nRepeats=1000)
desgn1

# Example 31 p. 430
library(daewr)
data(pastry); pastry
class(pastry$Block)
library(rsm)
blkrsm<-rsm(y~Block+SO(x1,x2,x3), data=pastry)
summary(blkrsm)

#Example 32 p. 432
# REML analysis using lmer
library(daewr)
data(cake)
library(lme4)
mmod<-lmer(y~x1+x2+x1:x2+x1sq+x2sq+(1|Ovenrun),data=cake)
summary(mmod)
# Least Squares Analysis using rsm
library(rsm)
mmodls<-rsm(y ~ SO(x1,x2), data=cake)
summary(mmodls)

# Example 33 EESPRS page 434
library("AlgDesign")
# uses gen.factorial function from AlgDesign to create Factorial portion of the design
sp<-gen.factorial(2,2,varNames=c("P","Q")); sp
Wp<-gen.factorial(2,2,varNames=c("A","B")); sp
A<-Wp$A; A
#stacks whole plots
wp<-c(rep((1:4),each=4)); wp
A<-c(rep((Wp$A),each=4)); A
B<-c(rep((Wp$B),each=4)); B
Fac<-cbind(wp,A,B,rbind(sp,sp,sp,sp)); Fac
# Subplot Axial Portion 
A<-c(rep(0,4))
B<-c(rep(0,4))
P<-c(-2,2,0,0)
Q<-c(0,0,-2,2)
wp<-c(rep(5,4))
spa<-cbind(wp,A,B,P,Q)
# Whole Plot Axial Portion
wp<-c(rep((6:9),each=4))
A<-c(rep(c(-2,2,0,0),each=4))
B<-c(rep(c(0,0,-2,2),each=4))
P<-c(rep(0,16))
Q<-c(rep(0,16))
wpa<-cbind(wp,A,B,P,Q)
# center points
wp<-c(rep((10:11),each=4)); wp
A<-c(rep(0,8)); A
B<-c(rep(0,8)); B
P<-c(rep(0,8)); P
Q<-c(rep(0,8)); Q
wpc<-cbind(wp,A,B,P,Q); wpc
SPDs<-rbind(Fac,spa,wpa,wpc); SPDs
############################################

# Example 34 p. 438
library(daewr)
EEw2s3( )
EEw2s3('EE21R7WP')
