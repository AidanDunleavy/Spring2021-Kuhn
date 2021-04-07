# Lecture Notes 11
# Chapter 10.7-10.9 Response Surface Designs

# 10.7 Determining Optimum Operating Conditions

# 10.7.1 visual analysis
# CCD design: 3 factor, 19 run
#   2^3 = 8 factorial, 3 in-center-factorial, 
#   6 axial, 2 in-center-axial
library(rsm) 
ccd.act <- ccd(defects~x1+x2+x3,
               n0=c(3,2), 
               alpha="rotatable", 
               coding=list(x1~(temp-15)/15, 
                           x2~(click-90)/30, 
                           x3~(hard-0.65)/0.15), 
               randomize=FALSE)
ccd.act
act <- as.coded.data(ccd.act); act # save as coded data
act$defects <- c(61, 69, 60, 45, 64, 66, 44, 57,
                 60, 61, 52, 41, 56, 61, 51, 63,
                 60, 60, 52); act # add defects response to coded data

# use of 2D contour plots and 3D plots to identify optimal conditions
library(rsm)
sp.quad <- rsm(defects ~ SO(x1, x2, x3), data = act) # fit ccd model to data
par (mfrow=c(2,2)) # various 2D contour plots for response # defects
contour(sp.quad, ~ x1+x2+x3 )

par (mfrow=c(2,2)) # various 3D contour plots for response # defects
persp(sp.quad, ~ x1+x2+x3, zlab="defects", contours=list(z="bottom") )

par (mfrow=c(1,2)) # 2D and 3D plots at x1 = -1 (low temp 0^o)
contour(sp.quad, x2~x3, at=list(x1=-1))
persp(sp.quad, x2~x3, at=list(x1=-1), zlab="defects", 
      contours=list(z="bottom"))

par (mfrow=c(1,2)) # 2D and 3D plots at x2 = 1 (high click 120)
contour(sp.quad, x1~x3, at=list(x2=1))
persp(sp.quad, x1~x3, at=list(x2=1), zlab="defects", 
      contours=list(z="bottom"))

par (mfrow=c(1,2)) # 2D and 3D plots at x3 = 1 (high hard 0.8)
contour(sp.quad, x1~x2, at=list(x3=1))
persp(sp.quad, x1~x2, at=list(x3=1), zlab="defects", 
      contours=list(z="bottom"))

# 10.7.2 canonical analysis
# CCD design: 3 factor, 19 run
#   2^3 = 8 factorial, 3 in-center-factorial, 
#   6 axial, 2 in-center-axial
library(rsm) 
ccd.act <- ccd(defects~x1+x2+x3,
               n0=c(3,2), 
               alpha="rotatable", 
               coding=list(x1~(temp-15)/15, 
                           x2~(click-90)/30, 
                           x3~(hard-0.65)/0.15), 
               randomize=FALSE)
ccd.act
act <- as.coded.data(ccd.act); act # save as coded data
act$defects <- c(61, 69, 60, 45, 64, 66, 44, 57,
                 60, 61, 52, 41, 56, 61, 51, 63,
                 60, 60, 52); act # add defects response to coded data
sp.quad <- rsm(defects ~ Block + SO(x1, x2, x3), data = act) # fit ccd model to data
summary(sp.quad) # determine stationary point

par (mfrow=c(2,2))
contour(sp.quad, ~ x1+x2+x3, at = xs(sp.quad) )

# 10.7.3 ridge analysis
#  using radius distance between center and factorial corner
ridge <- steepest(sp.quad, dist=seq(0, 1.73, by=.1), descent=FALSE) # maximum
ridge
# graphs of ridge analysis 
par (mfrow=c(1,2))
leg.txt<-c("temp","click","hard")
plot(ridge$dist,ridge$yhat, type="l",xlab="radius",ylab="Max. Predicted")
plot(ridge$dist,seq(0,85,by=5), type="n", xlab="radius", ylab="Factors")
lines(ridge$dist,ridge$temp,lty=1)
lines(ridge$dist,ridge$click,lty=2)
lines(ridge$dist,ridge$hard,lty=3)
legend(0.1,60,leg.txt,lty=c(1,2,3))
par (mfrow=c(1,1))

ridge <-steepest(sp.quad, dist=seq(0, 1.73, by=.1), descent=TRUE) # minimum
ridge
# graphs of ridge analysis 
par (mfrow=c(1,2))
leg.txt<-c("temp","click","hard")
plot(ridge$dist,ridge$yhat, type="l",xlab="radius",ylab="Min. Predicted")
plot(ridge$dist,seq(0,102,by=6), type="n", xlab="radius", ylab="Factors")
lines(ridge$dist,ridge$temp,lty=1)
lines(ridge$dist,ridge$click,lty=2)
lines(ridge$dist,ridge$hard,lty=3)
legend(0.2,85,leg.txt,lty=c(1,2,3))
par (mfrow=c(1,1))

# 10.7.4 nonlinear analysis
start <- c(20, 50, 0.5); start
prod <- function(x) {
  temp <- x[1]
  click <- x[2]
  hard <- x[3]
  f <- exp( -temp/15 * (click - 60)) + exp( -hard * (click - 60))
}
ui <- matrix(c(1, -1, 0, 0, 0, 0, 0, 0, 1, -1, 0, 0, 0, 0, 0, 0, 1, -1), 6, 3); ui
ci <- c(-10, -50, 40, -120, 0.1, -0.9); ci
constrOptim(start, prod, NULL, ui, ci)

par (mfrow=c(1,1)) # temp x click contour plot of nonlinear analysis
temp <- seq(12.2, 12.4, length=20); temp
click <- seq(80.0, 80.2, length=20); click
defect <- matrix(rep(0, 400), nrow=20); defect
for (i in 1:20)    {
  for (j in 1:20) {
    defect[i,j] <- exp( -temp[i]/15 * (click[j] - 60)) + exp( -0.84 * (click[j] - 60))
  }
}
library(graphics)
contour(temp, click, defect, xlab="temperature, degrees", ylab="clicks, number")

par (mfrow=c(1,1)) # temp x hard contour plot of nonlinear analysis
temp <- seq(12.2, 12.4, length=20); temp
hard <- seq(0.83, 0.84, length=20); hard
defect <- matrix(rep(0, 400), nrow=20); defect
for (i in 1:20)    {
  for (j in 1:20) {
    defect[i,j] <- exp( -temp[i]/15 * (80.1 - 60)) + exp( hard[j] * (80.1 - 60))
  }
}
library(graphics)
contour(temp, hard, defect, xlab="temperature, degrees", ylab="hardness")

par (mfrow=c(1,1)) # click x hard contour plot of nonlinear analysis
click <- seq(80.0, 80.2, length=20); click
hard <- seq(0.83, 0.85, length=20); hard
defect <- matrix(rep(0, 400), nrow=20); defect
for (i in 1:20)    {
  for (j in 1:20) {
    defect[i,j] <- exp( -12.3/15 * (click[i] - 60)) + exp( hard[j] * (click[i] - 60))
  }
}
library(graphics)
contour(click, hard, defect, xlab="clicks, number", ylab="hardness")


# 10.7.5 multi-response optimization (constrOptm)
#   two responses: number of defects and cost
#   cost function is nonlnear, cannot be included in linear ux - c > 0 constraints 
start <- c(20, 80, 0.5); start
prod <- function(x) {
  temp <- x[1]
  click <- x[2]
  hard <- x[3]
  ndefect <- exp( -temp/15 * (click - 60)) + exp( -hard * (click - 60)) # number defects
  cost <- max((8 + 3 * click + 6 * temp + 0.3 * temp*hard), 0) # cost > 0, cannot be negative 
  f <- ndefect + cost # minimize number of defects and cost
}
ui <- matrix(c(1, -1, 0, 0, 0, 0, 0, 0, 1, -1, 0, 0, 0, 0, 0, 0, 1, -1), 6, 3); ui
ci <- c(-10, -50, 40, -120, 0.1, -0.9); ci
constrOptim(start, prod, NULL, ui, ci)

# 10.8 Blocked Response Surface (BRS) Designs

# 10.8.1 orthogonal blocked three-factor CCD design
library(rsm)
des <- ccd(3, n0=2, alpha="orthogonal", randomize=FALSE, blocks=Block~( x1*x2*x3))
des 
X <- as.matrix(des[,-c(1,2,6)]); X # create matrix of ccd design
t(X)%*%X # matrix is diagonal showing overall design is orthogonal
XB1 <- X[1:6,]; XB1; t(XB1)%*%XB1 # each block matrix is also diagonal
XB2 <- X[7:12,]; XB2; t(XB2)%*%XB2
XB3 <- X[13:20,]; XB3; t(XB3)%*%XB3

library(rsm)
des <- ccd(3, n0=2, alpha="orthogonal", randomize=FALSE, blocks=Block~( x1*x2*x3))
des 
cand <- as.matrix(des[,-c(1,2,6)]); dim(cand); cand

# 10.8.2 blocked center-faced design which is D-s optimal
library(AlgDesign)
fact <- gen.factorial(levels=2, nVars=3); fact # 2^3 = 8 factorial design
fact <- rbind(fact, fact); fact
center <- data.frame(matrix(rep(c(0,0,0),6), ncol=3)); center # 6 center points
star <- data.frame(rbind(diag(3), -diag(3))); star # face-centered axial points
cand <- rbind(fact, center, star); cand
bdesign <- optBlock( ~ quad(.), # general quadratic model
                  cand,
                  blocksizes = c(7,7,7,7), # 4 blocks each of size 7
                  criterion = "Dp", # D_s design
                  nRepeats = 1000) # run numerical algorithm 1000 times
bdesign
des <- as.data.frame(bdesign$design); des; dim(des)

# 10.8.3 analyze blocked center-faced design 
library(rsm)
block <- c(rep(1, 7), rep(2, 7), rep(3, 7), rep(4, 7)); block; 
block <- as.data.frame(block); dim(block); class(block) # defects data
defects <- c(61, 69, 60, 55, 64, 66, 44, 
             57, 60, 61, 52, 51, 56, 61, 
             51, 63, 60, 60, 52, 56, 61,
             51, 53, 60, 64, 62, 52, 56); defects; 
defects <- as.data.frame(defects); dim(defects); class(defects) # defects data
d.defects <- cbind(block, des, defects); d.defects
a.defects <- rsm(defects ~ block + SO(x1,x2,x3), data=d.defects)
summary(a.defects)

# 10.9 Response Surface Split-Plot (RSSP) Designs

# 10.9.1 reorganize standard CCD to split-plot design
#   2^3 = 8 factorial, 3 in-center-factorial, 
#   6 axial, 2 in-center-axial
library(rsm) 
ccd.des <- ccd(defects~x1+x2+x3, n0=c(3,2), alpha="rotatable", randomize=FALSE)
ccd.des
des <- as.data.frame(ccd.des); des # convert to data.frame
des$defects <- c(61, 69, 60, 45, 64, 66, 44, 57,
                 60, 61, 52, 41, 56, 61, 51, 63,
                 60, 60, 52); des # add defects response to data.frame
des.ccd <- as.data.frame(des[,-c(1,2,7)]); des.ccd; dim(des.ccd) # edit to CCD
des.sp <- des.ccd[order(des.ccd$x1),]; des.sp # sort on x1 into split-plot design

# Least Squares Analysis using rsm
library(rsm)
mmodls<-rsm(defects ~ SO(x1,x2,x3), data=des.sp)
summary(mmodls)
# REML analysis using lmer
des.sp$x1sq <- des.sp$x1^2; # add squared terms to data
des.sp$x2sq <- des.sp$x2^2; 
des.sp$x3sq <- des.sp$x3^2; des.sp
des.sp$wholeplot <- c(1, rep(2,4), rep(3,9), rep(4,4), 5); des.sp # add whole plots to design
library(lme4)
mmod<-lmer(defects ~ x1 + x2 + x3 
           + x1:x2 + x1:x3 + x2:x3
           + x1sq + x2sq + x3sq + (1|wholeplot),
           data = des.sp)
summary(mmod)

# 10.9.2 EESRPS CCD design
#   k1 = 1 whole plot (WP) factor, with m = 7 whole plots 
#   k2 = 2 sub plot (SP) factors, each with 2^2 = 4 subplots 
#   2 WPs for W axials, 1 WP for SP axials, 2 factorial WPs, 2 center WPs 
library("AlgDesign")
# uses gen.factorial function from AlgDesign to create Factorial portion of the design
wp <- c(rep(c(2,6),each=4)); wp
A <- c(-1,-1,-1,-1,1,1,1,1)
sp <- gen.factorial(2,2,varNames=c("P","Q")); sp
Fac <- cbind(wp,A,rbind(sp,sp)); Fac
# subplot axial portion 
wp <- c(rep(4,4))
A <- c(rep(0,4))
P <- c(-1.68,1.68,0,0)
Q <- c(0,0,-1.68,1.68)
spa <- cbind(wp,A,P,Q); spa
# whole plot axial portion
wp<-c(rep(c(1,7),each=4))
A<-c(rep(c(-1.68,1.68),each=4))
P<-c(rep(0,8))
Q<-c(rep(0,8))
wpa<-cbind(wp,A,P,Q); wpa
# center points
wp <- c(rep(c(3,5),each=4)); wp
A <- c(rep(0,8)); A
P <- c(rep(0,8)); P
Q <- c(rep(0,8)); Q
wpc<-cbind(wp,A,P,Q); wpc
SPDs<-rbind(Fac,spa,wpa,wpc); SPDs
des <- SPDs[order(SPDs$wp),]; des # sort on A into EESRPS CCD design
des$defects <- c(56, 54, 61, 61, 60, 52, 64, 
                 66, 65, 43, 47, 50, 41, 61, 
                 60, 64, 44, 61, 51, 63, 60,
                 60, 61, 52, 60, 52, 69, 45); des # add defects to data.frame

# Least Squares Analysis using rsm
library(rsm)
mmodls<-rsm(defects ~ SO(A,P,Q), data=des)
summary(mmodls)
# REML analysis using lmer
des$Asq <- des$A^2; # add squared terms to data
des$Psq <- des$P^2; 
des$Qsq <- des$Q^2; des
library(lme4)
mmod<-lmer(defects ~ A + P + Q 
           + A:P + A:Q + P:Q
           + Asq + Psq + Qsq + (1|wp),
           data = des)
summary(mmod)

# 10.9.3 D-efficient EESPRS design
#   8 factorial, 4 axial, 2 center points = 14 run
library(daewr)
des <- EEw1s2('EE14R7WP')  
names(des)[names(des) == "w1"] <- "A"
names(des)[names(des) == "s1"] <- "P"
names(des)[names(des) == "s2"] <- "Q"
des$defects <- c(41, 61, 60, 64, 44, 61, 51, 63, 60,
                 60, 61, 52, 60, 52); des # add defects to data.frame
# Least Squares Analysis using rsm
library(rsm)
mmodls<-rsm(defects ~ SO(A,P,Q), data=des)
summary(mmodls)
# REML analysis using lmer
des$Asq <- des$A^2; # add squared terms to data
des$Psq <- des$P^2; 
des$Qsq <- des$Q^2; des
library(lme4)
mmod<-lmer(defects ~ A + P + Q 
           + A:P + A:Q + P:Q
           + Asq + Psq + Qsq + (1|WP),
           data = des)
summary(mmod)

###################################################
# Homework 11 and Test 6 questions
###################################################

# 10.8.a.1-2 pp 441-446
# cement example table 10.1 p 388 text
library(daewr) # retrieve cement data and design from daewr
data(cement)

# 10.8.a.3 pp 441-446
# cement example table 10.1 p 388 text
library(daewr) # retrieve cement data and design from daewr
data(cement)
# anova on only block 1 of data
cement.m <- rsm(y ~ SO(x1, x2, x3), data = cement, subset = (Block == 1))
anova(cement.m)

# 10.8.a.4 pp 441-446
library(rsm) # fit general quadratic model using rsm function
cement.m <- rsm(y ~ SO(x1, x2, x3), data = cement)
summary(cement.m)

# 10.8.b.1-2 pp 441-446
library(rsm) # canonical analysis using rsm
cement.m <- rsm(y ~ SO(x1, x2, x3), data = cement)
summary(cement.m)

# 10.8.b.3 pp 441-446
library(rsm) # ridge analysis using steepest in rsm
ridge <- steepest(cement.m, dist=seq(0, 1.73, by=.1), descent=FALSE) # maximum
ridge
# graphs of ridge analysis 
par (mfrow=c(1,2))
leg.txt<-c("WatCem","BlackL","SNF")
plot(ridge$dist,ridge$yhat, type="l",xlab="radius",ylab="Max. Predicted")
plot(ridge$dist,seq(0.04,0.38,by=0.02), type="n", xlab="radius", ylab="Factors")
lines(ridge$dist,ridge$WatCem,lty=1)
lines(ridge$dist,ridge$BlackL,lty=2)
lines(ridge$dist,ridge$SNF,lty=3)
legend(0.05,0.30,leg.txt,lty=c(1,2,3))
par (mfrow=c(1,1))

# 10.8.b.3 pp 441-446
library(rsm) # ridge analysis using steepest in rsm
ridge <- steepest(cement.m, dist=seq(0, 1.73, by=.1), descent=TRUE) # minimum
ridge
# graphs of ridge analysis 
par (mfrow=c(1,2))
leg.txt<-c("WatCem","BlackL","SNF")
plot(ridge$dist,ridge$yhat, type="l",xlab="radius",ylab="Max. Predicted")
plot(ridge$dist,seq(0.04,0.38,by=0.02), type="n", xlab="radius", ylab="Factors")
lines(ridge$dist,ridge$WatCem,lty=1)
lines(ridge$dist,ridge$BlackL,lty=2)
lines(ridge$dist,ridge$SNF,lty=3)
legend(0.05,0.30,leg.txt,lty=c(1,2,3))
par (mfrow=c(1,1))

# 10.10.a.1-2 pp 441-446
# orthogonal blocked CCD with 4 factors
library(rsm)
des <- ccd(4, n0=2, alpha="orthogonal", randomize=FALSE, blocks=Block~( x1*x2*x3*x4))
des 

# 10.10.b.1-2 pp 441-446
# orthogonal blocked BBD with 4 factors
library(rsm)
des <- bbd(y~ x1+x2+x3+x4, randomize=FALSE, n0=1)
des 

# 10.10.c pp 441-446 compare two blocked CCD designs
library(rsm)
des <- ccd(4, n0=2, alpha="orthogonal", randomize=FALSE, blocks=Block~( x1*x2*x3*x4))
des 
library(AlgDesign)
des.ccd <- as.data.frame(des[,-c(1,2,7)]); des.ccd; dim(des.ccd) # edit to data.frame
des.2 <- optBlock( ~ quad(.), # general quadratic model
                     des.ccd,
                     blocksizes = c(6,6,6,6,6), # 5 blocks each of size 6
                     criterion = "Dp", # D_s design
                     nRepeats = 1000) # run numerical algorithm 1000 times
des.2
library(Vdgraph) # variance dispersion graph
par(mfrow = c(1,2))
Vdgraph(des[ , 3:6]) 
Vdgraph(des.2$design[ , 1:4])
par(mfrow = c(1,1))

# 10.10.d pp 441-446 compare two blocked BBD designs
library(rsm)
des <- bbd(y~ x1+x2+x3+x4, randomize=FALSE, n0=1)
des 
library(AlgDesign)
des.bbd <- as.data.frame(des[,-c(1,2,3,8)]); des.bbd; dim(des.bbd) # edit to data.frame
center <- c(rep(0,4)); center
des.bbd.c <- rbind(des.bbd,center,center,center); des.bbd.c
des.2 <- optBlock( ~ quad(.), # general quadratic model
                   des.bbd.c,
                   blocksizes = c(6,6,6,6,6), # 5 blocks each of size 6
                   criterion = "Dp", # D_s design
                   nRepeats = 1000) # run numerical algorithm 1000 times
des.2
library(Vdgraph) # variance dispersion graph
par(mfrow = c(1,2))
Vdgraph(des[ , 4:7]) # graph 4 factors from bbd design 1
Vdgraph(des.2$design[ , 1:4]) # graph 4 factors from bbd + 3 center design 2
par(mfrow = c(1,1))

# 10.12.a pp 441-446 
# average daily gain (ADG) with factors x1 and x2 for cattle experiment
x1 <- c(6.723, 6.281, 5.788, 5.729, 6.411, 6.181,
        6.723, 6.281, 5.788, 5.729, 6.411, 6.181,
        6.723, 6.281, 5.788, 5.729, 6.411, 6.181,
        6.723, 6.281, 5.788, 5.729, 6.411, 6.181,
        6.723, 6.281, 5.788, 5.729, 6.411, 6.181,
        6.723, 6.281, 5.788, 5.729, 6.411, 6.181,
        6.723, 6.281,        5.729, 6.411, 6.181,
        6.723, 6.281,        5.729,        6.181)
x2 <- c(0.1095, 0.1273, 0.1157, 0.0669, 0.0729, 0.1004,
        0.1095, 0.1273, 0.1157, 0.0669, 0.0729, 0.1004,
        0.1095, 0.1273, 0.1157, 0.0669, 0.0729, 0.1004,
        0.1095, 0.1273, 0.1157, 0.0669, 0.0729, 0.1004,
        0.1095, 0.1273, 0.1157, 0.0669, 0.0729, 0.1004,
        0.1095, 0.1273, 0.1157, 0.0669, 0.0729, 0.1004,
        0.1095, 0.1273,         0.0669, 0.0729, 0.1004,
        0.1095, 0.1273,         0.0669,         0.1004)
block <-c(1, 1, 1, 1, 1, 1,
          2, 2, 2, 2, 2, 2,
          3, 3, 3, 3, 3, 3,
          4, 4, 4, 4, 4, 4,
          5, 5, 5, 5, 5, 5,
          6, 6, 6, 6, 6, 6,
          7, 7,    7, 7, 7,
          8, 8,    8,    8)
ADG <- c(2.45, 2.08, 2.97, 1.92, 4.17, 3.86,
          1.63, 2.86, 2.58, 2.63, 3.09, 2.39,
          1.28, 3.37, 2.51, 2.08, 1.37, 2.93,
          1.97, 2.37, 2.51, 2.08, 1.37, 2.93,
          1.80, 2.59, 2.57, 2.64, 2.38, 2.59,
          2.36, 2.16, 2.02, 2.17, 2.40, 3.23,
          1.55, 4.14,       1.35, 1.80, 2.19,
          1.89, 2.12,       4.21,       3.04)
cattle <- as.data.frame(cbind(x1,x2,block,ADG)); cattle
# compare first order (FO) with second order (SO) model to data, include block
library(rsm)
cattle.m <- rsm(ADG ~ block + FO(x1, x2), data = cattle) 
cattle.m2 = rsm(ADG ~ block + SO(x1, x2), data = cattle)
anova(cattle.m, cattle.m2) # investigate lack of fit
summary(cattle.m2) 

# 10.12.b pp 441-446 
# maximum average daily gain (ADG) for cattle experiment
cattle <- as.data.frame(cbind(x1,x2,block,ADG)); cattle
library(rsm)
cattle.m2 = rsm(ADG ~ block + SO(x1, x2), data = cattle)
summary(cattle.m2) # for eigenvalues and optimal point
ridge <- steepest(cattle.m2, dist=seq(0, 6.2, by=.1), descent=FALSE) # maximum point
ridge

# 10.12.c pp 441-446 
# check rotatability (ADG)
cattle <- as.data.frame(cbind(x1,x2,block,ADG)); cattle
Vdgraph(cattle[ , 1:2]) # graph 2 factors from irregular design 

# 10.12.d pp 441-446 
# identify maximum ADG using contour plot
cattle <- as.data.frame(cbind(x1,x2,block,ADG)); cattle
library(rsm)
cattle.m2 = rsm(ADG ~ block + SO(x1, x2), data = cattle)
contour(cattle.m2, ~ x1+x2 ) # 2D contour plot for ADG

# 10.14.a pp 441-446 
# WP factors A,B and SP factors P,Q for ceramic pipe experiment
library(daewr)
data(cpipe); cpipe

# 10.14.b pp 441-446 
# WP factors A,B and SP factors P,Q for ceramic pipe experiment
library(daewr)
data(cpipe); cpipe
# REML analysis using lmer
cpipe$Asq <- cpipe$A^2 # add squared terms to data
cpipe$Bsq <- cpipe$A^2
cpipe$Psq <- cpipe$P^2 
cpipe$Qsq <- cpipe$Q^2; cpipe
library(lme4)
cpipe.reml.m <- lmer(y ~ A + B + P + Q 
           + A:B + A:P + A:Q 
           + B:P + B:Q + P:Q
           + Asq + Bsq + Psq + Qsq + (1|WP),
           data = cpipe)
summary(cpipe.reml.m)

# 10.14.c pp 441-446 
# WP factors A,B and SP factors P,Q for ceramic pipe experiment
library(daewr)
data(cpipe); cpipe
# Least Squares Analysis using rsm
library(rsm)
cpipe.ls.m <- rsm(y ~ SO(A,B,P,Q), data = cpipe)
summary(cpipe.ls.m)
# REML analysis using lmer
cpipe$Asq <- cpipe$A^2 # add squared terms to data
cpipe$Bsq <- cpipe$A^2
cpipe$Psq <- cpipe$P^2 
cpipe$Qsq <- cpipe$Q^2; cpipe
library(lme4)
cpipe.reml.m <- lmer(y ~ A + B + P + Q 
                     + A:B + A:P + A:Q 
                     + B:P + B:Q + P:Q
                     + Asq + Bsq + Psq + Qsq + (1|WP),
                     data = cpipe)
summary(cpipe.reml.m)

# 10.16.abc pp 441-446 
library(daewr)
# D-efficient EESPRS design: EEw1w2; EE14R7WP
#   1 WP factor and 2 SP factors
#   7 whole plots, 2 sub-plots
library(daewr)
alumina <- EEw1s2('EE14R7WP'); alumina  
alumina$s.area <- c(186.8782, 131.9686, 210.127, 187.2568, 140.336,  
                 163.783, 171.3387, 171.0459, 170.2234, 217.1703, 
                 202.6545, 191.1687,  137.4807, 143.1752); des # add s.area to data.frame
# Least Squares Analysis using rsm
library(rsm)
alumina.ls.m <- rsm(s.area ~ SO(w1,s1,s2), data = alumina)
summary(alumina.ls.m)
# REML analysis using lmer
alumina$w1sq <- alumina$w1^2 # add squared terms to data
alumina$s1sq <- alumina$s1^2 
alumina$s2sq <- alumina$s2^2; alumina
library(lme4)
alumina.reml.m <- lmer(s.area ~ w1 + s1 + s2 
                     + w1:s1 + w1:s2 + s1:s2
                     + w1sq + s1sq + s2sq + (1|WP),
                     data = alumina)
summary(alumina.reml.m)

