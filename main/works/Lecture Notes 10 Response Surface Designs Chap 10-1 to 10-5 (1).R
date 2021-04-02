# Lecture Notes 10
# Chapter 10.1-10.5 Response Surface Designs

# 10.1 Introduction

# CCD design: 3 factor, 19 run
#   2^3 = 8 factorial, 3 in-center-factorial, 
#   3 x 2 = 6 axial, 2 in-center-axial
library(rsm) 
ccd.des <- ccd(defects~x1+x2+x3, n0=c(3,2), alpha="rotatable", randomize=FALSE)
ccd.des

des <- as.data.frame(ccd.des); des # convert to data.frame
des$defects <- c(61, 69, 60, 45, 64, 66, 44, 57,
  60, 61, 52, 41, 56, 61, 51, 63,
  60, 60, 52); des # add defects response to data.frame

# shows actual factor levels
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

# variance dispersion graph
library(Vdgraph)
Vdgraph(ccd.des[ , 3:5])
# fraction of design space plot
FDSPlot(ccd.des[ , 3:5], mod=2)

library(rsm) # lists all possible standard CCDs for 3 factors
ccd.pick(k=3)

# 10.2 Fundamentals of Response Surface Methodology

# BBD design: 3 factor, 17 run
#   2^2 factorial in each pair + other center point, 3 x 2^2 = 12
#   (arbitrary) 5 center points 
library(rsm) 
bbd.des <- bbd(defects ~ x1 + x2 + x3, n0 = 5, randomize = FALSE)
bbd.des

bbd.act <- bbd(defects ~ x1 + x2 + x3, 
    n0 = 5,  # 5 center points
    coding=list(x1~(temp-15)/15, # actual levels 
                x2~(click-90)/30, 
                x3~(hard-0.65)/0.15), 
    randomize = FALSE)
bbd.act
act <- as.coded.data(bbd.act); act # save as coded data
act$defects <- c(61, 69, 60, 45, 64, 66, 44, 57,
                 60, 61, 52, 41, 56, 61, 51, 63,
                 60); act # add defects response to coded data

# variance dispersion graph
library(Vdgraph)
Vdgraph(bbd.des[ , 3:5])
# fraction of design space plot
FDSPlot(bbd.des[ , 3:5], mod=2)

library(rsm) # compare relative prediction variance BBD with CCD
ccd2.des <- ccd.des[ ,3:5]
bbd2.des <- bbd.des[ ,3:5]
Compare2FDS(bbd2.des, ccd2.des, "BBD", "CCD", mod=2)

# 10.3 Standard Designs for Second Order Models

# SCD design: 3 factor, 11 run
#   4 factor points, 6 axial points, 1 center point
library(Vdgraph)
data(SCDH3)
scd.des <- SCDH3 

defects <- c(61, 69, 60, 45, 64, 
  66, 44, 57, 60, 61, 52); act # add defects response to coded data
act <- cbind(scd.des,defects); act

Vdgraph(SCDH3) # variance dispersion graph of 3-factor SCD

library(rsm) # compare relative prediction variance of BBD with CCD
ccd2.des <- ccd.des[ ,3:4] 
Compare2FDS(SCDH3, ccd2.des, "SCD", "CCD", mod=2)

# 10.4 Creating Standard Response Surface Designs in R

# hybrid design: 3 factor, 11 run
#   central composite design in k-1 = 3-1 = 2 factors
#   kth = 3rd factor makes design near rotatable
library(Vdgraph)
data(D311A)
hybrid.des <- D311A; D311A 

defects <- c(61, 69, 60, 45, 64, 
             66, 44, 57, 60, 61, 52); act # add defects response to coded data
act <- cbind(hybrid.des,defects); act

library(Vdgraph) # variance dispersion graph of 3-factor hybrid
data(D311A)
Vdgraph(D311A)

library(rsm) # compare relative prediction variance of hybrid with CCD
ccd2.des <- ccd.des[ ,3:4] 
Compare2FDS(D311A, ccd2.des, "D311A", "CCD", mod=2)

# 10.5 Non-Standard Response Surface Designs

# 3-parameter nonlinear model d = f(c) = exp(-T(c - c0)/15) + exp(-k(c - c0)) where
#   d = defects of sump pumps, response
#   c = number of clicks on/off of sump pump
#   T = (initial) temperature of water
#   k = hardness of water
#   c0 = initial number of clicks on/off of sump pump
T <- 15; k <- 0.50; c0 <- 60 # initial values for linear approx to model parameters
c <- seq(1:60) + 60; c # number of clicks ranges from 60 to 120
dfdT <- c(rep(0, 60)); dfdT # derivative wrt T, repeated 25 times
dfdk <- c(rep(0, 60)) # derivative wrt k, repeated 25 times
dfdc0 <- c(rep(0, 60)) # derivative wrt c0, repeated 25 times
for (i in 1:60) {
  dfdT[i] <- -1 * (c - c0)/15 * exp(-1 * T * (c[i] - c0)/15) 
  dfdk[i] <- -1 * (c - c0) * exp(-1 * k * (c[i] - c0)) 
  dfdc0[i] <- T * exp(-1 * T/15 * (c[i] - c0)/15) + k * exp(-1 * k * (c[i] - c0))
}
grid <- data.frame(c, dfdT, dfdk, dfdc0); grid
library(matrixcalc)
X <- as.matrix(grid); I <- t(X) %*% X; I
is.singular.matrix(I, tol = 1e-08)
det(I)

# collect I-optimal subset of grid of linear approx to model
library(AlgDesign)
nonlin.desI <- optFederov(~ -1 + dfdT + dfdk + dfdc0,
 data=grid,
 nTrials=3,
 center=TRUE, 
 criterion="I",
 nRepeats=20)
nonlin.desI$design

# collect D-optimal subset of grid of linear approx to model
library(AlgDesign)
nonlin.desD <- optFederov(~ -1 + dfdT + dfdk + dfdc0,
  data=grid,
  nTrials=3,
  center=TRUE, 
  criterion="D",
  nRepeats=20)
nonlin.desD$design

# 10.6 Fitting the Response Surface Model with R

# CCD design: 3 factor, 19 run
#   2^3 = 8 factorial, 3 in-center-factorial, 
#   3 x 2 = 6 axial, 2 in-center-axial
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

# ANOVA using Block 1 (factorial + center points) ONLY of sump pump runs
sp.lin <- rsm(defects ~ SO(x1, x2, x3), data = act, 
              subset = (Block == 1))
anova(sp.lin)

# ANOVA using both Block 1 and 2 runs of sump pump data
library(rsm)
sp.quad <- rsm(defects ~ SO(x1, x2, x3), data = act)
summary(sp.quad)

# ANOVA using both Block 1 and 2 runs of sump pump data
#  and also adding in Block term
library(rsm)
sp.quad <- rsm(defects ~ Block + SO(x1, x2, x3), data = act)
summary(sp.quad)

# nonlinear model fitting
clicks <- c(61, 62, 63, 70, 72, 73, 74, 75, 77, 79, 80, 83)
defects <- c(61, 69, 60, 45, 64, 66, 44, 57, 60, 61, 52, 41)
d <- as.data.frame(cbind(clicks, defects)); d
m.nls <- nls( defects  ~ exp( -Temp/15 * (clicks - c0)) + exp( -k * (clicks - c0)), 
              data = d,
              algorithm = "plinear",
              start = list(Temp = 15, k = 0.1, c0 = 60),
              trace = TRUE)
summary(m.nls) # does not give anova because nls does not converge

###################################################
# Homework 10 and Quiz 6 questions
###################################################

# 10.2.a,b,c pp 441-446
# CCD for 4 factors
library(rsm)
ccd.pick(k=4) # list standard CCDs for 4 factors
ccd.des <- ccd(defects ~ x1+x2+x3+x4,  # arbitrarily choose design 2
               n0=c(4,2),              # 2^4= 16 factorial + 4 center
               alpha="rotatable",      # 2^3= 8 axial runs + 2 center 
               randomize=FALSE)        
ccd.des
library(Vdgraph) # variance dispersion graph
par(mfrow = c(1,2))
Vdgraph(ccd.des[ , 3:6]) 
FDSPlot(ccd.des[ , 3:6], mod=2) # fraction of design space plot
par(mfrow = c(1,1))

# BBD with 4 factors
bbd.des <- bbd(defects ~ x1+x2+x3+x4, # 3 blocks of 10 runs
               n0 = 2,                # 2 center points per block 
               randomize = FALSE)   
bbd.des
library(Vdgraph) # variance dispersion graph
par(mfrow = c(1,2))
Vdgraph(bbd.des[ , 4:7]) 
FDSPlot(bbd.des[ , 4:7], mod=2) # fraction of BBD space plot
par(mfrow = c(1,1))

# SCD with 4 factors
library(Vdgraph)
data(SCDH4)
scd.des <- SCDH4; scd.des 
par(mfrow = c(1,2))
Vdgraph(scd.des) # variance dispersion graph of 4-factor SCD
FDSPlot(scd.des[ , 1:4], mod=2) # fraction of SCD space plot
par(mfrow = c(1,1))

# hybrid design
library(Vdgraph)
data(D416A)
hybrid.des <- D416A; hybrid.des 
par(mfrow = c(1,2))
Vdgraph(hybrid.des) # variance dispersion graph of 4-factor hybrid design
FDSPlot(hybrid.des[ , 1:4], mod=2) # fraction of hybrid design space plot
par(mfrow = c(1,1))

library(rsm) # compare relative prediction variance BBD with CCD
ccd.des <- ccd.des[ ,3:6]
bbd.des <- bbd.des[ ,4:7]
scd.des <- scd.des[ ,1:4]
hybrid.des <- hybrid.des[ ,1:4]
par(mfrow = c(1,3))
Compare2FDS(bbd.des, ccd.des, "BBD", "CCD", mod=2)
Compare2FDS(scd.des, ccd.des, "SCD", "CCD", mod=2)
Compare2FDS(hybrid.des, ccd.des, "Hybrid", "CCD", mod=2)
par(mfrow = c(1,1))

# 10.3 pp 441-446
RB <- c(10, 30, 10, 30, 10, 30,
        10, 30, 10, 30, 20, 20, 
        20, 20, 20, 20, 20)
AS <- c(0, 0, 2, 2, 0, 0,
        2, 2, 1, 1, 0, 2,
        1, 1, 1, 1, 1)
FT <- c(72, 72, 72, 72, 96, 96,
        96, 96, 84, 84, 84, 84,
        72, 96, 84, 84, 84)
Biomass <- c(3.83, 5.71, 6.74, 5.13, 5.55, 7.76, 
             12.45, 12.47, 11.54, 9.79, 7.13, 10.1, 
             6.29, 13.02, 10.66, 10.15, 10.97)

# 10.4 pp 441-446
# paper helicopter: CCD for 2 factors
library(rsm)
ccd.pick(k=2) # list standard CCDs for 2 factors
ccd.des <- ccd(flight.time ~ x1+x2,  # (arbitrarily) choose design 5
               n0=c(5,5),              # 2^2= 4 factorial + 5 center
               alpha="rotatable",      # 4 axial runs + 5 center 
               coding=list(x1~(wing.length-5.25)/1.75, 
                           x2~(heli.width-3.875)/1.625), 
               randomize=FALSE)        
ccd.des
des <- as.coded.data(ccd.des); des # save as coded data
des$flight.time <- c(
  5.2, 5.1, 5.3, 5.3, 5.2, 5.2, 5.0, 5.2, 5.2, 
  5.1, 5.1, 5.0, 5.2, 5.1, 5.1, 5.2, 5.1, 5.1); des
library(rsm)
des.m <- rsm(flight.time ~ Block + SO(x1, x2), data = des)
summary(des.m)

# 10.6 pp 441-446
# create grid x1 in (-1,1) and x2 in (-1,1)
# steps of 0.1, so 21 x 21 = 441 points
x1 <- seq(-1, 1, by=0.1); x1
x2 <- seq(-1, 1, by=0.1); x2
X1 <- c(rep(0,441)); X1
X2 <- c(rep(0,441)); X2
c <- 0
for (i in 1:21) {
  for (j in 1:21) {
    if ((x2[j] >= -2*x1[i] - 2) & (x2[j] <= -2*x1[i] + 1)) {
      c <- c + 1
      X1[c] <- x1[i]  
      X2[c] <- x2[j]
    }
  }
}
grid <- data.frame(cbind(X1,X2))
grid.constr <- grid[-(c+1):-nrow(grid),]; grid.constr # candidate points

# creating a D-optimal design
# for irregular experimental region
library(AlgDesign)
d <- data.frame(var=paste("x",1:2,sep=""),  
                low=-1,
                high=1,
                center=0,
                nLevels=21,
                round=1,
                factor=FALSE); d
constFcn <- function(x){x[2] >= -2*x[1] - 2 && x[2] <= -2*x[1] + 1}
desCon <- optMonteCarlo(~quad(.), 
                        d, 
                        constraint=constFcn,
                        criterion = "D"); desCon

# finer grid for irregular experimental region
library(AlgDesign)
d <- data.frame(var=paste("x",1:2,sep=""),  
                low=-1,
                high=1,
                center=0,
                nLevels=201,
                round=2,
                factor=FALSE); d
constFcn <- function(x){x[2] >= -2*x[1] - 2 && x[2] <= -2*x[1] + 1}
desCon <- optMonteCarlo(~quad(.), 
                        d, 
                        constraint=constFcn,
                        criterion = "D"); desCon

