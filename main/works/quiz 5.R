
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
ccd.des

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
bbd.des

# SCD with 4 factors
library(Vdgraph)
data(SCDH4)
scd.des <- SCDH4; scd.des 
par(mfrow = c(1,2))
Vdgraph(scd.des) # variance dispersion graph of 4-factor SCD
FDSPlot(scd.des[ , 1:4], mod=2) # fraction of SCD space plot
par(mfrow = c(1,1))
scd.des
# hybrid design
library(Vdgraph)
data(D416A)
hybrid.des <- D416A; hybrid.des 
par(mfrow = c(1,2))
Vdgraph(hybrid.des) # variance dispersion graph of 4-factor hybrid design
FDSPlot(hybrid.des[ , 1:4], mod=2) # fraction of hybrid design space plot
par(mfrow = c(1,1))
hybrid.des

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
contour(des.m, ~ x1 + x2)
des$x2
 des$x1
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
 
