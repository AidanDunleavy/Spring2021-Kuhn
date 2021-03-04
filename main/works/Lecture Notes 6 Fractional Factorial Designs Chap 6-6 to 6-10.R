# Lecture Notes 6
# Chapter 6.6-6.10 Fractional Factorial Designs

library(EMSaov)
library(AlgDesign)
library(BsMD)
library(car)
library(partitions)
library(DoE.base)
library(FrF2)
library(GAD)
library(gmodels)
library(leaps)
library(lme4)
library(lsmeans)
library(mixexp)
library(multcomp)
library(Vdgraph)
library(pracma)
library(zip)
library(later)
library(htmltools)
library(digest)
library(promises)
library(classInt)
library(units)
library(agricolae)
library(daewr)

# 6.6 Plackett-Burman (PB) and Model Robust Screening Designs

# 12-run PB design 1
library(FrF2)
pb( nruns = 12, randomize=FALSE)

# another 12-run PB design 2, for sump.pump data
library(BsMD)
data( PB12Des, package = "BsMD" )
colnames(PB12Des) <- c("e", "d", "c", "b", "a", "F", "E", "D", "C", "B", "A")
PB2 <- PB12Des[c(11,10,9,8,7,6,5,4,3,2,1)]; PB2
defects <- c(61, 69, 60, 45, 64, 66,
             44, 57, 60, 61, 52, 42)
sp.df <- cbind(PB2, defects); sp.df

# optimal subset of effects for 12-run, 6 factor PB design 3
library(daewr)
opt.sump.pump <- OptPB(12, 6, randomize = FALSE)
opt.sump.pump

# half-normal for sump.pump data, using PB 2
modpb <- lm( defects ~ (.), data = PB2 ); modpb
library(daewr)
cfs <- coef(modpb)[2:12]; cfs
names<-names(cfs)
halfnorm(cfs, names, alpha = .35, refline=FALSE)

# optimal PB: highest chance of nonsingularity 
library(daewr)
opt.model.nonsig <- OptPB(12, 6, randomize = FALSE)
opt.model.nonsig

# color map of PB design
library(daewr)
sp.df.map <- sp.df[ , c(1:6)] # pick out first six factors A to F
colormap(sp.df.map,mod=2) # color map of main and 2-factor interactions

# Altscreen, alternative to PB: optimal 16-run 6-factor 
library(daewr)
sp.model.alt <- Altscreen(6, randomize= FALSE)
head(sp.model.alt)
colormap(sp.model.alt,mod=2) # color map of main and 2-factor interactions

# ModelRobust, alternative to PB: optimal 16-run 7-factor 
library(daewr)
ModelRobust()
sp.model.robust <- ModelRobust('MR16m7g5')
sp.model.robust
colormap(sp.model.robust,mod=2) # color map of main and 2-factor interactions

# exhaustive regression subset selection 
#  of main and 2-factor interactions
sp.model <- sp.df[ , c(1:6, 12)]; sp.model # pick main factors and defects response
library(leaps)
modpbr <- regsubsets(
  defects ~ (.)^2, 
  data = sp.model, 
  method = "exhaustive", 
  nvmax = 7, # include at most 8 effects in model  
  nbest = 3) # keep best 3 models of each size
rs <- summary(modpbr); rs; rs$adjr2
plot(c(rep(1:8,each=3)), # plot of models and r^2_abj
     rs$adjr2, 
     xlab="Number of Effects in Model", 
     ylab="Adjusted R-square")
plot(modpbr,scale="r2") # r^2_adj for effects in various models

# forward stepwise to identify model
sp.model <- sp.df[ , c(1:6, 12)]; sp.model # pick main factors and defects response
null <- lm( defects ~ 1, data = sp.model ); null # only intercept
up <- lm( defects ~ (.)^2, data = sp.model ); up # main + 2-factor
step( null, 
      scope = list(lower = null, upper = up), 
      direction="forward", # add effects 
      steps=4) # 4 steps only

# forward stepwise with heredity
sp.des <- sp.df[ , c(1:6)]; sp.des # main factors design
defects <- sp.df[ , 12]; defects # defects response
library(daewr)
trm <- ihstep( defects, sp.des, m=0, c=6 ) # initial step, 6 1-factors
trm <- fhstep( defects, sp.des, trm, m=0, c=6) # forward step
trm <- fhstep( defects, sp.des, trm, m=0, c=6) # forward step
# bstep does NOT work

# 6.7 Mixed Level Factorials and Orthogonal Arrays (OAs)

# runs necessary for orthogonality
#  required for 4 x 3-level and 2 x 2-level factorial
# find orthogonal array
des <- oa.design(nlevels=c(3,3,3,3,2,2), 
  nruns=36,
  columns="min3", # minimize aliasing main with 2-factor interact
  randomize=TRUE,
  seed=104); des

# D-optimal subset of full factorial (orthogonal) design 
library(DoE.base)
cand <- oa.design(nlevels=c(3,3,3,3,2,2),
  nruns=36, # number of candidate runs
  columns="min3",
  seed=104) # generate candidate full ortho design runs
# select D-optimal 18-run subset from candidates
# to allow estimation of the main effects
library(AlgDesign)
optim <- optFederov(~A+B+C+D+E+F,
  cand,
  nRepeats=10,
  nTrials=18, # number of D-optimal runs
  criterion="D"); optim

# D-optimal subset of orthogonal fractional factorial: sump pump
defects <- c(61, 69, 60, 45, 64, 66,
             44, 57, 60, 61, 52, 42)
sp.df <- cbind(PB2, defects); sp.df
library(DoE.base)
cand <- oa.design(nlevels=c(3,3,3,3,2,2),
                  nruns=36, # number of candidate runs
                  columns="min3",
                  seed=104) # generate candidate full ortho design runs
# select D-optimal 18-run subset from candidates
# to allow estimation of the main effects
library(AlgDesign)
optim <- optFederov(~A+B+C+D+E+F,
                    cand,
                    nRepeats=10,
                    nTrials=18, # number of D-optimal runs
                    criterion="D"); optim

# library of orthogonal arrays
# often gives more factors than specified to ensure orthogonality
library("DoE.base") 
show.oas(factors=list(nlevels=c(3,2),number=c(4,2)),show= "all")
# 4 x 3-level and 2 x 2-level factors requested but
# 4 x 3-level and 16 x 2-levels needed for full orthogonality
oa.design(L36.2.16.3.4) 

# 6.8 Definitive Screening Designs

# sump pump example
library(daewr)
des <- DefScreen( 
  m = 4, # four 3-level main factors 
  c = 2, # two 2-level main factors
  randomize = FALSE )
des

# Figure 6.17 p. 249
library(daewr)
defs <- data.frame(cbind(test[, 2:21 ])); defs
colormap(defs, mod=1)


# 6.9 Review of Important Concepts

###################################################
# Homework 6 and Quiz 3 questions
###################################################

# 6.10.a.1 pp 252-259
# estimate from pb design
r1 <- c(1,-1, 1,-1,-1,-1, 1, 1, 1,-1, 1, 56)
r2 <- c(1, 1,-1, 1,-1,-1,-1, 1, 1, 1,-1, 93)
r3 <- c(-1, 1, 1,-1, 1,-1,-1,-1, 1, 1, 1, 67)
r4 <- c(1,-1, 1, 1,-1, 1,-1,-1,-1, 1, 1, 60)
r5 <- c(1, 1,-1, 1, 1,-1, 1,-1,-1,-1, 1, 77)
r6 <- c(1, 1, 1,-1, 1, 1,-1, 1,-1,-1,-1, 65)
r7 <- c(-1, 1, 1, 1,-1, 1, 1,-1, 1,-1,-1, 95)
r8 <- c(-1,-1, 1, 1, 1,-1, 1, 1,-1, 1,-1, 49)
r9 <- c(-1,-1,-1, 1, 1, 1,-1, 1, 1,-1, 1, 44)
r10 <- c(1,-1,-1,-1, 1, 1, 1,-1, 1, 1,-1, 63)
r11 <- c(-1, 1,-1,-1,-1, 1, 1, 1,-1, 1, 1, 63)
r12 <- c(-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, 61)
pb.cre.df <- as.data.frame(rbind(r1,r2,r3,r4,r5,r6,
                                 r7,r8,r9,r10,r11,r12)); pb.cre.df
colnames(pb.cre.df) <- c("A","B","C","D","E",
                         "c6","c7","c8","c9","c10","c11","reacted"); pb.cre.df
m.cre <- lm( reacted ~ A + B + C + D + E +
               c6 + c7 + c8 + c9 + c10 + c11, 
             data = pb.cre.df)
summary(m.cre) 

# 6.10.a.2 pp 252-259
# significant estimates from half-normal plot from pb design
library(daewr)
cfs <- coef(m.cre)[2:12]; cfs
names<-names(cfs)
halfnorm(cfs, names, alpha = .15, refline=FALSE)

# 6.10.bc.1 pp 252-259
# best main and 2-factor-interactions subset from pb design
cre.m <- pb.cre.df[ , c(1:5, 12)]; cre.m # pick main factors and response
library(leaps)
modpbr <- regsubsets(
  reacted ~ (.)^2, # main and 2-factor interactions
  data = cre.m, # main effects from model
  method = "exhaustive", 
  nvmax = 4, # include at most 4 effects in model  
  nbest = 4) # keep best 4 models of each size
rs <- summary(modpbr); rs; rs$adjr2

# 6.10.d pp 252-259
# use BD interaction plot to 
# determine factor settings to maximize reacted
B <- cre.m$B; D <- cre.m$D 
reacted <- cre.m$reacted; B; D; reacted
interaction.plot(B, D, 
                 reacted, 
                 type = "b", pch=c(18,24,22), 
                 leg.bty="o",
                 main="Interaction Plot for B and D",
                 xlab = "B", ylab = "% reacted")

# 6.12.ab pp 252-259
# 18-run orthogonal fractional factorial of 8 factors
r1 <- c(1,1,1,1,1,1,1,1)
r2 <- c(1,1,2,2,2,2,2,2)
r3 <- c(1,1,3,3,3,3,3,3)
r4 <- c(1,2,1,1,2,2,3,3)
r5 <- c(1,2,2,2,3,3,1,1)
r6 <- c(1,2,3,3,1,1,2,2)
r7 <- c(1,3,1,2,1,3,2,3)
r8 <- c(1,3,2,3,2,1,3,1)
r9 <- c(1,3,3,1,3,2,1,2)
r10 <- c(2,1,1,3,3,2,2,1)
r11 <- c(2,1,2,1,1,3,3,2)
r12 <- c(2,1,3,2,2,1,1,3)
r13 <- c(2,2,1,2,3,1,3,2)
r14 <- c(2,2,2,3,1,2,1,3)
r15 <- c(2,2,3,1,2,3,2,1)
r16 <- c(2,3,1,3,2,3,1,2)
r17 <- c(2,3,2,1,3,1,2,3)
r18 <- c(2,3,3,2,1,2,3,1)
df <- as.data.frame(
  rbind(r1,r2,r3,r4,r5,r6,
        r7,r8,r9,r10,r11,r12,
        r13,r14,r15,r16,r17,r18)
  ); df
df[] <- lapply( df, factor)
SN <- c(-0.66, -0.66, -0.69, -0.67, -0.82, -0.69, 
        -0.75, -0.65, -0.66, -0.80, -0.67, -0.65, 
        -1.04, -0.69, -0.67, -0.67, -0.68, -0.71)
d <- cbind(df,SN); d; class(d)
colnames(d) <- c("A","B","C","D","E",
      "F","G","H","SN"); d
attach(d)
contrasts(d$A) = contr.poly(2) # orthogonal contrast for 2-level factor
contrasts(d$B) = contr.poly(3) # orthogonal contrast for 3-level factor
contrasts(d$C) = contr.poly(3)
contrasts(d$D) = contr.poly(3)
contrasts(d$E) = contr.poly(3)
contrasts(d$F) = contr.poly(3)
contrasts(d$G) = contr.poly(3)
contrasts(d$H) = contr.poly(3)
m.LQ <- lm(SN ~ A + B + C + D + E + F + G + H, d)
summary(m.LQ)

# 6.12.c pp 252-259
# fit 5 most significant variables from last model
BQ <- contr.poly(3)[d$B,".Q"]
CL <- contr.poly(3)[d$C,".L"]
DQ <- contr.poly(3)[d$D,".Q"]
EL <- contr.poly(3)[d$E,".L"]
EQ <- contr.poly(3)[d$E,".Q"]
m.LQ <- lm(SN ~ BQ + CL + DQ + EL + EQ, d)
summary(m.LQ)

# 6.12.d pp 252-259
# fit AL, CL, DQ, BQ*DQ, FL*HQ
AL <- contr.poly(3)[d$A,".L"]
BQ <- contr.poly(3)[d$B,".Q"]
CL <- contr.poly(3)[d$C,".L"]
DQ <- contr.poly(3)[d$D,".Q"]
FL <- contr.poly(3)[d$F,".L"]
HQ <- contr.poly(3)[d$H,".Q"]
BQDQ <- BQ*DQ; FLHQ <- FL*HQ
m.LQ <- lm(SN ~ AL + CL + DQ + BQDQ + FLHQ, d)
summary(m.LQ)

# 6.12.e pp 252-259
# color map of A, B, C, D, AB, AC, AD, BC
a <- as.numeric(as.vector(d$A)); b <- as.numeric(as.vector(d$B)); a; b
c <- as.numeric(as.vector(d$C)); d <- as.numeric(as.vector(d$D)); c; d
ab <- a * b; ab
ac <- a * c; ac
ad <- a * d; ad
bc <- b * c; bc
ds <- as.data.frame(cbind(a,b,c,d,ab,ac,ad,bc)); ds; class(ds)
colnames(ds) <- c("A","B","C","D","AB",
                 "AC","AD","BC"); ds
library(daewr)
colormap(ds, mod=1)

