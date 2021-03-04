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

