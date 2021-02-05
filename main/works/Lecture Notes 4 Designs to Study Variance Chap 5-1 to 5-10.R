# Lecture Notes 4
# Chapter 5 Designs to Study Variances

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

# 5.1 Introduction
# 5.2 Random Factors and Random Sampling Experiments

# BALANCED
D <- expand.grid( lab = c('A', 'B', 'C', 'D'), # labs
                  patient = c(1, 2, 3, 4,
                              5, 6, 7, 8) ) # combinations of factors
D[] <- lapply(D, factor) # convert to factors
PTSD <- c(
  10.3, 6.8,  8.2,  2.4,
  9.1, 12.1,  6.5,  6.1,
  6.1, 5.1,  7.2,  8.1,
  7.2, 9.8,  8.1,  5.1,
  5.4, 4.2,  4.1,  7.2,
  7.0, 9.8,  8.1,  4.1,
  5.4, 4.3,  4.1,  7.2,
  6.2, 6.2,  4.9,  9.1); 
PTSDdata <- data.frame(cbind(D,PTSD)); PTSDdata # data.frame
m1 <- aov( PTSD ~ lab, data = PTSDdata ) # model 1, m1
s1 <- summary(m1); s1

# UNBALANCED
PTSD.unbal <- c(10.3,9.1,6.1,7.2,5.4,7.0,5.4,
                6.8,12.1,5.1,9.8,4.2,9.8,4.3,6.2,  
                8.2,6.5,7.2,8.1,4.1,8.1,4.1, 
                2.4,6.1,8.1,5.1,7.2,4.1,7.2,9.1); 
lab <- factor( c( rep ('A', 7), rep ('B', 8),
                  rep ('C', 7), rep ('D', 8) ) ); lab
PTSDdata.unbal <- data.frame(lab, PTSD.unbal); PTSDdata.unbal

# use UNbalanced PTSDdata
m1.unbal <- aov( PTSD.unbal ~ lab, data = PTSDdata.unbal ) # model 1, m1
s1.unbal <- summary(m1.unbal); s1.unbal

# 5.3 One-Factor Sampling Designs
# 5.3.1 method of moments estimation of variance components
# 5.3.1 equate expected MS equated with observed MS
# 5.3.1 balanced case
c <- r <- 8; c # coefficient = # replications 
# then get expected value of msE in ANOVA
sigma.2 <- as.matrix( s1[[1]][2,3] )
ms.lab <- as.matrix( s1[[1]][1,3] )
cat(" Mean Square for Lab = ",ms.lab,"\n",
    " Mean Square for Error = ", sigma.2,"\n",
    " Expected Mean Square for Lab","\n", 
    "Var(error)+",c,"Var(Lab)","\n")
# compare expected to observed: 
sigma.2t <- (ms.lab - sigma.2) / c
cat("Method of Moments Variance Component Estimates","\n", 
    "Var(error)=",sigma.2,"\n",
    "Var(Lab)=",sigma.2t,"\n")

# 5.3.1 MM UNbalanced case
Xmtx <- as.matrix( model.matrix(m1.unbal) ); Xmtx # extract lab indicators
Blab  <- Xmtx[,2]
Clab <- Xmtx[,3]
Dlab <- Xmtx[,4]
Alab <- Xmtx[,1]-(Blab+Clab+Dlab)
sum1 <- summary(aov (Alab ~ PTSDdata.unbal$lab )); sum1
ms1 <- as.matrix( sum1[[1]][1,3] ); ms1
sum2 <- summary(aov( Blab ~ PTSDdata.unbal$lab )); sum2
ms2 <- as.matrix( sum2[[1]][1,3] ); ms2
sum3 <- summary(aov(Clab ~ PTSDdata.unbal$lab )); sum3
ms3 <- as.matrix( sum3[[1]][1,3] )
sum4 <- summary(aov(Dlab ~ PTSDdata.unbal$lab )); sum4
ms4 <- as.matrix( sum4[[1]][1,3] )
c <- ms1+ms2+ms3+ms4; c # coefficient 
# then get expected value of msE in ANOVA
sigma.2 <- as.matrix( s1.unbal[[1]][2,3] )
ms.lab <- as.matrix( s1.unbal[[1]][1,3] )
cat(" Mean Square for Lab = ",ms.lab,"\n",
    " Mean Square for Error = ", sigma.2,"\n",
    " Expected Mean Square for Lab","\n", 
    "Var(error)+",c,"Var(Lab)","\n")
# compare expected to observed: 
sigma.2t <- (ms.lab - sigma.2) / c
cat("Method of Moments Variance Component Estimates","\n", 
    "Var(error)=",sigma.2,"\n",
    "Var(Lab)=",sigma.2t,"\n")

# 5.3.1 REML Balanced case
library(lme4)
rm <- lmer( PTSD ~ 1 + (1|lab), data = PTSDdata) # random model rm
summary(rm) # # splits into sigma^2_residual, sigma^2_lab

# 5.3.1 REML UNbalanced case
library(lme4)
rm.unbal <- lmer( PTSD.unbal ~ 1 + (1|lab), data = PTSDdata.unbal) 
summary(rm.unbal) # # splits into sigma^2_residual, sigma^2_lab

# 5.3.1 likelihood profile approximate CIs: balanced case 
likel.profile <- profile( # using likelihood profile method
  rm.likel <- lmer( PTSD ~ 1 + (1| lab), # rm.likel can be anything,
                    data = PTSDdata,     #  is just a dummary variable 
                    REML = FALSE))
confint (likel.profile) # 95% CI

# 5.3.1 likelihood profile approximate CIs: unbalanced case 
likel.profile <- profile( # using likelihood profile method
  rm.likel <- lmer( PTSD.unbal ~ 1 + (1| lab), # rm.likel can be anything,
                    data = PTSDdata.unbal,     #  is just a dummary variable 
                    REML = FALSE))
confint (likel.profile) # 95% CI

# 5.4 Estimating Variance Components

# 5.4.1 interval estimates
D <- expand.grid( lab = c('A', 'B', 'C', 'D'), # labs
                  patient = c(1, 2, 3, 4,
                              5, 6, 7, 8) ) # combinations of factors
D[] <- lapply(D, factor) # convert to factors
PTSD <- c(
  10.3, 6.8,  8.2,  2.4,
  9.1, 12.1,  6.5,  6.1,
  6.1, 5.1,  7.2,  8.1,
  7.2, 9.8,  8.1,  5.1,
  5.4, 4.2,  4.1,  7.2,
  7.0, 9.8,  8.1,  4.1,
  5.4, 4.3,  4.1,  7.2,
  6.2, 6.2,  4.9,  9.1); 
d <- PTSDdata <- data.frame(cbind(D,PTSD)); PTSDdata # data.frame
summary(aov(PTSD ~ lab, PTSDdata))

conf.coeff <- 0.95

var.component.1fac.eqrep.interval <- function(parameter, conf.coeff, data.frame) {
  columnNames <- names(data.frame); columnNames # simplify data frame variable names
  y <- get (columnNames[3], d); y 
  t <- get (columnNames[1], d); t
  
  sumANOVA <- summary(aov(y ~ t, d)); sumANOVA # ANOVA summary
  ssE <- as.matrix( sumANOVA[[1]][2,2] ); ssE # extract ssE, dfE, ssT, dfT from ANOVA
  dfE <- as.matrix( sumANOVA[[1]][2,1] ); dfE
  ssT <- as.matrix( sumANOVA[[1]][1,2] ); ssT
  dfT <- as.matrix( sumANOVA[[1]][1,1] ); dfT
  Tr <- dfT + 1; r <- (dfT + dfE + 1)/Tr; Tr; r # number of treatments, replications
  Fstat <- as.matrix( sumANOVA[[1]][1,4] ); Fstat
  F.crit.lower <- qf(alpha/2, dfT, dfE, lower.tail=TRUE); F.crit.lower
  F.crit.upper <- qf(alpha/2, dfT, dfE, lower.tail=FALSE); F.crit.upper

  if(parameter=="s^2") {
    alpha <- 1 - conf.coeff
    chi2.crit.lower <- qchisq(alpha/2, Tr*(r-1), lower.tail=TRUE); chi2.crit.lower
    chi2.crit.upper <- qchisq(alpha/2, Tr*(r-1), lower.tail=FALSE); chi2.crit.upper
    ci.lower <- ssE/chi2.crit.upper; ci.lower
    ci.upper <- ssE/chi2.crit.lower; ci.upper
    dat <- c(parameter, ci.lower, ci.upper)    
    names(dat) <- c(" ", "CI lower", "CI upper")
  }
  if(parameter=="s_t^2") {
    alpha <- (1 - conf.coeff)/2
    chi2.crit.lower <- qchisq(alpha/2, Tr-1, lower.tail=TRUE); chi2.crit.lower
    chi2.crit.upper <- qchisq(alpha/2, Tr-1, lower.tail=FALSE); chi2.crit.upper
    ci.lower <- ssT*(1-F.crit.upper/Fstat)/(r*chi2.crit.upper); ci.lower
    ci.upper <- ssT*(1-F.crit.lower/Fstat)/(r*chi2.crit.lower); ci.upper
    dat <- c(parameter, ci.lower, ci.upper)    
    names(dat) <- c(" ", "CI lower", "CI upper")
  }
  if(parameter=="s_t^2/(s_t^2+s^2)") {
    alpha <- 1 - conf.coeff
    ci.lower <- (Fstat/F.crit.upper - 1)/(r + Fstat/F.crit.upper - 1); ci.lower
    ci.upper <- (Fstat/F.crit.lower - 1)/(r + Fstat/F.crit.lower - 1); ci.upper
    dat <- c(parameter, ci.lower, ci.upper)    
    names(dat) <- c(" ", "CI lower", "CI upper")
  }
  if(parameter=="s^2/(s_t^2+s^2)") {
    alpha <- 1 - conf.coeff
    ci.lower <- r/(r + Fstat/F.crit.lower - 1); ci.lower
    ci.upper <- r/(r + Fstat/F.crit.upper - 1); ci.upper
    dat <- c(parameter, ci.lower, ci.upper)    
    names(dat) <- c(" ", "CI lower", "CI upper")
  }
  if(parameter=="s_t^2/s^2") {
    alpha <- 1 - conf.coeff
    ci.lower <- (Fstat/F.crit.upper - 1)/r; ci.lower
    ci.upper <- (Fstat/F.crit.lower - 1)/r; ci.upper
    dat <- c(parameter, ci.lower, ci.upper)    
    names(dat) <- c(" ", "CI lower", "CI upper")
  }
  return(dat) 
}
var.component.1fac.eqrep.interval("s^2", 0.95, d) # CI for sigma^2
var.component.1fac.eqrep.interval("s_t^2", 0.90, d) # CI for sigma_t^2
var.component.1fac.eqrep.interval("s_t^2/(s_t^2+s^2)", 0.95, d) 
var.component.1fac.eqrep.interval("s^2/(s_t^2+s^2)", 0.95, d) 
var.component.1fac.eqrep.interval("s_t^2/s^2", 0.95, d) 

# 5.4.2.a determine t and r; specifically, nu2 = t(r - 1) given CI s
alpha = 0.10
nu2 <- 36:44
chiu <- qchisq(1 - alpha/2, nu2)
chil <- qchisq(alpha/2, nu2)
width <- nu2 * (chiu - chil) / (chil * chiu)
halfw <- width/2
data.frame(nu2, width, halfw)

# 5.4.2.b rule of thumb if s^2_t > s^2, t = nu2, r = 2 
PTSDdata <- data.frame(cbind(D,PTSD)); PTSDdata # data.frame
m1 <- aov( PTSD ~ lab, data = PTSDdata ) # model 1, m1
s1 <- summary(m1); s1

# 5.4.2.c determine t and r; specifically, nu2 = t(r - 1) 
# given rho = sigma_t/sigma, power 
alpha <- .05
rho <- 2.0
t <- rep(5:7, each = 3)
r <- rep(2:4, 3)
nu_1 <- t-1
nu_2 <- t * (r - 1)
fcrit <- qf( 1 - alpha, nu_1, nu_2 )
factor <- 1 / ( 1 + r * rho )
plimit <- factor * fcrit
power <- 1 - pf( plimit, nu_1, nu_2 )
data.frame( r, t, power)

# 5.5. Two-Factor Sampling Designs

# 5.5.1 estimating variance components, balanced
D <- expand.grid( lab = c('A', 'B', 'C', 'D'), # design
                  patient = c(1, 2, 3, 4) ) # combinations of factors
D[] <- lapply(D, factor) # convert to factors
D <- rbind(D, D) # row-bind twice, for 2 replicates
PTSD2 <- c(
  10.3, 6.8, 8.2,  2.4,
  6.1, 5.1,  7.2,  8.1,
  5.4, 4.2,  4.1,  7.2,
  5.4, 4.3,  4.1,  7.2,
  9.1, 12.1, 6.5,  6.1,
  7.2, 9.8,  8.1,  5.1,
  7.0, 9.8,  8.1,  4.1,
  6.2, 6.2,  4.9,  9.1
);
PTSD2data <- data.frame(cbind(D,PTSD2)); PTSD2data # data.frame

# 5.5.1 test both main factors and interaction
m <- aov( PTSD2 ~ lab + patient + lab:patient, data = PTSD2data)
sm <- summary(m); sm # incorrect ANOVA with incorrect F ratios
msE <- as.matrix( sm[[1]][4,3] ); msE; dfE <- as.matrix( sm[[1]][4,1] ); dfE
msAB <- as.matrix( sm[[1]][3,3] ); msAB; dfAB <- as.matrix( sm[[1]][4,1] ); dfAB
msB <- as.matrix( sm[[1]][2,3] ); msB; dfB <- as.matrix( sm[[1]][4,1] ); dfB
msA <- as.matrix( sm[[1]][1,3] ); msA; dfA <- as.matrix( sm[[1]][4,1] ); dfA
F.AB <- msAB/msE; F.C; probF.AB <- 1 - pf(F.AB, dfAB, dfE) 
F.B <- msB/msAB; F.B; probF.B <- 1 - pf(F.B, dfB, dfAB) 
F.A <- msA/msAB; F.A; probF.A <- 1 - pf(F.A, dfA, dfAB) 
cat("corrected ANOVA F statistics and p-values for NSE design","\n",
    "F.A =",F.A, " P(>F.A) =", probF.A,"\n",
    "F.B =",F.B, " P(>F.B) =", probF.B,"\n",
    "F.AB =",F.AB, " P(>F.AB) =", probF.AB,"\n")

# 5.5.1 method of moments
m.2bal <- aov( PTSD2 ~ patient + lab + patient:lab, data = PTSD2data)
sum.m.2bal <- summary(m.2bal) # ANOVA summary for method of moments analysis
# using Bennett and Franklin
r <- 2; a <- 4; b <- 4 # 2 replications per cell; number of a and b
sigma2 <- msE <- as.matrix( sum.m.2bal[[1]][4,3] ); msE # extract from ANOVA
msAB <- as.matrix( sum.m.2bal[[1]][3,3] ); msAB 
msB <- as.matrix( sum.m.2bal[[1]][2,3] ); msB 
msA <- as.matrix( sum.m.2bal[[1]][1,3] ); msA 
sigma2patientlab <- (msAB - sigma2) / r; sigma2patientlab # method of moments
sigma2lab <- (msB - sigma2 - r * sigma2patientlab ) / (r*a); sigma2lab
sigma2patient <- (msA - sigma2 - r * sigma2patientlab ) / (r*b)
cat("Method of Moments Variance Component Estimates","\n",
    "Var(error)=",sigma2,"\n",
    "Var(patient x lab)=",sigma2patientlab,"\n",
    "Var(lab)=",sigma2lab,"\n",
    "Var(patient)=",sigma2patient,"\n")

# balanced, REML
library(lme4)
m.reml.2bal <- lmer(PTSD2 ~ 1 
                    + (1|patient) + (1|lab) 
                    + (1|patient:lab), data = PTSD2data)
summary(m.reml.2bal)

# 5.5.5 balanced, approximate CIs 2-factor sampling experiment
m.2bal <- aov( PTSD2 ~ patient + lab            # collect msb, msab
               + patient:lab, data = PTSD2data) # from ANOVA
summary(m.2bal)
library(daewr)
options(digits = 3)
# confidence level = 90% for sigma^2_b, and from ANOVA table,
# c1 = 1/8, ms1 = msb = 2.319, nu1 = 3,
# c2 = 1/8, ms2 = msab = 5.530, nu2 = 9
vci(confl = .90, 
    c1 = 1/8, ms1 = 2.319, nu1 = 3, 
    c2 = 1/8, ms2 = 5.530, nu2 = 9)
library(daewr)
options(digits = 3)
# confidence level = 90% for sigma^2_ab, and from ANOVA table,
# c1 = 1/2, ms1 = msab = 5.53, nu1 = 9,
# c2 = 1/2, ms2 = mse = 4.602, nu2 = 16
vci(confl = .90, 
    c1 = 1/2, ms1 = 5.53, nu1 = 9, 
    c2 = 1/2, ms2 = 4.602, nu2 = 16)

# 5.5.6 determine t, a, b; specifically, nu2 = (r - 1)ab given CI s
alpha = 0.10
nu2 <- 36:44
chiu <- qchisq(1 - alpha/2, nu2)
chil <- qchisq(alpha/2, nu2)
width <- nu2 * (chiu - chil) / (chil * chiu)
halfw <- width/2
data.frame(nu2, width, halfw)

# 5.5.8 2-factor, unbalanced data
PTSD2.unbal <- c(
  10.3, 6.8, 8.2,  2.4,
  12.1, 6.5,  6.1,
  6.1, 5.1,  7.2,  8.1,
  7.2, 9.8,  8.1,
  5.4, 4.2,  4.1,  7.2,
  7.0, 9.8,        4.1,
  5.4, 4.3,  4.1,  7.2,
  6.2,       4.9,  9.1
); 
patient.unbal <- factor( c( rep (1, 7), rep (2, 7),
                            rep (3, 7), rep (4, 7) ) ); patient.unbal
lab.unbal <- factor( c(
  'A','B','C','D',
  'B','C','D',
  'A','B','C','D',
  'A','B','C',
  'A','B','C','D',
  'A','B',    'D',
  'A','B','C','D',
  'A',    'C','D') ); lab.unbal
PTSD2data.unbal <- data.frame(lab.unbal, patient.unbal, PTSD2.unbal); PTSD2data.unbal

# 5.5.7 2-factor, unequal replicates, REML
library(lme4)
m.2unbal <- lmer( # ANOVA summary for REML analysis
  PTSD2.unbal ~ 1 + (1|lab.unbal) + (1|patient.unbal) 
  + (1|lab.unbal:patient.unbal), 
  data = PTSD2data.unbal)
sum.m.2unbal <- summary(m.2unbal); sum.m.2unbal

# 5.5.8 2-factor, unequal replicates, CI
cellmeans <- tapply( 
  PTSD2data.unbal$PTSD2.unbal, 
  list(PTSD2data.unbal$lab.unbal, PTSD2data.unbal$patient.unbal), 
  mean); cellmeans # collapse unbalanced replicated to balanced (1 ave per cell) unreplicated
dim(cellmeans) <- NULL
Lab <- factor( rep(c("A", "B", "C", "D"), 4))
Patient <- factor(rep( c( 1, 2, 3, 4), each = 4))
m.2unbal.CI <- aov( cellmeans ~ Lab + Patient + Lab:Patient )
sum.m.2unbal.CI <- summary(m.2unbal.CI); sum.m.2unbal.CI

options(digits = 3)
# c-bar = ab/[sum sum(1/r_ij)]; get a, b from ANOVA table
a <- as.matrix( sum.m.2unbal.CI[[1]][1,1] ) + 1; a 
b <- as.matrix( sum.m.2unbal.CI[[1]][2,1] ) + 1; b 
cbar <- a*b/sum(1/tapply( # a x b then divide by (count replicates, reciprocate, sum) 
  PTSD2data.unbal$PTSD2.unbal, 
  list(PTSD2data.unbal$lab.unbal, PTSD2data.unbal$patient.unbal), 
  FUN = length)); cbar 
c1 <- c2 <- 1/(b*cbar); c1; c2
# c1 = 1/(b*cbar) = 0.15625, ms1 = msa = 1.435, nu1 = df.a = 3,
# c2 = 1/(b*cbar) = 0.15625, ms2 = msab = 3.609, nu2 = df.ab = 6
# get msa, msab, nu1, nu2 from ANOVA table
msa <- ms1 <- as.matrix( sum.m.2unbal.CI[[1]][1,3] ); ms1 
msab <- ms2 <- as.matrix( sum.m.2unbal.CI[[1]][3,3] ); ms2
nu1 <- df.a <- as.matrix( sum.m.2unbal.CI[[1]][1,1] ); df.a 
nu2 <- df.ab <- as.matrix( sum.m.2unbal.CI[[1]][3,1] ); df.ab 
# confidence level conf1 = 0.90
confl <- 0.90
library(daewr)
vci(confl, c1, ms1, nu1, c2, ms2, nu2)

# 5.6 Nested Sampling Experiments

# 5.6.1 nested data
D <- expand.grid( lab = c('A', 'B', 'C'), # design
                  patient = c(1, 2, 3, 4) ) # combinations of factors
D[] <- lapply(D, factor) # convert to factors
D <- rbind(D, D) # row-bind twice, for 2 replicates
MRS.nest <- c(
  6.8, 8.2,  2.4,
  5.1,  7.2,  8.1,
  4.2,  4.1,  7.2,
  4.3,  4.1,  7.2,
  12.1, 6.5,  6.1,
  9.8,  8.1,  5.1,
  9.8,  8.1,  4.1,
  6.2,  4.9,  9.1
);
MRS.nestdata <- data.frame(cbind(D,MRS.nest)); MRS.nestdata # data.frame

# 5.6.2 patient nested in laboratory
m.nest1 <- aov( MRS.nest ~ lab + lab:patient, data = MRS.nestdata)
sm <- summary(m.nest1); sm # incorrect ANOVA with incorrect F ratios
msE <- as.matrix( sm[[1]][3,3] ); msE; dfE <- as.matrix( sm[[1]][3,1] ); dfE
msB <- as.matrix( sm[[1]][2,3] ); msB; dfB <- as.matrix( sm[[1]][2,1] ); dfB
msA <- as.matrix( sm[[1]][1,3] ); msA; dfA <- as.matrix( sm[[1]][1,1] ); dfA
F.B <- msB/msE; F.B; probF.B <- 1 - pf(F.B, dfB, dfE) 
F.A <- msA/msB; F.A; probF.A <- 1 - pf(F.A, dfA, dfB) 
cat("corrected ANOVA F statistics and p-values for NSE design","\n",
    "F.A =",F.A, " P(>F.A) =", probF.A,"\n",
    "F.B =",F.B, " P(>F.B) =", probF.B,"\n")
# 5.6.2 REML in estimates of variance components
library(lme4)
m.nest2 <- lmer( MRS.nest ~ 1 + (1|lab)  + (1|lab:patient), data = MRS.nestdata)
summary(m.nest2)

# 5.6.3 laboratory nested in patient
m.nest3 <- aov( MRS.nest ~ patient + patient:lab, data = MRS.nestdata)
sm <- summary(m.nest3); sm
msE <- as.matrix( sm[[1]][3,3] ); msE; dfE <- as.matrix( sm[[1]][3,1] ); dfE
msB <- as.matrix( sm[[1]][2,3] ); msB; dfB <- as.matrix( sm[[1]][2,1] ); dfB
msA <- as.matrix( sm[[1]][1,3] ); msA; dfA <- as.matrix( sm[[1]][1,1] ); dfA
F.B <- msB/msE; F.B; probF.B <- 1 - pf(F.B, dfB, dfE) 
F.A <- msA/msB; F.A; probF.A <- 1 - pf(F.A, dfA, dfB) 
cat("corrected ANOVA F statistics and p-values for NSE design","\n",
    "F.A =",F.A, " P(>F.A) =", probF.A,"\n",
    "F.B =",F.B, " P(>F.B) =", probF.B,"\n")
# REML in estimates of variance components
library(lme4)
m.nest4 <- lmer( MRS.nest ~ 1 + (1|patient)  + (1|patient:lab), data = MRS.nestdata)
summary(m.nest4)

# 5.6.2 three-stage random design
# mouse nested in litter nested in lab
D <- expand.grid( 
  repl = factor(c(1, 2, 3)), # replicates
  mouse = factor(c(1, 2, 3, 4)), # mouse factor 
  litter = factor(c(1,2,3,4)), # litter factor
  lab = factor(c("A","B")) # lab factor
) # combinations of factors
D[] <- lapply(D, factor); D # convert to factors
ROC <- c(
  1.4, 6.1, 3.4, 9.4, 8.1, 5.4, 7.4, 9.1, 8.4, 10.4, 16.1, 17.4,
  1.8, 5.1, 4.2, 4.8, 5.1, 7.2, 8.8, 5.1, 7.2, 14.8, 16.1, 19.2,
  1.2, 1.2, 4.1, 1.2, 4.2, 4.1, 1.2, 1.2, 6.1, 11.2, 15.2, 18.1,
  3.4, 5.1, 6.2, 12.4,8.1, 7.2,10.4, 8.1,17.2, 13.4, 18.1, 14.2,
  3.1, 4.2, 3.0, 4.1, 7.2, 7.0, 9.1, 7.7, 9.0, 11.1, 17.2, 16.0,
  2.1, 2.8, 3.8, 2.1, 9.8, 9.8,12.1, 9.8, 9.8, 12.1, 19.8, 17.8,
  1.5, 1.1, 1.1, 6.5, 4.1, 7.1, 6.5, 8.1, 7.1, 10.5, 10.1, 12.1,
  2.1, 5.1, 3.1, 3.1, 5.1, 7.1, 9.1, 8.1, 7.1, 16.1, 15.1, 14.1
); 
ROCdata <- data.frame(cbind(D,ROC)); ROCdata # data.frame

# 5.6.2 ANOVA F statistics and p-values for NSE design
m.nest3 <- aov( ROC ~ lab 
                + lab:litter 
                + lab:litter:mouse, data = ROCdata)
sm <- summary(m.nest3); sm
msE <- as.matrix( sm[[1]][4,3] ); msE; dfE <- as.matrix( sm[[1]][4,1] ); dfE
msC <- as.matrix( sm[[1]][3,3] ); msC; dfC <- as.matrix( sm[[1]][4,1] ); dfC
msB <- as.matrix( sm[[1]][2,3] ); msB; dfB <- as.matrix( sm[[1]][4,1] ); dfB
msA <- as.matrix( sm[[1]][1,3] ); msA; dfA <- as.matrix( sm[[1]][4,1] ); dfA
F.C <- msC/msE; F.C; probF.C <- 1 - pf(F.C, dfC, dfE) 
F.B <- msB/msC; F.B; probF.B <- 1 - pf(F.B, dfB, dfC) 
F.A <- msA/msB; F.A; probF.A <- 1 - pf(F.A, dfA, dfB) 
cat("corrected ANOVA F statistics and p-values for NSE design","\n",
    "F.A =",F.A, " P(>F.A) =", probF.A,"\n",
    "F.B =",F.B, " P(>F.B) =", probF.B,"\n",
    "F.C =",F.C, " P(>F.C) =", probF.C,"\n")

# 5.6.2 method of moments estimation of variance components
a <- 2; b <- 4; c <- 4; r <- 3
sigma2 <- msE <- as.matrix( sm[[1]][4,3] ); msE # extract from ANOVA
msC <- as.matrix( sm[[1]][3,3] ); msC
msB <- as.matrix( sm[[1]][2,3] ); msB
msA <- as.matrix( sm[[1]][1,3] ); msA
sigma2.C <- (msC - sigma2) / r # method of moments
sigma2.B <- (msB - sigma2 - r * sigma2.C ) / (r*c)
sigma2.A <- (msA - sigma2 - r * sigma2.C - r*c * sigma2.B) / (b*c*r)
cat("Method of Moments Variance Component Estimates","\n",
    "Var(lab)=",sigma2.A,"\n",
    "Var(lab:litter)=",sigma2.B,"\n",
    "Var(lab:litter:mouse)=",sigma2.C,"\n",
    "Var(error)=",sigma2,"\n")

# 5.6.2 REML estimation of variance components
library(lme4)
m.nest4 <- lmer( ROC ~ 1 + (1|lab)  
                 + (1|lab:litter) 
                 + (1|lab:litter:mouse), data = ROCdata)
summary(m.nest4)

# 5.7 Staggered Nested Designs

# 5.7.1 staggered nested data
patient <- factor(rep(1:12,4)); patient
lab <- factor(c(rep("A",36),rep("B",12))); lab
prep <- factor(c(rep(1,24),rep(2,12),rep(1,12))); prep
test <- factor(c(rep("alpha",12),rep("beta",12),rep("alpha",24))); test
PTSD3 <- c(
  6.8, 12.1, 5.1,  9.8,  4.2,  9.8,  4.3,  6.2, 14.2,  9.3, 10.3, 16.2,
  8.2,  6.5, 7.2,  8.1,  4.1,  8.1,  4.1,  4.9,  3.1,  4.1,  4.6,  3.9,
  2.4,  6.1, 8.1,  5.1,  7.2,  4.1,  7.2,  9.1,  8.2,  7.1,  5.2,  8.1,
  5.6,  4.5, 5.1,  6.3,  8.1,  9.0,  8.3,  8.3,  7.7,  8.6,  4.4,  6.5
);
PTSD3 <- as.numeric(PTSD3); PTSD3
PTSD3data <- data.frame(patient,lab,prep,test,PTSD3); PTSD3data # data.frame

# 5.7.2 staggered nested, method of moments
m.stag.nest <- aov(PTSD3 ~ patient + patient:lab + patient:lab:prep, data = PTSD3data)
summary.stag.nest <- summary(m.stag.nest); summary.stag.nest
sigma2 <- msE <- as.matrix( summary.stag.nest[[1]][4,3] ); msE # extract from ANOVA
msC <- as.matrix( summary.stag.nest[[1]][3,3] ); msC 
msB <- as.matrix( summary.stag.nest[[1]][2,3] ); msB 
msA <- as.matrix( summary.stag.nest[[1]][1,3] ); msA 
sigma2.C <- (msC - sigma2) / (4/3) # method of moments
sigma2.B <- (msB - sigma2 - (7/6) * sigma2.C ) / (3/2)
sigma2.A <- (msA - sigma2 - (3/2) * sigma2.C - (5/2) * sigma2.B) / 4
cat("Method of Moments Variance Component Estimates","\n",
    "Var(patient)=",sigma2.A,"\n",
    "Var(patient:lab)=",sigma2.B,"\n",
    "Var(patient:lab:prep)=",sigma2.C,"\n",
    "Var(error)=",sigma2,"\n")

# 5.7.3 staggered nested, REML
library(lme4)
m.REML <- lmer(PTSD3 ~ 1 + (1|patient) 
               + (1|patient:lab) 
               + (1|patient:lab:prep), data = PTSD3data)
summary(m.REML)

# 5.8 Designs with Fixed and Random Effects

# 5.8.1 random nested inside two fixed, data
Des <- expand.grid( patient = c(1, 2), # 2^3 design
                    prep = c(1, 2), # 2 x 2 x 2  = 8
                    test = c("alpha", "beta")) # 8 x 2 replicates = 16
Des[] <- lapply(Des, factor) # convert to factors
PTSD <- c(1.8613, 0.8753, 1.8240, 1.0818, 
       0.6843, 1.3462, 3.1200, 1.0620,
       1.6754, 1.1283, 1.4486, 4.1304, 
       0.8323, 1.5433, 4.0100, 1.1111); PTSD
PTSDdata <- data.frame(cbind(Des,PTSD)); PTSDdata # data.frame
# 5.8.1 fit model
library(lme4)
c1 <- c(-0.5,0.5) # pairwise contrast
m.rand.nested.2fix <- lmer(
  PTSD ~ 1 + prep + test + prep:test + (1|patient:prep:test),
  contrasts = list(prep = c1, test = c1), data = PTSDdata)
summary(m.rand.nested.2fix)
# pairwise comparisons
library(lmsmeans)
lsmeans(m.rand.nested.2fix, pairwise ~ test, adjust = c("tukey"))
lsmeans(m.rand.nested.2fix, pairwise ~ prep, adjust = c("tukey"))
# F-tests for fixed factors
anova(m.rand.nested.2fix) # F tests for fixed factors
1 - pf(1.8880,1,4); 1 - pf(0.0031,1,4); 1 - pf(0.1171,1,4) # p-values for F tests
# fit (incorrect F statistics) ANOVA model
library(daewr)
m.anova <- aov(
  PTSD ~ prep + test + prep:test + patient:prep:test,
  data = PTSDdata)
summary(m.anova)

# 5.8.2 mixed: random and fixed (versus factor, block)
D <- expand.grid( 
  repl = c(1, 2, 3), # replicates
  temp = c("0", "10", "20", "30"), # temperature factor 
  mouse = c(1, 2, 3, 4, 5, 6, 7, 8) # block
) # combinations of factors
D[] <- lapply(D, factor); D # convert to factors
ROC <- c(
  1.4, 6.1, 3.4, 9.4, 8.1, 5.4, 7.4, 9.1, 8.4, 10.4, 16.1, 17.4,
  1.8, 5.1, 4.2, 4.8, 5.1, 7.2, 8.8, 5.1, 7.2, 14.8, 16.1, 19.2,
  1.2, 1.2, 4.1, 1.2, 4.2, 4.1, 1.2, 1.2, 6.1, 11.2, 15.2, 18.1,
  3.4, 5.1, 6.2, 12.4,8.1, 7.2,10.4, 8.1,17.2, 13.4, 18.1, 14.2,
  3.1, 4.2, 3.0, 4.1, 7.2, 7.0, 9.1, 7.7, 9.0, 11.1, 17.2, 16.0,
  2.1, 2.8, 3.8, 2.1, 9.8, 9.8,12.1, 9.8, 9.8, 12.1, 19.8, 17.8,
  1.5, 1.1, 1.1, 6.5, 4.1, 7.1, 6.5, 8.1, 7.1, 10.5, 10.1, 12.1,
  2.1, 5.1, 3.1, 3.1, 5.1, 7.1, 9.1, 8.1, 7.1, 16.1, 15.1, 14.1
); 
ROCdata3 <- data.frame(cbind(D,ROC)); ROCdata3 # data.frame
# analysis of variance
library(daewr) # RCB design with replicates
m.fix.rand <- aov(ROC ~ temp + Error(mouse/temp), data = ROCdata3)
summary(m.fix.rand)


# 5.9 Graphical Methods to Check Model Assumptions

# 5.9.1 residual plots, one-factor
D <- expand.grid( lab = c('A', 'B', 'C', 'D'), # labs
                  patient = c(1, 2, 3, 4,
                              5, 6, 7, 8) ) # combinations of factors
D[] <- lapply(D, factor) # convert to factors
PTSD <- c(
  10.3, 6.8,  8.2,  2.4,
  9.1, 12.1,  6.5,  6.1,
  6.1, 5.1,  7.2,  8.1,
  7.2, 9.8,  8.1,  5.1,
  5.4, 4.2,  4.1,  7.2,
  7.0, 9.8,  8.1,  4.1,
  5.4, 4.3,  4.1,  7.2,
  6.2, 6.2,  4.9,  9.1); 
d <- PTSDdata <- data.frame(cbind(D,PTSD)); PTSDdata # data.frame
m1 <- aov(PTSD ~ lab, PTSDdata); m1

par( mfrow = c(2,2) )
plot( m1, which = 5 ) # should be constant stnd resid for different factor levels
plot( m1, which = 1 ) # should be constant for different fitted values
plot( m1, which = 2 ) # normal prob plot for normality should be linear
plot( residuals(m1) ~ patient, 
      main = "Residuals vs Exp. Unit", 
      font.main = 1, 
      data = PTSDdata)
abline( h=0, lty = 2 ) # should be constant
par(mfrow = c(1,1))

# 5.9.2 gamma and half-normal plot checks constant variance
# one-factor 
D <- expand.grid( lab = c('A', 'B', 'C', 'D'), # labs
                  patient = c(1, 2, 3, 4,
                              5, 6, 7, 8) ) # combinations of factors
D[] <- lapply(D, factor) # convert to factors
PTSD <- c(
  10.3, 6.8,  8.2,  2.4,
  9.1, 12.1,  6.5,  6.1,
  6.1, 5.1,  7.2,  8.1,
  7.2, 9.8,  8.1,  5.1,
  5.4, 4.2,  4.1,  7.2,
  7.0, 9.8,  8.1,  4.1,
  5.4, 4.3,  4.1,  7.2,
  6.2, 6.2,  4.9,  9.1); 
d <- PTSDdata <- data.frame(cbind(D,PTSD)); PTSDdata # data.frame
m1 <- aov(PTSD ~ lab, PTSDdata); m1
# gamma plot of lab variances shows two large SDs, but not which patients
s2 <- tapply( PTSDdata$PTSD, PTSDdata$lab, var ); s2
os2 <- sort(s2); r <- c( 1:length(s2) ) # sort, count s2
gscore <- qgamma( (r - .5 ) / length (s2), 2)
plot(gscore, 
     os2, 
     main = "Gamma plot of within lab variances", 
     xlab = "Gamma score", 
     ylab = "Lab Variance")

m.2bal <- aov( PTSD2 ~ patient + lab + patient:lab, data = PTSD2data)
sum.m.2bal <- summary(m.2bal) 

# 5.9.3 interaction plot, two-factors, 
# reducing factor levels to reduce interactions
D <- expand.grid( lab = c('A', 'B', 'C', 'D'), # design
                  patient = c(1, 2, 3, 4) ) # combinations of factors
D[] <- lapply(D, factor) # convert to factors
D <- rbind(D, D) # row-bind twice, for 2 replicates
PTSD2 <- c(
  10.3, 6.8, 8.2,  2.4,
  6.1, 5.1,  7.2,  8.1,
  5.4, 4.2,  4.1,  7.2,
  5.4, 4.3,  4.1,  7.2,
  9.1, 12.1, 6.5,  6.1,
  7.2, 9.8,  8.1,  5.1,
  7.0, 9.8,  8.1,  4.1,
  6.2, 6.2,  4.9,  9.1
);
PTSD2data <- data.frame(cbind(D,PTSD2)); PTSD2data # data.frame

# 5.9.3 VAR plot of two-factor design
library(VCA) 
varPlot(form=PTSD2 ~ patient + lab + patient:lab, Data = PTSD2data)

# 5.9.3 interaction plots for two-factor design
par( mfrow = c(1,2) )
with(PTSD2data, # 
     (interaction.plot(
       patient, 
       lab, 
       PTSD2, 
       type = "b", 
       pch = c(18,24,22,21), 
       leg.bty = "o",
       main = "GABA for each lab and patient",
       xlab = "Patient",
       ylab = "GABA"))
)
with(PTSD2data, # 
     (interaction.plot(
       lab, 
       patient, 
       PTSD2, 
       type = "b", 
       pch = c(18,24,22,21), 
       leg.bty = "o",
       main = "GABA for each patient and lab",
       xlab = "Lab",
       ylab = "GABA"))
)
par( mfrow = c(1,1) )

# 5.9.3 remove patient 4: reduce a = 4 to a = 3
D <- expand.grid( lab = c('A', 'B', 'C', 'D'), # design
                  patient = c(1, 2, 3) ) # combinations of factors
D[] <- lapply(D, factor) # convert to factors
D <- rbind(D, D) # row-bind twice, for 2 replicates
PTSD2 <- c(
  10.3, 6.8, 8.2,  2.4,
  6.1, 5.1,  7.2,  8.1,
  5.4, 4.2,  4.1,  7.2,
  5.4, 4.3,  4.1,  7.2,
  9.1, 12.1, 6.5,  6.1,
  7.2, 9.8,  8.1,  5.1
);
PTSD2data <- data.frame(cbind(D,PTSD2)); PTSD2data # data.frame
m.2bal <- aov( PTSD2 ~ patient + lab + patient:lab, data = PTSD2data)
sum.m.2bal <- summary(m.2bal) 
# ANOVA summary for method of moments analysis
r <- 2; a <- 3; b <- 4 # 2 replications per cell; number of a and b
sigma2 <- msE <- as.matrix( sum.m.2bal[[1]][4,3] ); msE # extract from ANOVA
msAB <- as.matrix( sum.m.2bal[[1]][3,3] ); msAB 
msB <- as.matrix( sum.m.2bal[[1]][2,3] ); msB 
msA <- as.matrix( sum.m.2bal[[1]][1,3] ); msA 
sigma2patientlab <- (msAB - sigma2) / r # method of moments
sigma2lab <- (msB - sigma2 - r * msAB ) / (r*b)
sigma2patient <- (msA - sigma2 - r * msAB ) / (r*a)
cat("Method of Moments Variance Component Estimates","\n",
    "Var(patient)=",sigma2patient,"\n",
    "Var(lab)=",sigma2lab,"\n",
    "Var(patient x lab)=",sigma2patientlab,"\n",
    "Var(error)=",sigma2,"\n")
# ANOVA summary for REML analysis 
library(lme4)
m.reml.2bal <- lmer(PTSD2 ~ 1 
                    + (1|patient) + (1|lab) 
                    + (1|patient:lab), data = PTSD2data)
summary(m.reml.2bal)

# 5.9.3  remove lab D: reduce b = 4 to b = 3
D <- expand.grid( lab = c('A', 'B', 'C'), # design
                  patient = c(1, 2, 3, 4) ) # combinations of factors
D[] <- lapply(D, factor) # convert to factors
D <- rbind(D, D) # row-bind twice, for 2 replicates
PTSD2 <- c(
  10.3, 6.8, 8.2,
  6.1, 5.1,  7.2,
  5.4, 4.2,  4.1,
  5.4, 4.3,  4.1,
  9.1, 12.1, 6.5,
  7.2, 9.8,  8.1,
  7.0, 9.8,  8.1,
  6.2, 6.2,  4.9
);
PTSD2data <- data.frame(cbind(D,PTSD2)); PTSD2data # data.frame
m.2bal <- aov( PTSD2 ~ patient + lab + patient:lab, data = PTSD2data)
sum.m.2bal <- summary(m.2bal) 
# ANOVA summary for method of moments analysis
r <- 2; a <- 4; b <- 3 # 2 replications per cell; number of a and b
sigma2 <- msE <- as.matrix( sum.m.2bal[[1]][4,3] ); msE # extract from ANOVA
msAB <- as.matrix( sum.m.2bal[[1]][3,3] ); msAB 
msB <- as.matrix( sum.m.2bal[[1]][2,3] ); msB 
msA <- as.matrix( sum.m.2bal[[1]][1,3] ); msA 
sigma2patientlab <- (msAB - sigma2) / r # method of moments
sigma2lab <- (msB - sigma2 - r * msAB ) / (r*b)
sigma2patient <- (msA - sigma2 - r * msAB ) / (r*a)
cat("Method of Moments Variance Component Estimates","\n",
    "Var(patient)=",sigma2patient,"\n",
    "Var(lab)=",sigma2lab,"\n",
    "Var(patient x lab)=",sigma2patientlab,"\n",
    "Var(error)=",sigma2,"\n")
# ANOVA summary for REML analysis 
library(lme4)
m.reml.2bal <- lmer(PTSD2 ~ 1 
                    + (1|patient) + (1|lab) 
                    + (1|patient:lab), data = PTSD2data)
summary(m.reml.2bal)

# 5.9.4 half-normal plot to identify patients with large lab variances
# for staggered nested design
patient <- factor(rep(1:12,4)); patient
lab <- factor(c(rep("A",36),rep("B",12))); lab
prep <- factor(c(rep(1,24),rep(2,12),rep(1,12))); prep
test <- factor(c(rep("alpha",12),rep("beta",12),rep("alpha",24))); test
PTSD3 <- c(
  6.8, 12.1, 5.1,  9.8,  4.2,  9.8,  4.3,  6.2, 14.2,  9.3, 10.3, 16.2,
  8.2,  6.5, 7.2,  8.1,  4.1,  8.1,  4.1,  4.9,  3.1,  4.1,  4.6,  3.9,
  2.4,  6.1, 8.1,  5.1,  7.2,  4.1,  7.2,  9.1,  8.2,  7.1,  5.2,  8.1,
  5.6,  4.5, 5.1,  6.3,  8.1,  9.0,  8.3,  8.3,  7.7,  8.6,  4.4,  6.5
);
PTSD3 <- as.numeric(PTSD3); PTSD3
PTSD3data <- data.frame(patient,lab,prep,test,PTSD3); PTSD3data # data.frame

# 5.9.4 VAR plot of staggered nested design
library(VCA) 
varPlot(form=PTSD3 ~ patient/lab/prep, Data = PTSD3data)

# order data.frame by patient
PTSD3data.ordered <- PTSD3data[order(PTSD3data$patient),]; PTSD3data.ordered 
# convert into 4 x 12 array of PTSD3 values
y <- array( PTSD3data.ordered$PTSD3, c(4,12) ); y 
# calculate SDs of each source
sd1 <- sqrt( (y[2,] - y[1,])**2 / 2) # SDs test(prep)
sd2 <- sqrt( (2/3) * ( y[3,] - (y[1,] + y[2,]) / 2)**2 ) # SDs prep(lab)
sd3 <- sqrt( (3/4) * (y[4,] - (y[1,] + y[2,] + y[3,] )/3 )**2) # SDs lab
# ordered SDs prep(lab) compared to z-score of ordered SDs prep(lab)
osd2 <- sort(sd2); osd2
r <- c( 1: length(sd2))
zscore <- qnorm( ( ( r - .5 ) / length(sd2) +1 )/ 2); zscore
plot( zscore, osd2, 
      main = "Half-normal plot of prep(lab) standard deviations", 
      xlab = "Half Normal Score", 
      ylab ="std. due to prep within lab")
# pretty table of data with SDs
library(formattable)
patient <- 1:12; patient
ty <- t(y); ty
colnames(ty) <- c('Y1', 'Y2', 'Y3', 'Y4') # naming columns
data.plus.sds <- data.frame(patient,ty,sd1,sd2,sd3); data.plus.sds # data.frame
formattable(data.plus.sds) # patient 1 and 6 have large SDs

# 5.9.4 remove patients 1 and 6, redo staggered nested design
patient <- factor(rep(c(2:5,7:12),4)); patient
lab <- factor(c(rep("A",30),rep("B",10))); lab
prep <- factor(c(rep(1,20),rep(2,10),rep(1,10))); prep
test <- factor(c(rep("alpha",10),rep("beta",10),rep("alpha",20))); test
PTSD3 <- c(
 12.1, 5.1,  9.8,  4.2,  4.3,  6.2, 14.2,  9.3, 10.3, 16.2,
  6.5, 7.2,  8.1,  4.1,  4.1,  4.9,  3.1,  4.1,  4.6,  3.9,
  6.1, 8.1,  5.1,  7.2,  7.2,  9.1,  8.2,  7.1,  5.2,  8.1,
  4.5, 5.1,  6.3,  8.1,  8.3,  8.3,  7.7,  8.6,  4.4,  6.5
);
PTSD3 <- as.numeric(PTSD3); PTSD3
PTSD3data <- data.frame(patient,lab,prep,test,PTSD3); PTSD3data # data.frame
m.stag.nest <- aov(PTSD3 ~ patient + patient:lab + patient:lab:prep, data = PTSD3data)
summary.stag.nest <- summary(m.stag.nest); summary.stag.nest

# 5.9.4 method of moments
sigma2 <- msE <- as.matrix( sum.m.2bal[[1]][4,3] ); msE # extract from ANOVA
msC <- as.matrix( sum.m.2bal[[1]][3,3] ); msC 
msB <- as.matrix( sum.m.2bal[[1]][2,3] ); msB 
msA <- as.matrix( sum.m.2bal[[1]][1,3] ); msA 
sigma2.C <- (msC - sigma2) / (4/3) # method of moments
sigma2.B <- (msB - sigma2 - (7/6) * sigma2.C ) / (3/2)
sigma2.A <- (msA - sigma2 - (3/2) * sigma2.C - (5/2) * sigma2.B) / 4
cat("Method of Moments Variance Component Estimates","\n",
    "Var(patient)=",sigma2.A,"\n",
    "Var(patient:lab)=",sigma2.B,"\n",
    "Var(patient:lab:prep)=",sigma2.C,"\n",
    "Var(error)=",sigma2,"\n")

# 5.9.4 REML
library(lme4)
m.REML <- lmer(PTSD3 ~ 1 + (1|patient) 
               + (1|patient:lab) 
               + (1|patient:lab:prep), data = PTSD3data)
summary(m.REML)

# 5.9.5 probability plots of BLUE random effects
# for staggered nested design
patient <- factor(rep(1:12,4)); patient
lab <- factor(c(rep("A",36),rep("B",12))); lab
prep <- factor(c(rep(1,24),rep(2,12),rep(1,12))); prep
test <- factor(c(rep("alpha",12),rep("beta",12),rep("alpha",24))); test
PTSD3 <- c(
  6.8, 12.1, 5.1,  9.8,  4.2,  9.8,  4.3,  6.2, 14.2,  9.3, 10.3, 16.2,
  8.2,  6.5, 7.2,  8.1,  4.1,  8.1,  4.1,  4.9,  3.1,  4.1,  4.6,  3.9,
  2.4,  6.1, 8.1,  5.1,  7.2,  4.1,  7.2,  9.1,  8.2,  7.1,  5.2,  8.1,
  5.6,  4.5, 5.1,  6.3,  8.1,  9.0,  8.3,  8.3,  7.7,  8.6,  4.4,  6.5
);
PTSD3 <- as.numeric(PTSD3); PTSD3
PTSD3data <- data.frame(patient,lab,prep,test,PTSD3); PTSD3data # data.frame
library(lme4)
m.stag.nest <- lmer(PTSD3 ~ 1 + (1|patient) 
                    + (1|patient:lab) 
                    + (1|patient:lab:prep), data = PTSD3data)
summary(m.stag.nest)
qqnorm( ranef(m.stag.nest)$"patient:lab:prep"[[1]], 
        main = "prep within lab and patient", 
        ylab="EBLUP",
        xlab = "Normal Score" )

###################################################
# Homework 4 and Quiz 2 questions
###################################################

# 5.1.a pp 190-192
y <- flight.time <- c(
  5.2, 5.1, 5.3, 5.3, 5.2, 5.2, 5.0, 5.2,
  5.2, 5.1, 5.1, 5.0, 5.2, 5.1, 5.1, 5.2,
  5.1, 5.1, 5.0, 5.0, 5.1, 5.0, 4.9, 5.1,
  4.8, 4.9, 4.8, 5.0, 4.8, 4.9, 4.9, 4.9
); y; class(y)
x <- wing.length <- factor(c(rep(4,8), rep(4.75, 8), rep(5.5,8), rep(6,8))); x; class(x)
helicopter <- as.numeric(1:32); helicopter; class(helicopter)
paper.heli <- data.frame(helicopter, wing.length, flight.time); paper.heli

m1 <- lm( flight.time ~ wing.length, data = paper.heli ) # model 1, m1
summary1 <- anova(m1); s1

# 5.2.a,b,c pp 190-192
D <- expand.grid( sample = c(1, 2, 3, 4, 5, 6), # samples
                  prep = c(1, 2, 3, 4, 5) ) # combinations of factors
D[] <- lapply(D, factor) # convert to factors
yield <- c(
1440, 1490, 1510, 1440, 1515, 1445,
1440, 1495, 1550, 1445, 1595, 1450,
1520, 1540, 1560, 1465, 1625, 1455,
1545, 1555, 1595, 1545, 1630, 1480,
1580, 1560, 1605, 1595, 1635, 1520); 
yielddata <- data.frame(cbind(D,yield)); yielddata # data.frame
m1 <- aov( yield ~ sample, data = yielddata ) # model 1, m1
s1 <- summary(m1); s1

# 5.2.d,e,f pp 190-192
c <- r <- 5; c # coefficient = # replications 
sigma.2 <- as.matrix( s1[[1]][2,3] )
ms.lab <- as.matrix( s1[[1]][1,3] )
cat(" Mean Square for Lab = ",ms.lab,"\n",
    " Mean Square for Error = ", sigma.2,"\n",
    " Expected Mean Square for Lab","\n", 
    "Var(error)+",c,"Var(Lab)","\n")
# compare expected to observed: 
sigma.2t <- (ms.lab - sigma.2) / c
cat("Method of Moments Variance Component Estimates","\n", 
    "Var(error)=",sigma.2,"\n",
    "Var(Lab)=",sigma.2t,"\n")

# 5.2.g pp 190-192
library(lme4)
rm <- lmer( yield ~ 1 + (1| sample), data = yielddata) # random model rm
summary(rm) # # splits into sigma^2_residual, sigma^2_lab

# 5.2.h pp 190-192
yield <- c(
  1440, 1440, 1520, 1545, 1580,
  1490, 1495, 1540, 1555,
  1510, 1550, 1560, 1595,
  1440, 1445, 1465, 1545, 1595,
  1515, 1595, 1625, 1630,
  1445, 1450, 1455, 1480, 1520); 
sample <- factor( c( rep ('1', 5), rep ('2', 4), rep('3',4),
                  rep ('4', 5), rep ('5', 4), rep('6', 5) ) ); sample
yielddata <- data.frame(yield, sample); yielddata
library(lme4)
rm <- lmer( yield ~ 1 + (1| sample), data = yielddata) # random model rm
summary(rm) # # splits into sigma^2_residual, sigma^2_lab

# 5.2.i pp 190-192 
yield <- c(
  1440, 1440, 1520, 1545, 1580,
  1490, 1495, 1540, 1555,
  1510, 1550, 1560, 1595,
  1440, 1445, 1465, 1545, 1595,
  1515, 1595, 1625, 1630,
  1445, 1450, 1455, 1480, 1520); 
sample <- factor( c( rep ('1', 5), rep ('2', 4), rep('3', 4),
                     rep ('4', 5), rep ('5', 4), rep('6', 5) ) ); sample
yielddata <- data.frame(yield, sample); yielddata
library(lme4)
likel.profile <- profile( # using likelihood profile method
  rm.likel <- lmer( yield ~ 1 + (1| sample), # rm.likel can be anything,
                    data = yielddata,     #  is just a dummary variable 
                    REML = FALSE))
confint (likel.profile) # 95% CI

# 5.4.c pp 190-192 
# determine nu2 = (r - 1)ab given CI = 75% of CI of sigma^2
alpha = 0.05
nu2 <- 50:65
chiu <- qchisq(1 - alpha/2, nu2)
chil <- qchisq(alpha/2, nu2)
width <- nu2 * (chiu - chil) / (chil * chiu)
halfw <- width/2
data.frame(nu2, chiu, chil, width)

# 5.5 pp 190-192 
# three-stage nested design
# sample nested in batch nested in supllier
library(daewr)
data(rubber); rubber # rubber data from daewr
m <- aov( elasticity ~ supplier 
                + batch:supplier 
                + sample:batch:supplier, data = rubber)
summary(m)
library(lme4)
m <- lmer( elasticity ~ 1 + (1|supplier)  
                 + (1|batch:supplier) 
                 + (1|sample:batch:supplier), data = rubber)
summary(m)
library(VCA)
varPlot(form=elasticity~supplier/batch/sample, Data=rubber)

# 5.6.a.1 pp 190-192
D <- expand.grid( 
  town = c('1', '2'),
  state = c('1', '2', '3'), # design
  household = c('1', '2', '3', '4') ) # combos
D[] <- lapply(D, factor); D # convert to factors
health.nest <- c(
  10,  7, 6,  6, 15,  12,
  13, 12, 5, 12, 18,  15,
  16, 11, 9,  7, 20,  18,
  12,  9, 3, 10, 19,  16
);
health.nestdata <- data.frame(cbind(D, health.nest)); health.nestdata # data.frame
library(lme4)
m.REML <- lmer(health.nest ~ 1 + (1|state) 
               + (1|state:town), data = health.nestdata)
summary(m.REML)

# 5.6.a.2 pp 190-192
D <- expand.grid( 
  town = c('1', '2'),
  state = c('1', '2', '3'), # design
  household = c('1', '2', '3', '4') ) # combos
D[] <- lapply(D, factor); D # convert to factors
health.nest <- c(
  10,  7, 6,  6, 15,  12,
  13, 12, 5, 12, 18,  15,
  16, 11, 9,  7, 20,  18,
  12,  9, 3, 10, 19,  16
);
health.nestdata <- data.frame(cbind(D, health.nest)); health.nestdata # data.frame
m.nest1 <- aov( 
  health.nest ~ state + state:town, 
  data = health.nestdata); m.nest1
sum.m.nest1 <- summary(m.nest1); sum.m.nest1
a <- 3; b <- 2; r <- 4
sigma2 <- msE <- as.matrix( sum.m.nest1[[1]][3,3] ); msE # extract from ANOVA
msB <- as.matrix( sum.m.nest1[[1]][2,3] ); msB 
msA <- as.matrix( sum.m.nest1[[1]][1,3] ); msA 
sigma2.B <- (msB - sigma2 ) / r
sigma2.A <- (msA - sigma2 - r * sigma2.B) / (b*r)
cat("Method of Moments Variance Component Estimates","\n",
    "Var(state)=",sigma2.A,"\n",
    "Var(state:town)=",sigma2.B,"\n",
    "Var(error)=",sigma2,"\n")

# 5.6.b.1 pp 190-192
# use results from either 5.6.a.1 or 5.6.a.2

# 5.6.b.2 pp 190-192
# compare 5.6.a.1 to 5.6.a.2

# 5.6.c pp 190-192
D <- expand.grid( 
  town = c('1', '2'),
  state = c('1', '2', '3'), # design
  household = c('1', '2', '3', '4') ) # combos
D[] <- lapply(D, factor); D # convert to factors
health.nest <- c(
  10,  7, 6,  6, 15,  12,
  13, 12, 5, 12, 18,  15,
  16, 11, 9,  7, 20,  18,
  12,  9, 3, 10, 19,  16
);
health.nestdata <- data.frame(cbind(D, health.nest)); health.nestdata # data.frame
m.nest1 <- aov( 
  health.nest ~ state + state:town, 
  data = health.nestdata); m.nest1
sum.m.nest1 <- summary(m.nest1); sum.m.nest1
a <- 3; b <- 2; r <- 4
library(daewr)
options(digits = 3)
# confidence level = 90% for sigma^2_b, and from ANOVA table,
# c1 = 1/r, ms1 = msb = 17.04, nu1 = 3,
# c2 = 1/r, ms2 = mse = 5.99, nu2 = 18
vci(confl = .90, 
    c1 = 1/r, ms1 = 17.04, nu1 = 3, 
    c2 = 1/r, ms2 = 5.99, nu2 = 18)

# 5.8.a,b
library(daewr)
data(rcb)
library(lme4)
modg<-lmer(cdistance~ 1 +(1|id) + teehgt +(1|id:teehgt), data = rcb)
summary(modg)
anova(modg)
library(lsmeans)
require(pbkrtest)
lsmeans(modg, pairwise ~ teehgt, adjust=("tukey"))


