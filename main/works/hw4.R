###################################################
# Homework 4 and Quiz 2 questions
###################################################

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
summary1 <- anova(m1); summary1



library(lme4)
heli <- lmer(flight.time ~ 1 
                    + (1|wing.length) , data = paper.heli)
summary(heli)



# 5.3.1 likelihood profile approximate CIs: balanced case 
likel.profile <- profile( # using likelihood profile method
  rm.likel <- lmer( flight.time ~ 1 + (1| wing.length), # rm.likel can be anything,
                    data = paper.heli,     #  is just a dummary variable 
                    REML = FALSE))
confint (likel.profile, level = 0.90) # 95% CI


################


# 5.6.2 method of moments estimation of variance components
sm <- summary(aov(m1))
sm
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


##################


# 5.5.1 method of moments
PTSD2data
m.2bal <- aov( PTSD2 ~ patient + lab + patient:lab, data = PTSD2data)
sum.m.2bal <- summary(m.2bal) # ANOVA summary for method of moments analysis
sum.m.2bal
# using Bennett and Franklin
r <- 2; a <- 4; b <- 4 # 2 replications per cell; number of a and b
sigma2 <- msE <- as.matrix( sum.m.2bal[[1]][4,3] ); msE # extract from ANOVA
msAB <- as.matrix( sum.m.2bal[[1]][3,3] ); msAB 
msB <- as.matrix( sum.m.2bal[[1]][2,3] ); msB 
msA <- as.matrix( sum.m.2bal[[1]][1,3] ); msA 
sigma2patientlab <- (msAB - sigma2) / r; sigma2patientlab # method of moments
sigma2patient <- (msB - sigma2 - r * sigma2patientlab ) / (r*a); sigma2lab
sigma2lab <- (msA - sigma2 - r * sigma2patientlab ) / (r*b)
cat("Method of Moments Variance Component Estimates","\n",
    "Var(lab)=          ",sigma2lab,"\n",
    "Var(patient)=      ",sigma2patient,"\n",
    "Var(lab x patient)=",sigma2patientlab,"\n",
    "Var(error)=        ",sigma2,"\n")

# balanced, REML
library(lme4)
m.reml.2bal <- lmer(PTSD2 ~ 1 
                    + (1|lab) + (1|patient)  
                    + (1|lab:patient), data = PTSD2data)
summary(m.reml.2bal)


#############################
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



# 5.6.2 method of moments estimation of variance components

m <- aov( elasticity ~ supplier 
          + batch:supplier 
          + sample:batch:supplier, data = rubber)

sm = summary(m)
a <- 4; b <- 4; c <- 2; r <- 3
sigma2 <- msE <- as.matrix( sm[[1]][4,3] ); msE # extract from ANOVA
msC <- as.matrix( sm[[1]][3,3] ); msC
msB <- as.matrix( sm[[1]][2,3] ); msB
msA <- as.matrix( sm[[1]][1,3] ); msA
sigma2.C <- (msC - sigma2) / r # method of moments
sigma2.B <- (msB - sigma2 - r * sigma2.C ) / (r*c)
sigma2.A <- (msA - sigma2 - r * sigma2.C - r*c * sigma2.B) / (b*c*r)
cat("Method of Moments Variance Component Estimates","\n",
    "Var(supplier)=",sigma2.A,"\n",
    "Var(batch:supplier)=",sigma2.B,"\n",
    "Var(sample:batch:supplier)=",sigma2.C,"\n",
    "Var(error)=",sigma2,"\n")


################################################
# QUIZ 2 code
####################################################

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


library(lme4)
naphthalene <- lmer(yield ~ 1 
             + (1|sample) , data = yielddata)
summary(naphthalene)



# 5.4.c pp 190-192 
# determine nu2 = (r - 1)ab given CI = 75% of CI of sigma^2
alpha = 0.05
nu2 <- 30:80
chiu <- qchisq(1 - alpha/2, nu2)
chil <- qchisq(alpha/2, nu2)
width <- nu2 * (chiu - chil) / (chil * chiu)
halfw <- width/2
data.frame(nu2, chiu, chil, width)



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
